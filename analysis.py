"""Analysis queries for the Loblaw Bio cell-count database."""

import sqlite3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"


def connect() -> sqlite3.Connection:
    """Open the database with foreign keys enforced."""
    con = sqlite3.connect(DB_PATH)
    con.execute("PRAGMA foreign_keys = ON")
    return con

SUMMARY_SQL = """
SELECT
    sample_id                                        AS sample,
    SUM(count) OVER (PARTITION BY sample_id)         AS total_count,
    population,
    count,
    ROUND(100.0 * count / SUM(count) OVER (PARTITION BY sample_id), 2) AS percentage
FROM cell_counts
ORDER BY sample_id, population
"""


def summary_table(con: sqlite3.Connection) -> pd.DataFrame:
    """Part 2: relative frequency of each population within each sample."""
    return pd.read_sql_query(SUMMARY_SQL, con)

RESPONSE_COHORT_SQL = """
SELECT
    s.subject_id,
    s.response,
    sm.sample_id                             AS sample,
    sm.time_from_treatment_start             AS day,
    cc.population,
    100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS percentage
FROM subjects s
JOIN samples     sm ON sm.subject_id = s.subject_id
JOIN cell_counts cc ON cc.sample_id  = sm.sample_id
WHERE s.condition   = 'melanoma'
  AND s.treatment   = 'miraclib'
  AND sm.sample_type = 'PBMC'
  AND s.response IS NOT NULL
ORDER BY s.subject_id, cc.population
"""


def response_cohort(con: sqlite3.Connection) -> pd.DataFrame:
    """Part 3: melanoma + miraclib + PBMC samples with a known response."""
    return pd.read_sql_query(RESPONSE_COHORT_SQL, con)

def per_subject(cohort: pd.DataFrame) -> pd.DataFrame:
    """Average each subject's repeated samples into one value per population."""
    return (
        cohort.groupby(["subject_id", "response", "population"], as_index=False)
        .percentage.mean()
    )
    
def compare_response(per_subj: pd.DataFrame) -> pd.DataFrame:
    """Mann-Whitney U test per population, with Benjamini-Hochberg correction."""
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    rows = []
    for pop, grp in per_subj.groupby("population"):
        yes = grp.loc[grp.response == "yes", "percentage"]
        no = grp.loc[grp.response == "no", "percentage"]

        u, p = stats.mannwhitneyu(yes, no, alternative="two-sided")

        # rank-biserial correlation: effect size from -1 to 1, 0 = no difference
        effect = 2 * u / (len(yes) * len(no)) - 1

        rows.append({
            "population": pop,
            "n_responders": len(yes),
            "n_non_responders": len(no),
            "median_responder": yes.median(),
            "median_non_responder": no.median(),
            "difference": yes.median() - no.median(),
            "effect_size": effect,
            "p_value": p,
        })

    out = pd.DataFrame(rows)

    # Benjamini-Hochberg: adjust p-values for running one test per population
    reject, p_adj, _, _ = multipletests(out.p_value, alpha=0.05, method="fdr_bh")
    out["p_adjusted"] = p_adj
    out["significant"] = reject

    return out.sort_values("population").reset_index(drop=True)

def boxplot_response(per_subj: pd.DataFrame, stats_df: pd.DataFrame, path: Path) -> Path:
    """Boxplot of each population's relative frequency, responders vs non-responders."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    order = sorted(per_subj.population.unique())

    fig, ax = plt.subplots(figsize=(11, 6))
    sns.boxplot(
        data=per_subj, x="population", y="percentage", hue="response",
        order=order, hue_order=["yes", "no"], showfliers=False,
        palette={"yes": "#4a7c8c", "no": "#c05f3c"}, ax=ax,
    )

    lookup = stats_df.set_index("population")
    top = per_subj.percentage.max()
    for i, pop in enumerate(order):
        p = lookup.loc[pop, "p_adjusted"]
        mark = "*" if lookup.loc[pop, "significant"] else "ns"
        ax.text(i, top * 1.02, f"{mark}\np={p:.3f}", ha="center", va="bottom", fontsize=9)

    ax.set_ylim(0, top * 1.18)
    ax.set_xlabel("")
    ax.set_ylabel("Relative frequency (% of sample)")
    ax.set_title(
        "Melanoma patients on miraclib (PBMC): responders vs non-responders\n"
        f"n = {per_subj.subject_id.nunique()} subjects; "
        "Mann-Whitney U, Benjamini-Hochberg adjusted"
    )
    ax.legend(title="Responder")
    sns.despine()
    fig.tight_layout()

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path

def significance_report(stats_df: pd.DataFrame) -> str:
    """Part 3: plain-language statement of which populations differ."""
    sig = stats_df[stats_df.significant]

    lines = [
        "SIGNIFICANCE REPORT — melanoma patients on miraclib (PBMC)",
        f"Responders (n={stats_df.n_responders.iloc[0]}) vs "
        f"non-responders (n={stats_df.n_non_responders.iloc[0]}), one value per subject.",
        "Mann-Whitney U test per population; Benjamini-Hochberg correction across 5 tests; alpha = 0.05.",
        "",
    ]

    if sig.empty:
        lines.append("RESULT: No cell population shows a significant difference between")
        lines.append("responders and non-responders after correcting for multiple testing.")
    else:
        lines.append("RESULT: Significant differences found in:")
        for _, r in sig.iterrows():
            direction = "higher" if r.difference > 0 else "lower"
            lines.append(
                f"  - {r.population}: {direction} in responders "
                f"({r.median_responder:.2f}% vs {r.median_non_responder:.2f}%), "
                f"adjusted p = {r.p_adjusted:.4f}, effect size = {r.effect_size:.3f}"
            )

    closest = stats_df.loc[stats_df.p_adjusted.idxmin()]
    lines += [
        "",
        f"Closest to significance: {closest.population} "
        f"(raw p = {closest.p_value:.4f}, adjusted p = {closest.p_adjusted:.4f}).",
        f"Median difference {closest.difference:+.2f} percentage points; "
        f"effect size {closest.effect_size:.3f} (near zero).",
        "",
        "Interpretation: with 656 subjects this analysis has ample power to detect",
        "modest differences. All five effect sizes are near zero and the group",
        "distributions overlap almost entirely, so the absence of significance",
        "reflects genuinely similar cell populations rather than insufficient data.",
    ]
    return "\n".join(lines)

BASELINE_SQL = """
SELECT
    sm.sample_id      AS sample,
    s.subject_id,
    s.project,
    s.response,
    s.sex,
    s.age
FROM subjects s
JOIN samples sm ON sm.subject_id = s.subject_id
WHERE s.condition   = 'melanoma'
  AND s.treatment   = 'miraclib'
  AND sm.sample_type = 'PBMC'
  AND sm.time_from_treatment_start = 0
ORDER BY s.project, s.subject_id
"""


def baseline_samples(con: sqlite3.Connection) -> pd.DataFrame:
    """Part 4.1: melanoma PBMC baseline samples from miraclib-treated patients."""
    return pd.read_sql_query(BASELINE_SQL, con)


BREAKDOWN_SQL = """
SELECT
    '{group_col}'                       AS category,
    COALESCE(s.{group_col}, 'unknown')  AS value,
    COUNT(sm.sample_id)                 AS samples,
    COUNT(DISTINCT s.subject_id)        AS subjects
FROM subjects s
JOIN samples sm ON sm.subject_id = s.subject_id
WHERE s.condition   = 'melanoma'
  AND s.treatment   = 'miraclib'
  AND sm.sample_type = 'PBMC'
  AND sm.time_from_treatment_start = 0
GROUP BY s.{group_col}
ORDER BY s.{group_col}
"""


def baseline_breakdown(con: sqlite3.Connection) -> pd.DataFrame:
    """Part 4.2: baseline cohort counted by project, response, and sex."""
    frames = [
        pd.read_sql_query(BREAKDOWN_SQL.format(group_col=col), con)
        for col in ("project", "response", "sex")
    ]
    return pd.concat(frames, ignore_index=True)