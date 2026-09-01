"""Loblaw Bio cell-count dashboard.

Run locally:  streamlit run dashboard.py
"""

from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

import analysis

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"

st.set_page_config(page_title="Loblaw Bio — Cell Counts", layout="wide")
st.title("Loblaw Bio — immune cell population dashboard")

if not DB_PATH.exists():
    st.warning("No database found. Run `make pipeline` to build it.")
    st.stop()


@st.cache_data
def load_all():
    con = analysis.connect()
    try:
        summary = analysis.summary_table(con)
        cohort = analysis.response_cohort(con)
        per_subj = analysis.per_subject(cohort)
        stats_df = analysis.compare_response(per_subj)
        report = analysis.significance_report(stats_df)
        baseline = analysis.baseline_samples(con)
        breakdown = analysis.baseline_breakdown(con)
    finally:
        con.close()
    return summary, per_subj, stats_df, report, baseline, breakdown


summary, per_subj, stats_df, report, baseline, breakdown = load_all()

tab2, tab3, tab4 = st.tabs(
    ["Part 2 — Frequencies", "Part 3 — Response", "Part 4 — Baseline"]
)

with tab2:
    st.subheader("Relative frequency of each population, per sample")
    st.caption(f"{len(summary):,} rows — one per population per sample.")

    samples = sorted(summary["sample"].unique())
    pick = st.multiselect("Filter to specific samples (leave empty for all)", samples)
    view = summary[summary["sample"].isin(pick)] if pick else summary

    st.dataframe(view, use_container_width=True, hide_index=True, height=400)
    st.download_button(
        "Download summary_table.csv",
        summary.to_csv(index=False),
        "summary_table.csv",
        "text/csv",
    )
    
with tab3:
    st.subheader("Responders vs non-responders — melanoma, miraclib, PBMC")

    c1, c2, c3 = st.columns(3)
    c1.metric("Subjects", per_subj.subject_id.nunique())
    c2.metric("Responders", int(stats_df.n_responders.iloc[0]))
    c3.metric("Non-responders", int(stats_df.n_non_responders.iloc[0]))

    pops = sorted(per_subj.population.unique())
    chosen = st.multiselect("Populations to show", pops, default=pops, key="p3_pops")
    plot_df = per_subj[per_subj.population.isin(chosen)] if chosen else per_subj

    labels = {
        r.population: f"{r.population}<br>p={r.p_adjusted:.3f}"
        for _, r in stats_df.iterrows()
    }
    plot_df = plot_df.assign(label=plot_df.population.map(labels))

    fig = px.box(
        plot_df,
        x="label",
        y="percentage",
        color="response",
        points=False,
        category_orders={
            "label": [labels[p] for p in pops if p in chosen],
            "response": ["yes", "no"],
        },
        color_discrete_map={"yes": "#4a7c8c", "no": "#c05f3c"},
        labels={
            "percentage": "Relative frequency (%)",
            "label": "",
            "response": "Responder",
        },
    )
    fig.update_layout(
        boxmode="group",
        height=500,
        margin=dict(t=30, b=10),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
    )
    st.plotly_chart(fig, use_container_width=True)
    st.caption(
        "Hover a box for quartiles; click a legend entry to hide that group. "
        "Benjamini-Hochberg adjusted p-values shown under each population. "
        "A static version is saved to outputs/response_boxplot.png."
    )

    st.markdown("**Statistical results**")
    st.dataframe(
        stats_df.round(4), use_container_width=True, hide_index=True
    )

    st.markdown("**Conclusion**")
    st.code(report, language=None)
with tab4:
    st.subheader("Melanoma PBMC baseline samples (day 0), miraclib-treated")
    st.metric("Baseline samples", len(baseline))

    for cat, label in [
        ("project", "Samples per project"),
        ("response", "Subjects by response"),
        ("sex", "Subjects by sex"),
    ]:
        st.markdown(f"**{label}**")
        st.dataframe(
            breakdown[breakdown.category == cat][["value", "samples", "subjects"]],
            use_container_width=True, hide_index=True,
        )

    with st.expander("All 656 baseline samples"):
        st.dataframe(baseline, use_container_width=True, hide_index=True)