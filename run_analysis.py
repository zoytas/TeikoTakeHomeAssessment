"""Run the full analysis and write outputs to outputs/.

Run with:  python run_analysis.py
"""

from pathlib import Path

import analysis

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def main() -> None:
    OUT.mkdir(exist_ok=True)
    con = analysis.connect()

    try:
        # Part 2
        summary = analysis.summary_table(con)
        summary.to_csv(OUT / "summary_table.csv", index=False)
        print(f"Part 2  summary_table.csv        {len(summary):,} rows")

        # Part 3
        cohort = analysis.response_cohort(con)
        per_subj = analysis.per_subject(cohort)
        stats_df = analysis.compare_response(per_subj)

        per_subj.to_csv(OUT / "response_cohort.csv", index=False)
        stats_df.to_csv(OUT / "response_statistics.csv", index=False)
        print(f"Part 3  response_cohort.csv      {len(per_subj):,} rows")
        print(f"Part 3  response_statistics.csv  {len(stats_df):,} rows")

        analysis.boxplot_response(per_subj, stats_df, OUT / "response_boxplot.png")
        print("Part 3  response_boxplot.png")

        report = analysis.significance_report(stats_df)
        (OUT / "significance_report.txt").write_text(report + "\n")
        print("Part 3  significance_report.txt")

        # Part 4
        baseline = analysis.baseline_samples(con)
        breakdown = analysis.baseline_breakdown(con)
        baseline.to_csv(OUT / "baseline_samples.csv", index=False)
        breakdown.to_csv(OUT / "baseline_breakdown.csv", index=False)
        print(f"Part 4  baseline_samples.csv     {len(baseline):,} rows")
        print(f"Part 4  baseline_breakdown.csv   {len(breakdown):,} rows")

    finally:
        con.close()

    print(f"\nOutputs written to {OUT}")
    print()
    print(report)


if __name__ == "__main__":
    main()