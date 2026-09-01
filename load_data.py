"""Initialize the SQLite database and load cell-count.csv into it.

Run with:  python load_data.py
"""

import sqlite3
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = (ROOT / "schema.sql").read_text()

def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        sys.exit(f"ERROR: could not find {path}")
        
    df = pd.read_csv(path)
    
    expected = {"project", "subject", "condition", "age", "sex", "treatment",
                "response", "sample", "sample_type", "time_from_treatment_start",
                *POPULATIONS,}
    
    missing = expected - set(df.columns)
    if missing:
        sys.exit(f"ERROR: CSV is missing columns: {sorted(missing)}")
    return df

def build_subjects(df: pd.DataFrame) -> pd.DataFrame:
    cols = ["subject", "project", "condition", "age", "sex", "treatment", "response"]
    
    """ collapses the identical rows into one, results in 3500 rows"""
    subjects = df[cols].drop_duplicates()
    
    clashes = subjects[subjects.duplicated("subject", keep=False)]
    if not clashes.empty:
        sys.exit(
             "ERROR: these subjects have conflicting attributes across rows:\n"
            f"{clashes.sort_values('subject')}"
        )
    return subjects.rename(columns={"subject": "subject_id"})
def build_samples(df: pd.DataFrame) -> pd.DataFrame:
    """One row per sample."""
    cols = ["sample", "subject", "sample_type", "time_from_treatment_start"]
    samples = df[cols].rename(
        columns={"sample": "sample_id", "subject": "subject_id"}
    )
    
    dupes = samples[samples.duplicated("sample_id", keep=False)]
    if not dupes.empty:
        sys.exit(f"ERROR: duplicate sample_id values founds:\n{dupes}")
    return samples

def build_cell_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Reshape the five population columns from wide to long format."""
    long = df.melt(
        id_vars=["sample"],
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    )
    return long.rename(columns={"sample": "sample_id"})

def to_rows(df: pd.DataFrame) -> list[tuple]:
    """Convert a DataFrame to tuples, turning pandas NaN into real SQL NULLs."""
    return list(
        df.astype(object).where(pd.notna(df), None).itertuples(index=False, name=None)
    )
    
def main() -> None:
    df = read_csv(CSV_PATH)
    print(f"Read {len(df):,} rows from {CSV_PATH.name}")
    subjects = build_subjects(df)
    samples = build_samples(df)
    cell_counts = build_cell_counts(df)
    
    if DB_PATH.exists():
        """makes sure we can build repeatedly"""
        DB_PATH.unlink()
    con = sqlite3.connect(DB_PATH)
    try:
        con.execute("PRAGMA foreign_keys = ON")
        con.executescript(SCHEMA)
        
        with con:
            con.executemany(
                "INSERT INTO subjects "
                "(subject_id, project, condition, age, sex, treatment, response) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                to_rows(subjects),
            )
            con.executemany(
                "INSERT INTO samples "
                "(sample_id, subject_id, sample_type, time_from_treatment_start) "
                "VALUES (?, ?, ?, ?)",
                to_rows(samples),
            )
            con.executemany(
                "INSERT INTO cell_counts (sample_id, population, count) "
                "VALUES (?, ?, ?)",
                to_rows(cell_counts),
            )
        for table in ("subjects", "samples", "cell_counts"):
            n = con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            print(f"  {table:<12} {n:>7,} rows")

        violations = con.execute("PRAGMA foreign_key_check").fetchall()
        if violations:
            sys.exit(f"ERROR: foreign key violations after load: {violations}")

    finally:
        con.close()

    print(f"\nDatabase written to {DB_PATH}")

if __name__ == "__main__":
    main()
