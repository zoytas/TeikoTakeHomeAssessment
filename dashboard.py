"""Loblaw Bio cell-count dashboard.

Run locally:  streamlit run dashboard.py
"""

import sqlite3
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"

st.set_page_config(page_title="Loblaw Bio — Cell Counts", layout="wide")

st.title("Loblaw Bio — immune cell population dashboard")

if not DB_PATH.exists():
    st.warning(
        f"No database found at `{DB_PATH.name}`. "
        "Run `python load_data.py` (or `make pipeline`) to build it."
    )
    st.stop()


@st.cache_data
def load_overview() -> pd.DataFrame:
    con = sqlite3.connect(DB_PATH)
    try:
        return pd.read_sql_query(
            """
            SELECT s.condition, s.treatment, sm.sample_type,
                   COUNT(DISTINCT s.subject_id) AS subjects,
                   COUNT(DISTINCT sm.sample_id) AS samples
            FROM subjects s
            JOIN samples sm ON sm.subject_id = s.subject_id
            GROUP BY s.condition, s.treatment, sm.sample_type
            ORDER BY s.condition, s.treatment, sm.sample_type
            """,
            con,
        )
    finally:
        con.close()


overview = load_overview()

col1, col2 = st.columns(2)
col1.metric("Subjects", f"{overview.subjects.sum():,}")
col2.metric("Samples", f"{overview.samples.sum():,}")

st.subheader("Cohort breakdown")
st.dataframe(overview, use_container_width=True, hide_index=True)

st.caption("Parts 2–4 to follow.")
