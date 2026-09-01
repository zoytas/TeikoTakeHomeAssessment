# Loblaw Bio — immune cell population analysis

Analysis of immune cell counts from a clinical trial of miraclib, covering database
design, relative frequency calculation, response prediction statistics, and baseline
cohort queries.

**Live dashboard:** [loblaw-bio-dashboard.onrender.com](https://loblaw-bio-dashboard.onrender.com/)

> Note: the dashboard is hosted on Render's free tier, which sleeps after 15 minutes
> of inactivity. The first load may take 30–60 seconds to wake.

---

## Running the code

Requires Python 3.10+. In GitHub Codespaces or any terminal:

```bash
make setup      # install dependencies from requirements.txt
make pipeline   # build the database and generate all outputs
make dashboard  # start the dashboard at localhost:8501
```

`make pipeline` runs `load_data.py` then `run_analysis.py`. It is idempotent — run it
as many times as you like; each run rebuilds `cell_counts.db` and `outputs/` from
`cell-count.csv`.

`make clean` removes the database and outputs if you want to verify a build from scratch.

### Generated outputs

| File | Contents |
|---|---|
| `cell_counts.db` | SQLite database (Part 1) |
| `outputs/summary_table.csv` | Relative frequencies, 52,500 rows (Part 2) |
| `outputs/response_cohort.csv` | Per-subject percentages for the response cohort (Part 3) |
| `outputs/response_statistics.csv` | Test results per population (Part 3) |
| `outputs/response_boxplot.png` | Responders vs non-responders boxplot (Part 3) |
| `outputs/significance_report.txt` | Plain-language conclusion (Part 3) |
| `outputs/baseline_samples.csv` | 656 baseline samples (Part 4) |
| `outputs/baseline_breakdown.csv` | Counts by project, response, sex (Part 4) |

---

## Part 1: Database schema

Three tables in third normal form:

```
subjects  (3,500 rows)      one row per person
    |  subject_id
samples   (10,500 rows)     one row per blood draw
    |  sample_id
cell_counts (52,500 rows)   one row per measurement
```

```sql
CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY NOT NULL,
    project    TEXT NOT NULL,
    condition  TEXT NOT NULL,
    age        INTEGER NOT NULL,
    sex        TEXT NOT NULL CHECK (sex IN ('M', 'F')),
    treatment  TEXT NOT NULL,
    response   TEXT CHECK (response IN ('yes', 'no'))
);

CREATE TABLE samples (
    sample_id                  TEXT PRIMARY KEY NOT NULL,
    subject_id                 TEXT NOT NULL,
    sample_type                TEXT NOT NULL,
    time_from_treatment_start  INTEGER,
    FOREIGN KEY (subject_id) REFERENCES subjects(subject_id)
);

CREATE TABLE cell_counts (
    id         INTEGER PRIMARY KEY,
    sample_id  TEXT NOT NULL,
    population TEXT NOT NULL,
    count      INTEGER NOT NULL CHECK (count >= 0),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    UNIQUE (sample_id, population)
);
```

### Why three tables

The CSV has 10,500 rows but only 3,500 subjects — each person contributes three
samples (days 0, 7, 14). Seven columns (`project`, `condition`, `age`, `sex`,
`treatment`, `response`, plus `subject`) are identical across a subject's three rows.
That is 49,000 duplicate values out of 73,500 stored.

The problem with duplication is not disk space, it is that copies can disagree.
Correcting a subject's age in one row and not the others leaves the file asserting
two different ages with nothing to flag the contradiction. Splitting subject-level
facts into their own table means each fact is stored once, so it cannot be
inconsistent with itself. The same split removes the insert anomaly (a subject cannot
be recorded before their first sample exists) and the delete anomaly (removing bad
samples would otherwise erase the subject entirely).

`sample_type` is a deliberate exception. It happens to be constant per subject in
this dataset, but PBMC versus whole blood describes how a specimen was processed,
not a property of the person — a future trial could vary it within a subject. It
therefore lives on `samples`.

### Why cell counts are stored long, not wide

`cell_counts` has one row per (sample, population) pair rather than five columns.
This means:

- Population names are data, not schema. `SUM(count)` totals a sample regardless of
  how many populations exist. With five columns, every query hard-codes the list
  `b_cell + cd8_t_cell + ...`, and adding a sixth population silently produces wrong
  totals in every query that was not updated — no error, just quietly incomplete sums.
- `UNIQUE (sample_id, population)` makes duplicate measurements structurally impossible.
- The shape matches what Part 2 asks for and what plotting and statistics libraries expect.

The cost is 5× the rows and a pivot to recover the CSV's original layout. Measured on
this dataset the query-time difference is 79 ms vs 82 ms, so it did not influence the
decision.

### Constraints

`CHECK` constraints on `sex`, `response`, and `count` reject invalid values at insert
time. `response` is deliberately nullable — the 474 healthy control subjects received
no treatment, so `NULL` correctly means "not applicable" rather than "missing". A
`CHECK` passes when the value is NULL, so the constraint permits blanks while still
rejecting values outside `yes`/`no`.

Note that SQLite ignores foreign keys unless `PRAGMA foreign_keys = ON` is set, and
the setting is per-connection rather than stored in the file. Every module that opens
the database goes through `analysis.connect()`, which sets it.

### Scaling to hundreds of projects and thousands of samples

At the target scale described — hundreds of projects, thousands of samples — this
schema is unchanged and SQLite remains adequate. The current database is 6 MB; a
100× larger trial would be roughly 600 MB, well within SQLite's practical range, and
the existing indexes on `samples(subject_id)`, `cell_counts(sample_id)` and
`cell_counts(population)` keep lookups at index-search rather than full-scan
(verified with `EXPLAIN QUERY PLAN`; indexed lookups measured ~178× faster than
unindexed on the current data).

Three things would change further out:

**Concurrency before size.** SQLite allows one writer at a time. The first real
pressure is not row count but multiple analysts or an ingestion process writing
concurrently — that is the point to move to PostgreSQL, and the schema ports directly.

**Columnar storage for analytics.** At hundreds of millions of rows, the row-oriented
layout becomes the bottleneck for aggregate queries. Parquet or DuckDB stores each
column contiguously, and the `population` column — millions of values drawn from five
strings — dictionary-compresses to almost nothing. Measured on this data, the long
format is 5× the rows of wide but only 11% larger as Parquet, so the storage argument
against long format disappears entirely at scale.

**Selective denormalization.** Normalization protects against update anomalies, which
arise from editing. A trial data warehouse is largely append-only, so at very large
scale it is common to copy `condition` and `treatment` onto the measurement rows to
avoid join cost — accepting redundancy precisely because nothing is edited in place.
The general rule: normalize where data is written and edited, denormalize where it is
read and aggregated. This dataset is firmly in the first category, since records are
corrected as the trial runs.

For "various types of analytics", the long format is what makes new questions cheap:
adding a population, a new timepoint, or a new grouping variable requires inserting
rows, not altering tables or rewriting queries.

---

## Part 2: Relative frequencies

For each sample, the five populations are summed and each expressed as a percentage
of that total. Output columns: `sample`, `total_count`, `population`, `count`,
`percentage` — 52,500 rows, one per population per sample.

Raw counts are not comparable across samples because total cells recovered varies by
draw. Percentages remove that.

The per-sample total is computed with a window function:

```sql
SUM(count) OVER (PARTITION BY sample_id) AS total_count
```

`GROUP BY` would collapse the five rows of a sample into one, losing the individual
counts. `OVER (PARTITION BY ...)` computes the same sum while keeping every row, so
each population row carries its sample's total and the percentage is a plain division.

A join against a subquery of totals produces identical results (verified) and runs
marginally faster; the window function was chosen for readability.

---

## Part 3: Response prediction

**Cohort:** melanoma patients treated with miraclib, PBMC samples only —
1,968 samples from 656 subjects (331 responders, 325 non-responders).

### Analytical decisions

**One value per subject.** Each subject contributes three samples. Treating them as
independent observations would overstate the evidence (pseudoreplication), so each
subject's three timepoints are averaged into one value per population, giving 656
independent observations. Mean percentages barely move across days 0/7/14
(e.g. b_cell 9.93 / 9.93 / 9.83), so this loses little signal.

**Mann-Whitney U test.** Non-parametric, so it does not assume normally distributed
percentages. Two-sided, since there is no prior expectation of direction.

**Benjamini-Hochberg correction.** Five populations means five tests, and five tests
at α=0.05 carry roughly a 23% chance of at least one false positive. BH controls the
false discovery rate across the family. Applied via `statsmodels.stats.multitest`.

**Effect sizes reported alongside p-values.** With 656 subjects the analysis is
powered to detect differences too small to be clinically meaningful, so rank-biserial
correlation is reported to separate "detectable" from "important".

### Result

No cell population differs significantly between responders and non-responders after
correction.

| Population | Median (resp) | Median (non-resp) | Effect size | Raw p | Adjusted p | Significant |
|---|---|---|---|---|---|---|
| b_cell | 9.67 | 9.84 | −0.042 | 0.346 | 0.432 | No |
| cd4_t_cell | 30.21 | 29.82 | 0.113 | 0.012 | 0.062 | No |
| cd8_t_cell | 24.90 | 25.01 | −0.022 | 0.622 | 0.622 | No |
| monocyte | 19.79 | 20.28 | −0.050 | 0.265 | 0.432 | No |
| nk_cell | 14.74 | 14.96 | −0.069 | 0.127 | 0.317 | No |

cd4_t_cell is closest, with a raw p of 0.012 that would appear significant if tested
in isolation. Adjusted for the five comparisons it reaches 0.062. Its median
difference is 0.39 percentage points with an effect size of 0.11 — near zero.

This is a negative result rather than an underpowered one. With 656 subjects the study
has ample power to detect modest differences; all five effect sizes are near zero and
the distributions overlap almost entirely (see `outputs/response_boxplot.png`).
On this data, none of the five populations predicts miraclib response.

---

## Part 4: Baseline subset

Melanoma PBMC samples at day 0 from miraclib-treated patients: **656 samples from
656 subjects** (exactly one baseline sample each).

| Category | Value | Samples | Subjects |
|---|---|---|---|
| project | prj1 | 384 | 384 |
| project | prj3 | 272 | 272 |
| response | no | 325 | 325 |
| response | yes | 331 | 331 |
| sex | F | 312 | 312 |
| sex | M | 344 | 344 |

`prj2` does not appear because it contributes no melanoma/miraclib/PBMC baseline
samples — `GROUP BY` reports only groups that exist, so an absent project yields no
row rather than a zero.

Samples and subjects are equal throughout because each subject has exactly one
baseline sample; both are reported because the requirement asks for samples per
project but subjects per response and sex, and the distinction matters if the trial
design changes.

---

## Code structure

```
├── Makefile              setup / pipeline / dashboard / clean
├── schema.sql            table definitions
├── load_data.py          Part 1: CSV -> SQLite
├── analysis.py           Parts 2-4: queries and statistics
├── run_analysis.py       batch runner, writes outputs/
├── dashboard.py          Streamlit interface
├── requirements.txt
├── cell-count.csv        input
├── cell_counts.db        generated
└── outputs/              generated
```

**`schema.sql` is separate from `load_data.py`** so the SQL is syntax-highlighted and
reviewable as SQL rather than buried in a Python string.

**`analysis.py` holds queries but writes nothing.** Every function takes a connection
and returns a DataFrame. This is what lets `run_analysis.py` and `dashboard.py` share
identical logic — the dashboard is not a reimplementation, so the numbers on screen
cannot drift from the numbers in `outputs/`.

**`run_analysis.py` handles all file output**, keeping the side effects in one place.

**Paths are anchored to the script location** (`Path(__file__).resolve().parent`)
rather than the working directory, so `make pipeline` behaves the same regardless of
where it is invoked from.

**Both scripts fail loudly and early.** `load_data.py` exits with a clear message if
the CSV is missing or has unexpected columns, and — importantly — if a subject's
attributes disagree across their rows, which would otherwise surface much later as an
opaque constraint violation.

## Verification

Results were checked against an independent path that reads `cell-count.csv` directly
into pandas and computes the statistics without touching SQLite. Both routes produce
identical p-values (cd4_t_cell: 0.0124), confirming the schema, loader, joins, and
window function are correct. Percentages sum to 100 per sample (±0.02 from rounding),
and cohort sizes match counts taken directly from the CSV.