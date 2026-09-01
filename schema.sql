-- Schema for the Loblaw Bio cell-count database.
-- Three tables: one row per subject, one per sample, one per measurement.

DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;

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

CREATE INDEX idx_samples_subject        ON samples(subject_id);
CREATE INDEX idx_cell_counts_sample     ON cell_counts(sample_id);
CREATE INDEX idx_cell_counts_population ON cell_counts(population);