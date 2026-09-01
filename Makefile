.PHONY: setup pipeline dashboard clean

PYTHON := $(shell command -v python3 || command -v python)

setup:
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) run_analysis.py

dashboard:
	$(PYTHON) -m streamlit run dashboard.py

clean:
	rm -f cell_counts.db
	rm -rf outputs