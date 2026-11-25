# ICFIFTY

ICFIFTY combines two IC50 determination workflows using the command-line and
Shiny application:

- **FRA** — focus-reduction assay processing and IC50/IC90 analysis
- **IRINA** — IC50/IC90 determination from reader plate data

Choose the method that matches the plate layout and its sample/dilution
workbooks. The input validator checks the plate format and required workbook
columns before analysis. If a format does not match the selected method, the
app reports that the user should check the files or select the correct method.

## Requirements and installation

Install [Miniforge](https://github.com/conda-forge/miniforge). From the repository root, create and activate
the environment:

```bash
conda env create -f environment.yml
conda activate icfifty
```

The environment includes Snakemake, R, and the R packages used by the workflows
and Shiny app.

## Run from the command line

From the repository root, select the workflow with the Snakemake `method`
configuration argument:

```bash
snakemake --cores 1 --config method=FRA
```

For IC50 determination with IRINA:

```bash
snakemake --cores 1 --config method=IRINA
```

Preview the selected workflow without running it:

```bash
snakemake --dry-run --config method=FRA
```

Snakemake reads the selected workflow's bundled inputs under `FRA/` or `IRINA/`.
To analyse a different experiment, replace the files in that method's `data/`
folder and its `sample_sheet.xlsx` and `Drug_Dilution_Series.xlsx` files. Keep
the documented workbook layouts and filenames. Results and logs are separated
by method:

- `result/FRA/` and `logs/FRA/`
- `result/IRINA/` and `logs/IRINA/`

The processed measurements, dose-response PDF, and IC50/IC90 summary CSV are
written to the selected method's result folder. On an input or analysis error,
Snakemake prints the relevant log and advises checking the selected method and
input files.

## Run the Shiny app

```bash
Rscript run_app.R
```

Open <http://127.0.0.1:3838>. Use the full-page **Input workbooks** tab to
select **FRA** or **IRINA**, upload plate-reader files, a sample sheet, and a
dilution-series workbook, then select **Run analysis**. All three input boxes
are required; the dashboard lists selected filenames and reports an error if
any input is missing. The **Summary**,
**Dose-response plots**, and **Downloads** tabs contain the results. The app
validates the selected method against the plate layout and displays guidance
when files do not match. Uploaded data and generated files are kept in a
session-specific temporary directory and are removed when the browser session
ends.

The plate files and sample sheet must cover the same plate numbers. A mismatch is reported before
processing, with the missing plate numbers and what to change.

The app binds to `127.0.0.1:3838` by default. Configure the bind address or port
with `ICFIFTY_SHINY_HOST` and `ICFIFTY_SHINY_PORT`, respectively.

## Input formats

### FRA

- Plate files named `*Plate_<number>.xls` or `.xlsx`.
- Spot counts at `Results!A55:N62` (or the same range on the first sheet),
  with row labels A-H.
- A sample sheet with `Plate`, `Row`, `Sample`, `Drug`, and
  `Drug_Dilution_Start_nM` columns.
- A dilution workbook with a `Dilution` column and concentration columns
  matching `Drug_Dilution_Start_nM`.

FRA reads eight drug dilutions from wells 2-9, virus controls from well 10, and
background controls from well 11. It reports replicate and combined-sample
estimates.

### IRINA

- Plate files named `*Plate_<number>.xlsx`.
- A `Microplate End point` sheet with the reader results in `A15:M23`, row
  labels A-H, and well numbers 1-12.
- The same sample-sheet and dilution-workbook column names listed for FRA.

IRINA uses the selected sample/dilution mappings and the plate's control
readings to produce per-sample dose-response curves and IC50/IC90 estimates.

## Repository layout

```text
.
├── app.R
├── run_app.R
├── Snakefile
├── validate_workflow_inputs.R
├── environment.yml
├── FRA/
│   ├── data/
│   ├── data_processing.R
│   └── dose_response_analysis.R
└── IRINA/
    ├── data/
    ├── data_processing.R
    └── dose_response_analysis.R
```

Analysis outputs, logs, Snakemake state, and local editor files are excluded by
[.gitignore](.gitignore).
