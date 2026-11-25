from pathlib import Path
import re


method = config.get("method")
if method not in {"FRA", "IRINA"}:
    raise ValueError(
        "Choose an analysis method with --config method=FRA or --config method=IRINA."
    )

workflow_dir = Path(method)
data_dir = workflow_dir / "data"
plate_pattern = r"Plate_[0-9]+\.xls[x]?$" if method == "FRA" else r"Plate_[0-9]+\.xlsx$"
plate_files = sorted(
    str(path)
    for path in data_dir.glob("*Plate_*")
    if re.search(plate_pattern, path.name, re.IGNORECASE)
)
if not plate_files:
    raise ValueError(
        f"No plate files for {method} were found in {data_dir}. "
        "Check the selected method and its input files."
    )

processed = f"result/{method}/processed_data.csv"
curves = f"result/{method}/dose_response_curves.pdf"
summary = f"result/{method}/ic50_ic90_summary_table.csv"
processing_script = str(workflow_dir / "data_processing.R")
analysis_script = str(workflow_dir / "dose_response_analysis.R")
sample_sheet = str(workflow_dir / "sample_sheet.xlsx")
dilution_series = str(workflow_dir / "Drug_Dilution_Series.xlsx")
validator_script = "validate_workflow_inputs.R"


rule all:
    input:
        curves,
        summary,


rule process_plate_data:
    input:
        plates=plate_files,
        sample_sheet=sample_sheet,
        dilution_series=dilution_series,
        processing_script=processing_script,
        validator_script=validator_script,
    output:
        processed,
    log:
        f"logs/{method}/data_processing.log",
    conda:
        "environment.yml",
    shell:
        """
        Rscript {input.validator_script:q} \
          --method {method:q} \
          --data-dir {data_dir:q} \
          --sample-sheet {input.sample_sheet:q} \
          --dilution-series {input.dilution_series:q} \
          > {log:q} 2>&1 || {{ cat {log:q}; exit 1; }}
        Rscript {input.processing_script:q} \
          --data-dir {data_dir:q} \
          --sample-sheet {input.sample_sheet:q} \
          --dilution-series {input.dilution_series:q} \
          --output {output:q} \
          >> {log:q} 2>&1 || {{
            cat {log:q}
            echo "Check that the selected method matches your input data and that you chose the correct files." >&2
            exit 1
          }}
        """


rule fit_dose_response:
    input:
        processed=processed,
        analysis_script=analysis_script,
    output:
        curves=curves,
        summary=summary,
    log:
        f"logs/{method}/dose_response_analysis.log",
    conda:
        "environment.yml",
    shell:
        """
        Rscript {input.analysis_script:q} \
          --input-csv {input.processed:q} \
          --output-pdf {output.curves:q} \
          --output-csv {output.summary:q} \
          > {log:q} 2>&1 || {{
            cat {log:q}
            exit 1
          }}
        """
