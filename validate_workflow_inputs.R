library(readxl)

parse_cli <- function(args) {
  options <- list(
    method = NULL,
    `data-dir` = NULL,
    `sample-sheet` = NULL,
    `dilution-series` = NULL
  )

  if (length(args) %% 2 != 0L) {
    stop("Options must be supplied as --name value pairs.", call. = FALSE)
  }
  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[[i]])
    if (!startsWith(args[[i]], "--") || !key %in% names(options)) {
      stop(sprintf("Unknown option: %s", args[[i]]), call. = FALSE)
    }
    options[[key]] <- args[[i + 1L]]
  }
  if (any(vapply(options, is.null, logical(1)))) {
    stop(
      "Supply --method, --data-dir, --sample-sheet, and --dilution-series.",
      call. = FALSE
    )
  }
  options
}

validate_workflow_inputs <- function(method, data_dir, sample_sheet, dilution_series) {
  if (!method %in% c("FRA", "IRINA")) {
    stop("Method must be FRA or IRINA.", call. = FALSE)
  }

  mismatch <- function(details) {
    stop(
      sprintf(
        paste(
          "Input data do not match the selected method '%s'.",
          "%s Check that you selected the correct files or choose the matching",
          "method (FRA or IRINA)."
        ),
        method,
        details
      ),
      call. = FALSE
    )
  }

  extension_pattern <- if (method == "FRA") {
    "Plate_[0-9]+\\.xls[x]?$"
  } else {
    "Plate_[0-9]+\\.xlsx$"
  }
  all_files <- list.files(data_dir, full.names = TRUE)
  plates <- all_files[grepl(extension_pattern, basename(all_files), ignore.case = TRUE)]
  if (length(plates) == 0L) {
    mismatch(sprintf("No expected %s plate workbooks were found in '%s'.", method, data_dir))
  }
  if (length(plates) != length(all_files)) {
    invalid_files <- setdiff(basename(all_files), basename(plates))
    mismatch(sprintf(
      "These plate input files do not use the expected %s filename/extension: %s.",
      method,
      paste(invalid_files, collapse = ", ")
    ))
  }
  plate_ids <- sub(
    ".*Plate_([0-9]+)\\.[^.]+$",
    "\\1",
    basename(plates),
    ignore.case = TRUE
  )
  if (anyDuplicated(plate_ids)) {
    mismatch("More than one plate workbook has the same Plate_<number> identifier.")
  }

  for (path in plates) {
    tryCatch({
      sheets <- excel_sheets(path)
      if (method == "FRA") {
        sheet <- if ("Results" %in% sheets) "Results" else 1
        plate <- read_excel(
          path,
          sheet = sheet,
          range = "A55:N62",
          col_names = FALSE,
          .name_repair = "minimal"
        )
        if (nrow(plate) != 8L || ncol(plate) != 14L ||
            !identical(
              toupper(trimws(as.character(unlist(plate[1:8, 2], use.names = FALSE)))),
              LETTERS[1:8]
            )) {
          mismatch(sprintf(
            "'%s' is not in the FRA plate layout (Results rows 55-62).",
            basename(path)
          ))
        }
      } else {
        if (!"Microplate End point" %in% sheets) {
          mismatch(sprintf(
            "'%s' has no 'Microplate End point' sheet required by IRINA.",
            basename(path)
          ))
        }
        plate <- read_excel(
          path,
          sheet = "Microplate End point",
          range = "A15:M23",
          col_names = FALSE,
          .name_repair = "minimal"
        )
        well_headers <- suppressWarnings(as.numeric(unlist(plate[1, 2:13], use.names = FALSE)))
        row_labels <- toupper(trimws(as.character(unlist(plate[2:9, 1], use.names = FALSE))))
        if (nrow(plate) != 9L || ncol(plate) != 13L ||
            !identical(well_headers, as.numeric(1:12)) ||
            !identical(row_labels, LETTERS[1:8])) {
          mismatch(sprintf(
            "'%s' is not in the IRINA plate layout (Microplate End point rows 15-23).",
            basename(path)
          ))
        }
      }
    }, error = function(error) {
      message_text <- conditionMessage(error)
      if (grepl("Input data do not match the selected method", message_text, fixed = TRUE)) {
        stop(message_text, call. = FALSE)
      }
      mismatch(sprintf(
        "Could not read '%s' using the expected %s plate layout (%s).",
        basename(path),
        method,
        message_text
      ))
    })
  }

  if (!file.exists(sample_sheet) || !file.exists(dilution_series)) {
    mismatch("The sample sheet or drug dilution-series workbook is missing.")
  }
  sample_columns <- tryCatch(
    names(read_excel(sample_sheet, n_max = 0, .name_repair = "minimal")),
    error = function(error) mismatch(sprintf("Could not read the sample sheet (%s).",
                                             conditionMessage(error)))
  )
  required_sample_columns <- c(
    "Plate", "Row", "Sample", "Drug", "Drug_Dilution_Start_nM"
  )
  missing_columns <- setdiff(required_sample_columns, sample_columns)
  if (length(missing_columns) > 0L) {
    mismatch(sprintf(
      "The sample sheet is missing required columns: %s.",
      paste(missing_columns, collapse = ", ")
    ))
  }
  dilution_columns <- tryCatch(
    names(read_excel(dilution_series, n_max = 0, .name_repair = "minimal")),
    error = function(error) mismatch(sprintf("Could not read the dilution series (%s).",
                                             conditionMessage(error)))
  )
  if (!"Dilution" %in% dilution_columns || length(dilution_columns) < 2L) {
    mismatch("The dilution-series workbook must contain a 'Dilution' column and concentration columns.")
  }

  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  options <- parse_cli(commandArgs(trailingOnly = TRUE))
  validate_workflow_inputs(
    method = options$method,
    data_dir = options[["data-dir"]],
    sample_sheet = options[["sample-sheet"]],
    dilution_series = options[["dilution-series"]]
  )
  message(sprintf("Input files match the selected %s method.", options$method))
}
