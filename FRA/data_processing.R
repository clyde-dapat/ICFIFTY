library(tidyverse)
library(readxl)

parse_cli <- function(args) {
  options <- list(
    `data-dir` = "data",
    `sample-sheet` = "sample_sheet.xlsx",
    `dilution-series` = "Drug_Dilution_Series.xlsx",
    output = "result/processed_data.csv"
  )

  if (length(args) %% 2 != 0) {
    stop("Options must be supplied as --name value pairs.", call. = FALSE)
  }

  if (length(args) > 0L) {
    for (i in seq(1, length(args), by = 2)) {
      key <- sub("^--", "", args[[i]])
      if (!startsWith(args[[i]], "--") || !key %in% names(options)) {
        stop(sprintf("Unknown option: %s", args[[i]]), call. = FALSE)
      }
      options[[key]] <- args[[i + 1]]
    }
  }

  options
}

read_plate <- function(path) {
  match <- regmatches(
    basename(path),
    regexec("Plate_([0-9]+)\\.xls[x]?$", basename(path), ignore.case = TRUE)
  )[[1]]
  if (length(match) != 2L) {
    stop(sprintf("Plate filename must end in Plate_<number>.xls or .xlsx: %s", path),
         call. = FALSE)
  }

  plate_number <- as.integer(match[[2]])
  plate <- read_excel(
    path,
    sheet = if ("Results" %in% excel_sheets(path)) "Results" else 1,
    range = "A55:N62",
    col_names = FALSE,
    .name_repair = "minimal"
  )
  if (nrow(plate) != 8L || ncol(plate) != 14L) {
    stop(sprintf("Unexpected plate layout in %s; expected A55:N62 with rows A-H and wells 1-12.",
                 path), call. = FALSE)
  }

  values <- plate[, 2:14]
  names(values) <- c("Row", paste0("Well_", 1:12))
  values <- values %>%
    mutate(
      Row = toupper(trimws(as.character(Row))),
      Plate = plate_number,
      across(starts_with("Well_"), ~ suppressWarnings(as.numeric(.x)))
    )

  if (!identical(values$Row, LETTERS[1:8])) {
    stop(sprintf("Invalid row labels in %s; expected A-H at Excel row 55-62.", path),
         call. = FALSE)
  }

  values
}

read_sample_map <- function(path) {
  samples <- read_excel(path, .name_repair = "minimal")
  required_columns <- c("Plate", "Row", "Sample", "Drug", "Drug_Dilution_Start_nM")
  missing_columns <- setdiff(required_columns, names(samples))
  if (length(missing_columns) > 0) {
    stop(sprintf("Sample sheet is missing columns: %s",
                 paste(missing_columns, collapse = ", ")), call. = FALSE)
  }

  samples <- samples %>%
    transmute(
      Plate = as.numeric(Plate),
      Row = toupper(trimws(as.character(Row))),
      Sample = trimws(as.character(Sample)),
      Drug = trimws(as.character(Drug)),
      Drug_Dilution_Start_nM = as.numeric(Drug_Dilution_Start_nM)
    ) %>%
    filter(!is.na(Sample), Sample != "")

  if (nrow(samples) == 0L ||
      anyNA(samples) ||
      any(!is.finite(samples$Plate)) ||
      any(samples$Plate != floor(samples$Plate)) ||
      any(!is.finite(samples$Drug_Dilution_Start_nM)) ||
      any(!samples$Row %in% LETTERS[1:8]) ||
      any(samples$Drug == "") ||
      anyDuplicated(samples[c("Plate", "Row")])) {
    stop(paste(
      "Sample-sheet assignments must have finite numeric plate/start values,",
      "A-H rows, nonempty sample/drug names, and unique Plate/Row pairs."
    ), call. = FALSE)
  }

  samples %>%
    mutate(
      Plate = as.integer(Plate),
      Sample_Group = sub("([_ -]?[Rr]ep[0-9]+)$", "", Sample),
      Replicate = if_else(
        grepl("[_ -]?[Rr]ep[0-9]+$", Sample),
        stringr::str_extract(Sample, "[0-9]+$"),
        as.character(row_number())
      )
    )
}

process_assay <- function(data_dir = "data",
                          sample_sheet = "sample_sheet.xlsx",
                          dilution_series = "Drug_Dilution_Series.xlsx") {
  plate_files <- list.files(
    data_dir,
    pattern = "Plate_[0-9]+\\.xls[x]?$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(plate_files) == 0L) {
    stop(sprintf("No *Plate_<number>.xls or .xlsx files found in %s.", data_dir),
         call. = FALSE)
  }

  samples <- read_sample_map(sample_sheet)
  plate_numbers <- vapply(plate_files, function(path) {
    as.integer(regmatches(
      basename(path),
      regexec("Plate_([0-9]+)\\.xls[x]?$", basename(path), ignore.case = TRUE)
    )[[1]][2])
  }, integer(1))
  if (anyDuplicated(plate_numbers)) {
    stop("Multiple input workbooks have the same Plate_<number> identifier.", call. = FALSE)
  }
  missing_plates <- setdiff(unique(samples$Plate), plate_numbers)
  if (length(missing_plates) > 0L) {
    stop(sprintf("Sample sheet refers to plate(s) with no matching input file: %s",
                 paste(missing_plates, collapse = ", ")), call. = FALSE)
  }

  raw_data <- purrr::map_dfr(plate_files, read_plate) %>%
    left_join(samples, by = c("Plate", "Row"))
  unmatched_samples <- anti_join(
    samples,
    raw_data %>% distinct(Plate, Row),
    by = c("Plate", "Row")
  )
  if (nrow(unmatched_samples) > 0L) {
    stop("One or more sample-sheet Plate/Row assignments are absent from the reader workbooks.",
         call. = FALSE)
  }
  raw_data <- raw_data %>% filter(!is.na(Sample))

  replicate_counts <- raw_data %>%
    distinct(Plate, Sample_Group, Sample) %>%
    count(Plate, Sample_Group, name = "Replicate_Count")
  if (any(replicate_counts$Replicate_Count != 3L)) {
    invalid <- replicate_counts %>%
      filter(Replicate_Count != 3L) %>%
      transmute(group = paste0("plate ", Plate, " / ", Sample_Group,
                               " (", Replicate_Count, " replicates)")) %>%
      pull(group)
    stop(sprintf("Each sample must have three replicate rows: %s",
                 paste(invalid, collapse = "; ")), call. = FALSE)
  }

  controls <- raw_data %>%
    group_by(Plate, Sample_Group) %>%
    summarise(
      Avg_Virus_Control = mean(Well_10, na.rm = TRUE),
      Avg_Background = mean(Well_11, na.rm = TRUE),
      .groups = "drop"
    )
  if (any(!is.finite(controls$Avg_Virus_Control)) ||
      any(controls$Avg_Virus_Control <= 0) ||
      any(!is.finite(controls$Avg_Background))) {
    stop("Every sample needs finite well 10 virus controls and well 11 background controls.",
         call. = FALSE)
  }

  dilution <- read_excel(dilution_series, .name_repair = "minimal")
  if (!"Dilution" %in% names(dilution)) {
    stop("Dilution workbook must contain a 'Dilution' column.", call. = FALSE)
  }
  dilution_long <- dilution %>%
    pivot_longer(
      cols = -Dilution,
      names_to = "Drug_Dilution_Start_nM",
      values_to = "Concentration_nM"
    ) %>%
    mutate(
      Dilution = as.character(Dilution),
      Drug_Dilution_Start_nM = as.numeric(Drug_Dilution_Start_nM),
      Concentration_nM = as.numeric(Concentration_nM)
    ) %>%
    filter(!is.na(Dilution))

  processed <- raw_data %>%
    left_join(controls, by = c("Plate", "Sample_Group")) %>%
    pivot_longer(
      cols = all_of(paste0("Well_", 2:9)),
      names_to = "Plate_Well",
      names_pattern = "Well_([0-9]+)",
      values_to = "raw_RFU"
    ) %>%
    mutate(
      Plate_Well = as.integer(Plate_Well),
      Well = Plate_Well,
      Dilution = paste("Dil", Plate_Well - 1L)
    ) %>%
    left_join(
      dilution_long,
      by = c("Dilution", "Drug_Dilution_Start_nM")
    ) %>%
    mutate(
      RFU = raw_RFU - Avg_Background,
      Percent_Virus_Control = raw_RFU / Avg_Virus_Control * 100,
      Dataset = Sample
    )

  if (anyNA(processed$Concentration_nM) ||
      any(!is.finite(processed$Concentration_nM)) ||
      any(!is.finite(processed$raw_RFU))) {
    stop("Some drug wells have no finite count or matching concentration in the dilution workbook.",
         call. = FALSE)
  }

  processed
}

if (sys.nframe() == 0L) {
  options <- parse_cli(commandArgs(trailingOnly = TRUE))
  processed <- process_assay(
    data_dir = options[["data-dir"]],
    sample_sheet = options[["sample-sheet"]],
    dilution_series = options[["dilution-series"]]
  )
  dir.create(dirname(options$output), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(processed, options$output)
  message(sprintf(
    "Processed %d plate(s), %d sample group(s), and %d measurements into %s.",
    dplyr::n_distinct(processed$Plate),
    dplyr::n_distinct(processed$Sample_Group),
    nrow(processed),
    options$output
  ))
}
