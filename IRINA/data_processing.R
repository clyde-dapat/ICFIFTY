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

  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[[i]])
    if (!startsWith(args[[i]], "--") || !key %in% names(options)) {
      stop(sprintf("Unknown option: %s", args[[i]]), call. = FALSE)
    }
    options[[key]] <- args[[i + 1]]
  }

  options
}

read_plate <- function(path) {
  plate_number <- as.integer(sub("^.*Plate_([0-9]+)\\.xlsx$", "\\1", basename(path)))
  if (is.na(plate_number)) {
    stop(sprintf("Plate filename must end in Plate_<number>.xlsx: %s", path), call. = FALSE)
  }

  plate <- read_excel(
    path,
    sheet = "Microplate End point",
    range = "A15:M23",
    col_names = FALSE,
    .name_repair = "minimal"
  )
  if (nrow(plate) != 9 || ncol(plate) != 13 ||
      !identical(as.integer(unlist(plate[1, -1], use.names = FALSE)), 1:12)) {
    stop(sprintf("Unexpected plate layout in %s; expected A-H rows and 12 numbered wells.", path),
         call. = FALSE)
  }

  values <- plate[-1, ]
  names(values) <- c("Row", paste0("Well_", 1:12))
  values <- values %>%
    mutate(
      Row = as.character(Row),
      Plate = plate_number,
      across(starts_with("Well_"), as.numeric)
    )

  if (!all(values$Row %in% LETTERS[1:8]) || anyDuplicated(values$Row)) {
    stop(sprintf("Invalid or duplicate row labels in %s.", path), call. = FALSE)
  }

  values
}

options <- parse_cli(commandArgs(trailingOnly = TRUE))
plate_files <- list.files(
  options[["data-dir"]],
  pattern = "Plate_[0-9]+\\.xlsx$",
  full.names = TRUE
)
if (length(plate_files) == 0) {
  stop(sprintf("No files matching *Plate_<number>.xlsx found in %s.", options[["data-dir"]]),
       call. = FALSE)
}

samples <- read_excel(options[["sample-sheet"]], .name_repair = "minimal")
required_sample_columns <- c("Plate", "Row", "Sample", "Drug", "Drug_Dilution_Start_nM")
missing_sample_columns <- setdiff(required_sample_columns, names(samples))
if (length(missing_sample_columns) > 0) {
  stop(sprintf("Sample sheet is missing columns: %s",
               paste(missing_sample_columns, collapse = ", ")), call. = FALSE)
}

samples <- samples %>%
  transmute(
    Plate = as.numeric(Plate),
    Row = toupper(trimws(as.character(Row))),
    Sample = trimws(as.character(Sample)),
    Drug = trimws(as.character(Drug)),
    Drug_Dilution_Start_nM = as.numeric(Drug_Dilution_Start_nM)
  )
if (anyNA(samples) ||
    any(!is.finite(samples$Plate)) ||
    any(samples$Plate != floor(samples$Plate)) ||
    any(!is.finite(samples$Drug_Dilution_Start_nM)) ||
    any(!samples$Row %in% LETTERS[1:8]) ||
    any(samples$Sample == "") ||
    anyDuplicated(samples[c("Plate", "Row")])) {
  stop(paste(
    "Sample sheet values must be finite, Plate values must be whole numbers,",
    "rows must be A-H, sample names must be nonempty, and Plate/Row assignments unique."
  ), call. = FALSE)
}
samples$Plate <- as.integer(samples$Plate)

plate_numbers <- as.integer(sub("^.*Plate_([0-9]+)\\.xlsx$", "\\1", basename(plate_files)))
if (anyDuplicated(plate_numbers)) {
  stop("Multiple input workbooks have the same Plate_<number> identifier.", call. = FALSE)
}
missing_plates <- setdiff(unique(samples$Plate), plate_numbers)
if (length(missing_plates) > 0) {
  stop(sprintf("Sample sheet refers to plate(s) with no matching input file: %s",
               paste(missing_plates, collapse = ", ")), call. = FALSE)
}

raw_data <- purrr::map_dfr(plate_files, read_plate) %>%
  left_join(samples, by = c("Plate", "Row"))
unmatched_samples <- anti_join(samples, raw_data %>% distinct(Plate, Row),
                               by = c("Plate", "Row"))
if (nrow(unmatched_samples) > 0) {
  stop("One or more sample-sheet Plate/Row assignments are absent from the reader workbooks.",
       call. = FALSE)
}
if (anyNA(raw_data$Sample)) {
  unmatched <- raw_data %>%
    filter(is.na(Sample)) %>%
    distinct(Plate, Row) %>%
    transmute(location = paste0("plate ", Plate, " row ", Row)) %>%
    pull(location)
  stop(sprintf("Plate rows have no sample-sheet assignment: %s",
               paste(unmatched, collapse = ", ")), call. = FALSE)
}

dilution <- read_excel(options[["dilution-series"]], .name_repair = "minimal")
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
  )

raw_data <- raw_data %>%
  group_by(Plate) %>%
  mutate(Avg_CC = mean(Well_12[Row %in% LETTERS[1:6]], na.rm = TRUE)) %>%
  ungroup()

if (any(!is.finite(raw_data$Avg_CC))) {
  stop("Could not calculate a finite no-virus control mean from column 12, rows A-F, on every plate.",
       call. = FALSE)
}

processed <- raw_data %>%
  pivot_longer(
    cols = all_of(paste0("Well_", 1:11)),
    names_to = "Well",
    values_to = "raw_RFU"
  ) %>%
  mutate(
    Well = as.integer(sub("^Well_", "", Well)),
    Dilution = paste("Dil", Well)
  ) %>%
  left_join(dilution_long, by = c("Dilution", "Drug_Dilution_Start_nM")) %>%
  mutate(
    RFU = raw_RFU - Avg_CC,
    Dataset = Sample
  )

if (anyNA(processed$Concentration_nM) || any(!is.finite(processed$Concentration_nM))) {
  stop("Some sample dilution steps have no finite matching concentration in the dilution workbook.",
       call. = FALSE)
}
if (any(!is.finite(processed$Avg_CC))) {
  stop("Could not determine a finite cell-control value for every plate.", call. = FALSE)
}

output_dir <- dirname(options$output)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(processed, options$output)
message(sprintf(
  "Processed %d plate(s), %d sample(s), and %d measurements into %s.",
  length(plate_files), n_distinct(processed$Dataset), nrow(processed), options$output
))
