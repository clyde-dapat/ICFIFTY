library(drc)
library(tidyverse)
library(gridExtra)

parse_cli <- function(args) {
  options <- list(
    `input-csv` = "result/processed_data.csv",
    `output-pdf` = "result/dose_response_curves.pdf",
    `output-csv` = "result/ic50_ic90_summary_table.csv"
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

cli_options <- parse_cli(commandArgs(trailingOnly = TRUE))
input_csv <- cli_options[["input-csv"]]
group_col   <- "Dataset"
sample_col  <- "Sample"
conc_col   <- "Concentration_nM"
resp_col   <- "RFU"
resp_floor  <- 1                   
output_pdf  <- cli_options[["output-pdf"]]
output_csv  <- cli_options[["output-csv"]]


# 1. Fit the 4-parameter log-logistic model for one dataset

fit_dose_response_model <- function(df, conc_col = "concentration", resp_col = "viability") {
  df <- df %>%
    filter(
      is.finite(.data[[conc_col]]),
      .data[[conc_col]] > 0,
      is.finite(.data[[resp_col]])
    )
  if (nrow(df) < 4) {
    stop("At least four finite, positive-concentration observations are required.")
  }
  formula_str <- as.formula(paste(resp_col, "~", conc_col))
  model <- drm(formula_str, data = df, fct = LL.4())
  model
}


clean_response_values <- function(df, resp_col = "viability", floor_value = 1) {
  df[[resp_col]] <- ifelse(df[[resp_col]] < 0, floor_value, df[[resp_col]])
  df
}


# 2. Calculate IC50 / IC90 point estimates + 95% delta-method CIs

calculate_ic_values <- function(model) {
  coefs <- coefficients(model)
  c_val <- coefs["c:(Intercept)"]   # lower asymptote
  d_val <- coefs["d:(Intercept)"]   # upper asymptote
  
  estimates <- ED(model, c(50, 90), interval = "delta", display = FALSE)
  
  tibble(
    Label = c("IC50", "IC90"),
    X = estimates[c(1, 2), "Estimate"],
    X_min = estimates[c(1, 2), "Lower"],
    X_max = estimates[c(1, 2), "Upper"],
    Y = c(d_val - 0.50 * (d_val - c_val),
          d_val - 0.90 * (d_val - c_val))
  )
}


# 3. Generate a prediction curve with 95% confidence band

generate_curve_with_ci <- function(model, conc_range, conc_col = "concentration", n = 500) {
  # Log-spacing requires strictly positive values
  positive_conc <- conc_range[is.finite(conc_range) & conc_range > 0]
  if (length(positive_conc) == 0) {
    stop("No positive, finite concentration values available to build curve range.")
  }
  curve_conc <- exp(seq(log(min(positive_conc)), log(max(positive_conc)), length.out = n))

  prediction_data <- data.frame(curve_conc)
  names(prediction_data) <- conc_col
  pred <- predict(model, newdata = prediction_data, interval = "confidence", level = 0.95)
  curve_grid <- data.frame(concentration = curve_conc)
  
  curve_grid$viability <- pred[, "Prediction"]
  curve_grid$Lower     <- pred[, "Lower"]
  curve_grid$Upper     <- pred[, "Upper"]
  curve_grid
}


# 4. Build the ggplot for a single dataset

plot_dose_response <- function(data, curve_grid, ic_points, strain_name,
                               conc_col = "concentration", resp_col = "viability",
                               y_label = "Response") {
  
  
  exp_range <- -3:3
  x_axis_breaks <- 10^exp_range
  x_axis_labels <- parse(text = paste0("10^", exp_range))
 
  y_values <- c(data[[resp_col]], curve_grid$viability, curve_grid$Lower, curve_grid$Upper, ic_points$Y)
  y_range <- range(y_values, na.rm = TRUE)
  y_pad <- diff(y_range) * 0.05
  y_limits <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  
  # IC50/IC90 summary text block in the upper-right corner
  ic_label_text <- ic_points %>%
    mutate(line = sprintf("%s: %.3g nM (%.3g\u2013%.3g)", Label, X, X_min, X_max)) %>%
    pull(line) %>%
    paste(collapse = "\n")
  
  p <- ggplot() +
    # 95% confidence band of the fitted sigmoid curve
    geom_ribbon(data = curve_grid, aes(x = concentration, ymin = Lower, ymax = Upper),
                fill = "darkgray", alpha = 0.25) +
    
    # Raw experimental data points
    geom_point(data = data, aes(x = .data[[conc_col]], y = .data[[resp_col]]),
               size = 3, alpha = 0.8, color = "black") +
    
    # Fitted model curve
    geom_line(data = curve_grid, aes(x = concentration, y = viability),
              color = "darkgray", linewidth = 1) +
    
    # Guide lines to IC50 / IC90
    geom_segment(data = ic_points, aes(x = 0, xend = X, y = Y, yend = Y),
                 lty = "dashed", color = "red", alpha = 0.4) +
    geom_segment(data = ic_points, aes(x = X, xend = X, y = y_limits[1], yend = Y),
                 lty = "dashed", color = "red", alpha = 0.4) +
    
    # Horizontal 95% CI error bars for IC50 / IC90
    geom_errorbarh(data = ic_points, aes(y = Y, xmin = X_min, xmax = X_max),
                   color = "blue", height = diff(y_limits) * 0.03, linewidth = 0.8) +
    
    # IC50 / IC90 point markers
    geom_point(data = ic_points, aes(x = X, y = Y), color = "blue", size = 4) +
    
    # IC50 / IC90 labels on the curve itself
    geom_text(data = ic_points, aes(x = X, y = Y, label = Label),
              vjust = -1.2, color = "blue", fontface = "bold") +
    
    # IC50 / IC90 numeric summary
    annotate("text", x = Inf, y = Inf, label = ic_label_text,
             hjust = 1.05, vjust = 1.3, size = 3, color = "blue4",
             fontface = "bold", lineheight = 1.1) +
    
    scale_x_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 0.001),
      breaks = x_axis_breaks,
      labels = x_axis_labels,
      limits = c(0, max(x_axis_breaks))
    ) +
    scale_y_continuous(limits = y_limits, labels = scales::label_number(big.mark = ",")) +
    
    labs(
#      title = paste("Dose-Response Curve:", strain_name),
      title = strain_name,
#      subtitle = "LL.4 model (drc) \u2014 IC50/IC90 with delta-method 95% CI; gray band = 95% CI of fitted curve",
      x = "Drug concentration [nM]",
      y = y_label
    ) +
    
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0.5),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
  
  p
}


# 5. Main driver: read data, split by Strain_Name, fit/plot/summarize each

if (!file.exists(input_csv)) {
  stop(sprintf("Processed input CSV does not exist: %s", input_csv), call. = FALSE)
}
raw_data <- read_csv(input_csv, show_col_types = FALSE)
required_columns <- c(group_col, sample_col, "Plate", "Row", "Well", conc_col, resp_col)
missing_columns <- setdiff(required_columns, names(raw_data))
if (length(missing_columns) > 0) {
  stop(sprintf("Input CSV is missing required columns: %s",
               paste(missing_columns, collapse = ", ")), call. = FALSE)
}
if (nrow(raw_data) == 0 || anyNA(raw_data[[group_col]])) {
  stop("Input CSV must contain at least one row with a dataset name.", call. = FALSE)
}

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

raw_data <- raw_data %>%
  mutate(Curve_ID = paste(.data[["Plate"]], .data[["Row"]], sep = "_"))
grouped_data <- raw_data %>% group_by(.data[["Curve_ID"]])
datasets <- group_split(grouped_data)
names(datasets) <- group_keys(grouped_data) %>% pull(1) %>% as.character()

all_ic_tables <- list()
all_plots <- list()

for (strain_name in names(datasets)) {
  
  df <- datasets[[strain_name]]
  sample_label <- unique(df[[sample_col]])
  if (length(sample_label) != 1L || is.na(sample_label) || sample_label == "") {
    stop(sprintf("Plate/row group '%s' must map to exactly one nonempty sample-sheet label.",
                 strain_name), call. = FALSE)
  }
  df <- df %>%
    filter(Well != 11) %>%
    clean_response_values(resp_col, floor_value = resp_floor)
  stage <- "model fitting"
  
  result <- tryCatch({
    model <- fit_dose_response_model(df, conc_col, resp_col)
    stage <- "IC50/IC90 estimation"
    ic_points <- calculate_ic_values(model)
    stage <- "curve prediction"
    curve_grid <- generate_curve_with_ci(model, df[[conc_col]], conc_col)
    stage <- "plot construction"
    p <- plot_dose_response(df, curve_grid, ic_points, sample_label,
                            conc_col, resp_col, y_label = resp_col)
    
    ic_points <- ic_points %>%
      mutate(!!group_col := sample_label) %>%
      relocate(!!group_col)
    
    list(plot = p, ic_table = ic_points)
  }, error = function(e) {
    message(sprintf("Skipping '%s' during %s (%s)",
                    strain_name, stage, conditionMessage(e)))
    NULL
  })
  
  if (!is.null(result)) {
    all_plots[[strain_name]] <- result$plot
    all_ic_tables[[strain_name]] <- result$ic_table
  } else {
    all_ic_tables[[strain_name]] <- tibble(
      Label = c("IC50", "IC90"),
      X = NA_real_,
      X_min = NA_real_,
      X_max = NA_real_,
      !!group_col := sample_label
    ) %>%
      relocate(all_of(group_col))
  }
}

# Fit each sample again using all of its plate/row replicates together.
combined_data <- raw_data %>%
  mutate(Combined_Sample = sub("_Rep[0-9]+$", "", as.character(.data[[group_col]])))
grouped_combined_data <- combined_data %>% group_by(.data[["Combined_Sample"]])
combined_datasets <- group_split(grouped_combined_data)
names(combined_datasets) <- group_keys(grouped_combined_data) %>%
  pull(1) %>%
  as.character()

for (sample_label in names(combined_datasets)) {
  df <- combined_datasets[[sample_label]] %>%
    filter(Well != 11) %>%
    clean_response_values(resp_col, floor_value = resp_floor)
  combined_label <- paste0(sample_label, " (combined replicates)")
  combined_plot_label <- paste0(sample_label, "\n(combined replicates)")
  stage <- "model fitting"

  result <- tryCatch({
    model <- fit_dose_response_model(df, conc_col, resp_col)
    stage <- "IC50/IC90 estimation"
    ic_points <- calculate_ic_values(model)
    stage <- "curve prediction"
    curve_grid <- generate_curve_with_ci(model, df[[conc_col]], conc_col)
    stage <- "plot construction"
    p <- plot_dose_response(df, curve_grid, ic_points, combined_plot_label,
                            conc_col, resp_col, y_label = resp_col)

    ic_points <- ic_points %>%
      mutate(!!group_col := combined_label) %>%
      relocate(!!group_col)

    list(plot = p, ic_table = ic_points)
  }, error = function(e) {
    message(sprintf("Skipping combined replicates for '%s' during %s (%s)",
                    sample_label, stage, conditionMessage(e)))
    NULL
  })

  if (!is.null(result)) {
    all_plots[[paste0("combined_", sample_label)]] <- result$plot
    all_ic_tables[[paste0("combined_", sample_label)]] <- result$ic_table
  } else {
    all_ic_tables[[paste0("combined_", sample_label)]] <- tibble(
      Label = c("IC50", "IC90"),
      X = NA_real_,
      X_min = NA_real_,
      X_max = NA_real_,
      !!group_col := combined_label
    ) %>%
      relocate(all_of(group_col))
  }
}


# 6. Combine plots into a 2-column x 4-row grid per page, A4 portrait

if (length(all_plots) > 0) {
  n_row <- 4
  n_col <- 2
  
  ml <- marrangeGrob(
    grobs = all_plots,
    nrow = n_row, ncol = n_col,
    layout_matrix = matrix(
      seq_len(n_row * n_col),
      nrow = n_row, ncol = n_col,
      byrow = TRUE 
    ),
    top = NULL  
  )
  
  # A4 portrait = 8.27 x 11.69 inches
  ggsave(output_pdf, ml, width = 8.27, height = 11.69, units = "in", dpi = 1200)
}


# 7. Combine and export the IC50/IC90 summary table across all datasets

if (length(all_plots) == 0) {
  
  message(sprintf(
    "No datasets fit successfully. No PDF or summary CSV written.\n%s",
    "Check the messages above for the reason each dataset was skipped."
  ))
  stop("No dose-response curves could be fitted.", call. = FALSE)
  
} else {
  
  ic_summary <- bind_rows(all_ic_tables) %>%
    select(all_of(group_col), Label, X, X_min, X_max) %>%
    rename(
      Estimate = X,
      CI_Lower = X_min,
      CI_Upper = X_max
    )
  
  write_csv(ic_summary, output_csv)
  
  message(sprintf("Done: %d dose-response curve(s) fit successfully. Plots written to '%s', summary written to '%s'.",
                  length(all_plots), output_pdf, output_csv))
}
