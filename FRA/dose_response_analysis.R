library(tidyverse)

has_drc <- requireNamespace("drc", quietly = TRUE)

parse_cli <- function(args) {
  options <- list(
    `input-csv` = "result/processed_data.csv",
    `output-pdf` = "result/dose_response_curves.pdf",
    `output-csv` = "result/ic50_ic90_summary_table.csv"
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

fit_dose_response_model <- function(data,
                                    conc_col = "Concentration_nM",
                                    response_col = "Percent_Virus_Control") {
  fit_data <- data %>%
    filter(
      is.finite(.data[[conc_col]]),
      .data[[conc_col]] > 0,
      is.finite(.data[[response_col]])
    )
  if (nrow(fit_data) < 4L || n_distinct(fit_data[[conc_col]]) < 4L) {
    stop("At least four finite observations at distinct positive concentrations are required.")
  }
  if (has_drc) {
    return(drc::drm(
      stats::reformulate(conc_col, response = response_col),
      data = fit_data,
      fct = drc::LL.4()
    ))
  }

  x <- fit_data[[conc_col]]
  y <- fit_data[[response_col]]
  response_range <- max(y) - min(y)
  if (!is.finite(response_range) || response_range <= 0) {
    stop("The response values do not vary enough to fit a dose-response curve.")
  }
  log_x <- log(x)
  objective <- function(parameters) {
    lower <- parameters[[1]]
    span <- exp(parameters[[2]])
    log_ic50 <- parameters[[3]]
    slope <- exp(parameters[[4]])
    prediction <- lower + span / (1 + exp(slope * (log_x - log_ic50)))
    sum((y - prediction)^2)
  }
  fitted <- stats::optim(
    par = c(min(y), log(response_range), median(log_x), 0),
    fn = objective,
    method = "L-BFGS-B",
    lower = c(min(y) - 2 * response_range, log(response_range * 1e-4),
              min(log_x) - 8, log(0.05)),
    upper = c(max(y), log(response_range * 10), max(log_x) + 8, log(20)),
    control = list(maxit = 20000)
  )
  if (fitted$convergence != 0L || !is.finite(fitted$value)) {
    stop(sprintf("4-parameter logistic optimization failed (code %d).", fitted$convergence))
  }

  structure(
    list(
      parameters = fitted$par,
      objective = fitted$value,
      observations = length(y),
      x = x,
      y = y,
      conc_col = conc_col,
      response_col = response_col
    ),
    class = "frafifty_ll4"
  )
}

predict.frafifty_ll4 <- function(object, newdata, ...) {
  parameters <- object$parameters
  lower <- parameters[[1]]
  span <- exp(parameters[[2]])
  ic50 <- parameters[[3]]
  slope <- exp(parameters[[4]])
  concentration <- newdata[[object$conc_col]]
  lower + span / (1 + exp(slope * (log(concentration) - ic50)))
}

calculate_ic_values <- function(model) {
  if (inherits(model, "frafifty_ll4")) {
    parameters <- model$parameters
    ic50_log <- parameters[[3]]
    slope <- exp(parameters[[4]])
    ic_logs <- c(ic50_log, ic50_log + log(9) / slope)
    hessian <- tryCatch(
      stats::optimHess(parameters, function(value) {
        lower <- value[[1]]
        span <- exp(value[[2]])
        log_ic50 <- value[[3]]
        hill <- exp(value[[4]])
        prediction <- lower + span /
          (1 + exp(hill * (log(model$x) - log_ic50)))
        sum((model$y - prediction)^2)
      }),
      error = function(error) NULL
    )
    if (is.null(hessian)) {
      ci <- matrix(NA_real_, nrow = 2L, ncol = 2L)
    } else {
      covariance <- tryCatch(
        2 * model$objective / (model$observations - 4) * solve(hessian),
        error = function(error) NULL
      )
      if (is.null(covariance) || any(!is.finite(covariance))) {
        ci <- matrix(NA_real_, nrow = 2L, ncol = 2L)
      } else {
        gradients <- rbind(
          c(0, 0, 1, 0),
          c(0, 0, 1, -log(9) / slope)
        )
        standard_errors <- sqrt(pmax(
          0,
          rowSums((gradients %*% covariance) * gradients)
        ))
        ci <- cbind(
          exp(ic_logs - 1.96 * standard_errors),
          exp(ic_logs + 1.96 * standard_errors)
        )
      }
    }
    return(tibble(
      IC50_nM = exp(ic_logs[[1]]),
      IC50_CI_Lower_nM = ci[1, 1],
      IC50_CI_Upper_nM = ci[1, 2],
      IC90_nM = exp(ic_logs[[2]]),
      IC90_CI_Lower_nM = ci[2, 1],
      IC90_CI_Upper_nM = ci[2, 2]
    ))
  }

  estimates <- drc::ED(model, c(50, 90), interval = "delta", display = FALSE)
  estimate_data <- as.data.frame(estimates)
  if (!all(c("Estimate", "Lower", "Upper") %in% names(estimate_data))) {
    stop("The dose-response model did not return IC estimates with confidence intervals.")
  }
  tibble(
    IC50_nM = as.numeric(estimate_data$Estimate[[1]]),
    IC50_CI_Lower_nM = as.numeric(estimate_data$Lower[[1]]),
    IC50_CI_Upper_nM = as.numeric(estimate_data$Upper[[1]]),
    IC90_nM = as.numeric(estimate_data$Estimate[[2]]),
    IC90_CI_Lower_nM = as.numeric(estimate_data$Lower[[2]]),
    IC90_CI_Upper_nM = as.numeric(estimate_data$Upper[[2]])
  )
}

predict_curve_with_ci <- function(model, concentrations) {
  prediction_data <- data.frame(Concentration_nM = concentrations)
  if (!inherits(model, "frafifty_ll4")) {
    prediction <- tryCatch(
      stats::predict(
        model,
        newdata = prediction_data,
        interval = "confidence",
        level = 0.95
      ),
      error = function(error) NULL
    )
    if (!is.null(prediction) && is.matrix(prediction) &&
        all(c("Prediction", "Lower", "Upper") %in% colnames(prediction))) {
      return(tibble(
        Concentration_nM = concentrations,
        Response = prediction[, "Prediction"],
        Lower = prediction[, "Lower"],
        Upper = prediction[, "Upper"]
      ))
    }
  }

  if (inherits(model, "frafifty_ll4")) {
    parameters <- model$parameters
    lower <- parameters[[1]]
    span <- exp(parameters[[2]])
    log_ic50 <- parameters[[3]]
    slope <- exp(parameters[[4]])
    log_x <- log(concentrations)
    z <- slope * (log_x - log_ic50)
    exp_z <- exp(z)
    denominator <- 1 + exp_z
    response <- lower + span / denominator
    gradients <- cbind(
      1,
      span / denominator,
      span * slope * exp_z / denominator^2,
      -span * exp_z * z / denominator^2
    )
    hessian <- stats::optimHess(parameters, function(value) {
      fit <- value[[1]] + exp(value[[2]]) /
        (1 + exp(exp(value[[4]]) * (log(model$x) - value[[3]])))
      sum((model$y - fit)^2)
    })
    covariance <- tryCatch(
      2 * model$objective / (model$observations - 4) * solve(hessian),
      error = function(error) NULL
    )
    if (!is.null(covariance) && all(is.finite(covariance))) {
      standard_error <- sqrt(pmax(
        0,
        rowSums((gradients %*% covariance) * gradients)
      ))
      return(tibble(
        Concentration_nM = concentrations,
        Response = response,
        Lower = response - 1.96 * standard_error,
        Upper = response + 1.96 * standard_error
      ))
    }
    return(tibble(
      Concentration_nM = concentrations,
      Response = response,
      Lower = NA_real_,
      Upper = NA_real_
    ))
  }

  prediction <- as.numeric(stats::predict(model, newdata = prediction_data))
  tibble(
    Concentration_nM = concentrations,
    Response = prediction,
    Lower = NA_real_,
    Upper = NA_real_
  )
}

make_curve_plot <- function(data, model, ic_values, title,
                            show_interval_band = FALSE) {
  concentrations <- data$Concentration_nM[
    is.finite(data$Concentration_nM) & data$Concentration_nM > 0
  ]
  plot_x_min <- 10^-4
  plot_x_max <- max(10^4, max(concentrations) * 10)
  curve_concentrations <- exp(seq(
    log(min(concentrations)),
    log(max(concentrations) * 10),
    length.out = 500
  ))
  curve_data <- predict_curve_with_ci(model, curve_concentrations)
  curve_data <- curve_data %>%
    mutate(Lower = pmax(Lower, -20))
  ic_concentrations <- c(ic_values$IC50_nM, ic_values$IC90_nM)
  ic_responses <- as.numeric(stats::predict(
    model,
    newdata = data.frame(Concentration_nM = ic_concentrations)
  ))
  ic_markers <- tibble(
    Label = c("IC50", "IC90"),
    Concentration_nM = ic_concentrations,
    Response = ic_responses,
    CI_Lower_nM = c(
      ic_values$IC50_CI_Lower_nM,
      ic_values$IC90_CI_Lower_nM
    ),
    CI_Upper_nM = c(
      ic_values$IC50_CI_Upper_nM,
      ic_values$IC90_CI_Upper_nM
    )
  )
  ic_markers <- ic_markers %>%
    mutate(
      CI_Lower_nM = if_else(
        is.finite(CI_Lower_nM) & CI_Lower_nM > 0,
        CI_Lower_nM,
        NA_real_
      ),
      CI_Upper_nM = if_else(
        is.finite(CI_Upper_nM) & CI_Upper_nM > 0,
        CI_Upper_nM,
        NA_real_
      )
    ) %>%
    mutate(
      CI_Lower_nM = pmax(CI_Lower_nM, plot_x_min),
      CI_Upper_nM = pmin(CI_Upper_nM, plot_x_max)
    ) %>%
    filter(
      is.na(CI_Lower_nM) | is.na(CI_Upper_nM) |
        CI_Lower_nM <= CI_Upper_nM
    )

  y_values <- c(data$Percent_Virus_Control, curve_data$Response)
  if (show_interval_band) {
    y_values <- c(y_values, curve_data$Lower, curve_data$Upper)
  }
  y_min <- -20
  y_max <- ceiling(max(c(110, y_values), na.rm = TRUE) / 10) * 10
  x_min <- min(concentrations)
  ic_text <- sprintf(
    "IC50: %.3g nM (%.3g-%.3g)\nIC90: %.3g nM (%.3g-%.3g)",
    ic_values$IC50_nM,
    ic_values$IC50_CI_Lower_nM,
    ic_values$IC50_CI_Upper_nM,
    ic_values$IC90_nM,
    ic_values$IC90_CI_Lower_nM,
    ic_values$IC90_CI_Upper_nM
  )

  plot <- ggplot(data, aes(
    x = Concentration_nM,
    y = Percent_Virus_Control
  ))
  if (show_interval_band) {
    plot <- plot +
      geom_ribbon(
      data = curve_data,
      aes(x = Concentration_nM, ymin = Lower, ymax = Upper),
      inherit.aes = FALSE,
      fill = "grey75",
      alpha = 0.45
      )
  }

  plot +
    geom_line(
      data = curve_data,
      aes(x = Concentration_nM, y = Response),
      inherit.aes = FALSE,
      color = "grey60",
      linewidth = 0.8
    ) +
    geom_point(color = "grey15", size = 1.4, alpha = 0.9) +
    geom_segment(
      data = ic_markers,
      aes(x = x_min, xend = Concentration_nM, y = Response, yend = Response),
      inherit.aes = FALSE,
      color = "#ff6666",
      linetype = "dashed",
      linewidth = 0.35
    ) +
    geom_segment(
      data = ic_markers,
      aes(x = Concentration_nM, xend = Concentration_nM, y = y_min, yend = Response),
      inherit.aes = FALSE,
      color = "#ff6666",
      linetype = "dashed",
      linewidth = 0.35
    ) +
    geom_point(
      data = ic_markers,
      aes(x = Concentration_nM, y = Response),
      inherit.aes = FALSE,
      color = "blue",
      size = 2
    ) +
    geom_errorbarh(
      data = ic_markers %>% filter(
        is.finite(CI_Lower_nM),
        is.finite(CI_Upper_nM)
      ),
      aes(
        y = Response,
        xmin = CI_Lower_nM,
        xmax = CI_Upper_nM
      ),
      inherit.aes = FALSE,
      color = "blue",
      linewidth = 0.55,
      height = 1.8
    ) +
    geom_text(
      data = ic_markers,
      aes(x = Concentration_nM, y = Response, label = Label),
      inherit.aes = FALSE,
      color = "blue",
      fontface = "bold",
      size = 2.7,
      vjust = -0.65
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = ic_text,
      hjust = 1.02,
      vjust = 1.15,
      size = 2.5,
      color = "blue4",
      fontface = "bold",
      lineheight = 1.05
    ) +
    scale_x_log10(
      breaks = 10^seq(-4, 4),
      limits = c(
        plot_x_min,
        plot_x_max
      ),
      labels = function(x) {
        parse(text = sprintf("10^%d", as.integer(round(log10(x)))))
      }
    ) +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = scales::breaks_pretty(n = 5),
      expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
      title = title,
      x = "Drug concentration [nM]",
      y = "% spot count vs control"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
      plot.margin = margin(5, 6, 12, 5),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.35)
    )
}

make_failed_plot <- function(data, title, reason) {
  ggplot(data, aes(
    x = Concentration_nM,
    y = Percent_Virus_Control,
    color = Replicate
  )) +
    geom_point(size = 2.7, alpha = 0.85) +
    scale_x_log10(labels = scales::label_number()) +
    scale_color_brewer(palette = "Dark2", name = "Replicate") +
    labs(
      title = title,
      subtitle = paste("Model could not be fitted:", reason),
      x = "Drug concentration [nM]",
      y = "% spot count vs control"
    ) +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

analyze_dose_response <- function(data,
                                  output_pdf = NULL,
                                  output_csv = NULL) {
  required_columns <- c(
    "Dataset", "Sample_Group", "Replicate", "Drug", "Plate", "Row",
    "Concentration_nM", "Percent_Virus_Control"
  )
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(sprintf("Processed data is missing columns: %s",
                 paste(missing_columns, collapse = ", ")), call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("Processed data contains no measurements.", call. = FALSE)
  }
  data <- data %>%
    mutate(
      Sample = as.character(Dataset),
      Replicate = as.character(Replicate),
      Row = as.character(Row)
    )

  replicate_keys <- data %>%
    distinct(Plate, Row, Sample_Group, Sample = Dataset, Drug, Replicate) %>%
    mutate(
      .row_order = match(Row, LETTERS),
      .replicate_order = suppressWarnings(as.numeric(Replicate))
    ) %>%
    arrange(Plate, .row_order, .replicate_order, Sample) %>%
    select(-.row_order, -.replicate_order)
  replicate_data <- lapply(seq_len(nrow(replicate_keys)), function(i) {
    key <- replicate_keys[i, ]
    data %>%
      filter(
        Plate == key$Plate[[1]],
        Row == key$Row[[1]],
        Sample_Group == key$Sample_Group[[1]],
        Dataset == key$Sample[[1]],
        Drug == key$Drug[[1]],
        Replicate == key$Replicate[[1]]
      )
  })

  sample_order <- data %>%
    distinct(Plate, Row, Sample_Group, Drug) %>%
    mutate(.row_order = match(Row, LETTERS)) %>%
    group_by(Plate, Sample_Group, Drug) %>%
    summarise(.row_order = min(.row_order), .groups = "drop")
  combined_keys <- data %>%
    distinct(Plate, Sample_Group, Drug) %>%
    left_join(sample_order, by = c("Plate", "Sample_Group", "Drug")) %>%
    arrange(Plate, .row_order, Sample_Group) %>%
    select(-.row_order)
  combined_data <- lapply(seq_len(nrow(combined_keys)), function(i) {
    key <- combined_keys[i, ]
    data %>%
      filter(
        Plate == key$Plate[[1]],
        Sample_Group == key$Sample_Group[[1]],
        Drug == key$Drug[[1]]
      )
  })

  curve_data <- c(replicate_data, combined_data)
  curve_keys <- bind_rows(
    replicate_keys %>% mutate(Curve_Type = "Replicate"),
    combined_keys %>%
      mutate(
        Row = NA_character_,
        Sample = Sample_Group,
        Replicate = "Combined",
        Curve_Type = "Combined"
      ) %>%
      select(all_of(names(replicate_keys)), Curve_Type)
  ) %>%
    mutate(Curve_ID = paste(Curve_Type, Plate, Sample_Group, Replicate, sep = "_"))

  plots <- vector("list", length(curve_data))
  summaries <- vector("list", length(curve_data))

  for (i in seq_along(curve_data)) {
    curve <- curve_data[[i]]
    key <- curve_keys[i, ]
    curve_type <- key$Curve_Type
    title <- if (curve_type == "Combined") {
      paste0(key$Sample_Group, " - combined replicates")
    } else {
      as.character(key$Sample)
    }

    fit_result <- tryCatch({
      model <- fit_dose_response_model(curve)
      list(model = model, ic_values = calculate_ic_values(model), error = NULL)
    }, error = function(error) {
      list(model = NULL, ic_values = NULL, error = conditionMessage(error))
    })

    if (is.null(fit_result$error)) {
      plots[[i]] <- make_curve_plot(
        curve,
        fit_result$model,
        fit_result$ic_values,
        title,
        show_interval_band = curve_type == "Combined"
      )
      summaries[[i]] <- bind_cols(
        key,
        fit_result$ic_values,
        tibble(Fit_Status = "Success", Failure_Reason = NA_character_)
      )
    } else {
      message(sprintf("Could not fit %s: %s", title, fit_result$error))
      plots[[i]] <- make_failed_plot(curve, title, fit_result$error)
      summaries[[i]] <- bind_cols(
        key,
        tibble(
          IC50_nM = NA_real_,
          IC50_CI_Lower_nM = NA_real_,
          IC50_CI_Upper_nM = NA_real_,
          IC90_nM = NA_real_,
          IC90_CI_Lower_nM = NA_real_,
          IC90_CI_Upper_nM = NA_real_,
          Fit_Status = "Failed",
          Failure_Reason = fit_result$error
        )
      )
    }
  }

  summary_wide <- bind_rows(summaries) %>%
    select(
      Curve_Type, Plate, Row, Sample_Group, Sample, Drug, Replicate,
      IC50_nM, IC50_CI_Lower_nM, IC50_CI_Upper_nM,
      IC90_nM, IC90_CI_Lower_nM, IC90_CI_Upper_nM,
      Fit_Status, Failure_Reason
    )
  summary_metadata <- summary_wide %>%
    select(
      Curve_Type, Plate, Row, Sample_Group, Sample, Drug, Replicate,
      Fit_Status, Failure_Reason
    )
  summary <- bind_rows(
    summary_metadata %>%
      bind_cols(summary_wide %>% transmute(
        IC_Type = "IC50",
        Estimate_nM = IC50_nM,
        CI_Lower_nM = IC50_CI_Lower_nM,
        CI_Upper_nM = IC50_CI_Upper_nM
      )),
    summary_metadata %>%
      bind_cols(summary_wide %>% transmute(
        IC_Type = "IC90",
        Estimate_nM = IC90_nM,
        CI_Lower_nM = IC90_CI_Lower_nM,
        CI_Upper_nM = IC90_CI_Upper_nM
      ))
  ) %>%
    select(
      Curve_Type, Plate, Row, Sample_Group, Sample, Drug, Replicate,
      IC_Type, Estimate_nM, CI_Lower_nM, CI_Upper_nM,
      Fit_Status, Failure_Reason
    ) %>%
    arrange(Plate, Row, Replicate, IC_Type)

  if (!is.null(output_pdf)) {
    dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
    plot_pages <- gridExtra::marrangeGrob(
      grobs = plots,
      nrow = 4,
      ncol = 2,
      layout_matrix = matrix(seq_len(8), nrow = 4, ncol = 2, byrow = TRUE),
      top = NULL
    )
    ggplot2::ggsave(
      output_pdf,
      plot_pages,
      width = 8.27,
      height = 11.69,
      units = "in",
      limitsize = FALSE
    )
  }
  if (!is.null(output_csv)) {
    dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(
      select(summary, -Fit_Status, -Failure_Reason),
      output_csv
    )
  }

  list(summary = summary, plots = plots, plot_labels = curve_keys$Curve_ID)
}

if (sys.nframe() == 0L) {
  options <- parse_cli(commandArgs(trailingOnly = TRUE))
  input_csv <- options[["input-csv"]]
  if (!file.exists(input_csv)) {
    stop(sprintf("Processed input CSV does not exist: %s", input_csv), call. = FALSE)
  }
  processed <- readr::read_csv(input_csv, show_col_types = FALSE)
  result <- analyze_dose_response(
    processed,
    output_pdf = options[["output-pdf"]],
    output_csv = options[["output-csv"]]
  )
  message(sprintf(
    "Analyzed %d curves (%d successful); plots: %s; summary: %s.",
    length(result$plots),
    sum(result$summary$Fit_Status == "Success" & result$summary$IC_Type == "IC50"),
    options[["output-pdf"]],
    options[["output-csv"]]
  ))
}
