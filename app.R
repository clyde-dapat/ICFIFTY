library(shiny)
library(bslib)
library(readr)

app_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
methods <- c(
  "Focus Reduction Assay (FRA)" = "FRA",
  "Influenza Replication Inhibition Neuraminidase-based Assay (IRINA)" = "IRINA"
)

ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
  navbarPage(
    title = tags$span("ICFIFTY | Curve fitting analysis using the ", em("drc"), " R package"),
    windowTitle = "ICFIFTY | Curve fitting analysis using the drc R package",
    id = "main_nav",
    header = tags$head(tags$style(HTML("
      .run-button { margin: 8px 0 18px; }
      .table-wrap { overflow-x: auto; }
      .status { white-space: pre-wrap; }
      .upload-box {
        border: 2px dashed #7b8a97;
        border-radius: 8px;
        padding: 24px;
        margin: 14px 0 22px;
        min-height: 190px;
        background: #f8fafb;
        cursor: pointer;
        transition: border-color 0.15s ease, background-color 0.15s ease;
      }
      .upload-box:hover,
      .upload-box.dragover {
        border-color: #337ab7;
        background: #eef7ff;
      }
      .upload-box h4 { margin-top: 0; }
      .upload-list {
        border: 1px solid #d5dadd;
        border-radius: 4px;
        background: #f8fafb;
        padding: 12px 18px;
        margin: 8px 0 18px;
      }
      .upload-list ul { margin-bottom: 0; }
      .progress-text { margin-top: 8px; }
       #method label { width: auto; max-width: none; }
      .action-buttons {
        display: flex;
        flex-direction: column;
        align-items: flex-start;
        gap: 8px;
      }
    ")),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('reset-upload-inputs', function(ids) {
        ids.forEach(function(id) {
          var input = document.getElementById(id);
          if (!input) return;
          input.value = '';
          var container = input.closest('.shiny-input-container');
          if (container) {
            var fileLabel = container.querySelector('input[type=text]');
            if (fileLabel) fileLabel.value = '';
            var progress = container.querySelector('.shiny-file-input-progress');
            if (progress) {
              progress.classList.remove('active');
              var bar = progress.querySelector('.progress-bar');
              if (bar) {
                bar.style.width = '0%';
                bar.style.visibility = 'hidden';
              }
            }
          }
          Shiny.setInputValue(id, null, {priority: 'event'});
        });
      });
      $(document).on('click', '.upload-box', function(event) {
        if ($(event.target).closest('.shiny-input-container').length) return;
        var input = $(this).find('input[type=file]')[0];
        if (input) input.click();
      });
      $(document).on('dragenter dragover', '.upload-box', function(event) {
        event.preventDefault();
        event.stopPropagation();
        $(this).addClass('dragover');
      });
      $(document).on('dragleave drop', '.upload-box', function(event) {
        event.preventDefault();
        event.stopPropagation();
        $(this).removeClass('dragover');
      });
      $(document).on('drop', '.upload-box', function(event) {
        var input = $(this).find('input[type=file]')[0];
        var files = event.originalEvent.dataTransfer.files;
        if (!input || !files || files.length === 0) return;
        if (!input.multiple && files.length !== 1) {
          Shiny.setInputValue('upload_error',
            'Upload exactly one workbook to this box.', {priority: 'event'});
          return;
        }
        var acceptedExtensions = (input.accept || '')
          .split(',')
          .map(function(extension) { return extension.trim().toLowerCase(); })
          .filter(Boolean);
        var invalidFile = Array.from(files).find(function(file) {
          return acceptedExtensions.length &&
            !acceptedExtensions.some(function(extension) {
              return file.name.toLowerCase().endsWith(extension);
            });
        });
        if (invalidFile) {
          Shiny.setInputValue('upload_error',
            'One or more files have an unsupported file extension.',
            {priority: 'event'});
          return;
        }
        if (typeof DataTransfer === 'undefined') {
          Shiny.setInputValue('upload_error',
            'Drag-and-drop is not supported by this browser. Use Browse instead.',
            {priority: 'event'});
          return;
        }
        var transfer = new DataTransfer();
        Array.from(files).forEach(function(file) { transfer.items.add(file); });
        input.files = transfer.files;
        input.dispatchEvent(new Event('change', {bubbles: true}));
      });
    "))),
    tabPanel(
      "Input workbooks",
      value = "inputs",
      h4("Choose method and input files"),
      p("Upload the plate-reader spreadsheet(s), sample sheet, and dilution series for the selected method."),
      radioButtons(
        "method",
        "Analysis method",
        choices = methods,
        selected = "FRA",
        width = "100%"
      ),
      uiOutput("method_help"),
      hr(),
      h4("1. Plate results"),
      p(
        "Upload all plate spreadsheets for this run. Rename each file by adding,",
        tags$code("_Plate_n"), "at the end, where n is the plate number in the",
        "sample sheet. For example, ", tags$code("Experiment_Plate_1.xlsx"), " corresponds to Plate 1."
      ),
      div(
        class = "upload-box",
        fileInput(
          "plates",
          "Plate-reader workbooks",
          multiple = TRUE,
          accept = c(".xls", ".xlsx")
        ),
        uiOutput("plate_file_list")
      ),
      h4("2. Sample sheet"),
      p("Upload a sample sheet with Plate values correspond to the plate workbooks selected above."),
      div(
        class = "upload-box",
        fileInput("sample_sheet", "Sample sheet (.xlsx)", accept = ".xlsx"),
        uiOutput("sample_file_list")
      ),
      h4("3. Drug dilution series"),
      p("Upload the dilution workbook containing the dilution column and concentration columns."),
      div(
        class = "upload-box",
        fileInput("dilution_series", "Drug dilution series (.xlsx)", accept = ".xlsx"),
        uiOutput("dilution_file_list")
      ),
      div(
        class = "action-buttons",
        actionButton("run", "Run analysis", class = "btn-primary run-button"),
        actionButton("reset", "Reset", class = "btn-default")
      ),
      div(class = "status", textOutput("run_status"))
    ),
    tabPanel(
      "Summary",
      value = "summary",
      h2("IC50 / IC90 summary"),
      div(class = "table-wrap", tableOutput("summary_table"))
    ),
    tabPanel(
      "Dose-response plots",
      value = "curves",
      h2("Dose-response plots"),
      uiOutput("curve_view")
    ),
    tabPanel(
      "Downloads",
      value = "downloads",
      h2("Export results"),
      p("Download the plots, IC50/IC90 summary, or processed measurements."),
      downloadButton("download_pdf", "Download plots PDF"),
      downloadButton("download_csv", "Download IC50/IC90 CSV"),
      downloadButton("download_processed", "Download processed data CSV")
    )
  )
)

server <- function(input, output, session) {
  session_dir <- tempfile("icfifty-session-")
  if (!dir.create(session_dir, recursive = TRUE)) {
    stop("Unable to create a temporary analysis directory.")
  }

  state <- reactiveValues(result = NULL, status = "Choose a method and run the analysis.")
  resource_prefix <- NULL

  session$onSessionEnded(function() {
    if (!is.null(resource_prefix)) {
      removeResourcePath(resource_prefix)
    }
    unlink(session_dir, recursive = TRUE)
  })

  observeEvent(input$reset, {
    session$sendCustomMessage(
      "reset-upload-inputs",
      c("plates", "sample_sheet", "dilution_series")
    )
    updateRadioButtons(session, "method", selected = "FRA")
    state$result <- NULL
    state$status <- "Inputs reset. Choose a method and upload files to run an analysis."
    if (!is.null(resource_prefix)) {
      removeResourcePath(resource_prefix)
      resource_prefix <<- NULL
    }
    updateNavbarPage(session, "main_nav", selected = "inputs")
  }, ignoreInit = TRUE)

  observeEvent(input$upload_error, {
    showNotification(input$upload_error, type = "error")
  }, ignoreInit = TRUE)

  output$method_help <- renderUI({
    if (identical(input$method, "FRA")) {
      helpText(
        ""
      )
    } else {
      helpText(
        ""
      )
    }
  })

  file_list_ui <- function(files, empty_message) {
    if (is.null(files) || nrow(files) == 0L) {
      return(div(class = "upload-list text-muted", empty_message))
    }
    div(
      class = "upload-list",
      tags$strong(sprintf("%d file(s) selected:", nrow(files))),
      tags$ul(lapply(files$name, tags$li))
    )
  }
  output$plate_file_list <- renderUI({
    file_list_ui(
      input$plates,
      "No plate files uploaded."
    )
  })
  output$sample_file_list <- renderUI({
    file_list_ui(
      input$sample_sheet,
      "No sample sheet uploaded."
    )
  })
  output$dilution_file_list <- renderUI({
    file_list_ui(
      input$dilution_series,
      "No dilution series uploaded."
    )
  })

  copy_file <- function(source, destination, description) {
    if (!file.copy(source, destination, overwrite = TRUE)) {
      stop(sprintf("Could not copy %s into the session workspace.", description),
           call. = FALSE)
    }
  }

  copy_inputs <- function(run_dir) {
    data_dir <- file.path(run_dir, "data")
    dir.create(data_dir, recursive = TRUE)

    uploaded_names <- basename(input$plates$name)
    if (anyDuplicated(tolower(uploaded_names))) {
      stop("The plate upload contains duplicate filenames.", call. = FALSE)
    }
    for (i in seq_len(nrow(input$plates))) {
      copy_file(
        input$plates$datapath[i],
        file.path(data_dir, uploaded_names[i]),
        input$plates$name[i]
      )
    }

    sample_path <- file.path(run_dir, "sample_sheet.xlsx")
    copy_file(input$sample_sheet$datapath, sample_path, "the uploaded sample sheet")

    dilution_path <- file.path(run_dir, "Drug_Dilution_Series.xlsx")
    copy_file(
      input$dilution_series$datapath,
      dilution_path,
      "the uploaded dilution series"
    )

    list(data_dir = data_dir, sample_sheet = sample_path, dilution_series = dilution_path)
  }

  validate_sample_plate_mapping <- function(data_dir, sample_sheet, method) {
    samples <- readxl::read_excel(sample_sheet, .name_repair = "minimal")
    if (!"Plate" %in% names(samples)) {
      return(invisible(NULL))
    }

    sample_plate_ids <- suppressWarnings(as.numeric(samples$Plate))
    sample_plate_ids <- sort(unique(sample_plate_ids[is.finite(sample_plate_ids)]))
    filename_pattern <- if (method == "FRA") {
      "Plate_([0-9]+)\\.xls[x]?$"
    } else {
      "Plate_([0-9]+)\\.xlsx$"
    }
    plate_files <- list.files(
      data_dir,
      pattern = filename_pattern,
      full.names = TRUE,
      ignore.case = TRUE
    )
    plate_matches <- regmatches(
      basename(plate_files),
      regexec(filename_pattern, basename(plate_files), ignore.case = TRUE)
    )
    plate_ids <- sort(unique(as.numeric(vapply(
      plate_matches,
      function(match) if (length(match) >= 2L) match[[2]] else NA_character_,
      character(1)
    ))))
    plate_ids <- plate_ids[is.finite(plate_ids)]
    missing_ids <- setdiff(sample_plate_ids, plate_ids)

    if (length(missing_ids) > 0L) {
      stop(
        paste0(
          "The sample sheet refers to plate(s) ",
          paste(missing_ids, collapse = ", "),
          " that are not in the selected plate files (plate(s) ",
          paste(plate_ids, collapse = ", "),
          "). Upload all referenced plate files or upload a sample sheet ",
          "that matches the plates you selected. The FRA example sample sheet ",
          "maps to plates 1-6."
        ),
        call. = FALSE
      )
    }
    invisible(NULL)
  }

  check_analysis_packages <- function(method) {
    required_packages <- c("tidyverse", "readxl", "gridExtra")
    if (identical(method, "IRINA")) {
      required_packages <- c(required_packages, "drc")
    }
    missing_packages <- required_packages[
      !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_packages) > 0L) {
      stop(
        paste0(
          "Missing R package(s): ",
          paste(missing_packages, collapse = ", "),
          ". Activate the icfifty Conda environment and restart the app with ",
          "`Rscript run_app.R`."
        ),
        call. = FALSE
      )
    }
    invisible(NULL)
  }

  run_script <- function(script, arguments, log_path) {
    rscript <- file.path(R.home("bin"), "Rscript")
    output_lines <- system2(
      rscript,
      args = c(shQuote(script), vapply(arguments, shQuote, character(1))),
      stdout = TRUE,
      stderr = TRUE
    )
    exit_status <- attr(output_lines, "status")
    if (is.null(exit_status)) {
      exit_status <- 0L
    }
    writeLines(output_lines, log_path, useBytes = TRUE)
    if (exit_status != 0L) {
      details <- paste(tail(output_lines, 25L), collapse = "\n")
      stop(sprintf("Analysis script failed (exit status %s):\n%s", exit_status, details),
           call. = FALSE)
    }
  }

  results <- eventReactive(input$run, {
    state$status <- "Checking the selected method and input files..."
    state$result <- NULL
    withProgress(message = "Running analysis", value = 0, {
      tryCatch({
        method <- input$method
        missing_inputs <- c(
          if (is.null(input$plates) || nrow(input$plates) == 0L) {
            "plate-reader spreadsheet(s)"
          },
          if (is.null(input$sample_sheet) || nrow(input$sample_sheet) == 0L) {
            "sample sheet"
          },
          if (is.null(input$dilution_series) || nrow(input$dilution_series) == 0L) {
            "drug dilution series"
          }
        )
        if (length(missing_inputs) > 0L) {
          stop(
            paste0(
              "No files uploaded for: ",
              paste(missing_inputs, collapse = ", "),
              ". Upload a file in each required input box."
            ),
            call. = FALSE
          )
        }
        incProgress(0.05, detail = "Checking required packages")
        check_analysis_packages(method)
        run_id <- as.integer(input$run)
        run_dir <- file.path(session_dir, paste0("run-", run_id))
        result_dir <- file.path(run_dir, "result")
        dir.create(result_dir, recursive = TRUE)
        inputs <- copy_inputs(run_dir)
        incProgress(0.10, detail = "Checking selected input files")
        validate_sample_plate_mapping(inputs$data_dir, inputs$sample_sheet, method)

        processed_path <- file.path(result_dir, "processed_data.csv")
        pdf_path <- file.path(result_dir, "dose_response_curves.pdf")
        csv_path <- file.path(result_dir, "ic50_ic90_summary_table.csv")
        workflow_dir <- file.path(app_dir, method)

        incProgress(0.15, detail = "Validating workbook layouts")
        run_script(
          file.path(app_dir, "validate_workflow_inputs.R"),
          c(
            "--method", method,
            "--data-dir", inputs$data_dir,
            "--sample-sheet", inputs$sample_sheet,
            "--dilution-series", inputs$dilution_series
          ),
          file.path(run_dir, "input_validation.log")
        )
        incProgress(0.45, detail = "Processing plate measurements")
        run_script(
          file.path(workflow_dir, "data_processing.R"),
          c(
            "--data-dir", inputs$data_dir,
            "--sample-sheet", inputs$sample_sheet,
            "--dilution-series", inputs$dilution_series,
            "--output", processed_path
          ),
          file.path(run_dir, "data_processing.log")
        )
        incProgress(0.20, detail = "Fitting dose-response curves")
        run_script(
          file.path(workflow_dir, "dose_response_analysis.R"),
          c(
            "--input-csv", processed_path,
            "--output-pdf", pdf_path,
            "--output-csv", csv_path
          ),
          file.path(run_dir, "dose_response_analysis.log")
        )

        if (!file.exists(pdf_path) || !file.exists(csv_path) ||
            !file.exists(processed_path)) {
          stop("Analysis completed without producing all expected result files.",
               call. = FALSE)
        }

        if (!is.null(resource_prefix)) {
          removeResourcePath(resource_prefix)
        }
        resource_prefix <<- paste0(
          "icfifty_", gsub("[^A-Za-z0-9]", "", session$token), "_", run_id
        )
        addResourcePath(resource_prefix, result_dir)
        state$status <- sprintf(
          "%s analysis complete. Results are ready to view or download.",
          method
        )
        result <- list(
          method = method,
          result_dir = result_dir,
          processed_path = processed_path,
          pdf_path = pdf_path,
          csv_path = csv_path,
          summary = read_csv(csv_path, show_col_types = FALSE)
        )
        state$result <- result
        incProgress(0.05, detail = "Preparing results")
        updateNavbarPage(session, "main_nav", selected = "summary")
        result
      }, error = function(error) {
        state$result <- NULL
        state$status <- paste(
          "Analysis failed:",
          conditionMessage(error),
          "Check that the selected method matches your data and that you chose the correct files."
        )
        showNotification(state$status, type = "error", duration = NULL)
        NULL
      })
    })
  }, ignoreInit = TRUE)

  observeEvent(input$run, {
    results()
  }, ignoreInit = TRUE)

  output$run_status <- renderText({
    state$status
  })

  output$summary_table <- renderTable({
    result <- state$result
    req(result)
    result$summary
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$curve_view <- renderUI({
    result <- state$result
    if (is.null(result)) {
      return(p("Run an analysis to view the dose-response plots."))
    }
    tags$iframe(
      title = paste(result$method, "dose-response plots"),
      src = paste0("/", resource_prefix, "/dose_response_curves.pdf"),
      style = "width: 100%; height: 75vh; border: 1px solid #d5dadd;"
    )
  })

  output$download_pdf <- downloadHandler(
    filename = function() "dose_response_curves.pdf",
    content = function(file) {
      req(state$result)
      if (!file.copy(state$result$pdf_path, file)) {
        stop("Could not prepare the PDF download.")
      }
    },
    contentType = "application/pdf"
  )
  output$download_csv <- downloadHandler(
    filename = function() "ic50_ic90_summary_table.csv",
    content = function(file) {
      req(state$result)
      if (!file.copy(state$result$csv_path, file)) {
        stop("Could not prepare the summary CSV download.")
      }
    },
    contentType = "text/csv"
  )
  output$download_processed <- downloadHandler(
    filename = function() "processed_data.csv",
    content = function(file) {
      req(state$result)
      if (!file.copy(state$result$processed_path, file)) {
        stop("Could not prepare the processed data download.")
      }
    },
    contentType = "text/csv"
  )
}

shinyApp(ui = ui, server = server)
