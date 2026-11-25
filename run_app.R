script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_args) != 1L) {
  stop("Could not determine the app launcher path.", call. = FALSE)
}

app_dir <- normalizePath(
  dirname(sub("^--file=", "", script_args)),
  winslash = "/",
  mustWork = TRUE
)
if (!file.exists(file.path(app_dir, "app.R"))) {
  stop(sprintf("No app.R found in %s.", app_dir), call. = FALSE)
}

host <- Sys.getenv("ICFIFTY_SHINY_HOST", unset = "127.0.0.1")
port_text <- Sys.getenv("ICFIFTY_SHINY_PORT", unset = "3838")
port <- suppressWarnings(as.integer(port_text))
if (is.na(port) || port < 1L || port > 65535L) {
  stop("ICFIFTY_SHINY_PORT must be an integer from 1 to 65535.", call. = FALSE)
}

setwd(app_dir)
shiny::runApp(
  appDir = app_dir,
  host = host,
  port = port,
  launch.browser = FALSE
)
