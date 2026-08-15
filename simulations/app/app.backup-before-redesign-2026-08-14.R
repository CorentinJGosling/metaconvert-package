## =============================================================================
## metaConvert simulation results
##
## Browses the ADEMP aggregates written by run_study() into
## simulations/data/aggregated/*.csv. Entirely generic: every study is discovered
## from disk and its condition columns, methods, targets and replication counts
## are read from the file, so a new study appears here as soon as it has been run.
##
## Visual identity follows metaconvert.org: rose accent, indigo for data, warm
## neutrals, Work Sans + Cascadia Code. The brand contains no green.
##
## Run from the simulations/ root:   shiny::runApp("app")
## =============================================================================

library(shiny)
library(bslib)
library(ggplot2)
library(DT)

## ---- brand -------------------------------------------------------------------
MC <- list(
  paper  = "#f7f7f3",
  white  = "#ffffff",
  ink    = "#211f38",
  ink2   = "#2c2a47",
  soft   = "#6b6a7d",
  rule   = "#e6e6e0",
  rose   = "#d66268",
  rose_d = "#a23c42",
  amber  = "#8a6712"
)

## Categorical scale for methods. Anchored on the brand's rose and indigo and
## extended with amber, plum, blue and brown; deliberately green-free.
##
## Ordered so that neighbouring positions are far apart in hue. Methods are taken
## alphabetically, so the two rose tones must not land on adjacent slots -- with
## the deep rose in position 4, study 09's "2x2_tetrachoric" and "lipsey_cooper"
## came out as two near-identical reds.
MC_METHODS <- c("#d66268", "#2c2a47", "#8a6712", "#6b4370",
                "#4a6fa5", "#8a5a3c", "#a23c42", "#c98a5b",
                "#3f4a6b", "#b07d4f", "#7d5e8c", "#c4525f")

## ---- locate the aggregates ---------------------------------------------------
.find_agg_dir <- function() {
  cands <- c(file.path("..", "data", "aggregated"),
             file.path("data", "aggregated"),
             file.path(dirname(getwd()), "data", "aggregated"))
  for (p in cands) if (dir.exists(p)) return(normalizePath(p, winslash = "/"))
  stop("could not locate simulations/data/aggregated")
}
AGG_DIR <- .find_agg_dir()

## Human-readable titles. A study missing from this list still works; it simply
## shows its file stem.
STUDY_LABELS <- c(
  "01a_smd_to_cor_r"      = "SMD to correlation (r)",
  "01b_smd_to_cor_z"      = "SMD to correlation (z)",
  "02a_cor_to_smd_GROUPS" = "Correlation to SMD (grouped data)",
  "02b_cor_to_smd_CONT"   = "Correlation to SMD (continuous data)",
  "03a_2x2_to_cor_CAT"    = "2x2 table to correlation (categorical latent)",
  "03b_2x2_to_cor_CONT"   = "2x2 table to correlation (continuous latent)",
  "04_or_to_rr"           = "Odds ratio to risk ratio",
  "05_rr_to_or"           = "Risk ratio to odds ratio",
  "06_or_se_imputation"   = "Odds ratio standard error imputation",
  "07a_ancova_to_smd_d"   = "ANCOVA to SMD (Cohen's d)",
  "07b_ancova_to_smd_g"   = "ANCOVA to SMD (Hedges' g)",
  "08a_pre_post_to_smd_d" = "Pre-post to SMD (Cohen's d)",
  "08b_pre_post_to_smd_g" = "Pre-post to SMD (Hedges' g)",
  "09a_or_to_cor_CONT"    = "Odds ratio to correlation (continuous latent)",
  "09b_or_to_cor_CAT"     = "Odds ratio to correlation (categorical latent)"
)

## The package argument each study evaluates. Shown beside the title so the app
## answers "which setting am I looking at?" without opening the study file.
STUDY_ARG <- c(
  "01a_smd_to_cor_r" = "smd_to_cor", "01b_smd_to_cor_z" = "smd_to_cor",
  "02a_cor_to_smd_GROUPS" = "cor_to_smd", "02b_cor_to_smd_CONT" = "cor_to_smd",
  "03a_2x2_to_cor_CAT" = "table_2x2_to_cor", "03b_2x2_to_cor_CONT" = "table_2x2_to_cor",
  "04_or_to_rr" = "or_to_rr", "05_rr_to_or" = "rr_to_or",
  "06_or_se_imputation" = "es_from_or", "07a_ancova_to_smd_d" = "cov_outcome_r",
  "07b_ancova_to_smd_g" = "cov_outcome_r", "08a_pre_post_to_smd_d" = "pre_post_to_smd",
  "08b_pre_post_to_smd_g" = "pre_post_to_smd",
  "09a_or_to_cor_CONT" = "or_to_cor", "09b_or_to_cor_CAT" = "or_to_cor"
)

## ---- performance measures ----------------------------------------------------
METRICS <- list(
  bias        = list(label = "Bias",                        mcse = "bias_mcse",     ref = 0,    better = "closer to 0"),
  rmse        = list(label = "Root mean squared error",     mcse = NA,              ref = 0,    better = "smaller"),
  coverage    = list(label = "Confidence interval coverage", mcse = "coverage_mcse", ref = 0.95, better = "closer to 0.95"),
  se_ratio    = list(label = "SE ratio (model / empirical)", mcse = "se_ratio_mcse", ref = 1,    better = "closer to 1"),
  emp_se      = list(label = "Empirical standard error",    mcse = "emp_se_mcse",   ref = NA,   better = "smaller"),
  ci_width    = list(label = "Mean interval width",         mcse = NA,              ref = NA,   better = "smaller"),
  nonest_rate = list(label = "Non-estimable rate",          mcse = NA,              ref = 0,    better = "smaller")
)

## ---- discovery ---------------------------------------------------------------
scan_files <- function() {
  f <- list.files(AGG_DIR, pattern = "_nrep[0-9]+[.]csv$", full.names = TRUE)
  if (!length(f)) return(NULL)
  stem <- sub("_nrep[0-9]+[.]csv$", "", basename(f))
  nrep <- as.integer(sub(".*_nrep([0-9]+)[.]csv$", "\\1", basename(f)))
  d <- data.frame(path = f, stem = stem, nrep = nrep, stringsAsFactors = FALSE)
  d <- d[order(d$stem, -d$nrep), ]
  d$label <- ifelse(d$stem %in% names(STUDY_LABELS), STUDY_LABELS[d$stem], d$stem)
  ## the leading number of the study file, used as the sidebar's index column
  d$idx <- toupper(sub("^([0-9]+[a-z]?)_.*$", "\\1", d$stem))
  d
}
FILES <- scan_files()

condition_cols <- function(df) {
  i <- match("method", names(df))
  if (is.na(i) || i < 2) return(character(0))
  names(df)[seq_len(i - 1)]
}

## ---- plot theme --------------------------------------------------------------
theme_mc <- function() {
  theme_minimal(base_size = 13, base_family = "Work Sans") +
    theme(
      text             = element_text(colour = MC$ink),
      axis.title       = element_text(size = 11, colour = MC$soft),
      axis.text        = element_text(size = 10, colour = MC$soft, family = "Cascadia Code"),
      panel.grid.major = element_line(colour = MC$rule, linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(1.1, "lines"),
      strip.background = element_rect(fill = MC$paper, colour = NA),
      strip.text       = element_text(size = 10.5, colour = MC$ink, face = "bold",
                                      margin = margin(5, 6, 5, 6)),
      legend.position  = "top",
      legend.justification = "left",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 10.5),
      legend.key.width = unit(1.4, "lines"),
      legend.key.height = unit(0.9, "lines"),
      legend.margin    = margin(0, 0, 2, 0),
      legend.box.spacing = unit(2, "pt"),
      plot.margin      = margin(2, 10, 2, 2)
    )
}

## ---- UI ----------------------------------------------------------------------
## The stylesheet is a PLAIN string, never a sprintf() format. CSS is full of
## literal per-cent signs (width:100%, background-size:98% 62%, @media rules),
## and inside a format string each one is read as a conversion specification and
## aborts the app at start-up with "too few arguments". The custom properties are
## built from the MC palette by pasting, so the two still cannot drift.
CSS_VARS <- list(paper = MC$paper, ink = MC$ink, ink2 = MC$ink2, soft = MC$soft,
                 rule = MC$rule, rose = MC$rose, "rose-d" = MC$rose_d,
                 amber = MC$amber)

app_css <- paste0(
":root{",
paste0("--", names(CSS_VARS), ":", unlist(CSS_VARS), ";", collapse = ""),
'--mono:"Cascadia Code","Cascadia Mono",ui-monospace,SFMono-Regular,Menlo,monospace;}
body{background:var(--paper);}

/* ---- header ---- */
.navbar{
  background:#fff !important; border-bottom:1px solid var(--rule);
  box-shadow:none; min-height:58px;
}
.navbar .navbar-brand{
  font-weight:700; font-size:1.05rem; letter-spacing:-.01em; color:var(--ink) !important;
}
.navbar .navbar-brand .sub{
  font-weight:400; color:var(--soft); margin-left:.55rem;
  padding-left:.55rem; border-left:1px solid var(--rule);
}

/* ---- sidebar ---- */
.sidebar, .bslib-sidebar-layout>.sidebar{ background:#fff; border-right:1px solid var(--rule); }
.sidebar .form-label, .sidebar label{ font-weight:600; font-size:.82rem; color:var(--ink); }
.grp{
  font-size:.68rem; letter-spacing:.1em; text-transform:uppercase;
  color:var(--soft); font-weight:700; margin:1.35rem 0 .5rem;
  padding-bottom:.35rem; border-bottom:1px solid var(--rule);
}
.grp:first-of-type{ margin-top:.15rem; }
.form-select, .form-control{
  border-color:var(--rule); border-radius:7px; font-size:.88rem; color:var(--ink);
}
.form-select:focus,.form-control:focus{
  border-color:var(--rose); box-shadow:0 0 0 3px rgba(214,98,104,.16);
}
.form-check-input:checked{ background-color:var(--rose); border-color:var(--rose); }
.form-check-input:focus{ box-shadow:0 0 0 3px rgba(214,98,104,.16); border-color:var(--rose); }
.form-check-label{ font-weight:400 !important; font-size:.86rem; }

/* the target help is reference material, not a permanent fixture */
details.note{ margin:.4rem 0 .2rem; }
details.note>summary{
  cursor:pointer; font-size:.78rem; color:var(--rose-d); font-weight:600;
  list-style:none; padding:.15rem 0;
}
details.note>summary::-webkit-details-marker{ display:none; }
details.note>summary::before{ content:"+ "; font-family:var(--mono); }
details.note[open]>summary::before{ content:"- "; }
details.note .body{
  font-size:.79rem; line-height:1.5; color:var(--soft);
  border-left:2px solid var(--rule); padding:.1rem 0 .1rem .7rem; margin-top:.3rem;
}
details.note .body b{ color:var(--ink); }
.hint{ font-size:.75rem; line-height:1.45; color:var(--soft); margin:-.35rem 0 .2rem; }

/* ---- tabs ---- */
.nav-tabs{ border-bottom:1px solid var(--rule); gap:.15rem; }
.nav-tabs .nav-link{
  border:0; color:var(--soft); font-weight:600; font-size:.9rem;
  padding:.65rem .95rem; border-bottom:2px solid transparent; border-radius:0;
}
.nav-tabs .nav-link:hover{ color:var(--ink); background:transparent; }
.nav-tabs .nav-link.active{
  color:var(--ink) !important; background:transparent !important;
  border-bottom:2px solid var(--rose) !important;
}
.card{ border:1px solid var(--rule); border-radius:12px; box-shadow:0 1px 2px rgba(33,31,56,.04); }

/* ---- reading strip above the plot ---- */
.strip{
  display:flex; flex-wrap:wrap; align-items:baseline; gap:.5rem .85rem;
  padding:0 0 .85rem; margin-bottom:.9rem; border-bottom:1px solid var(--rule);
}
.strip h2{
  font-size:1.12rem; font-weight:700; color:var(--ink); margin:0; letter-spacing:-.01em;
}
.strip .arg{
  font-family:var(--mono); font-size:.8rem; color:var(--rose-d);
  background:rgba(214,98,104,.09); padding:.14rem .45rem; border-radius:5px;
}
.strip .facts{
  margin-left:auto; font-family:var(--mono); font-size:.76rem; color:var(--soft);
}
.strip .facts b{ color:var(--ink2); font-weight:600; }
.lede{ font-size:.85rem; color:var(--soft); margin:.1rem 0 1rem; max-width:78ch; line-height:1.55; }
.lede code{ font-family:var(--mono); font-size:.8rem; color:var(--ink2); background:var(--paper);
  padding:.05rem .3rem; border-radius:4px; }

/* one compact line under the title: what is plotted, and the band key */
.subline{
  display:flex; flex-wrap:wrap; align-items:center; gap:.5rem 1.1rem;
  margin:-.35rem 0 .55rem;
}
.subline .l{ font-size:.83rem; color:var(--soft); }
.band-key{
  display:inline-flex; align-items:center; gap:.4rem;
  font-size:.76rem; color:var(--soft); font-family:var(--mono);
  cursor:help;
}
.band-key i{
  width:26px; height:11px; display:inline-block; border-radius:2px;
  background:rgba(44,42,71,.10); border-top:1px dashed rgba(44,42,71,.45);
  border-bottom:1px dashed rgba(44,42,71,.45);
}

/* ---- tables ---- */
table.dataTable{ font-size:.84rem; }
table.dataTable thead th{
  font-size:.72rem; text-transform:uppercase; letter-spacing:.05em;
  color:var(--soft); border-bottom:1px solid var(--rule) !important;
}
table.dataTable td{ font-family:var(--mono); color:var(--ink2); }
table.dataTable td:first-child{ font-family:"Work Sans"; font-weight:600; color:var(--ink); }
.dataTables_wrapper .dataTables_filter input{ border:1px solid var(--rule); border-radius:6px; }
.btn{ border-radius:7px; font-weight:600; font-size:.84rem; }
.btn-default,.btn-outline-secondary{ border-color:var(--rule); color:var(--ink); }
/* shiny gives downloadButton width:100% via .btn-block-ish defaults */
#dl{ display:inline-block; width:auto; }
')

TARGET_HELP <- HTML(
  '<div class="body">
   <b>Shared benchmarks</b> — one quantity, applied to every method. Use these to
   compare methods against each other.<br><br>
   &bull; <b>Population parameter</b> — the value the data were generated from. This
   is the quantity you actually want, so it is the default, and coverage and
   standard-error calibration are meaningful only here.<br>
   &bull; <b>Same-sample statistic</b> — the same statistic recomputed on that
   replication&rsquo;s own sample: how faithful is the conversion? Being random, it
   has no 0.95 reference, so coverage cannot be read against it.<br><br>
   <b>Per-method benchmark</b><br>
   &bull; <b>Each method&rsquo;s own estimand</b> — every method scored against a
   <i>different</i> quantity, the one it is designed to estimate. This is a
   diagnostic, not a ranking: it answers &ldquo;does this method compute its own
   value correctly?&rdquo;, never &ldquo;is that the value I want?&rdquo;. Two methods
   scored here are being held to two different standards.<br><br>
   The gap between the two is the estimand mismatch, and it is not bias — the
   method is computing its own quantity correctly, but that quantity is not the one
   being asked for. The <b>Estimand check</b> tab shows both at once, which is the
   only place the comparison is safe.</div>')

ui <- page_sidebar(
  title = span(class = "navbar-brand", "metaConvert",
               span(class = "sub", "simulation results")),
  window_title = "metaConvert simulation results",
  theme = bs_theme(
    version = 5,
    bg = MC$white, fg = MC$ink, primary = MC$rose,
    base_font = font_google("Work Sans"),
    code_font = font_collection("Cascadia Code", "Cascadia Mono", font_google("IBM Plex Mono")),
    "border-color" = MC$rule
  ),
  tags$head(tags$style(HTML(app_css))),

  sidebar = sidebar(
    width = 320, bg = MC$white,
    if (is.null(FILES)) {
      div(class = "text-danger", "No aggregate files found in data/aggregated.")
    } else {
      tagList(
        div(class = "grp", "Study"),
        selectInput("study", NULL,
                    choices = stats::setNames(unique(FILES$stem),
                                              sprintf("%s — %s",
                                                      FILES$idx[!duplicated(FILES$stem)],
                                                      FILES$label[!duplicated(FILES$stem)]))),
        uiOutput("nrep_ui"),

        div(class = "grp", "Scored against"),
        uiOutput("target_ui"),
        tags$details(class = "note",
                     tags$summary("What the targets mean"), TARGET_HELP),

        div(class = "grp", "Measure"),
        selectInput("metric", NULL,
                    choices = stats::setNames(names(METRICS),
                                              vapply(METRICS, `[[`, "", "label")),
                    selected = "coverage"),

        div(class = "grp", "Panels"),
        uiOutput("xvar_ui"),
        uiOutput("facet_ui"),
        checkboxInput("free_y", "Free vertical scale per panel", FALSE),
        div(class = "hint", "A shared scale is the default so panels can be compared."),

        uiOutput("scale_ui"),

        div(class = "grp", "Methods"),
        uiOutput("methods_ui"),

        uiOutput("filters_ui")
      )
    }
  ),

  navset_card_tab(
    id = "tabs",
    nav_panel("Plot",       card_body(uiOutput("strip"), plotOutput("main_plot", height = "700px"))),
    nav_panel("Ranking",    card_body(uiOutput("rank_lede"), DTOutput("rank_tbl"))),
    nav_panel("Estimand check", card_body(uiOutput("estimand_lede"), plotOutput("estimand_plot", height = "520px"))),
    nav_panel("Data",       card_body(uiOutput("data_lede"),
                                      downloadButton("dl", "Download these rows", class = "btn-outline-secondary btn-sm"),
                                      br(), br(), DTOutput("data_tbl")))
  )
)

## ---- server ------------------------------------------------------------------
server <- function(input, output, session) {

  output$nrep_ui <- renderUI({
    req(input$study)
    n <- sort(unique(FILES$nrep[FILES$stem == input$study]), decreasing = TRUE)
    selectInput("nrep", "Replications per condition", choices = n, selected = n[1])
  })

  raw <- reactive({
    req(input$study, input$nrep)
    p <- FILES$path[FILES$stem == input$study & FILES$nrep == as.integer(input$nrep)]
    req(length(p) == 1)
    utils::read.csv(p, stringsAsFactors = FALSE)
  })

  cond_cols <- reactive(condition_cols(raw()))

  ## Target names are study-specific: most studies record own/population/sample,
  ## but study 01 records biserial_* and pointbiserial_* separately and study 08
  ## adds own_true_r. Read them from the file rather than assuming.
  ## Targets are of two different kinds, and presenting them as one flat list
  ## invites the mistake of ranking methods against `own`. A SHARED target is one
  ## benchmark applied to every method, and is what you need to choose between
  ## them. `own` is a DIFFERENT benchmark per method -- it measures whether each
  ## method computes its own quantity correctly, which is a diagnostic, not a
  ## basis for comparison. Shared targets are listed first and one of them is the
  ## default; `own` is separated and labelled for what it is.
  output$target_ui <- renderUI({
    tg <- unique(raw()$target)
    own_like <- grep("^own", tg, value = TRUE)
    shared <- setdiff(tg, own_like)
    pop <- grep("population", shared, value = TRUE)
    ord <- c(sort(pop), sort(setdiff(shared, pop)), sort(own_like))
    pretty <- vapply(ord, function(z) {
      if (z == "own") "Each method's own estimand"
      else if (z == "own_true_r") "Own estimand, at the true r"
      else if (z == "population") "Population parameter"
      else if (z == "sample") "Same-sample statistic"
      else gsub("_", " ", z)
    }, character(1))
    tagList(
      radioButtons("target", NULL,
                   choices = stats::setNames(ord, pretty),
                   selected = if (length(pop)) sort(pop)[1] else ord[1]),
      div(class = "hint",
          if (length(own_like))
            HTML(paste0("The first option", if (length(shared) > 1) "s are" else " is",
                        " a <b>shared</b> benchmark, the same quantity for every method — ",
                        "use that to compare methods. <b>Each method's own estimand</b> ",
                        "scores every method against a different quantity, so it says ",
                        "whether a method computes <i>its</i> value correctly, not whether ",
                        "that value is the one you want. See the Estimand check tab."))
          else NULL)
    )
  })

  ## Several studies encode the output scale in the method name -- "bonett (r)"
  ## and "bonett (z)" are the same formula reported on the correlation and on the
  ## Fisher's z scale. Plotting both against one vertical axis compares numbers
  ## that are not on the same scale, and the z values, being larger, flatten the
  ## r values. Where such suffixes exist, offer them as a scale switch and show
  ## one scale at a time.
  method_scales <- reactive({
    m <- unique(raw()$method)
    sfx <- ifelse(grepl("\\s\\([^()]+\\)$", m), sub(".*\\s\\(([^()]+)\\)$", "\\1", m), NA)
    if (length(unique(stats::na.omit(sfx))) < 2) return(NULL)
    stats::setNames(sfx, m)
  })

  output$scale_ui <- renderUI({
    s <- method_scales()
    if (is.null(s)) return(NULL)
    lv <- sort(unique(stats::na.omit(s)))
    tagList(
      div(class = "grp", "Reported scale"),
      radioButtons("scale", NULL, choices = lv, selected = lv[1], inline = TRUE),
      div(class = "hint", "These are the same formulas on different scales, so only one is shown at a time.")
    )
  })

  output$methods_ui <- renderUI({
    s <- method_scales()
    m <- sort(unique(raw()$method))
    if (!is.null(s) && !is.null(input$scale)) m <- sort(names(s)[!is.na(s) & s == input$scale])
    checkboxGroupInput("methods", NULL, choices = m, selected = m)
  })

  output$xvar_ui <- renderUI({
    cc <- cond_cols(); req(length(cc) > 0)
    nlev <- vapply(cc, function(k) length(unique(raw()[[k]])), integer(1))
    selectInput("xvar", "Horizontal axis", choices = cc, selected = cc[which.max(nlev)])
  })

  output$facet_ui <- renderUI({
    cc <- cond_cols(); req(length(cc) > 0)
    tagList(
      selectInput("facet_row", "Split panels by", choices = c("(none)", cc),
                  selected = if (length(cc) >= 2) cc[2] else "(none)"),
      selectInput("facet_col", "and by", choices = c("(none)", cc), selected = "(none)")
    )
  })

  output$filters_ui <- renderUI({
    cc <- cond_cols()
    rest <- setdiff(cc, c(input$xvar, input$facet_row, input$facet_col))
    if (!length(rest)) return(NULL)
    tagList(
      div(class = "grp", "Held constant"),
      lapply(rest, function(k) {
        v <- sort(unique(raw()[[k]]))
        selectInput(paste0("f_", k), k, choices = v, selected = v[1])
      })
    )
  })

  dat <- reactive({
    d <- raw(); req(input$target, input$methods)
    d <- d[d$target == input$target & d$method %in% input$methods, , drop = FALSE]
    for (k in setdiff(cond_cols(), c(input$xvar, input$facet_row, input$facet_col))) {
      v <- input[[paste0("f_", k)]]
      if (!is.null(v)) d <- d[as.character(d[[k]]) == as.character(v), , drop = FALSE]
    }
    d
  })

  ## ---- reading strip ---------------------------------------------------------
  output$strip <- renderUI({
    d <- dat(); req(nrow(d) > 0)
    m <- METRICS[[input$metric]]
    lab <- if (input$study %in% names(STUDY_LABELS)) STUDY_LABELS[[input$study]] else input$study
    arg <- STUDY_ARG[[input$study]]
    tagList(
      div(class = "strip",
          h2(lab),
          if (!is.null(arg)) span(class = "arg", arg),
          span(class = "facts",
               sprintf("%d conditions · %s reps · %d methods",
                       nrow(d) / max(1, length(unique(d$method))),
                       format(as.integer(input$nrep), big.mark = ","),
                       length(unique(d$method))))),
      div(class = "subline",
          span(class = "l",
               sprintf("%s against the %s target", m$label, input$target),
               if (!is.na(m$ref)) sprintf(" · dashed line at %s is no error", m$ref) else ""),
          if (!is.na(m$mcse))
            span(class = "band-key", tags$i(),
                 tags$span(title = paste("Half-width is twice the median Monte Carlo SE",
                                         "of the cells in that panel. Differences inside",
                                         "the band are smaller than the simulation's own",
                                         "noise and are not results."),
                           "resolution limit, per panel"))
          else NULL)
    )
  })

  ## ---- main plot -------------------------------------------------------------
  output$main_plot <- renderPlot({
    d <- dat()
    validate(need(nrow(d) > 0, "No rows match these filters."))
    m <- METRICS[[input$metric]]
    req(input$xvar)

    d$.x <- d[[input$xvar]]
    d$.y <- suppressWarnings(as.numeric(d[[input$metric]]))
    validate(need(any(is.finite(d$.y)),
                  paste0(m$label, " was not recorded for this study.")))

    fr <- if (identical(input$facet_row, "(none)")) NULL else input$facet_row
    fc <- if (identical(input$facet_col, "(none)")) NULL else input$facet_col
    keys <- c(fr, fc)

    p <- ggplot(d, aes(.x, .y, colour = method, group = method))

    ## SIGNATURE: the Monte Carlo resolution band. Morris, White & Crowther (2019)
    ## report every measure with its MCSE precisely so a difference smaller than
    ## the simulation's own noise is not mistaken for a finding; this draws that
    ## threshold instead of leaving the reader to compute it. Half-width is 2x the
    ## median MCSE, the k = 2 convention of mc_resolvable() in R/03_performance.R.
    ##
    ## Computed PER PANEL, not once for the whole plot: MCSE shrinks with the
    ## sample size and the replication count, so a single global band would be far
    ## too narrow in the small-n panels and too wide in the large-n ones -- i.e. it
    ## would misstate resolvability in exactly the panels where the question
    ## matters most. Per panel, the band visibly narrows as precision improves.
    if (!is.na(m$mcse) && m$mcse %in% names(d) && !is.na(m$ref)) {
      e <- suppressWarnings(as.numeric(d[[m$mcse]]))
      grp <- if (length(keys)) interaction(d[keys], drop = TRUE) else factor(rep(1, nrow(d)))
      band <- do.call(rbind, lapply(split(seq_len(nrow(d)), grp), function(i) {
        h <- 2 * stats::median(e[i][is.finite(e[i])], na.rm = TRUE)
        if (!is.finite(h) || h <= 0) return(NULL)
        cbind(d[i[1], keys, drop = FALSE],
              data.frame(ymin = m$ref - h, ymax = m$ref + h), row.names = NULL)
      }))
      if (!is.null(band) && nrow(band))
        p <- p + geom_rect(data = band, inherit.aes = FALSE,
                           aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
                           fill = MC$ink2, alpha = 0.07)
    }
    if (!is.na(m$ref))
      p <- p + geom_hline(yintercept = m$ref, linetype = "22",
                          colour = MC$ink2, linewidth = 0.45)

    if (!is.na(m$mcse) && m$mcse %in% names(d)) {
      e <- suppressWarnings(as.numeric(d[[m$mcse]]))
      p <- p + geom_linerange(aes(ymin = .y - e, ymax = .y + e),
                              alpha = 0.5, linewidth = 0.45)
    }

    p <- p +
      geom_line(alpha = 0.85, linewidth = 0.55) +
      geom_point(size = 2, stroke = 0) +
      scale_colour_manual(values = rep(MC_METHODS, length.out = 60)) +
      scale_y_continuous(labels = function(x) format(x, trim = TRUE)) +
      labs(x = input$xvar, y = m$label) +
      theme_mc()

    if (length(keys)) {
      ## facet_wrap keeps the strips on TOP and horizontal. facet_grid puts row
      ## strips on the right rotated 90 degrees, which is unreadable at these
      ## label lengths.
      p <- p + facet_wrap(keys, labeller = label_both,
                          ncol = if (length(keys) > 1) NULL else 2,
                          scales = if (isTRUE(input$free_y)) "free_y" else "fixed")
    }
    p
  }, res = 108)

  ## ---- ranking ---------------------------------------------------------------
  ## A scorecard, not a single ranking. No one measure is sufficient: a
  ## conversion can be almost unbiased in the point estimate and badly
  ## miscalibrated in its standard error, and the two orderings genuinely
  ## disagree in these data. In study 09a the best method on bias
  ## (2x2_tetrachoric) is the WORST on coverage, and across the shipped studies
  ## the rank correlation between |bias| and |coverage - 0.95| is weak and
  ## sometimes negative. So all four measures are shown side by side; the
  ## selected measure only decides the sort order.
  output$rank_lede <- renderUI({
    m <- METRICS[[input$metric]]
    tagList(
      div(class = "lede",
          sprintf("Every measure for every method, averaged over the conditions currently shown, sorted by %s (%s). ",
                  m$label, m$better),
          HTML(paste0(
            "Read them together: <b>bias</b> is the point estimate only, ",
            "<b>RMSE</b> adds the estimator's variability, <b>SE ratio</b> asks whether the ",
            "reported standard error matches the real one, and <b>coverage</b> is the only ",
            "measure that degrades under <i>either</i> a wrong estimate or a wrong standard ",
            "error — which is why it is the default. The sort column is the mean ",
            "deviation <i>per condition</i>, which is not the deviation of the mean beside ",
            "it: a method covering 0.90 half the time and 1.00 the rest averages to a ",
            "perfect 0.95 while being miscalibrated throughout. <code>worst</code> is the ",
            "least favourable single condition, which is what a review should plan for."))),
      if (grepl("sample", input$target))
        div(class = "hint", style = "max-width:78ch;",
            HTML(paste0("Coverage and the SE ratio are <b>not interpretable</b> against the ",
                        "<b>", input$target, "</b> target: a nominal 95% interval is not built to ",
                        "cover a quantity that is itself random. Switch to the population ",
                        "parameter to read them.")))
      else if (grepl("^own", input$target))
        div(class = "hint", style = "max-width:78ch;",
            HTML(paste0("Every method here is scored against a <b>different</b> quantity, so ",
                        "this table says how accurately each computes its own value — not ",
                        "which method to use. Switch to the population parameter to compare ",
                        "them, and see the Estimand check tab for the difference.")))
      else NULL
    )
  })

  rank_df <- reactive({
    d <- dat(); m <- METRICS[[input$metric]]
    validate(need(nrow(d) > 0, "No rows match these filters."))
    sortv <- suppressWarnings(as.numeric(d[[input$metric]]))
    sortdev <- if (is.na(m$ref)) abs(sortv) else abs(sortv - m$ref)
    num <- function(col) suppressWarnings(as.numeric(d[[col]]))
    s <- split(seq_along(sortv), d$method)
    out <- do.call(rbind, lapply(names(s), function(k) {
      i <- s[[k]]
      cov_i <- num("coverage")[i]
      data.frame(
        method       = k,
        conditions   = length(i),
        ## The sort key, shown because it is NOT the same as a deviation computed
        ## from the mean column below: a method covering 0.90 in half the
        ## conditions and 1.00 in the rest averages to a perfect 0.95 while being
        ## miscalibrated everywhere. Averaging the per-condition deviation catches
        ## that; deviating the average does not.
        sortkey      = round(mean(sortdev[i], na.rm = TRUE), 4),
        bias         = round(mean(num("bias")[i], na.rm = TRUE), 4),
        `abs bias`   = round(mean(abs(num("bias")[i]), na.rm = TRUE), 4),
        rmse         = round(mean(num("rmse")[i], na.rm = TRUE), 4),
        `SE ratio`   = round(mean(num("se_ratio")[i], na.rm = TRUE), 3),
        coverage     = round(mean(cov_i, na.rm = TRUE), 3),
        `worst coverage` = round(suppressWarnings(min(cov_i, na.rm = TRUE)), 3),
        `max non-estimable` = round(suppressWarnings(max(num("nonest_rate")[i], na.rm = TRUE)), 3),
        .sort        = mean(sortdev[i], na.rm = TRUE),
        check.names  = FALSE, stringsAsFactors = FALSE)
    }))
    out <- out[order(out$.sort), ]
    out$.sort <- NULL
    ## numeric columns only: is.finite() is FALSE for every element of a
    ## character vector, so applying this to `method` blanks the whole column.
    isnum <- vapply(out, is.numeric, logical(1))
    out[isnum] <- lapply(out[isnum], function(x) { x[!is.finite(x)] <- NA; x })
    names(out)[names(out) == "sortkey"] <-
      if (is.na(m$ref)) paste0("mean |", input$metric, "|")
      else paste0("mean |", input$metric, " - ", m$ref, "|")
    rownames(out) <- NULL
    out
  })

  output$rank_tbl <- renderDT({
    r <- rank_df()
    dt <- datatable(r, rownames = FALSE,
                    options = list(dom = "t", pageLength = 30, ordering = TRUE))
    ## bar on the absolute bias, and a colour cue on coverage: red as it falls
    ## away from nominal in either direction.
    dt <- formatStyle(dt, "abs bias",
                      background = styleColorBar(range(c(0, r[["abs bias"]]), na.rm = TRUE),
                                                 "rgba(214,98,104,.18)"),
                      backgroundSize = "98% 62%", backgroundRepeat = "no-repeat",
                      backgroundPosition = "center")
    dt <- formatStyle(dt, "coverage",
                      color = styleInterval(c(0.90, 0.93, 0.97),
                                            c("#9c2b2b", "#8a6712", "#211f38", "#8a6712")),
                      fontWeight = "bold")
    dt
  })

  ## ---- estimand check --------------------------------------------------------
  output$estimand_lede <- renderUI({
    tagList(
      div(class = "lede",
          HTML(paste0(
            "Two questions have to be kept apart. <b>Is the method computing its own ",
            "quantity correctly?</b> — that is the dark bar. <b>Is that quantity the one ",
            "you asked for?</b> — the gap between the two bars. A long bar with a short gap ",
            "is a method with a real accuracy problem. A short bar with a long gap is a ",
            "method working perfectly on a different estimand, and reporting that gap as ",
            "bias is a mistake."))),
      div(class = "lede", style = "margin-top:-.5rem;",
          "Bars are the mean absolute bias over the conditions currently shown."))
  })

  output$estimand_plot <- renderPlot({
    d0 <- raw()
    validate(need("own" %in% d0$target,
                  "This study did not record an 'own' target, so there is no estimand gap to show."))
    pop <- grep("population", unique(d0$target), value = TRUE)
    validate(need(length(pop) > 0, "This study records no population-scale target."))
    req(input$methods)
    a <- d0[d0$target == "own" & d0$method %in% input$methods, , drop = FALSE]
    b <- d0[d0$target == pop[1] & d0$method %in% input$methods, , drop = FALSE]
    key <- function(x) do.call(paste, c(x[c(cond_cols(), "method")], sep = "\r"))
    i <- match(key(a), key(b)); ok <- !is.na(i)
    validate(need(any(ok), "No matching own/population cells."))
    df <- data.frame(method = a$method[ok], bias_own = a$bias[ok], bias_pop = b$bias[i[ok]])

    ## One row per method: how far it misses its OWN estimand (computational
    ## error) and how far it misses the POPULATION parameter (the total). The
    ## distance between the two is the estimand gap. A dumbbell reads far more
    ## directly than the bias-vs-bias scatter this replaced, where the reader had
    ## to infer the gap from a point's distance off the identity line.
    s <- do.call(rbind, lapply(split(df, df$method), function(x) data.frame(
      method = x$method[1],
      own    = mean(abs(x$bias_own), na.rm = TRUE),
      pop    = mean(abs(x$bias_pop), na.rm = TRUE),
      stringsAsFactors = FALSE)))
    s <- s[is.finite(s$own) & is.finite(s$pop), , drop = FALSE]
    validate(need(nrow(s) > 0, "No comparable cells."))
    s$gap <- s$pop - s$own
    s <- s[order(s$pop), ]
    s$method <- factor(s$method, levels = s$method)

    lab <- function(v) formatC(v, format = "f", digits = 3)
    ggplot(s) +
      geom_segment(aes(y = method, yend = method, x = own, xend = pop),
                   colour = MC$rule, linewidth = 2.6, lineend = "round") +
      geom_point(aes(y = method, x = pop), colour = MC$rose, size = 3.6) +
      geom_point(aes(y = method, x = own), colour = MC$ink, size = 3.6) +
      ## own above the dot, population below: when a method has no estimand gap
      ## the two dots coincide, and labels on the same side would overlap into a
      ## dark dot carrying a rose number.
      geom_text(aes(y = method, x = own, label = lab(own)),
                colour = MC$ink, size = 3.1, vjust = -1.5, family = "Cascadia Code") +
      geom_text(aes(y = method, x = pop, label = lab(pop)),
                colour = MC$rose_d, size = 3.1, vjust = 2.3, family = "Cascadia Code") +
      ## start at zero so the dot positions and the gap length are both read
      ## against a true origin
      scale_x_continuous(limits = c(0, NA),
                         expand = expansion(mult = c(0.02, 0.12))) +
      labs(x = "Mean absolute bias", y = NULL,
           title = NULL,
           subtitle = "dark = against the method's own estimand   |   rose = against the population parameter") +
      theme_mc() +
      theme(legend.position = "none",
            panel.grid.major.y = element_blank(),
            axis.text.y = element_text(family = "Work Sans", size = 11,
                                       colour = MC$ink, face = "bold"),
            plot.subtitle = element_text(size = 10.5, colour = MC$soft,
                                         margin = margin(0, 0, 12, 0)))
  }, res = 108)

  ## ---- data ------------------------------------------------------------------
  output$data_lede <- renderUI({
    ## emp_se and rmse sit next to each other here, and on a random target they do
    ## not compose the way readers expect. Say so at the point of confusion.
    random_target <- !is.null(input$target) && !grepl("population", input$target)
    tagList(
      div(class = "lede",
          "The rows behind the current view, with every ADEMP measure and its Monte Carlo standard error."),
      if (random_target)
        div(class = "hint", style = "margin:-.5rem 0 .9rem;max-width:78ch;",
            HTML(paste0(
              "Note on <code>rmse</code> and <code>emp_se</code>: the <b>", input$target,
              "</b> target is recomputed on each replication's own sample, so it is not a ",
              "fixed constant and <code>RMSE&sup2; = bias&sup2; + EmpSE&sup2;</code> does not apply. ",
              "The estimate and the target are strongly correlated, so <code>rmse</code> can sit ",
              "below <code>emp_se</code> with nothing wrong. Compare <code>rmse</code> only with ",
              "other <code>rmse</code> values on the same target.")))
      else NULL
    )
  })

  output$data_tbl <- renderDT({
    d <- dat()
    num <- vapply(d, is.numeric, logical(1))
    d[num] <- lapply(d[num], round, 4)
    datatable(d, rownames = FALSE, filter = "top",
              options = list(pageLength = 20, scrollX = TRUE))
  })

  output$dl <- downloadHandler(
    filename = function() paste0(input$study, "_nrep", input$nrep, "_", input$target, ".csv"),
    content = function(f) utils::write.csv(dat(), f, row.names = FALSE)
  )
}

shinyApp(ui, server)
