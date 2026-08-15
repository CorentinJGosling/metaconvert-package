library(tidyverse)
library(shinycssloaders)
library(shinyWidgets)
library(shinyjs)
library(shinyBS)
ui <- fluidPage(
  tags$style("@import url(https://use.fontawesome.com/releases/v5.14.1/css/all.css);"),
  tags$style(
    HTML("
    # TABPANEL ======================================
    .tabbable > .nav > li > a:hover {
       background-color: #fff;
       border-bottom: 5px solid #FCC8CB;
       border-left: 1.5px solid #fff;
       border-top: 1.5px solid #fff;
       border-right: 1.5px solid #fff;
      }

   .tabbable > .nav > li[class=active] > a {
       color: #000;
       # background-color: #FCC8CB !important;
       border-left: 1.5px solid #fff;
       border-top: 1.5px solid #fff;
       border-right: 1.5px solid #fff;
       border-bottom: 6px solid #d66268;
   }

   .tabbable > .nav > li[class=active] > a:hover {
       border-bottom: 6px solid #d66268 !important;
   }
    .tabbable .nav-tabs li a {
        font-size: 1.7rem;
       color: darkgrey;
       border-left: 1.5px solid #fff;
       border-top: 1.5px solid #fff;
       border-right: 1.5px solid #fff;
       # border-bottom: 1px solid darkgrey;
    }
   .tabbable .nav-tabs li a:hover {
       color: #000;
       border-bottom: 6px solid #D99097;
       background-color: #fff !important;
    }
  .split {
     margin-top: 1rem;
     display: flex !important;
     gap: 2rem;
     # justify-content: space-between;
     max-height: 95vh !important;
     max-width: 99vw !important;

 }
 .left {
   background-color: #fff;
   min-height: 95vh !important;
   # width: 25vw !important;
   # overflow: hidden;
   padding: 1rem;
   padding-right: 3rem;
    box-shadow: rgba(14, 30, 37, 0.12) 0px 2px 4px 0px,
                rgba(14, 30, 37, 0.32) 0px 2px 16px 0px;  }
  .right {
     # border: 2px solid red;
     background-color: #fff;
     overflow-x: hidden;
     width: 100% !important;
     max-height: 95vh !important;
     # width: 75vw !important;

  }
  .tab_dashboard{
      padding: 1rem;
      # border: 1.5px solid #000;
      min-height: 80vh !important;
      # background-color: #FCC8CB;
  }
  .tab_fix{
    display: grid;
    gap: 2rem;
  }
  .title {
      font-size: 2.5rem !important;
      font-weight: bold  !important:
      font: 'Fira Sans', sans-serif;
  }
  .paragraph{
    display: flex;
  }
  .title_small,
  .title_bold,
  .title_ref {
      font-size: 1.5rem !important;
      font: 'Fira Sans', sans-serif;
  }
  .title_bold {
    paddingn-top: 1rem;
    font-weight: bold;
  }
  .title_ref{
  font-style: italic;
  }

  .boxes_group{
    display: grid;
    grid-template-columns: 1fr 1fr 1fr;
    gap: 3rem;
  }
  .box {
    display: grid;
    grid-auto-flow: column;
    justify-content: space-between;
    # grid-template-columns: 30% 70%;
    background-color: #FCC8CB;
    align-items: center;
    box-shadow: rgba(60, 64, 67, 0.3) 0px 1px 2px 0px, rgba(60, 64, 67, 0.15) 0px 1px 3px 1px;
  }
  .left_card{
    color: #fff;
    font-size: 6rem;
    background-color: #D99097;
    height: 100% !important;
    width: 100% !important;
    padding: 2rem;

  }
  .right_card{
    height: 100% !important;
    width: 100% !important;
    padding: 2rem;
    color: #474747;
  }

  .card_number {
    font-size: 6rem;
    text-align: center;
    font-weight: bold;
  }
  .card_text {
    padding-top: -6rem !important;
    font-weight: bold;
    font-size: 1.3rem;
  }

    #wrap_cons{
      animation: moving 3s forwards;
      animation-iteration-count: infinite;
    }
    @keyframes moving {
      0%   {margin-top: 0rem;}
      80%   {margin-top: 0rem;}
      90%  {margin-top: -1rem;}
      100% {margin-top: 0rem;}
    }

   .icon_cons{
     text-align: center;
     color: red;
     font-weight: bold;
     font-size: 2.5rem;
   }
   .text_cons{
     text-align: center;
     color: red;
     font-weight: bold;
     font-size: 1.5rem;
     max-width: 15ch;
     margin-left: auto;
     margin-right: auto;
   }
   .tab_smd_cor{
    margin-bottom: 3rem;
   }
   .warning{
    margin-top: 2rem;
    font-size: 1.5rem;
    color: red;
    padding: 2rem;
    border: 2px solid red;
    width: auto;
   }
   .plot_control_1{
      display: grid !important;
      grid-template-column: 1fr 1fr;
      display: inline-block;
      vertical-align:top

   }

  @media (max-width: 800px) {
   .left {
   min-height: auto;
   }

   .split {
     overflow-y: scroll !important;
     max-height: 200vh !important;
     flex-direction: column !important;
   }
  .right {
     width: 100% !important;
  }
  .left{
      width: 100% !important;
      display: flex !important;
      flex-direction: column !important;
      align-items: center !important;
      text-align : center !important;
      overflow-y: scroll !important;
  }

")),
  shinyjs::useShinyjs(),
  inlineCSS(list(.border_table = "border: 1px solid #000 !important;",
                 .disappear = "display: none !important;")),
  div(class = "split",
      # Selecting panel ========================
      div(class = "left",
          pickerInput("dataset", "1. Select a conversion:",
            width="auto",
            choices = list(
              'You want to convert SMD to' = c(
                "Correlation (R)" = "smd_r",
                "Correlation (Z)" = "smd_z"),
              'You want to convert MD to' = c(
                "Standardized mean difference (SMD)" = "md_smd"),
              'You want to convert adjusted means to' = c(
                "Standardized mean difference (SMD)" = "means_adj_smd"),
              'You want to convert t-value from ANCOVA to' = c(
                "Standardized mean difference (SMD)" = "f_ancova_smd"),
              'You want to convert OR to' = c(
                "Risk ratio (RR)" = "or_rr",
                "Odds ratio + standard error" = "se_or",
                "Correlation (R)" = "or_r",
                "Correlation (Z)" = "or_z"),
              'You want to convert RR to' = c(
                "Odds ratio (OR)" = "rr_or"),
              'You want to convert contingency (2x2) table to' = c(
                "Correlation (R)" = "contingency_r",
                "Correlation (Z)" = "contingency_z")),
            options = list(title = 'Nothing selected')
         ),
          div(id = "wrap_cons",
              div(class="icon_cons",icon("arrow-alt-circle-up")),
              div(class="text_cons", "Start by selecting a conversion")),
          div(class="wrap_dis",
              pickerInput(inputId = "grouper",
                      label = "2. Select splitting variables",
                      choices = c("Converting method" = "method"),
                      selected = "method",
                      multiple = TRUE)),

          div(class="wrap_dis", uiOutput("filtering_method"))

          ),
      div(class = "right",
          tabsetPanel(
            # First tab ========================
            tabPanel("Dashboard",
                     div(class="tab_dashboard",
                         div(class="tab_fix",
                           div(class = "title", "General information on the MC simulations."),
                           div(class = "title_small", "We conducted several Monte Carlo (MC) simulations exploring the approximations caused by some of the formulas available in the es.utils tools.
                           In this section, we give you the opportunity to interactively visualize the results of these simulations. Deviations led by the converting formulas were observed on 3 indicators (effect size value, standard error, and p-value).
                               "),
                           div(class = "boxes_group",
                             div(class="box box1",
                                 div(class="left_card",
                                     icon("square-root-alt")),
                                 div(class="right_card",
                                     div(class="card_number", "10"),
                                     div(class="card_text", "simulations"))),
                             div(class="box box2",
                                 div(class="left_card",
                                     icon("superscript")),
                                 div(class="right_card",
                                     div(class="card_number", "31"),
                                     div(class="card_text", "formulas"))),
                             div(class="box box3",
                                 div(class="left_card",
                                     icon("exclamation-circle")),
                                 div(class="right_card",
                                     div(class="card_number", "3"),
                                     div(class="card_text", "deviation indicators")))
                             ),
                           div(class = "title", "Overview of the methods."),
                           div(class = "title_small", "For each of these simulations, we started by generating
                           a large amount of fictitious datasets containing always 1) an effect size value and its variance,
                           plus 2) input data that can be converted to estimate the generated effect size.",
                           "We controlled various parameters
                           parameters in the datasets (such as the number of participants, the strength of an effect size and its variance,
                           the ratio of cases among participants, etc.), and we generated 5000 datasets for each parameter combination.
                           Hundreds of thousands of dataset were thus generated for each simulation."),

                           div(class = "title_small", "Then, we applied some of the formulas contained in es.utils to convert
                           back the input data into the generated effect size. The comparison of this effect size estimated from the
                           input data with the effect size directly generated in the simulation allowed to quantify the approximation
                           led by these formulas."),
                           tags$figure(
                             class = "centerFigure",
                             tags$img(
                               src = "fig_back.jpg",
                               width = 600,
                               alt = ""
                             ),
                             # tags$figcaption("Image of Astragalus by Yaan, 2007")
                           )
                           ),
                       # textOutput("checkinput"),
                       # div(class="tab_smd_cor",
                       #     div(class = "title", "Conversion of SMD to COR"),
                       #     div(class = "title_small", "To convert SMD to COR (R or Z), we relied on two approaches."),
                       #     div(class = "paragraph",
                       #         div(class = "title_bold", "1. Viechtbauer."),
                       #         div(class = "title_small", " In this first approach, we converted the SMD to a 'biserial' correlation.")),
                       #     div(class = "title_small", "Jacobs, P., & Viechtbauer, W. (2017). Estimation of the biserial correlation and its sampling variance for use in meta-analysis. Research synthesis methods, 8(2), 161–180. https://doi.org/10.1002/jrsm.1218"),
                       # ),
                       div(class="title title_table", "Results of the MC simulations"),
                       uiOutput('summary'))
            ),
            # Second tab ========================
            tabPanel("Plot (ES & SE)",
                     div(id = "plot_control_1",
                           uiOutput('plotOptions_1'),
                           uiOutput('plotSize_1')),
                     div(id = "plot_output_1",
                         uiOutput('plot1'))
                     ),

            # Third tab ========================
            tabPanel("Plot (p-value)",
                     div(id = "plot_box_2",
                         div(id = "plot_control_2",
                             uiOutput('plotOptions_2'),
                             uiOutput('plotSize_2')),
                         div(id = "plot_output_2",
                             uiOutput('plot2')))
            )

          )
      )

  )
)

server <- function(input, output, session) {

  # Reactive panel ----------
  observe({
    disable_grouper <- input$dataset == ""

    toggleClass(selector = ".wrap_dis", class="disappear",
                condition = disable_grouper)
    toggleClass(selector = "#wrap_cons", class="disappear",
                condition = !disable_grouper)

    toggleClass(selector = ".tab_fix", class="disappear",
                condition = !disable_grouper)
    toggleClass(selector = ".title_table", class="disappear",
                condition = disable_grouper)

    input_smd_cor <- input$dataset %in% c("smd_r", "smd_z")
    toggleClass(selector = ".tab_smd_cor", class="disappear",
                condition = !input_smd_cor)



  })
  # Dataset ----------
  mydata <- reactive({

    if (input$dataset == "smd_r") {
      res = read.delim("D:/simulations/data/SMD_to_COR_AGG.txt") %>%
        filter(grepl("(r)", method, fixed = TRUE))
      res$n = factor(res$n, ordered = TRUE,
                     levels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
      data.frame(row = 1:nrow(res), res)

    } else if (input$dataset == "smd_z") {
      res = read.delim("D:/simulations/data/SMD_to_COR_AGG.txt") %>%
        filter(grepl("(z)", method, fixed = TRUE))
      res$r = gsub("r", "z", res$r)
      res$n = factor(res$n, ordered = TRUE,
                     levels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
      data.frame(row = 1:nrow(res), res)

    } else if (input$dataset == "contingency_r") {
      res = read.delim("./data/2x2_to_COR_AGG.txt") %>%
        filter(grepl("(r)", method, fixed = TRUE))
      res$n = factor(res$n, ordered = TRUE,
                     levels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
      data.frame(row = 1:nrow(res), res)

    } else if (input$dataset == "contingency_z") {
      res = read.delim("./data/2x2_to_COR_AGG.txt") %>%
        filter(grepl("(z)", method, fixed = TRUE))
      res$r = gsub("r", "z", res$r)
      res$n = factor(res$n, ordered = TRUE,
                     levels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
      data.frame(row = 1:nrow(res), res)
    } else if (input$dataset == "or_r") {
      res = read.delim("./data/OR_to_COR_AGG.txt") %>%
        filter(grepl("(r)", method, fixed = TRUE))
      res$n = factor(res$n, ordered = TRUE,
                     levels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
      data.frame(row = 1:nrow(res), res)

    } else if (input$dataset == "or_z") {
      res = read.delim("./data/OR_to_COR_AGG.txt") %>%
        filter(grepl("(z)", method, fixed = TRUE))
      res$r = gsub("r", "z", res$r)
      res$n = factor(res$n, ordered = TRUE,
                     levels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
      data.frame(row = 1:nrow(res), res)
    } else if (input$dataset == "or_rr") {
      res = read.delim("D:/simulations/data/OR_to_RR_AGG.txt")
      res$n = factor(res$n, ordered = TRUE,
                     levels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
      data.frame(row = 1:nrow(res), res)
    } else if (input$dataset == "means_adj_smd") {
      res = read.delim("./data/MEANS_ADJ_to_SMD_AGG.txt")
      res$n = factor(res$n, ordered = TRUE,
                     levels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
      data.frame(row = 1:nrow(res), res)
    } else if (input$dataset == "md_smd") {
      res = read.delim("./data/MD_to_SMD_AGG.txt")
      res$n = factor(res$n, ordered = TRUE,
                     levels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
      data.frame(row = 1:nrow(res), res)
    } else if (input$dataset == "se_or") {
      res = read.delim("./data/OR_SE_AGG.txt")
      res$n = factor(res$n, ordered = TRUE,
                     levels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
      data.frame(row = 1:nrow(res), res)
    } else if (input$dataset == "f_ancova_smd") {
      res = read.delim("./data/F_ANCOVA_to_SMD_AGG.txt")
      res$n = factor(res$n, ordered = TRUE,
                     levels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
      data.frame(row = 1:nrow(res), res)
    } else if (input$dataset == "rr_or") {
      res = read.delim("./data/RR_to_OR_AGG.txt")
      res$n = factor(res$n, ordered = TRUE,
                     levels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
      data.frame(row = 1:nrow(res), res)
    }
  })



  # Choices input ------------
  choices_group <- reactive({
    # req(input$dataset)
    if (input$dataset == "smd_r") {
      c('Population correlation (R)' = "r",
        'Sample size (n)' = "n",
        'Percentage of cases (p)' = "p",
        'Converting method' = "method")
    } else if (input$dataset == "smd_z") {
      c('Population correlation (Z)' = "r",
        'Sample size (n)' = "n",
        'Percentage of cases (p)' = "p",
        'Converting method' = "method")
    } else if (input$dataset %in% c("or_r", "contingency_r")) {
      c('Population correlation (R)' = "r",
        'Sample size (n)' = "n",
        # 'Percentage of cases (p)' = "p",
        'Converting method' = "method")
    } else if (input$dataset %in% c("or_z", "contingency_z")) {
      c('Population correlation (Z)' = "r",
        'Sample size (n)' = "n",
        # 'Percentage of cases (p)' = "p",
        'Converting method' = "method")
    } else if (input$dataset == "or_rr") {
      c('Population RR' = "rr",
        'Sample size (n)' = "n",
        'TRUE baseline risk (br)' = "br",
        'GUESSED baseline risk (br)' = "br_guess",
        'Converting method' = "method")
    } else if (input$dataset %in% c("se_or")) {
      c('Population OR' = "or",
        'Sample size (n)' = "n",
        'Baseline risk (br)' = "br",
        'Converting method' = "method")
    } else if (input$dataset %in% c("rr_or")) {
      c('Population OR' = "or",
        'Sample size (n)' = "n",
        'TRUE baseline risk (br)' = "br",
        'GUESSED baseline risk (br)' = "br_guess",
        'Converting method' = "method")
    } else if (input$dataset %in% c("f_ancova_smd", "means_adj_smd")) {
      c('Population SMD' = "d",
        'Sample size (n)' = "n",
        'TRUE outcome-covariate cor (r)' = "r_cov",
        'GUESSED outcome-covariate cor (r)' = "r_cov_guess",
        'Converting method' = "method")
    } else if (input$dataset == "md_smd") {
      c('Population SMD' = "d",
        'Sample size (n)' = "n",
        'Converting method' = "method")
    }
  })
  # Selected choices =============
  rv <- reactiveValues()
  observeEvent(input$grouper_open, {
    if (!isTRUE(input$grouper_open)) {
      rv$grouper_delayed <- input$grouper
    }
  }, ignoreNULL = FALSE)

  observe({
    selected_group <- input$grouper
    Choices = choices_group()
    updatePickerInput(
      session = session,
      inputId = "grouper",
      choices = Choices,
      selected = selected_group
    )
  })

  observe({
    if (is.null(rv$grouper_delayed)) {
      updatePickerInput(
        session = session,
        inputId = "grouper",
        choices = choices_group(),
        selected = NULL
      )

    }
  })


  # Filtering -----------------------
  output$filtering_method <- renderUI({
    req(rv$grouper_delayed)
    if ("method" %in% rv$grouper_delayed) {
      selectizeInput(inputId = "filter_method", label = "3. Restrict to a specific converting method",
                     choices = unique(mydata()[, "method"]), multiple = TRUE,
                     options = list(plugins = list('remove_button'),
                                    placeholder = 'Please select an option below',
                                    onInitialize = I('function() { this.setValue(""); }')))
    }

  })


  # OUTPUT 1 - Table ------------------------------------
  summary_data <- reactive({
    if (!is.null(rv$grouper_delayed)) {
      if ("method" %in% rv$grouper_delayed & !is.null(input$filter_method)) {
        res = subset(mydata(), method %in% input$filter_method)
        res %>%
          dplyr::group_by(!!!rlang::syms(rv$grouper_delayed)) %>%
          dplyr::summarise(bias_es = round(mean(bias_es, na.rm = TRUE), 3),
                           bias_var = round(mean(bias_var, na.rm = TRUE), 3),
                           bias_ci = round(mean(bias_ci, na.rm = TRUE), 3))

      } else {
        mydata() %>%
          dplyr::group_by(!!!rlang::syms(rv$grouper_delayed)) %>%
          dplyr::summarise(bias_es = round(mean(bias_es, na.rm = TRUE), 3),
                           bias_var = round(mean(bias_var, na.rm = TRUE), 3),
                           bias_ci = round(mean(bias_ci, na.rm = TRUE), 3))
      }
    } else {
      mydata()
    }
  })

  output$summary <- renderUI({
    req(input$dataset)
    DT::renderDataTable({
      DT::datatable(summary_data(),
                    rownames = FALSE,
                    extensions = 'Buttons',
                    options = list(
                      dom = c('tB'),
                      scrollX = TRUE,
                      # scrollY = "600px",
                      # pageLength = 80,
                      buttons =
                        list("copy", list(
                          extend = "collection"
                          , buttons = c("csv", "excel")
                          , text = "Download"
                        ))
                    )
      )
    })
  })

  # OUTPUT 2 - Plot ES/SE ------------------------------------
  observeEvent(input$plotOptions_1, {

    output$plot_1 <- renderPlot({

      res_plot = summary_data()

      # title of x-axis
      if (input$dataset %in% c("smd_z", "contingency_z", "or_z")) {
        x_lab = "Deviation (Z value)"
        if ("r" %in% rv$grouper_delayed) {
          res_plot$es = res_plot$r
        } else {
          res_plot$es = ""
        }
      } else if (input$dataset %in% c("smd_r", "contingency_r", "or_r")) {
        x_lab = "Deviation (R value)"
        if ("r" %in% rv$grouper_delayed) {
          res_plot$es = res_plot$r
        } else {
          res_plot$es = ""
        }
      } else if (input$dataset %in% c("f_ancova_smd", "means_adj_smd", "md_smd")) {
        x_lab = "Deviation (SMD value)"
        if ("d" %in% rv$grouper_delayed) {
          res_plot$es = res_plot$d
        } else {
          res_plot$es = ""
        }
      } else if (input$dataset %in% c("or_rr")) {
        x_lab = "Deviation (RR value)"
        if ("rr" %in% rv$grouper_delayed) {
          res_plot$es = res_plot$rr
        } else {
          res_plot$es = ""
        }
      } else if (input$dataset %in% c("rr_or", "se_or")) {
        x_lab = "Deviation (OR value)"
        if ("or" %in% rv$grouper_delayed) {
          res_plot$es = res_plot$or
        } else {
          res_plot$es = ""
        }
      }

      if ("method" %in% rv$grouper_delayed) {
        Method = res_plot$method
      } else {
        Method = "Across methods"
      }

      if ("p" %in% rv$grouper_delayed) {
        if (Method[1] != "Across methods") {
          Method = paste0(Method, " [p=", res_plot$p, "]")
        } else {
          Method = paste0("Across methods: [p=", res_plot$p, "]")
        }
      }

      plot_es = ggplot(res_plot,
                       aes(x = bias_es, y = bias_se,
                           color = Method,
                           id = Method,
                           group = Method,
                           label = Method)) +
        geom_point(size = 5, alpha = 0.6) +
        ylab("Deviation (standard error)") +
        xlab(x_lab) +
        geom_hline(yintercept=0, linetype="dashed", color = "red") +
        geom_vline(xintercept=0, linetype="dashed", color = "red") +
        theme_bw() +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 12, face="bold"),
              legend.position="top") +
        ggrepel::geom_text_repel(max.overlaps = 200)

      # n ,
      # Filters: N
      if ("n" %in% rv$grouper_delayed) {
        if (res_plot$es[1] != "") {
          if ("br" %in% rv$grouper_delayed) {
            if ("br_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(es + n ~ br + br_guess)
            } else {
              plot_es + facet_grid(es + n ~ br)
            }
          } else if ("r_cov" %in% rv$grouper_delayed) {
            if ("r_cov_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(es + n ~ r_cov + r_cov_guess)
            } else {
              plot_es + facet_grid(es + n ~ r_cov)
            }
          } else {
            if ("br_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(es + n ~ br_guess)
            } else if ("r_cov_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(es + n ~ r_cov_guess)
            } else {
              plot_es + facet_grid(es ~ n)
            }
          }
        } else {
          if ("br" %in% rv$grouper_delayed) {
            if ("br_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(n ~ br + br_guess)
            } else {
              plot_es + facet_grid(n ~ br)
            }
          } else if ("r_cov" %in% rv$grouper_delayed) {
            if ("r_cov_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(n ~ r_cov + r_cov_guess)
            } else {
              plot_es + facet_grid(n ~ r_cov)
            }
          } else {
            if ("br_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(n ~ br_guess)
            } else if ("r_cov_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(n ~ r_cov_guess)
            } else {
              plot_es + facet_grid( ~ n)
            }
          }
        }
      } else {
          if (res_plot$es[1] != "") {
            if ("br" %in% rv$grouper_delayed) {
              if ("br_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid(es ~ br + br_guess)
              } else {
                plot_es + facet_grid(es ~ br)
              }
            } else if ("r_cov" %in% rv$grouper_delayed) {
              if ("r_cov_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid(es ~ r_cov + r_cov_guess)
              } else {
                plot_es + facet_grid(es ~ r_cov)
              }
            } else {
              if ("br_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid(es ~ br_guess)
              } else if ("r_cov_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid(es  ~ r_cov_guess)
              } else {
                plot_es + facet_grid( ~ es)
              }
            }
          } else {
            if ("br" %in% rv$grouper_delayed) {
              if ("br_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid(br ~ br_guess)
              } else {
                plot_es + facet_grid( ~ br)
              }
            } else if ("r_cov" %in% rv$grouper_delayed) {
              if ("r_cov_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid(r_cov ~ r_cov_guess)
              } else {
                plot_es + facet_grid( ~ r_cov)
              }
            } else {
              if ("br_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid( ~ br_guess)
              } else if ("r_cov_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid( ~ r_cov_guess)
              } else {
                plot_es
              }
            }
          }
      }

      # if ("r_cov" %in% rv$grouper_delayed) {
      #   if ("r_cov_guess" %in% rv$grouper_delayed) {
      #   }
      # }

    }, res = input$plotOptions_1)
  })
  ### options for plot 2 ------------------
  output$plotOptions_1 <- renderUI({
    if (input$dataset != "") {
      req(input$dataset)
      sliderInput('plotOptions_1', "Font size", min=20, max=200, value=75)
    } else {
      div(class="warning", "Before being able to visualize the results, you must select a specific conversion using the left menu.")
    }
  })
  output$plotSize_1 <- renderUI({
    req(input$dataset)
    sliderInput('plotSize_1', "Height", min=200, max=2000, value=600)
  })
  output$plot1 <- renderUI({
    req(input$dataset)
    plotOutput('plot_1', height = input$plotSize_1)  %>% withSpinner(color = "#333366", type = 8)
  })

  # OUTPUT 3 - Plot pval ------------------------------------
  observeEvent(input$plotOptions_2, {

    output$plot_2 <- renderPlot({

      res_plot = summary_data()

      if (input$dataset %in% c("smd_z", "contingency_z", "or_z")) {
        if ("r" %in% rv$grouper_delayed) {
          res_plot$es = res_plot$r
        } else {
          res_plot$es = ""
        }
      } else if (input$dataset %in% c("smd_r", "contingency_r", "or_r")) {
        if ("r" %in% rv$grouper_delayed) {
          res_plot$es = res_plot$r
        } else {
          res_plot$es = ""
        }
      } else if (input$dataset %in% c("f_ancova_smd", "means_adj_smd", "md_smd")) {
        if ("d" %in% rv$grouper_delayed) {
          res_plot$es = res_plot$d
        } else {
          res_plot$es = ""
        }
      } else if (input$dataset %in% c("or_rr")) {
        if ("rr" %in% rv$grouper_delayed) {
          res_plot$es = res_plot$rr
        } else {
          res_plot$es = ""
        }
      } else if (input$dataset %in% c("rr_or", "se_or")) {
        x_lab = "Deviation (OR value)"
        if ("or" %in% rv$grouper_delayed) {
          res_plot$es = res_plot$or
        } else {
          res_plot$es = ""
        }
      }

      if ("method" %in% rv$grouper_delayed) {
        Method = res_plot$method
      } else {
        Method = "Across methods"
      }

      if ("p" %in% rv$grouper_delayed) {
        if (Method[1] != "Across methods") {
          Method = paste0(Method, " [p=", res_plot$p, "]")
        } else {
          Method = paste0("Across methods: [p=", res_plot$p, "]")
        }
      }

      plot_es = res_plot %>%
        mutate(bias_pval = as.numeric(bias_pval)) %>%
        ggplot() +
        aes(x = Method, y = bias_pval, fill = Method) +
        geom_bar(stat="identity", alpha=.4) +
        ylab("Deviation (p-value)") +
        geom_hline(yintercept=0, linetype="dashed", color = "red") +
        geom_vline(xintercept=0, linetype="dashed", color = "red") +
        theme_bw() +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 12, face="bold"),
              legend.position="top") +
        coord_flip()

      if ("n" %in% rv$grouper_delayed) {
        if (res_plot$es[1] != "") {
          if ("br" %in% rv$grouper_delayed) {
            if ("br_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(es + n ~ br + br_guess)
            } else {
              plot_es + facet_grid(es + n ~ br)
            }
          } else if ("r_cov" %in% rv$grouper_delayed) {
            if ("r_cov_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(es + n ~ r_cov + r_cov_guess)
            } else {
              plot_es + facet_grid(es + n ~ r_cov)
            }
          } else {
            if ("br_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(es + n ~ br_guess)
            } else if ("r_cov_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(es + n ~ r_cov_guess)
            } else {
              plot_es + facet_grid(es ~ n)
            }
          }
        } else {
          if ("br" %in% rv$grouper_delayed) {
            if ("br_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(n ~ br + br_guess)
            } else {
              plot_es + facet_grid(n ~ br)
            }
          } else if ("r_cov" %in% rv$grouper_delayed) {
            if ("r_cov_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(n ~ r_cov + r_cov_guess)
            } else {
              plot_es + facet_grid(n ~ r_cov)
            }
          } else {
            if ("br_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(n ~ br_guess)
            } else if ("r_cov_guess" %in% rv$grouper_delayed) {
              plot_es + facet_grid(n ~ r_cov_guess)
            } else {
              plot_es + facet_grid( ~ n)
            }
          }
        }
      } else {
        {
          if (res_plot$es[1] != "") {
            if ("br" %in% rv$grouper_delayed) {
              if ("br_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid(es ~ br + br_guess)
              } else {
                plot_es + facet_grid(es ~ br)
              }
            } else if ("r_cov" %in% rv$grouper_delayed) {
              if ("r_cov_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid(es ~ r_cov + r_cov_guess)
              } else {
                plot_es + facet_grid(es ~ r_cov)
              }
            } else {
              if ("br_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid(es ~ br_guess)
              } else if ("r_cov_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid(es  ~ r_cov_guess)
              } else {
                plot_es + facet_grid( ~ es)
              }
            }
          } else {
            if ("br" %in% rv$grouper_delayed) {
              if ("br_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid(br ~ br_guess)
              } else {
                plot_es + facet_grid( ~ br)
              }
            } else if ("r_cov" %in% rv$grouper_delayed) {
              if ("r_cov_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid(r_cov ~ r_cov_guess)
              } else {
                plot_es + facet_grid( ~ r_cov)
              }
            } else {
              if ("br_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid( ~ br_guess)
              } else if ("r_cov_guess" %in% rv$grouper_delayed) {
                plot_es + facet_grid( ~ r_cov_guess)
              } else {
                plot_es
              }
            }
          }
        }
      }
    }, res = input$plotOptions_2)
  })

  ## options for plot 3 --------------
  output$plotOptions_2 <- renderUI({
    if (input$dataset != "") {
      req(input$dataset)
      sliderInput('plotOptions_2', "Font size", min=20, max=200, value=75)
    } else {
      div(class="warning", "Before being able to visualize the results, you must select a specific conversion using the left menu.")
    }
  })

  output$plotSize_2 <- renderUI({
    req(input$dataset)
    sliderInput('plotSize_2', "Height", min=200, max=2000, value=600)
  })

  output$plot2 <- renderUI({
    req(input$dataset)
    plotOutput('plot_2', height = input$plotSize_2)  %>% withSpinner(color = "#333366", type = 8)})

}

shinyApp(ui, server)
