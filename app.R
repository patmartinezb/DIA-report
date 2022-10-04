
# Load libraries ----------------------------------------------------------

# Bioconductor

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

library(BiocManager)


# Import libraries

lib <-  c("shiny", "readxl", "dplyr", "ggplot2", "tidyr", "ggsci", "stringr", "readr", "bslib",
          "purrr", "tools", "factoextra", "EnhancedVolcano", "limma", "officer", "ggplotify",
          "pheatmap", "imputeLCMD", "org.Hs.eg.db", "clusterProfiler", "venn", "ggpolypath")


# Checking missing packages from list

new.packages <- setdiff(lib, installed.packages())

if (length(new.packages) != 0) invisible(lapply(new.packages, BiocManager::install, update = FALSE))

# Load

invisible(lapply(lib, require, character.only = TRUE))

# The types of files that are going to be uploaded are bigger than the specified shiny limit

options(shiny.maxRequestSize = 30*1024^2)

# Define UI ---------------------------------------------------------------

ui <- fluidPage(
  titlePanel("DIA Tims-TOF data analysis report"),
  theme = bslib::bs_theme(bootswatch = "sandstone"),
  tabsetPanel(
    tabPanel("Report",
             br(),
             fileInput("file", "Upload data (.xlsx/.tsv):", buttonLabel = "Upload data", 
                       multiple = FALSE, accept = c(".tsv", ".xlsx")),
             hr(),
             fileInput("meta", "Upload metadata (.xlsx):", buttonLabel = "Upload metadata", 
                       multiple = FALSE, accept = ".xlsx"),
             hr(),
             fileInput("comp", "Upload comparison file (.xlsx):", buttonLabel = "Upload comparisons", 
                       multiple = FALSE, accept = ".xlsx"),
             hr(),
             uiOutput("cityControls"),
             hr(),
             sliderInput("slider", "Choose a threshold for NA filtering:", 0, 1, 0.7),
             hr(),
             paste0("Downloads:"),
             br(),
             div(style="display:inline-block", downloadButton("report", "Generate report")),
             div(style="display:inline-block", downloadButton("tables", "Generate tables")),
             div(style="display:inline-block", downloadButton("figs", "Generate figures")),
             hr()
             ),
  
  tabPanel("QC figures",
           fluidRow(
             column(6,
                    plotOutput(outputId = "hist")
                    ),
             
             column(6,
                    plotOutput(outputId = "misvals")
                    )
             )
           )
  )
)


# Define server -----------------------------------------------------------

  
  server <-  function(input, output, session) {
    
    # Raw data, metadata and comparisons
    
    data <- reactive({
      req(input$file)
      
      ext <- tools::file_ext(input$file$name)
      switch(ext,
             xlsx = read_excel(input$file$datapath) %>% readr::type_convert(),
             tsv = vroom::vroom(input$file$datapath, delim = "\t"),
             validate("Invalid file; Please upload a .tsv or .xlsx file")
      )
    })
    
    
    metadata <- reactive({
      
      req(input$meta)
      
      df <- read_excel(input$meta$datapath)
      
      return(df)
      
    })
    
    comp <- reactive({
      
      req(input$comp)
      
      df <- read_excel(input$comp$datapath)
      
      return(df)
      
    })
    
    
    # Processed data and metadata
    
    data2 <- reactive({
      
      df <- data()
      
      colnames(df) <- make.names(colnames(df))
      
      if (length(input$samples) == 0){
        
        df <- df
        
      } else {
        
        df <- df %>% 
          select(-all_of(input$samples))
        
      }
      
      df <- df %>%
        mutate(Protein.Ids = sub(";.*", "", Protein.Ids),
               Genes = sub(";.*", "", Genes),
               across(where(is.numeric), na_if, 0),
               across(where(is.numeric), na_if, Inf),
               across(where(is.numeric), na_if, -Inf),
               across(where(is.numeric), log2)) %>% 
        distinct(Protein.Ids, .keep_all = TRUE)
      
      return(df)
      
    })
    
    metadata2 <- reactive({
      
      df <- metadata()
      
      df <- as_tibble(apply(df, 2, make.names))
      
      if (length(input$samples) == 0){
        
        df <- df
        
      } else {
        
        df <- df %>% 
          filter(!BioReplicate %in% input$samples)
        
      }
      
      return(df)
      
    })
    
    # Optional removal of samples
    
    output$cityControls <- renderUI({
      
      req(input$meta)
      
      samples <- metadata() %>% 
        select(BioReplicate) %>% 
        pull() %>% 
        make.names(.)
      
      selectInput(inputId = "samples", choices = samples, 
                  label = "Do you want to remove one or multiple samples? If so, choose which one(s):", multiple = TRUE)
      
    })
    
    # Plot QC figures
    
    output$hist <- renderPlot({
      
      req(input$file)
      req(input$meta)
      
      plot_n_prots_sample(data2(), metadata2())
      
    })
    
    output$misvals <- renderPlot({
      
      req(input$file)
      req(input$meta)
      
      plot_na(data2(), metadata2())
      
    })
    
    
    # Download tables
    
    listy <- reactive({
      
      req(input$file)
      req(input$meta)
      req(input$comp)
      
      list <- all_analysis(data2(), metadata2(), comp(), thres = input$slider)
      return(list)
      
    })
    
    
    output$tables <- downloadHandler(
      filename = function() {
        paste0("DIA_data_", Sys.Date(), ".xlsx")},
      content = function(file) {
        openxlsx::write.xlsx(listy()$tables, file)
      }
    )
    
    # Download figures
    
    output$figs <- downloadHandler(
      filename = function() paste0("DIA_figures_", Sys.Date(), ".pptx"),

      content = function(file) {
        file_pptx <- tempfile(fileext = ".pptx")
        gen_pptx(listy()$figs, file_pptx)
        file.rename(from = file_pptx, to = file )
      }
    )
    
  
    # Generate report
    
    output$report <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "report.html",
      content = function(file) {
        
        withProgress(message = 'Rendering, please wait!', value = 0, {
        
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed)
        tempReport <- file.path(tempdir(), "dia_analysis_report.Rmd")
        file.copy("dia_analysis_report.Rmd", tempReport, overwrite = TRUE)

        # Set up parameters to pass to Rmd document
        params <- list(data = input$file$datapath,
                       data_name = input$file$name,
                       metadata = input$meta$datapath,
                       comp = input$comp$datapath,
                       filter = input$slider,
                       rm.samples = input$samples,
                       rendered_by_shiny = TRUE)

        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app)
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
                          )
        })
      }
    )
  }


# Run the app -------------------------------------------------------------

shinyApp(ui = ui, server = server)
