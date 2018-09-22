library(shiny)

# Define UI for allele binning app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("MicroBinR"),
  
  # Tabset layout ----
  tabsetPanel(
    
    # Define upload and binning tab ----
    tabPanel("Upload",
             
             # Tab title ----
             titlePanel("Upload"),
             
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               
               # Sidebar panel for inputs ----
               sidebarPanel(
                 
                 # Input: Select a file ----
                 fileInput("file1", "Choose CSV File",
                           multiple = TRUE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 
                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Select separator ----
                 radioButtons("sep", "Separator",
                              choices = c(Comma = ",",
                                          Semicolon = ";",
                                          Tab = "\t"),
                              selected = ","),
                 
                 # Input: Select quotes ----
                 radioButtons("quote", "Quote",
                              choices = c(None = "",
                                          "Double Quote" = '"',
                                          "Single Quote" = "'"),
                              selected = '"'),
                 
                 
                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Silder select ploidy ----
                 sliderInput("ploidy", "Ploidy Level", min = 1, max = 6, value = 2, step = 1),
                 
                 # Input: Action Button bin ----
                 actionButton(inputId = "binAll", label = "Go!"),
                 
                 # Loading indicator ----
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  tags$div("Loading...",id="loadmessage"))
                 
               ),
               
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 # Output: Data file ----
                 tableOutput("unbinned")
                 
               )
               
             )
    ),
    
    # Define Bin tab ----
    tabPanel("Bin",
             
             # Tab Title ----
             titlePanel("Visualize and Bin"),
             
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               
               # Sidebar panle for input ----
               sidebarPanel(
                 
                 # Input: Select a locus ----
                 htmlOutput("selectLocus"), 
                 
                 # Input: Sigma slider ----
                 sliderInput("sgma", "Select sigma", min = 0.1, max = 1, value = 0.4, step = 0.05),
                 
                 
                 # Input: Range slider ----
                 sliderInput("range", "Select Range", min = 5, max = 25, value = 20, step = 1),
                 
                 # Input: Action Button recalculate ----
                 actionButton("recalculate", "recalculate"),
                 
                 # Loading indicator ----
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  tags$div("Loading...",id="loadmessage"))
                 
               ),
               
               
               tabsetPanel(
                 tabPanel("plot", 
                          # Output: Density plot ----
                          plotOutput("plot")
                 ),
                 tabPanel("Stats",
                          # Output: bin statisitcs ----
                          tableOutput("stat")
                 )
                 
                 
               )
               
               
             )
    ),
    
    # Define Download tab ----
    tabPanel("Download",
             
             # Tab titel ----
             titlePanel("Download"),
             
             # Sidebar layout ----
             sidebarLayout(
               
               # Sidebar Pannel ----
               sidebarPanel(
                 downloadButton("download")
                 
               ),
               
               # Main Pannel ----
               mainPanel(
                 tableOutput("binned")
               )
             )
             
    )
    
  )
)

# Define server logic to read selected file ----
server <- function(input, output, session) {
  
  
  # First tab  ----------------------------------------------------------------------------------------------
  repo1 <- reactive({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, Database$repo will be generated.
    
    
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = TRUE,
                   sep = input$sep,
                   quote = input$quote)
    
  repo <- MicroBinner::sat_makeDatabase(df, input$ploidy)
  
  return(repo)
  
  })
  
  # generate a table of binned data ----
  output$unbinned <- renderTable({
    
    req(repo1())
    
    unbinned <- MicroBinner::sat_return(repo1(), "data")
    
    return(unbinned)
    
  })
  
  # Update locus selection ----
  output$selectLocus <- renderUI({
    
    if (input$binAll == FALSE) {
      x <- "No Loci Avalible"
    } else {
      x <- names(repo1())
    }
    
    selectInput("locus", "Select a locus", x)
    
  })
  
  # Second Tab -----------------------------------------------------------------------------------------
  
   #
  # Build reactive data repository ----
  data <- reactiveValues(repo = NULL)
  
  # Bin all loci with standard settings ----
 observeEvent(input$binAll, {
   
   data$repo <- MicroBinner::sat_assess(repo1())
   
 })

  # when settings are changed recalculate specific loci ----
  observeEvent(input$recalculate, {
    locus <- names(data$repo)
    s <- purrr::map_dbl(locus, ~data$repo[[.]]$info$sgma)
    r <- purrr::map_dbl(locus, ~data$repo[[.]]$info$range)
    
    n <- seq_along(locus)[locus == input$locus]
    
    s[n] <- input$sgma
    r[n] <- input$range
    
    data$repo <- MicroBinner::sat_assess(data$repo, h = s, range = r)
  })
  
  
  # generate plot ----
  output$plot <- renderPlot({
    
    req(input$binAll, data$repo)
    
    return(data$repo[[input$locus]]$plot)
    
  })
  
  # generate bin stats ----
  output$stat <- renderTable({
    
    req(input$binAll)
    
    return(data$repo[[input$locus]]$stat)
  })
  
  # Third tab -------------------------------------------------------------------------------------------------------
  
  # generate binned data ----
  output$binned <- renderTable({
    
    req(input$binAll)
    
    return(MicroBinner::sat_getBinned(MicroBinner::sat_assign(data$repo)))
  })
  
  # generate download object ----
  output$download <- downloadHandler("Binned.csv",
                                     content = function(file) {
                                       write.csv(MicroBinner::sat_getBinned(MicroBinner::sat_assign(data$repo)), file)
                                       }
                                     )
}

shinyApp(ui, server)
