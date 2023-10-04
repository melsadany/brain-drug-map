################################################################################
#                         make a shiny app for drug maps                       #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(shiny)
library(plotly)
library(ggseg3d)
library(ggsegDKT)
library(ggseg)
scene <- list(camera = list(eye = list(x = 1.5, y = 0, z = 0)))
#################################################################################
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
source("src/functions/brain_3d.R")
#################################################################################
#################################################################################
drug.correlations <- read_csv("/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map/data/all-drug-map-activity-averaged-by-anatomical-region-092723.csv")
file_names <- colnames(drug.correlations)[-1]

# Define UI
ui <- fluidPage(
  titlePanel("Predicted drug activity maps"),
  
  # Description text
  p("The map is built as a correlation between the drug transcriptomic signature from the CMAP dataset and the predicted brain gene expression from a deep learning model."),
  
  # Sidebar for selection and table
  sidebarLayout(
    sidebarPanel(
      selectInput("fileSelect", "Choose a drug to view its activity map.", 
                  choices = file_names, 
                  width = "200px",
                  selectize = TRUE, 
                  selected = "methylphenidate"),
      hr(),
      helpText("The drug list provided here is FDA-approved and has a probability of passing the BBB > 0.5"), 
      width = 3,
      style = "margin-bottom:2px;", 
      uiOutput("data_table_ui")  # Define the UI for the data table
    ),
    
    # Main panel for plot
    mainPanel(
      uiOutput("plotFrame")  # Define the UI for the plot
    )
  ),
  tags$head(
    tags$style(
      HTML(
        "
      .shiny-plotly-output-container {
        height: 100% !important;
      }, 
      .dataTables_filter input {
        width: 50px !important; /* Adjust the width as needed */
      }
      "
      )
    )
  )
)

# Define server
server <- function(input, output) {
  output$plotFrame <- renderUI({
    brain_3d(drug.correlations, drug = input$fileSelect)
  })
  
  # Define the UI for the data table
  output$data_table_ui <- renderUI({
    dataTableOutput("data_table")
    
  })
  
  # Render the data table
  output$data_table <- renderDataTable({
    data <- read_csv("/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map/data/all-drug-map-activity-averaged-by-anatomical-region-092723.csv")
    data[, c(1, which(colnames(data) == input$fileSelect))]
  }, options = list(pageLength = 10))
  
}

shinyApp(ui, server)
#################################################################################