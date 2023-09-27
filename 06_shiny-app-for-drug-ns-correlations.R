################################################################################
#           make a shiny app for drug correlations with neurosynth             #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(shiny)
library(plotly)
#################################################################################
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
#################################################################################
# read correlations
data <- pdsload(fname = "data/drug-features-correlations-by-xyz-in-predicted-MNI-wpval-beh.rds.pxz") %>% 
  as.data.frame() %>%
  rownames_to_column("c_drug") %>%
  pivot_longer(cols = colnames(.)[-1], names_to = "feature", values_to = "value")


# Define UI
ui <- fluidPage(
  titlePanel("drug correlations with Neurosynth terms"),
  
  fluidRow(
    column(
      width = 12,
      style = "padding-left: 15px;",  # Add padding to the left
      selectizeInput("x_var", "Select drugs of interest:", 
                     choices = unique(data$c_drug), 
                     width = "100%",
                     multiple = T, 
                     selected = c("methylphenidate", "sertraline", "venlafaxine", "fluoxetine",
                                  "caffeine", "bupropion", "trazodone", "zolpidem",
                                  "atomoxetine", "citalopram", "clonidine", "clonazepam",
                                  "clozapine", "diazepam", "escitalopram", "fluvoxamine",
                                  "haloperidol", "histamine", "ibuprofen", "ketamine",
                                  "lidocaine", "nicotine", "orlistat", "paroxetine",
                                  "risperidone", "sumatriptan", "topiramate")),
      style = "margin-bottom: 15px;",  # Add margin at the bottom
      selectizeInput("y_var", "Select terms/features:", 
                     choices = unique(data$feature), 
                     width = "100%",
                     multiple = T, 
                     selected = c("attention", "visual", "language", "depression"))
    )
  ),
  
  # Main panel with the heatmap
  fluidRow(
    plotlyOutput("heatmap_plot", height = "auto")
  )
)

# Define server
server <- function(input, output) {
  # Create a reactive subset of the data based on selected x and y variables
  filtered_data <- reactive({
    data %>%
      filter(c_drug %in% input$x_var,
             feature %in% input$y_var) %>%
      select(x = c_drug, y = feature, value)
  })
  # print(filtered_data)
  
  # Calculate the number of unique x and y items
  num_x_items <- reactive({
    length(unique(filtered_data()$x))
  })
  
  num_y_items <- reactive({
    length(unique(filtered_data()$y))
  })
  
  # # Automatically set height and width based on the number of items
  # auto_height <- reactive({
  #   # Calculate the height based on the number of y items
  #   height_factor <- max(1, num_y_items() / 5)  # Adjust the divisor as needed
  #   height_px <- 400 * height_factor  # Adjust the base height as needed
  #   paste0(height_px, "px")
  # })
  # 
  # auto_width <- reactive({
  #   # Calculate the width based on the number of x items
  #   width_factor <- max(1, num_x_items() / 5)  # Adjust the divisor as needed
  #   width_percent <- 100 / width_factor  # Adjust the base width as needed
  #   paste0(width_percent, "%")
  # })
  
  # Render the heatmap plot
  output$heatmap_plot <- renderPlotly({
    heatmap_data <- filtered_data()
    plot_ly(data = heatmap_data, x = ~x, y = ~y, z = ~value, type = "heatmap", colors = c(redblu.col[2], "white", redblu.col[1])) %>%
      layout(
        xaxis = list(title = "", tickangle = -90),
        yaxis = list(title = ""),
        title = "")
  })
  # Generate the plotlyOutput with auto-adjusted height and width
  output$heatmap_output <- renderUI({
    # plotlyOutput("heatmap_plot", height = auto_height(), width = auto_width())
    plotlyOutput("heatmap_plot", height = "auto", width = "auto")
  })
}

# Run the Shiny app
shinyApp(ui, server)

