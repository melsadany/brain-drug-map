################################################################################
#            make a shiny app for drug correlations with BrainMap              #
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
data <- pdsload(fname = "data/drug-correlations-w-BrainMap-mni-maps.rds.pxz") %>% 
  as.data.frame() %>%
  rename(c_drug = V2, feature_lab = V1, value = r)
beh.meta <- read_csv("/wdata/msmuhammad/data/BrainMap/behavioral-domain/metadata.csv")
data <- inner_join(data, beh.meta %>% rename(feature_lab = full_label, feature = sub_domain))

# Define UI
ui <- fluidPage(
  titlePanel("drug correlations with BrainMap behavior domains"),
  sidebarLayout(
    sidebarPanel(
      
      fluidRow(
        column(12,
          # style = "padding-left: 15px;",  # Add padding to the left
          selectizeInput("x_var", "Select drugs of interest:", 
                         choices = unique(data$c_drug), 
                         # width = "100%",
                         multiple = T, 
                         selected = c("methylphenidate", "sertraline", "venlafaxine", "fluoxetine",
                                      "caffeine", "bupropion", "trazodone", "zolpidem",
                                      "atomoxetine", "citalopram", "clonidine", "clonazepam",
                                      "clozapine", "diazepam", "escitalopram", "fluvoxamine",
                                      "haloperidol", "histamine", "ibuprofen", "ketamine",
                                      "lidocaine", "nicotine", "orlistat", "paroxetine",
                                      "risperidone", "sumatriptan", "topiramate")))),
          # fluidRow(style = "margin-bottom: 15px;"),  # Add margin at the bottom
          fluidRow(
            column(12,
                   # style = "padding-left: 15px;",  # Add padding to the left
                   selectizeInput("y_var", "Select terms/features:", 
                         choices = unique(data$feature), 
                         # width = "100%",
                         multiple = T, 
                         selected = beh.meta$sub_domain[1:10]),
                   ),
          ),
      ),
    # Main panel with the heatmap
    mainPanel(
      plotlyOutput("heatmap_plot")
    ))
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
  
  # # Calculate the number of unique x and y items
  # num_x_items <- reactive({
  #   length(unique(filtered_data()$x))
  # })
  # 
  # num_y_items <- reactive({
  #   length(unique(filtered_data()$y))
  # })
  # Render the heatmap plot
  output$heatmap_plot <- renderPlotly({
    heatmap_data <- filtered_data()
    plot_ly(data = heatmap_data, x = ~x, y = ~y, z = ~value, 
            type = "heatmap", colors = c(redblu.col[2], "white", redblu.col[1])) %>%
      layout(
        xaxis = list(title = "", tickangle = -90),
        yaxis = list(title = ""),
        title = "")
  })
  # Generate the plotlyOutput with auto-adjusted height and width
  output$heatmap_output <- renderUI({
    # plotlyOutput("heatmap_plot", height = auto_height(), width = auto_width())
    plotlyOutput("heatmap_plot", height = "100%", width = "auto")
  })
}

# Run the Shiny app
shinyApp(ui, server)

