################################################################################
#                         make a shiny app for gene maps                       #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(shiny)
library(oro.nifti)
library(plotly)
# library(ggseg3d)
# library(ggsegDKT)
# library(ggseg)
scene <- list(camera = list(eye = list(x = 1.5, y = 0, z = 0)))
#################################################################################
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
#################################################################################
# Get full file paths
gene_paths <- list.files("/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map/data/maps/model-092723/gene-exp/", 
                         full.names = F)
# Filter out files with underscores
genes <- gene_paths[!grepl("_", basename(gene_paths))]

# Define UI
ui <- fluidPage(
  titlePanel("Predicted gene expression map"),

  # Description text
  p("The map represents the predicted gene expression in brain regions using a deep learning model. The model was trained on the Allen Institue microarray gene expression data"),
  tags$head(
    tags$style(
      HTML(
        "
        .sidebar {
          width: 250px; /* Adjust the width as needed */
        }
        "
      )
    )
  ),
  # Sidebar for selection and table
  sidebarLayout(fluid = F,
    sidebarPanel(
      selectizeInput("fileSelect", "Choose a gene to view its expression map.",
                     choices = sub(".nii.gz","", genes),
                     width = "200px",
                     selected = "DCC"),
      # sliderInput("x", "X Coordinate", min = -91, max = 90, value = 0),
      # sliderInput("y", "Y Coordinate", min = -126, max = 91, value = 0),
      # sliderInput("z", "Z Coordinate", min = -72, max = 109, value = 0)
      fluidRow(
        column(3, numericInput("x", "X ", min = -91, max = 90, value = 0, width = "200px")),
        column(3, numericInput("y", "Y ", min = -126, max = 91, value = 0, width = "200px")),
        column(3, numericInput("z", "Z ", min = -72, max = 109, value = 0, width = "200px"))
      )
      # numericInput("x", "X Coordinate", min = -91, max = 90, value = 0, width = "200px"),
      # numericInput("y", "Y Coordinate", min = -126, max = 91, value = 0, width = "200px"),
      # numericInput("z", "Z Coordinate", min = -72, max = 109, value = 0, width = "200px")
    ),
    mainPanel(
      fluidRow(
        column(4, plotlyOutput("viewX", width = "400px", height = "450px")),
        column(4, plotlyOutput("viewY", width = "400px", height = "450px")),
        column(4, plotlyOutput("viewZ", width = "400px", height = "450px"))
      )
    )
    
  )
)

# Define server
server <- function(input, output) {
  observe({
    req(input$fileSelect)
    # print(input$fileSelect)
    img <- readNIfTI(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map/data/maps/model-092723/gene-exp/", input$fileSelect, ".nii.gz"))
    
    origin <- c(91,127,73)
    img.hit <- which(abs(img)>0, arr.ind = T) %>%
      as.data.frame() %>%
      mutate(mni_x = -1*(dim1-origin[1])) %>%
      mutate(mni_y = dim2-origin[2]) %>%
      mutate(mni_z = dim3-origin[3])
    img.hit.whole <- img.hit %>%
      mutate(exp = img[as.matrix(img.hit[1:3])])
    
    output$viewX <- renderPlotly({
      plot_ly(data = img.hit.whole %>%
                filter(mni_x == input$x),
                # filter(mni_x == 0),
              x = ~mni_x, 
              y = ~mni_y, 
              z = ~mni_z, 
              color = ~exp,
              mode = "markers", type = "scatter3d",
              marker = list(size = 4)) %>% 
        layout(scene =  list(camera = list(eye = list(x = 2, y = 0, z = 0))))
    })
    output$viewY <- renderPlotly({
      plot_ly(data = img.hit.whole %>%
                filter(mni_y == input$y),
                # filter(mni_y == 0),
              x = ~mni_x, 
              y = ~mni_y, 
              z = ~mni_z, 
              color = ~exp,
              mode = "markers", type = "scatter3d",
              marker = list(size = 4)) %>% 
        layout(scene =  list(camera = list(eye = list(x = 0, y = -2, z = 0))))
    })
    output$viewZ <- renderPlotly({
      plot_ly(data = img.hit.whole %>%
                filter(mni_z == input$z),
                # filter(mni_z == 0),
              x = ~mni_x, 
              y = ~mni_y, 
              z = ~mni_z, 
              color = ~exp,
              mode = "markers", type = "scatter3d",
              marker = list(size = 4)) %>% 
        layout(scene =  list(camera = list(eye = list(x = 0, y = 0, z = 2))))
      
    })
    # Add observeEvent for the plot clicks
    observeEvent(event_data("plotly_click"), {
      click <- event_data("plotly_click")
      if (!is.null(click)) {
        session <- getDefaultReactiveDomain()
        updateTextInput(session, "x", value = round(click$x))
        updateTextInput(session, "y", value = round(click$y))
        updateTextInput(session, "z", value = round(click$z))
      }
    })
  })
  }

shinyApp(ui, server)
#################################################################################
