################################################################################
#                   make a shiny app for gene maps in MRI theme                #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(shiny)
library(RNifti)
library(shinythemes)
library(shinyWidgets)
library(markdown)
#################################################################################
#################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map"
setwd(project.dir)
#################################################################################
#################################################################################
# needed functions
lazyr <- function(inputs, FUN) {
  obj = new.env(parent = globalenv())
  
  len = length(inputs)
  obj$FUN = FUN
  obj$inputs = inputs
  obj$checks = vector(mode = "logical", length = len)
  obj$returns = rep(NA, len)
  obj$length = len
  
  class(obj) = "lzr_obj"
  return(obj)
}
image_render <- function(x, h, v) {
  par(oma = rep(0, 4), mar = rep(0, 4), bg = "black")
  graphics::image(1:dim(x)[1], 1:dim(x)[2], x,
                  col = colorRampPalette(redblu.col[c(2,1)])(64), asp = 1,
                  xlab = "", ylab = "", axes = FALSE, 
                  useRaster = T)
  abline(h = h, v = v, col = "red")
}

raster3d_interactive_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(4, offset = 8, uiOutput(ns("time_slider_UI")))
    ),
    fluidRow(
      lapply(c("x", "y", "z"), function(x) {
        column(4, align = "center", uiOutput(ns(paste0(x, "_slider_UI"))))
      })
    ),
    fluidRow(
      lapply(c("x_plot", "y_plot", "z_plot"), function(x) {
        column(4, align = "center", 
               plotOutput(ns(x), click = ns(paste0(x, "_click"))))
      })
    )
  )
}

raster3d_interactive_Module <- function(input, output, session, im) {
  ns <- session$ns
  
  time_im <- reactive({
    req(im())
    if (length(dim(im())) == 3) return(im())
    return(im()[,,, input$time_slider])
  })
  
  output$x_plot <- renderPlot({
    req(im())
    image_render(time_im()[input$x_slider,,], input$z_slider, input$y_slider)
  })
  
  output$y_plot <- renderPlot({
    req(im())
    image_render(time_im()[, input$y_slider, ], input$z_slider, input$x_slider)
  })
  
  output$z_plot <- renderPlot({
    req(im())
    image_render(time_im()[,, input$z_slider], input$y_slider, input$x_slider)
  })
  
  output$x_slider_UI <- renderUI({
    req(im())
    sliderInput(
      ns("x_slider"), label = NULL, min = 1, max = dim(im())[1],
      value = ceiling(dim(im())[1] / 2)
    )
  })
  
  output$y_slider_UI <- renderUI({
    req(im())
    sliderInput(
      ns("y_slider"), label = NULL, min = 1, max = dim(im())[2],
      value = ceiling(dim(im())[2] / 2)
    )
  })
  
  output$z_slider_UI <- renderUI({
    req(im())
    sliderInput(
      ns("z_slider"), label = NULL, min = 1, max = dim(im())[3],
      value = ceiling(dim(im())[3] / 2)
    )
  })
  
  output$time_slider_UI <- renderUI({
    req(im())
    if (length(dim(im())) == 3) return(NULL)
    sliderInput(ns("time_slider"), label = "Time", min = 1, max = dim(im())[4],
                value = 1, step = 1)
  })
  
  observeEvent(input$y_plot_click, {
    updateSliderInput(session = session, "x_slider",
                      value = input$y_plot_click$x)
    updateSliderInput(session = session, "z_slider",
                      value = input$y_plot_click$y)
  })
  
  observeEvent(input$x_plot_click, {
    updateSliderInput(session = session, "y_slider",
                      value = input$x_plot_click$x)
    updateSliderInput(session = session, "z_slider",
                      value = input$x_plot_click$y)
  })
  
  observeEvent(input$z_plot_click, {
    updateSliderInput(session = session, "x_slider",
                      value = input$z_plot_click$x)
    updateSliderInput(session = session, "y_slider",
                      value = input$z_plot_click$y)
  })
}
print.lzr_obj <- function(obj) {
  print(list("inputs" = obj$inputs, "returns" = obj$returns))
}

`[.lzr_obj` <- function(obj, idx, ...) {
  if (any(idx > obj$length)) stop("Out of bounds")
  if (is.na(obj$checks[idx]) | !obj$checks[idx]) {
    obj$returns[idx] = obj$FUN(obj$inputs[idx], ...)
    obj$checks[idx] = T
  }
  return(obj$returns[idx])
}

# Cache invalidation
`[<-.lzr_obj` <- function(obj, idx, value) {
  obj$inputs[idx] <- value
  obj$checks[idx] <- F
  invisible(obj)
}
raster3d_animation_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(4, tagList(
        actionButton(ns("zoom_m"), "", icon = icon("minus-square-o"), width = "40px", 
                     style = "border-radius: 25px; padding: 0px;"),
        actionButton(ns("zoom_p"), "", icon = icon("plus-square-o"), width = "40px", 
                     style = "border-radius: 25px; padding: 0px;")
      )),
      column(4, offset = 4, uiOutput(ns("time_slider_UI")))
    ),
    fluidRow(
      lapply(c("x", "y", "z"), function(x) {
        column(4, align = "center", uiOutput(ns(paste0(x, "_slider_UI"))))
      })
    ),
    fluidRow(
      lapply(c("x_plot", "y_plot", "z_plot"), function(x) {
        column(4, align = "center", plotOutput(ns(x)))
      })
    )
  )
}

raster3d_animation_Module <- function(input, output, session, im) {
  temp_dir <- file.path(tempdir(), "shinyMRI")
  dir.create(temp_dir, showWarnings = FALSE)
  img_slice <- function(x) {
    temp_img <- tempfile("", temp_dir, fileext = ".png")
    png(temp_img)
    par(oma = rep(0, 4), mar = rep(0, 4), bg = "black")
    graphics::image(1:dim(x)[1], 1:dim(x)[2], x,
                    col = gray(0:64/64), asp = 1,
                    xlab = "", ylab = "", axes = FALSE, 
                    useRaster = T)
    dev.off()
    return(temp_img)
  }
  
  ns <- session$ns
  rv <- reactiveValues(zoom = 1)
  
  observeEvent(im(), {
    im_dim <- dim(im())
    if (length(im_dim) == 3) {
      rv$is_4d <- FALSE
      rv$x_p <- list(lazyr(seq(im_dim[1]), function(x) {img_slice(im()[x,,])}))
      rv$y_p <- list(lazyr(seq(im_dim[2]), function(x) {img_slice(im()[,x,])}))
      rv$z_p <- list(lazyr(seq(im_dim[3]), function(x) {img_slice(im()[,,x])}))
    } else {
      rv$is_4d <- TRUE
      rv$x_p <- lapply(seq(im_dim[4]), function(t) {
        lazyr(seq(im_dim[1]), function(x) {img_slice(im()[x,,,t])})
      })
      rv$y_p <- lapply(seq(im_dim[4]), function(t) {
        lazyr(seq(im_dim[2]), function(x) {img_slice(im()[,x,,t])})
      })
      rv$z_p <- lapply(seq(im_dim[4]), function(t) {
        lazyr(seq(im_dim[3]), function(x) {img_slice(im()[,,x,t])})
      })
    }
  })
  
  time_slot <- reactive({
    if (is.null(input$time_slider)) return(1)
    input$time_slider
  })
  
  output$x_plot <- renderImage({
    req(input$x_slider)
    index_buff <- min(input$x_slider, dim(im())[1])
    if (index_buff > dim(im())[1]) index_buff <- ceiling(dim(im())[1]/2)
    list(
      src = rv$x_p[[time_slot()]][index_buff],
      width = dim(im())[2] * rv$zoom,
      height = dim(im())[3] * rv$zoom
    )
  }, deleteFile = F)
  
  output$y_plot <- renderImage({
    req(input$y_slider)
    index_buff <- input$y_slider
    if (index_buff > dim(im())[2]) index_buff <- ceiling(dim(im())[2]/2)
    list(
      src = rv$y_p[[time_slot()]][index_buff],
      width = dim(im())[1] * rv$zoom,
      height = dim(im())[3] * rv$zoom
    )
  }, deleteFile = F)
  
  output$z_plot <- renderImage({
    req(input$z_slider)
    index_buff <- input$z_slider
    if (index_buff > dim(im())[3]) index_buff <- ceiling(dim(im())[3]/2)
    list(
      src = rv$z_p[[time_slot()]][input$z_slider],
      width = dim(im())[1] * rv$zoom,
      height = dim(im())[2] * rv$zoom
    )
  }, deleteFile = F)
  
  output$x_slider_UI <- renderUI({
    req(im())
    sliderInput(
      ns("x_slider"), label = NULL, min = 1, max = dim(im())[1], step = 1,
      value = ceiling(dim(im())[1] / 2), 
      animate = animationOptions(interval = 100, loop = TRUE)
    )
  })
  
  output$y_slider_UI <- renderUI({
    req(im())
    sliderInput(
      ns("y_slider"), label = NULL, min = 1, max = dim(im())[2], step = 1,
      value = ceiling(dim(im())[2] / 2), 
      animate = animationOptions(interval = 100, loop = TRUE)
    )
  })
  
  output$z_slider_UI <- renderUI({
    req(im())
    sliderInput(
      ns("z_slider"), label = NULL, min = 1, max = dim(im())[3], step = 1,
      value = ceiling(dim(im())[3] / 2), 
      animate = animationOptions(interval = 100, loop = TRUE)
    )
  })
  
  output$time_slider_UI <- renderUI({
    req(im())
    if (length(dim(im())) == 3) return(NULL)
    sliderInput(ns("time_slider"), NULL, min = 1, max = dim(im())[4],
                value = 1, step = 1, 
                animate = animationOptions(interval = 100, loop = TRUE))
  })
  
  observeEvent(input$zoom_p, {
    if (rv$zoom < 1.5) rv$zoom <- rv$zoom + 0.05
  })
  observeEvent(input$zoom_m, {
    if (rv$zoom > 0.5) rv$zoom <- rv$zoom - 0.05
  })
  
  session$onSessionEnded(function() {
    unlink(temp_dir, recursive = T)
  })
}
read_img_as_array <- function(path) {
  img_raw <- RNifti::readNifti(path)
  if (length(dim(img_raw)) == 3) return(img_raw[,,])
  return(img_raw[,,,])
}
#################################################################################
#################################################################################
#################################################################################
####

# Get full file paths
gene_paths <- list.files("/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map/data/maps/model-092723/gene-exp/", 
                         full.names = F)
# Filter out files with underscores
genes <- gene_paths[!grepl("_", basename(gene_paths))]

# Define UI
ui <- fluidPage(
  titlePanel("Predicted gene expression map"),
  theme = shinytheme("cyborg"),
  # Description text
  p("The map represents the predicted gene expression in brain regions using a deep learning model. The model was trained on the Allen Institue microarray gene expression data"),
  tabPanel(
    "Home",
    fluidRow(
      column(4, selectizeInput("fileSelect", "Choose a gene to view its expression map.",
                                                choices = sub(".nii.gz","", genes),
                                                width = "200px",
                                                selected = "DCC")),
    ),
    uiOutput("raster_panel")
  )
)

# Define server
server <- function(input, output) {
  options(shiny.maxRequestSize = 500*1024^2)
  
  app_dt <- reactive({
    datapath <- paste0("/Dedicated/jmichaelson-wdata/msmuhammad/projects/brain-drug-map/data/maps/model-092723/gene-exp/", input$fileSelect, ".nii.gz")
    out <- read_img_as_array(datapath)
    out[out == 0] <- NA
    return(out)
  })
  output$raster_panel <- renderUI({
    callModule(raster3d_interactive_Module, "mri_3d", im = app_dt)
    raster3d_interactive_UI("mri_3d")
  })
}

shinyApp(ui, server)
#################################################################################
