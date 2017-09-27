# User can select and correlate data time points vs another time point on another gene

library(shiny)
dat <- read.table("C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\Lab2\\spellman.txt",header=T,row.names=1)
dat <- as.data.frame(dat);

cdc15 <- dat[, 23:46]


shinyServer(function(input, output, session) {
  # Get the data from the variables declared on the ui.R file
  selectedData <- reactive({cdc15[, c(input$xcol, input$ycol)]})

  # Create the plot
  output$plot <- renderPlot({plot(selectedData(), pch = 20, cex = 1, col = "black",
                                  main = "CDC15 Coorelation")})
})

