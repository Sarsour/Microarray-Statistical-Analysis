# User can select and correlate data time points vs another time point on another gene

library(shiny)
dat <- read.table("C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\Lab2\\spellman.txt",header=T,row.names=1)
dat <- as.data.frame(dat);

cdc15 <- dat[, 23:46]

shinyUI(pageWithSidebar(
  # Application title
  titlePanel('Shiny Application - Lab 2'),
  
  # Sidebar with 2 select inputs and a numeric input
  sidebarPanel(
    selectInput('xcol', 'X Variable', names(cdc15)),
    selectInput('ycol', 'Y Variable', names(cdc15))),
  
  # Shows the plot
  mainPanel(plotOutput('plot'))
))