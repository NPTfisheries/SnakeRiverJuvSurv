library(shiny)

# UI
ui = fluidPage(
  titlePanel("Snake River Juvenile Marking"),
  mainPanel(style = 'width:100% !important; height:90vh;',
    helpText('Please be patient while the map loads.'),
    includeHTML(path = 'leaflet/sr_juv_marking_2025.html')
    )
)

# Server
server = function(input, output) {
}

# App 
shinyApp(ui = ui, server = server)
