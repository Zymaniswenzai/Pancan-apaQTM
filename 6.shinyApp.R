library(shiny)
library(DT)


ui <- fluidPage(
  titlePanel("Pancan apaQTM v0.0.1"),
  sidebarLayout(
    sidebarPanel(
      selectInput("cancer", "Choose Cancer Type", choices = c('',"BLCA", "BRCA",  "COAD")),
      selectInput("type", "Choose cis or trans", choices = c('',"cis", "trans")),
      textInput("cpg", "Input methylation site", placeholder = "e.g. cg08269402"),
      actionButton("search_button", "Search"),
      br(),
      br(),
      downloadButton("downloadZipBtn", "Download ZIP file")
    ),
    mainPanel(
      dataTableOutput("results_table")
      
    )
  ))

server <- function(input, output, session) {
  loaded_data <- reactiveVal(NULL)
  
  observe({
    # 替换以下路径为您自己的数据文件的路径
    data_file <- "data/cancer_aqtm.txt.gz"
    loaded_data(read.table(gzfile(data_file),sep='\t',header=T))
  })
  
  # 处理搜索按钮点击
  observeEvent(input$search_button, {
    # 检查数据是否已加载
    if (is.null(loaded_data())) {
      return(NULL)
    }
    
    filter_data <- loaded_data()
    if (input$cancer != "All") {
      filter_data <- filter_data[filter_data$cancer == input$cancer, ]
    }
    if (input$type != "All") {
      filter_data <- filter_data[filter_data$type == input$type, ]
    }
    if (!is.null(input$cpg) && input$cpg != "") {
      selected_sites <- as.character(input$cpg)
      filter_data <- filter_data[filter_data$cpg %in% selected_sites, ]
    }
    output$results_table <- renderDataTable({
      filter_data
    })
  })
  
  
  zip_file_path <- "data/aqtm_data.zip"  # 替换为您的ZIP文件的路径和名称
  
  output$downloadZipBtn <- downloadHandler(
    filename = function() {
      "downloaded_file.zip"  # 设置下载的文件名
    },
    content = function(file) {
      file.copy(zip_file_path, file)
    }
  )
  
  output$downloadResult <- downloadHandler(
    filename = function() {
      "search_results.csv"  # 设置下载的文件名
    },
    content = function(file) {
      write.csv(filtered_data, file)
    }
  )
  
}

shinyApp(ui, server)