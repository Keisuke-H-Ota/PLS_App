
rm( list=ls(all=TRUE) ) 
library(shiny)
library(tidytidbits) #isTruthy() に必要
library(ggplot2)
library(plotly)

# Define UI for dataset viewer app ----
ui <- fluidPage(
  navbarPage(
    title = "PLS ROG App",
    tabPanel("Analysis results",
      sidebarLayout(
        sidebarPanel(
          width = 4,
          uiOutput("tab"),
          fileInput("RawData", "Choose TSV File", accept = ".tsv"),
          selectInput(
            inputId = "normalization",
            label = "Normalization method",
            choices = c("none","composition"),
            selected = NULL
          ),
          selectInput(
            inputId = "select",
            label = "Analytics method",
            choices = c("naïve pls","pls-rog"),
            selected = NULL
          ),
          sliderInput(inputId = "level",
                        label = "Confidence interval",
                        min = 0.01,
                        max = 0.99,
                        value = 0.95),
          conditionalPanel(
            condition = "input.select == 'pls-rog'",
            sliderInput(inputId = "kappa",
                        label = "Kappa (smoothing parameter)",
                        min = 0.00,
                        max = 0.99,
                        value = 0.00)
          ),
          selectInput(inputId = "Axis1", 
                      label = "Axes1", 
                      choices = ""),
          selectInput(inputId = "Axis2", 
                      label = "Axes2", 
                      choices = ""),
        ),
        mainPanel(
          tabsetPanel(type = "tabs",
           tabPanel("score_X",
            plotlyOutput(
              outputId = "score_X_plots",
              width = "700px",
              height = "600px",
              inline = FALSE
              )
            ),
           tabPanel("score_Y",
            plotlyOutput(
              outputId = "score_Y_plots",
              width = "700px",
              height = "600px",
              inline = FALSE
              )
            ),
           tabPanel("Weight_X",
            plotlyOutput(
              outputId = "Wx_plots",
              width = "600px",
              height = "600px",
              inline = FALSE
              )
            ),
           tabPanel("R",
            plotlyOutput(
              outputId = "R_plots",
              width = "600px",
              height = "600px",
              inline = FALSE
              )
            ),
           tabPanel("Q",
            plotlyOutput(
              outputId = "Q_plots",
              width = "600px",
              height = "600px",
              inline = FALSE
              )
            )
          )
        )
      )
    ),
    tabPanel("Raw Data",
      sidebarLayout(
        sidebarPanel(
          width = 4,
          selectInput(inputId = "feature", 
                      label = "select feature", 
                      choices = ""
                      )
        ),
        mainPanel(
          plotlyOutput(
            outputId = "boxplot",
            width = "600px",
            height = "600px",
            inline = FALSE
            ),
        )
      )
    )
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {
  url <- a("here", href="https://github.com/Keisuke-H-Ota/PLS_App")
    output$tab <- renderUI({
      tagList("Refer input data format", url)
    })
  observeEvent(input$RawData, {
    # データの読み込み
    mbx <- read.table(file=input$RawData$datapath, sep="\t",header=F, quote='', comment.char="")
    nr <- nrow(mbx) # Number of rows
    nc <- ncol(mbx) # Numbers of columns
    sn <- mbx[2:nr,1] # Sample Name
    fn <- as.vector(t(mbx)[2:(nc-1),1]) # Feature Name
    gr <- mbx[2:nr,nc] # Group
    DF <- mbx[2:nr,2:(nc-1)]
    DF <- data.frame(sapply(DF, function(x) as.numeric(as.character(x))))
    
    # サンプルごとに欠損値を最小値の半分で置換
    for(i in 1:nrow(DF)){
      m <- min(DF[i,],na.rm=TRUE)
      na_index <- which(is.na(DF[i,]))
      DF[i,na_index] <- m/2
    }
    
    # 全サンプルで 0 以下の列を削除する。
    nanidx <- c()
    for(i in 1:ncol(DF)){
      if(sum(DF[,i]) <= 0){
        nanidx <- c(nanidx,i)
      }
    }
    if(length(nanidx)>0){
      DF <- DF[,-nanidx]
      fn <- fn[-nanidx]
    }

    # 全 feature で 0 の行を除く
    nanidx <- c()
    for(i in 1:nrow(DF)){
      if(sum(DF[i,])<=0){
        nanidx <- c(nanidx,i)
      }
    }
    if(length(nanidx)>0){
      DF <- DF[-nanidx,]
      sn <- sn[-nanidx]
      gr <- gr[-nanidx]
    }
   
    
    # 前処理
    df <- reactive({
      if(input$normalization=='none'){
        tmp = as.data.frame(DF)
        colnames(tmp) = fn
        return(tmp)
      }else if(input$normalization=='composition'){
        tmp2 = as.data.frame(DF/rowSums(DF))
        colnames(tmp2) = fn
        return(tmp2)}
    })
    
    # 関数の読み込み
    source('plsrog.R')
    
    # pls-rog の実行
    res <- reactive({
      if(input$select=='naïve pls'){
        return(pls(df(), gr))}else{
        return(plsrog(df(), gr, input$kappa))
        }
    })

    axis1 <- reactive({
      score_X <- cbind(as.data.frame(res()[[1]]))
      if(!isTruthy(input$Axis1)){
        return(colnames(score_X)[1])
      }else{
        return(input$Axis1)
        }
      })

    axis2 <- reactive({
      score_X <- cbind(as.data.frame(res()[[1]]))
      if(!isTruthy(input$Axis2)){
        return(colnames(score_X)[2])
      }else{
        return(input$Axis2)
        }
      })

    output$score_X_plots <- renderPlotly({
      score_X <- cbind(as.data.frame(res()[[1]]),sn,gr)
      ng <- ncol(score_X)-2
      observe({
        updateSelectInput(session = session, 
                          inputId = "Axis1", 
                          choices = colnames(score_X)[1:ng],
                          selected = input$Axis1)
      })
      observe({
        updateSelectInput(session = session, 
                          inputId = "Axis2", 
                          choices = colnames(score_X)[1:ng],
                          selected = input$Axis2)
      })
      data = data.frame(x = score_X[,axis1()], y = score_X[,axis2()])
      plot <- ggplot(data, aes(x = x, y = y, color = gr)) + 
           geom_point(aes(color = gr, text = paste("Sample Name:",sn)),alpha = 1.0, size = 1) +
           labs(x = 'Axis1', y = 'Axis2', color = "diagnosis") +
           stat_ellipse(aes(color = gr), type = "t", level = input$level)
      plot
    })
    
    output$score_Y_plots <- renderPlotly({
      score_Y <- cbind(as.data.frame(res()[[2]]),gr)
      data = data.frame(x = score_Y[,axis1()], y = score_Y[,axis2()], color = gr)
      plot <- ggplot(data, aes(x = x, y = y)) + 
           geom_point(aes(color = gr),alpha = 1, size = 2) +
           labs(x = axis1(), y = axis2())
      plot
    })

    output$Wx_plots <- renderPlotly({
      Wx <- cbind(as.data.frame(res()[[3]]),fn)
      data = data.frame(x = Wx[,axis1()], y = Wx[,axis2()], fn = fn)
      plot <- ggplot(data, aes(x = x, y = y, fn = fn)) + 
           geom_point(alpha = 1, size = 1) +
           labs(x = axis1(), y = axis2())
      plot
    })
    
    output$R_plots <- renderPlotly({
      R <- res()[[5]]
      colnames(R) <- NULL
      R <- cbind(as.data.frame(R),fn)
      data = data.frame(x = R[,axis1()], y = R[,axis2()], fn = fn)
      plot <- ggplot(data, aes(x = x, y = y, fn = fn)) + 
           geom_point(alpha = 1, size = 1) +
           labs(x = axis1(), y = axis2())
      plot
    })
    
    output$Q_plots <- renderPlotly({
      Q <- res()[[7]]
      colnames(Q) <- NULL
      Q <- cbind(as.data.frame(Q),fn)
      data = data.frame(x = Q[,axis1()], y = Q[,axis2()], fn = fn)
      plot <- ggplot(data, aes(x = x, y = y, fn = fn)) + 
           geom_point(alpha = 1, size = 1) +
           labs(x = axis1(), y = axis2())
      plot
      })

    feature <- reactive({
      if(!isTruthy(input$feature)){
        return(fn[1])
      }else{
        return(input$feature)
        }
    })
    
    output$boxplot <- renderPlotly({
      observe({
        updateSelectInput(session = session, 
                          inputId = "feature", 
                          choices = fn,
                          selected = input$feature)
      })
      data = cbind(data.frame(y = df()[,feature()]),gr)
      plot <- ggplot(data, aes(y = y, x = gr)) + 
           geom_boxplot() +
           labs(y = feature(), x = gr)
      plot
      })
    }
  )
}

# Create Shiny app ----
shinyApp(ui, server)