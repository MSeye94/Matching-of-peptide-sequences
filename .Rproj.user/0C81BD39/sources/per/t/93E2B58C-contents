#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# library(shiny)
# library(shinyjs)
# library(shinythemes)
# library(shinyWidgets)
# library(shinycssloaders)
# library(shinyalert)
# 
# library(ggplot2)
# library(gridExtra)
# library(dplyr)
# library(np)
# library(MASS)
# library(readr)

source("lib/lib.R")






# Define UI for application that draws a histogram
ui <- fluidPage(
  
  includeCSS("www/myStyle.css"),
  
  # Application title
  titlePanel(title = div(style="font-weight:bold;color: #760001;font-size: 0.6em"
                         ,"MAtching of Peptide SEQuences (MAPSEQ)"),
             
             windowTitle = "Matching of peptide sequences"),
  
  # Sidebar with a slider input for number of bins 
  
  navbarPage(
    id = "analysisNavbar",
    title = NULL,
    inverse = TRUE,
    collapsible = TRUE,
    
    
    tabPanel("Peaks Grouping",
             
             sidebarLayout(
               
               
               
               sidebarPanel(
                 h4("Import files", icon("question-circle",
                                         class = "myIcoInfo")),
                 
                 fileInput(
                   inputId = "Input_massif_list",
                   label = "Upload massif-list files (accept .csv)",
                   accept = c(".csv"),
                   multiple = TRUE,
                   width = "100%"
                 ),
                 
                 hr(),
                 h4("Grouping parameters"),
                 
                 
                 numericInput(
                   inputId="rt_tolGroupingRefMap",
                   label="CE-time tolerence (Second)",
                   value=3*60,
                   min = 30,
                   max = 20*60,
                   step = 5
                 ),
                 
                 hr(),
                 fluidRow(column(6,
                                 actionButton(inputId="GroupingButton",
                                              label="Start",
                                              icon=icon("rocket"),
                                              class="btn-primary"),
                                 ),
                          column(6,
                                 
                                 shinySaveButton(id =  "saveFeatures_list",
                                                 label = "Save data",
                                                 title = "Save file as...",
                                                 filename = "Features-list",
                                                 filetype = list(text = "csv"),
                                                 viewtype = "icon",
                                                 icon = icon("save", lib = "glyphicon"),
                                                 class="btn-primary")
                                 )),
                 
                 
                
                 
                 
               ),
               mainPanel(
                 
                 fluidRow(
                   column(12,
                          div( class = "selectSampleViewer",
                               pickerInput(
                                 inputId = "levelSample",
                                 label = "Select/deselect sample(s) (to viewer):",
                                 choices = NULL,
                                 selected = NULL,
                                 options = pickerOptions(
                                   actionsBox = TRUE,
                                   size = 10,
                                   selectedTextFormat = "count",
                                   #showTick = TRUE,
                                   "max-options" = 6,
                                   style = "btn-primary"
                                 ),
                                 multiple = TRUE
                               ),
                               
                          ))),
                 
                 hr(),
                 div(
                   style = "position:relative",
                   shinycssloaders::withSpinner(
                     plotOutput("FateauresList_Plot",
                                width = "100%",
                                height = "600px",
                                click = "FateauresList_Plot_click",
                                dblclick = "FateauresList_Plot_dblclick",
                                brush = brushOpts(
                                  id = "FateauresList_Plot_brush",
                                  resetOnNew = TRUE
                                ),
                                hover = hoverOpts("FateauresList_Plot_hover", delay = 100, delayType = "debounce")
                                
                     ), type = 1, size = 0.8),
                   uiOutput("FateauresList_Plot_hover_info"),
                   hidden(
                     div(
                       id = "FateauresList_Plot_id",
                       class="pull-right",
                       style="display:inline-block;font-weight:bold;color: #760001;",
                       #disabled(
                       verbatimTextOutput("FateauresList_Plot_info")

                       #)
                     )

                   )
                   
                 )
                 
               )
             )
             
    ),
    
    
    
    
    tabPanel("Matching",
             
             sidebarLayout(
               
               
               
               sidebarPanel(
                 
                 h4("Import files", icon("question-circle",
                                         class = "myIcoInfo")),
                 
                 fileInput(
                   inputId = "refData",
                   label = "Upload reference file (accept .csv)",
                   accept = c(".csv"),
                   multiple = FALSE,
                   width = "100%"
                 ),
                 
                 fileInput(
                   inputId = "newSample",
                   label = "Upload file to align (accept .csv)",
                   accept = c(".csv"),
                   multiple = FALSE,
                   width = "100%"
                 ),
                 
                 
                 
                 hr(),
                 h4("Matching mz~ rt"),
                 
                 numericInput(
                   inputId = "rtTolerence",
                   label = "CE-time tolerence (Second)",
                   min = 0,
                   step = 5,
                   value = 720  
                 ),
                 
                 fluidRow(column(5,
                                 actionButton(inputId="match",
                                              label="Match",
                                              icon=icon("rocket"),
                                              class="btn-primary")
                                 ),
                          column(7,
                                 shinySaveButton(id =  "export_match",
                                                 label = "Export data",
                                                 title = "Save file as...",
                                                 filename = "Matched",
                                                 filetype = list(text = "csv"),
                                                 viewtype = "icon",
                                                 icon = icon("save", lib = "glyphicon"),
                                                 class="btn-primary")
                                 ))
                
                 

                
                 
                 
               ),
               
               
               # Show a plot of the generated distribution
               mainPanel(
                 
                 div(
                   style = "position:relative",
                   shinycssloaders::withSpinner(plotOutput("matchPlot",
                                                           width = "100%",
                                                           height = "600px",
                                                           dblclick = "matchPlot_dblclick",
                                                           brush = brushOpts(
                                                             id = "matchPlot_brush",
                                                             resetOnNew = TRUE
                                                           ),
                                                           hover = hoverOpts("matchPlot_hover", 
                                                                             delay = 100, delayType = "debounce")
                   ), 
                   type = 1, size = 0.8,
                   id = "Id_matchPlot"),
                   #uiOutput("matchPlot_hover_info")
                 )

               )
             )
             
             
             ),
    
    tabPanel("CE-time correction",
             
             
             sidebarLayout(
               
               
               
               sidebarPanel(
                 h4("Filter data with kernel density"),
                 
                 h6("CE-time file to align (Second):"),
                 fluidRow(
                   column(6,
                          align = "center",
                          numericInput(
                            inputId="rt_min_cut",
                            label="Min",
                            value=NULL,
                            min = 0)
                   ),
                   column(6,
                          align = "center",
                          numericInput(
                            inputId="rt_max_cut",
                            label="Max",
                            value=NULL,
                            min = 0)
                   )),
                     actionButton(inputId="cutOff",
                                  label="Cut-off",
                                  class="btn-primary",
                                  icon = icon("scissors")),
                 hr(),
                 
                 sliderInput("bandwidth",
                             "Bandwidth",
                             min = 1,
                             max = 1000,
                             value = 50),
                 
                 fluidRow(column(6,
                                 numericInput(
                                   inputId = "gridSize",
                                   label = "Grid size (n)",
                                   value = 500
                                 )
                 ),
                 column(6,
                        numericInput(
                          inputId = "minDensity",
                          label = "Min density",
                          min = 0,
                          max = 1,
                          step = 0.01,
                          value = 0  
                        )
                 )),
                 
                 hr(),
                 h4("Fit the model"),
                 
                 awesomeRadio(
                   inputId="ckertype",
                   label = "Kernel type",
                   inline=TRUE,
                   checkbox = TRUE,
                   choices=list("gaussian" = "gaussian",
                                "truncated gaussian" = "truncated gaussian",
                                "epanechnikov" = "epanechnikov",
                                "uniform" = "uniform")
                 ),
                 
                 sliderInput("bandwidth_fit_Model",
                             "Bandwidth",
                             min = 1,
                             max = 1000,
                             value = 100),
                 
                 
                 actionButton(inputId="fitModel",
                              label="Fit model",
                              icon=icon("rocket"),
                              class="btn-primary"),
                 
                 shinySaveButton(id =  "saveData",
                                 label = "Export aligned file",
                                 title = "Save file as...",
                                 filename = "ALigned",
                                 filetype = list(text = "csv"),
                                 viewtype = "icon",
                                 icon = icon("save", lib = "glyphicon"),
                                 class="btn-primary"),
                 
               ),
               
               
               # Show a plot of the generated distribution
               mainPanel(
                 tabsetPanel(
                   tabPanel("Density plot 2D and and Fit model ",
                            
                            shinycssloaders::withSpinner(plotOutput("distributionPlot",
                                                                    height = "500px"), 
                                                         type = 1, size = 0.8),
                            
                            shinycssloaders::withSpinner(plotOutput("distPlot",
                                                                    height = "265px"), 
                                                         type = 1, size = 0.8)
                            
                   ),
                   
                   tabPanel("Plot of aligned samples",
                            
                            shinycssloaders::withSpinner(
                              plotOutput("AdjustedSamplePlot",
                                         height = "500px"),
                              type = 1, size = 0.8)
                            
                            
                   )
                 )
                 
               )
             )
             
             )
    
)

)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Init packages~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  initPackages(session)
  
  data_rt_filter<-reactiveValues(data = NULL)
  data_rt<-reactiveValues(data = NULL)
  data_match_toExport<-reactiveValues(data = NULL)
  
  model.np<-reactiveValues(val = NULL)
  
  ref_data<-reactiveValues(data = NULL)
  new_sample<-reactiveValues(data = NULL)
  new_sample_adjusted<-reactiveValues(data = NULL)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Matching~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$refData
    }
    ,{
      
      observe({
        ref_data$data<-read.csv(input$refData$datapath, sep = req(input$delim_ref))
      })
      
      
      output$refDataView<-renderDT({
        
        datatable(
          req(ref_data$data)
        )
      })
      
      showModal(modalDialog(
        title = "View Data",
        
        h4("Options"),
        awesomeRadio(
          inputId="delim_ref",
          label = "Delimiter :",
          inline=FALSE,
          checkbox = TRUE,
          choices=list("Comma" = ",",
                       "Semicolon" = ";",
                       "Tab" = "\t")
        ),
        
        DTOutput("refDataView"),
        
       
          
        
        size = "l", easyClose = FALSE, fade=TRUE, footer = modalButton("Close (Esc)")
      ))
      
     
      

    
    })
  
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$newSample
    }
    ,{
      
      
      observe({
        new_sample$data<-read.csv(input$newSample$datapath, sep = req(input$delim_newSample))
      })
      
      
      output$newSampleDataView<-renderDT({
        
        datatable(
          req(new_sample$data)
        )
      })
      
      showModal(modalDialog(
        title = "View Data",
        
        h4("Options"),
        awesomeRadio(
          inputId="delim_newSample",
          label = "Delimiter :",
          inline=FALSE,
          checkbox = TRUE,
          choices=list("Comma" = ",",
                       "Semicolon" = ";",
                       "Tab" = "\t")
        ),
        
        DTOutput("newSampleDataView"),
        
        
        
        
        size = "l", easyClose = FALSE, fade=TRUE, footer = modalButton("Close (Esc)")
      ))
      
      
    })
  
  
  
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$match
    }
    ,{
      
      if(!is.null(ref_data$data) & !is.null(new_sample$data)) {
        
        res_match<-matchRtMz(x =  ref_data$data,
                             table = new_sample$data,
                             rt_tolerance = input$rtTolerence,
                             session = session)
        
        if(!is.null(res_match)){
          
          data_match<-res_match$MatchTable
          
          data_match$Match<-"No match"
          
          data_rt$data<-data_match[,c("rt1","rt2")]
          data_rt$data<-data_rt$data[!is.na(data_rt$data$rt2),]
          
          if(is.element(TRUE, !is.na(data_match$rt2)))
            data_match[!is.na(data_match$rt2),]$Match<-"Matched"
          
          data_match$Match<-factor(data_match$Match, levels = c("No match", "Matched"))
          
          data_match_toExport$data<-data_match
          
          rangesZoommatchViewer <- reactiveValues(x = NULL, y = NULL)
          
          # When a double-click happens, check if there's a brush on the plot.
          # If so, zoom to the brush bounds; if not, reset the zoom.
          observeEvent({input$matchPlot_dblclick},{
            brush <- input$matchPlot_brush
            if (!is.null(brush)) {
              rangesZoommatchViewer$x <- c(brush$xmin, brush$xmax)
              rangesZoommatchViewer$y <- c(brush$ymin, brush$ymax)
              
            } else {
              rangesZoommatchViewer$x <- NULL
              rangesZoommatchViewer$y <- NULL
            }
          })
          
          
          output$matchPlot<-renderPlot({
            
            par(fig = c(0,1,0,1),mar=c(4,4,5,0))
            plot(x = data_match$rt1,
                 y = data_match$mz1,
                 axes = T, pch=20, col= "gray", cex= 1,
                 xlab = paste("Ce-time(Second)"),
                 ylab = "M+H(Da)",
                 xlim = rangesZoommatchViewer$x,
                 ylim = rangesZoommatchViewer$y,
                 main = paste("REFRENCE MAP :","\n",
                              "Number of features matched :",res_match$numberMatch,"\n",
                              "Percentage matched :", res_match$percentageMatch,"%"),
                 col.main = "#760001",
                 
                 font.main = 2,
                 cex.lab = 1.2,
                 font.lab = 2,
                 cex.axis = 1.1,
                 font.axis = 2,
                 cex.main = 1.2)
            
            
            
            
            points(x = data_match[data_match$Match == "Matched",]$rt1,
                   y = data_match[data_match$Match== "Matched",]$mz1,
                   pch=20, col= "blue"#, cex=input$sizePointsMatch,
            )
            
            # axis(1,col="black",col.axis="black",
            #      font = 2, cex.lab = 1.3, cex.axis = 1.2, lwd = 2, line = 0)
            # # mtext("Ce-time(second)",side=1,line=3,col="black", font = 2, cex=1.2)
            # #
            # axis(2,col="black",col.axis="black",
            #      font = 2, cex.lab = 1.3, cex.axis = 1.2, lwd = 2, line = 0)
            # # mtext("M+H(Da)",side=2,line=3,col="black", font = 2, cex=1.2)
            
            legend("topleft", pch = 20, pt.cex = 2, inset = 0, text.font = 2.3,
                   legend = levels(data_match$Match),  bty = "n", xpd = NA, cex = 1.2,
                   col = c("gray", "blue"), horiz=FALSE)
            grid(nx = 14)
          })
        }
        
        
      }
        
      
    })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save table match ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$export_match
    }
    ,
    handlerExpr = {
      
      
      shinyFileSave(input,
                    id = "export_match",
                    roots = volumes,
                    session = session)
      
      path_Save_Files_origine<-parseSavePath(volumes, input$export_match)$datapath
      
      
      if(length(path_Save_Files_origine)>0){
        
        
        path_Save_Files<-strsplit(path_Save_Files_origine, split = "/")[[1]]
        directory_Save_Files<-paste0(path_Save_Files[-length(path_Save_Files)], collapse = "/")
        
        if(!is.null(data_match_toExport$data)){
          
          data_match_toSave<-data_match_toExport$data
          
          #data_match_toSave<-data_match_toSave[,c("ID", "mz2", "rt2", "Sequence", "sample", "Match")]
          colnames_data<-colnames(data_match_toSave)
          colnames(data_match_toSave)[which(colnames_data %in% c("mz2", "rt2"))]<-c("mz","rt")
          
          
          #sample_name<-unique(data_match_toSave$sample[!is.na(data_match_toSave$sample)])
          
         
          sample_name<-sub(input$newSample$name, pattern = ".csv",
                           replacement = "",
                           fixed = TRUE)
          
          write.table(data_match_toSave, 
                      file = paste0(c(directory_Save_Files, 
                                      paste(
                                      sample_name,
                                      path_Save_Files[length(path_Save_Files)],
                                      collapse = "-", sep = "-")),
                                    collapse = "/"),
                      sep = ",", row.names = FALSE)
        }
        
      }
      
    })
  
  

  observe({
    if(is.null(ref_data$data) || is.null(new_sample$data)) {
      
      output$matchPlot<-renderPlot({})
    }
  })
  
  
  
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CE-time correction~~~~~~~~~~~~~~~~~~~~~~~~~~#
 
  
  observe({
    if(!is.null(data_rt)){
      
      if(is.na(input$rt_min_cut) & is.na(input$rt_max_cut)){
        data_rt_filter$data<-data_rt$data
      }
      }
  })
  
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$cutOff
    }
    ,{
      if(!is.null(data_rt)){
        
        if(is.na(input$rt_min_cut) & !is.na(input$rt_max_cut)){
          
          data_rt_filter$data<-data_rt$data %>%
          dplyr::filter(rt2<=as.numeric(input$rt_max_cut))
        } else if(!is.na(input$rt_min_cut) & is.na(input$rt_max_cut)){
          
          data_rt_filter$data<-data_rt$data %>%
            dplyr::filter(rt2>=as.numeric(input$rt_min_cut))
        } else if(!is.na(input$rt_min_cut) & !is.na(input$rt_max_cut)){
          
          data_rt_filter$data<-data_rt$data %>%
            dplyr::filter(rt2>=as.numeric(input$rt_min_cut), rt2<=as.numeric(input$rt_max_cut))
        }
      }
      
      
    })

                          
    
  data_dens_filter<-reactive(
    if(!is.null(data_rt_filter$data)){
      dens<-kde2d(data_rt_filter$data$rt2, data_rt_filter$data$rt1, h = input$bandwidth, 
            n = input$gridSize)
      nx <- nrow(data_rt_filter$data)
      df <- expand.grid(x = dens$x, y = dens$y)
      df$density <- (as.vector(dens$z))
      df$density<-df$density/max(df$density)
      #df$group <- data$group[1]
      df$ndensity <- df$density / max(df$density, na.rm = TRUE)
      df$count <- nx * df$density
      df$n <- nx
      df$level <- 1
      df$piece <- 1
      df_filter<-df %>%
        filter(density>=req(input$minDensity))
      
      colnames(df_filter)[1:2]<-c("rt2","rt1")
      
      df_filter
    }
  )
  
  
  output$distributionPlot <- renderPlot({
    if(!is.null(data_rt_filter$data) & !is.null(data_dens_filter())){
      
      p1<- ggplot(data_rt_filter$data, aes(x=rt2, y=rt1)) +
        geom_point() +
        geom_density_2d(h = input$bandwidth, n = input$gridSize)+
        scale_x_continuous(n.breaks = 14)+
        scale_y_continuous(n.breaks = 14)+
        
        coord_cartesian(xlim =c(min(range(data_rt_filter$data$rt2)[1], range(data_rt_filter$data$rt1)[1]),
                                max(range(data_rt_filter$data$rt2)[2], range(data_rt_filter$data$rt1)[2])),
                        ylim = c(min(range(data_rt_filter$data$rt2)[1], range(data_rt_filter$data$rt1)[1]),
                                 max(range(data_rt_filter$data$rt2)[2], range(data_rt_filter$data$rt1)[2]))) +
        
        
        xlab("CE-time (sample to align)")+
        ylab("CE-time (reference sample)")+
        theme_ben()
      
      p2<-ggplot(data_dens_filter(), aes(x=rt2, y=rt1)) +
        geom_point(aes(color = density)) +
        
        scale_x_continuous(n.breaks = 14)+
        scale_y_continuous(n.breaks = 14)+
        
        coord_cartesian(xlim =c(min(range(data_dens_filter()$rt2)[1], range(data_dens_filter()$rt1)[1]),
                                max(range(data_dens_filter()$rt2)[2], range(data_dens_filter()$rt1)[2])),
                        ylim = c(min(range(data_dens_filter()$rt2)[1], range(data_dens_filter()$rt1)[1]),
                                 max(range(data_dens_filter()$rt2)[2], range(data_dens_filter()$rt1)[2])))+
        
        xlab("CE-time (sample to align)")+
        ylab("CE-time (reference sample)")+
        
        theme_ben()+
        theme(legend.title = element_text(size = rel(0.95), face = "bold.italic", hjust = 0.5),
              legend.text = element_text(size = rel(0.85), face = "bold.italic"))
      
      grid.arrange(p1, p2, nrow = 2)
      
      }
    
    })

   
      
      
        
        observeEvent(
          ignoreNULL = TRUE,
          eventExpr = {
            input$fitModel
          }
          ,{
              
              if(!is.null(data_dens_filter())){
                
                withProgress(message = 'Fit model...', value = 0, {
                  
                  incProgress(1/3, detail = paste(round(1*100/3,0),"%..."))
                  model.np$val<-npreg(rt1 ~ rt2,
                                      bws = req(input$bandwidth_fit_Model),
                                      bwtype = c("fixed","generalized_nn","adaptive_nn")[1],
                                      regtype = "ll",
                                      ckertype = input$ckertype,
                                      #bwmethod = "cv.aic",
                                      gradients = TRUE,
                                      data = isolate(data_dens_filter()))
                  
                  if(!is.null(model.np$val)){
                    
                    incProgress(1/3, detail = paste(round(2*100/3,0),"%..."))
                    output$distPlot <- renderPlot({
                      
                      #predict_data_model.np<-data.frame(rt1 = fitted(req(model.np$val)), rt2 = isolate(data_dens_filter()$rt2))
                      predict_data_model.np<-data.frame(rt1 = predict(model.np$val, 
                                                                      newdata = data.frame(rt2 = seq(min(range(data_rt_filter$data$rt2)[1], 
                                                                                                         range(data_rt_filter$data$rt1)[1]),
                                                                                                     max(range(data_rt_filter$data$rt2)[2], 
                                                                                                         range(data_rt_filter$data$rt1)[2])), by = 50)), 
                                                        rt2 = seq(min(range(data_rt_filter$data$rt2)[1], 
                                                                      range(data_rt_filter$data$rt1)[1]),
                                                                  max(range(data_rt_filter$data$rt2)[2], 
                                                                      range(data_rt_filter$data$rt1)[2])), by = 50)
                      
                      
                      median_line<-seq(min(range(data_rt_filter$data$rt2)[1], 
                                           range(data_rt_filter$data$rt1)[1]),
                                       max(range(data_rt_filter$data$rt2)[2], 
                                           range(data_rt_filter$data$rt1)[2]), by = 50)
                      
                      median_line_data<-data.frame(x =  median_line, y  =  median_line)
                      
                      
                      
                      ggplot(data_rt_filter$data, aes(x=rt2, y=rt1)) +
                        geom_point() +
                        
                        coord_cartesian(xlim =c(min(range(data_rt_filter$data$rt2)[1], range(data_rt_filter$data$rt1)[1]),
                                                max(range(data_rt_filter$data$rt2)[2], range(data_rt_filter$data$rt1)[2])),
                                        ylim = c(min(range(data_rt_filter$data$rt2)[1], range(data_rt_filter$data$rt1)[1]),
                                                 max(range(data_rt_filter$data$rt2)[2], range(data_rt_filter$data$rt1)[2])))+
                        
                        scale_x_continuous(n.breaks = 14)+
                        scale_y_continuous(n.breaks = 14)+
                        
                        # xlab("CE-time (sample to align)")+
                        # ylab("CE-time (reference sample)")+
                        
                        
                        
                        geom_line(data = predict_data_model.np, aes(x = rt2, y = rt1, color = "Model"), lwd = 1, size = 1.5) +
                        #geom_smooth(method ="loess", color="blue", fill="#69b3a2", se=TRUE, span = input$span) +
                        # geom_abline(yintercept = 0, lwd = 1, color = "green")+
                        geom_line(data = median_line_data, aes(x=x,y=x, color = "Median"), lwd = 1, size = 1.5)+
                        
                        labs(x = "CE-time (sample to align)",
                             y = "CE-time (reference sample)",
                             color = "Legend") +
                        
                        scale_color_manual(values = c("Model" = "red",
                                                      "Median" = "green"))+
                        
                        
                        theme_ben()+
                        theme(legend.title = element_text(size = rel(0.95), face = "bold.italic", hjust = 0.5),
                              legend.text = element_text(size = rel(0.85), face = "bold.italic"))
                      
                    })
                    
                    incProgress(1/3, detail = paste(round(3*100/3,0),"%..."))
                    ### Adjust the CE-time new sample
                    
                    new_sample_adjusted$data<-new_sample$data
                    
                    new_sample_adjusted$data$rt<-predict(model.np$val, newdata = data.frame(rt2 = new_sample$data$rt))
                    
                    new_sample_adjusted$data<-new_sample_adjusted$data %>%
                      dplyr::filter(rt>0)
                    
                   
                    
                    
                    ###Plot data
                    output$AdjustedSamplePlot <- renderPlot({
                      
                      
                      ## Before
                      dataCEtime_Before_RefSample<-data.frame(ref_data$data[,c("mz","rt")],
                                                              sample = "Reference sample",
                                                              Adjusted = "Before")
                      
                      
                      dataCEtime_BeforeNewSample<-data.frame(new_sample$data[,c("mz","rt")],
                                                             sample = "Sample to align",
                                                             Adjusted = "Before")
                      
                     
                      
                      dataCEtime_Before<-rbind(dataCEtime_Before_RefSample,
                                               dataCEtime_BeforeNewSample)
                      
                      
                      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
                      ## After
                      
                      dataCEtime_After_RefSample<-data.frame(ref_data$data[,c("mz","rt")],
                                                             sample = "Reference sample",
                                                             Adjusted = "After")
                      
                      
                      dataCEtime_AfterNewSample<-data.frame(new_sample_adjusted$data[,c("mz","rt")],
                                                             sample = "Sample to align",
                                                             Adjusted = "After")
                      
                      
                      
                      dataCEtime_After<-rbind(dataCEtime_After_RefSample,
                                              dataCEtime_AfterNewSample)
                      
                      ## Data Before and After
                      
                      dataCEtime<-rbind(dataCEtime_Before,
                                        dataCEtime_After)
                      dataCEtime$Adjusted<-factor(dataCEtime$Adjusted, levels = c("Before", "After"))
                      
                      
                      dataCEtime %>%
                        ggplot()+
                        aes(x = rt, y = mz, colour = sample) +
                        geom_point(shape = "circle", size = 0.5) +
                        ylab("Mass (M+H) (Da)") +
                        xlab("CE-time (Second)") +
                        scale_x_continuous(n.breaks = 14)+
                        # scale_y_continuous(n.breaks = 14)+
                        
                        #geom_hline(yintercept = 0, color = "black", linetype = "dashed", lwd = 0.8)+
                        scale_color_hue(direction = 1) +
                        facet_grid(vars(sample), vars(Adjusted)) +
                        theme_ben()+
                        theme(#legend.position = "none",
                              
                              plot.title = element_text(size = rel(1),
                                                        face = "bold",
                                                        color = "#760001",
                                                        margin = margin(0,0,5,0), hjust = 0.5),
                              legend.title = element_text(size = rel(0.95), face = "bold.italic", hjust = 0.5),
                              legend.text = element_text(size = rel(0.90), face = "bold"),
                              panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                          color = "blue", 
                                                          linewidth = 0.1,
                                                          linetype = "dashed"))
                      
                    })
                    
                    
                    
                    
                    
                    
                  }
                  
                  
                  
                  
                })
                
                
                
              } 
          
            
          })
        
        
        
        observeEvent(
          ignoreNULL = TRUE,
          eventExpr = {
            input$saveData
          }
          ,
          handlerExpr = {
            
            
            shinyFileSave(input,
                          id = "saveData",
                          roots = volumes,
                          session = session)
            
            path_Save_Files_origine<-parseSavePath(volumes, input$saveData)$datapath
            
            
            if(length(path_Save_Files_origine)>0){
              
              
              path_Save_Files<-strsplit(path_Save_Files_origine, split = "/")[[1]]
              directory_Save_Files<-paste0(path_Save_Files[-length(path_Save_Files)], collapse = "/")
              
              if(!is.null(new_sample_adjusted$data)){
                
                new_sample_adjusted_toSave<-new_sample_adjusted$data
                
                # colnames_data<-colnames(new_sample_adjusted_toSave)
                # 
                # colnames(new_sample_adjusted_toSave)[which(colnames_data %in% c("mz","rt"))]<-c("M+H", "CE-time")
                
                write.table(new_sample_adjusted_toSave, 
                          file = paste0(c(directory_Save_Files,
                                          paste(
                                            sub(input$newSample$name, pattern = ".csv",
                                                replacement = "",
                                                fixed = TRUE),
                                        path_Save_Files[length(path_Save_Files)],
                                        collapse = "", sep = "-")),
                                        collapse = "/"),
                          sep = ",", row.names = FALSE)
              }
              
            }
            
          })
        
        observe({
            if(is.null(model.np$val)){
              
              output$distPlot <- renderPlot({})
              output$AdjustedSamplePlot <- renderPlot({})
            }
          })
        
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~ Grouping Peaks~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
        
        rValue<-reactiveValues(massifList_list = NULL,
                               pheno_Data = NULL,
                               Features_list = NULL,
                               sample_name = NULL)
        
        observeEvent(
          ignoreNULL = TRUE,
          eventExpr = {
            input$Input_massif_list
          }
          ,{
            
            files <- input$Input_massif_list
            ext <- tools::file_ext(files$datapath)
            reqCols<-c('M.H','M.H.min','M.H.max','CE.time','CE.time.min','CE.time.max',
                       'integrated.intensity','intensity', 'sn', 'sample')
            
            reqColsUsed<-c('mz','mzmin','mzmax','rt','rtmin','rtmax',
                           'into','maxo', 'sn', 'sample')
            
            reqColsMess<-" 'M+H', 'M+H.min', 'M+H.max', 'CE-time', 'CE-time.min', 'CE-time.max',
                   'integrated-intensity', 'intensity', 'sn', 'sample'"
            
            if(length(ext)>1){
              if(all.equal(ext,rep("csv",length(ext)))!=TRUE){
                message("Please upload a csv file...!")
                # hide("progressWrapper") 
                # hide("analysisProgress")
                # shinyjs::html("analysisProgress",html = NULL)
                
                sendSweetAlert(
                  session = session,
                  title = "Warning !",
                  text = "Please upload a csv file !",
                  type = "warning",
                  width = "80%"
                )
              } else {
                # req(files)
                # validate(need(ext == "csv", "Please upload a csv file"))
                
                a<-lapply(input$Input_massif_list$datapath,read.csv)
                tryCatch({
                  a<-do.call("rbind", a)
                },
                error = function(e) {
                  message("Please upload a csv file...!")
                  # hide("progressWrapper") 
                  # hide("analysisProgress")
                  # shinyjs::html("analysisProgress",html = NULL)
                  
                  sendSweetAlert(
                    session = session,
                    title = "Warning !",
                    text = paste("Required columns :",reqColsMess,"not found in files.", sep = " "),
                    type = "warning",
                    width = "80%"
                  )
                })
                
                if(!all(reqCols %in% colnames(a))){
                  message("Please upload a csv file...!")
                  # hide("progressWrapper") 
                  # hide("analysisProgress")
                  # shinyjs::html("analysisProgress",html = NULL)
                  
                  sendSweetAlert(
                    session = session,
                    title = "Warning !",
                    text = paste("Required columns :",reqColsMess,"not found in files.", sep = " "),
                    type = "warning",
                    width = "80%"
                  )
                } else {
                  
                  a<-a[,reqCols]
                  colnames(a)<-reqColsUsed
                  #a<-split(a, f = a$sample)
                  #names(a)<-1:length(a)
                  rValue$massifList_list<-a
                  #hide("runPeakPicking")
                  
                  # Lunch the Reset button
                  enable("RetentionTimeFiltering_id")
                  hide("id_ResetCuttingRun")
                  shinyjs::show("id_validCuttingRun")
                  disable("SaveCuttingSample")
                  shinyjs::reset("rt_min_cut")
                  shinyjs::reset("rt_max_cut")
                  disable("GoToSelectRefranceSample")
                  disable("id_SelectRefrenceSamplePanel")
                  
                  
               
                }
                
              }
            } else {
              if(ext != "csv"){
                message("Please upload a csv file...!")
                # hide("progressWrapper") 
                # hide("analysisProgress")
                # shinyjs::html("analysisProgress",html = NULL)
                
                sendSweetAlert(
                  session = session,
                  title = "Warning !",
                  text = "Please upload a csv file !",
                  type = "warning",
                  width = "80%"
                )
              } else {
                # req(files)
                # validate(need(ext == "csv", "Please upload a csv file"))
                a<- lapply(input$Input_massif_list$datapath,read.csv)
                tryCatch({
                  a<-do.call("rbind", a)
                },
                error = function(e) {
                  message("Please upload a csv file...!")
                  # hide("progressWrapper") 
                  # hide("analysisProgress")
                  # shinyjs::html("analysisProgress",html = NULL)
                  
                  sendSweetAlert(
                    session = session,
                    title = "Warning !",
                    text = paste("Required columns :",reqColsMess,"not found in files.", sep = " "),
                    type = "warning",
                    width = "80%"
                  )
                })
                
                if(!all(reqCols %in% colnames(a))){
                  message("Please upload a csv file...!")
                  # hide("progressWrapper") 
                  # hide("analysisProgress")
                  # shinyjs::html("analysisProgress",html = NULL)
                  
                  sendSweetAlert(
                    session = session,
                    title = "Warning !",
                    text = paste("Required columns :",reqColsMess,"not found in files.", sep = " "),
                    type = "warning",
                    width = "80%"
                  )
                } else {
                  a<-a[,reqCols]
                  colnames(a)<-reqColsUsed
                  #a<-split(a, f = a$sample)
                  #names(a)<-1:length(a)
                  rValue$massifList_list<-a
                  #hide(id = "runPeakPicking")
                  
                  
                }
              }
            }
            
            
          }
        )
        
        
        
        observeEvent(
          ignoreNULL = TRUE,
          eventExpr = {
            input$GroupingButton
          }
          ,
          handlerExpr = {
            
            if(!is.null(rValue$massifList_list)){
              
              rValue$pheno_Data<-data.frame(
                Filename=unique(rValue$massifList_list$sample),
                Class="X")
              
              
              
              massif_list_samples_list<-split(rValue$massifList_list, f = rValue$massifList_list$sample)
              
             
              
              transform_massif_list<-function(massif_list){
                
                
                massif_list$sample<-1L
               
                
                #Coversion de la colonne sample
                massif_list$sample<-as.double(massif_list$sample)
                
                ## On tris le tableau par numero d'echantillon et on change le numerotation des lignes
                massif_list<-
                  massif_list[order(massif_list[,"sample"]),]
                
                rownames(massif_list)<-1:nrow(massif_list)
                
                massif_list_Matrix<-as.matrix.data.frame(massif_list)
                
                return(massif_list_Matrix)
              }
              
             
              
              





              ### convert peak table on Matrix

              Massif_list_Sample_Matrix_list<-lapply(massif_list_samples_list, transform_massif_list)
              
             
              rValue$Features_list<-processing_Grouping_xcms(peaks_sample_list = Massif_list_Sample_Matrix_list,
                                                      pheno_Data = rValue$pheno_Data

              )
              
              
              ### View sample 
              observe({
                if(!is.null(rValue$Features_list)){
                  
                  rValue$sample_name<-names(rValue$Features_list)
                  
                  # updateSelectInput(session = session,
                  #                   inputId = "levelSample",
                  #                   choices = rValues$sample_name)
                  
                  
                  updatePickerInput(session = session,
                                    inputId = "levelSample",
                                    choices = rValue$sample_name,
                                    selected = rValue$sample_name[1])
                  
                  
                }
              })
              
              
              rangesZoomFateauresList_Plot <- reactiveValues(x = NULL, y = NULL)
              
              # When a double-click happens, check if there's a brush on the plot.
              # If so, zoom to the brush bounds; if not, reset the zoom.
              observeEvent({input$FateauresList_Plot_dblclick},{
                brush <- input$FateauresList_Plot_brush
                if (!is.null(brush)) {
                  rangesZoomFateauresList_Plot$x <- c(brush$xmin, brush$xmax)
                  rangesZoomFateauresList_Plot$y <- c(brush$ymin, brush$ymax)
                  
                } else {
                  ## reset coord_cartesian
                  rangesZoomFateauresList_Plot$x <- NULL
                  rangesZoomFateauresList_Plot$y <- NULL
                }
              })
              
              output$FateauresList_Plot<-renderPlot({
                
                if(!is.null(rValue$Features_list)){
                  
                  if(length(rValue$Features_list)>=1 & !is.null(input$levelSample)){
                    
                    data_toPlot<-rValue$Features_list[input$levelSample]
                    
                    data_toPlot<-do.call("rbind", data_toPlot)
                    
                    data_toPlot<-as.data.frame(data_toPlot)
                    rownames(data_toPlot)<-1:nrow(data_toPlot)
                    
                    data_toPlot %>%
                      ggplot() +
                      aes(x = rt, y = mz, colour = log2(intensity)) +
                      geom_point() +
                      scale_color_viridis_c(option = "inferno", direction = -1) +
                      theme_gray() +
                      facet_grid(vars(sample), vars())  +
                      # facet_wrap(~req(sample) , dir="v", ncol = 1)  +
                      coord_cartesian(xlim = rangesZoomFateauresList_Plot$x, 
                                      ylim = rangesZoomFateauresList_Plot$y, expand = TRUE)+
                      ylab("Mass (M+H) (Da)") +
                      xlab("CE-time (Second)") +
                      
                      #ggtitle(paste("Nubmer of features :", number_of_features()))+
                      
                      scale_x_continuous(n.breaks = 14)+
                      scale_y_continuous(n.breaks = 14)+
                      #scale_y_continuous(breaks=seq(input$mzRange[1],input$mzRange[2],by=200))+
                      
                      
                      #facet_grid(vars(sample), vars())+
                      #option_graphe.2 +
                      theme_ben() +
                      labs(colour ="log2 Intensity")+
                      theme(plot.title = element_text(size = rel(1),
                                                      face = "bold",
                                                      color = "blue",
                                                      margin = margin(0,0,5,0),
                                                      hjust = 0.5),
                            plot.subtitle = element_text(hjust = 0.5),
                            panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                        color = "blue",
                                                        linewidth = 0.5,
                                                        linetype = "dashed")
                      )
                    }
                  
                  
                }
              })
              
              
              ### Position mouse
              shinyjs::show(id = "FateauresList_Plot_id")
              output$FateauresList_Plot_info <- renderText({
                paste0("Mouse position : CE-time = ", round(req(input$FateauresList_Plot_click$x),0)," (Second) ",
                       " M+H = ", round(req(input$FateauresList_Plot_click$y),4)," (Da)")
              })
              
              
              
              #################
              output$FateauresList_Plot_hover_info <- renderUI({
                hover <- req(input$FateauresList_Plot_hover)
                if(!is.null(rValue$Features_list)){
                  
                  if(length(rValue$Features_list)>=1 & !is.null(input$levelSample)){
                    
                    data_toPlot<-rValue$Features_list[input$levelSample]
                    
                    data_toPlot<-do.call("rbind", data_toPlot)
                    
                    data_toPlot<-as.data.frame(data_toPlot)
                    rownames(data_toPlot)<-1:nrow(data_toPlot)
                    
                    
                    
                    point <- nearPoints(req(data_toPlot), hover, threshold = 5, maxpoints = 1, addDist = TRUE)
                    if (nrow(point) == 0) return(NULL)
                    
                    colnames_data<-colnames(point)
                    colnames(point)[which(colnames_data %in% c("mz","rt"))]<-c("M.H","CE.time")
                    
                    # calculate point position INSIDE the image as percent of total dimensions
                    # from left (horizontal) and from top (vertical)
                    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
                    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
                    
                    # calculate distance from left and bottom side of the picture in pixels
                    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
                    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
                    
                    # create style property fot tooltip
                    # background color is set so tooltip is a bit transparent
                    # z-index is set so we are sure are tooltip will be on top
                    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                                    "left:", left_px + 2, "px; top:", top_px + 2, "px;")
                    
                    # actual tooltip created as wellPanel
                    div(
                      class = "well well-sm",
                      style = style,
                      p(HTML(paste0(
                        "<span class='bullText'> CE-time: </span>", round(point$CE.time,2), "(Second)<br/>",
                        "<span class='bullText'> M+H: </span>", round(point$M.H,4), "(Da)<br/>",
                        "<span class='bullText'> log2 Intensity: </span>", round(log2(point$intensity),2), "<br/>",
                        "<span class='bullText'> Intensity: </span>", round((point$intensity),0), "<br/>"
                      )))
                    )
                    
                    }
                  }
                
                
              })
             


              
            }
            
            
          })
        
        #~~~~~~~~~~~~~~~~~~~~~Save files ~~~~~~~~~~~~~~~~#
        
        observeEvent(
          ignoreNULL = TRUE,
          eventExpr = {
            input$saveFeatures_list
          }
          ,
          handlerExpr = {
            
            
            
            shinyFileSave(input,
                          id = "saveFeatures_list",
                          roots = volumes,
                          session = session)
            
            
            path_Save_Files_origine<-parseSavePath(volumes, input$saveFeatures_list)$datapath
            
            if(length(path_Save_Files_origine)>0){
              
              path_Save_Files<-strsplit(path_Save_Files_origine, split = "/")[[1]]
              directory_Save_Files<-paste0(path_Save_Files[-length(path_Save_Files)], collapse = "/")
              
              if(length(rValue$Features_list)>0 & !is.null(rValue$Features_list)){
                
                if(length(rValue$Features_list)==1){
                  write.table(rValue$Features_list[[1]], #file = path_Save_Files_origine,
                              file = paste0(c(directory_Save_Files,
                                              paste(names(rValue$Features_list)[1],
                                                    path_Save_Files[length(path_Save_Files)],
                                                    collapse = "", sep = "-")), collapse = "/"),
                              sep = ",", row.names = FALSE)
                  print("Save ok")
                } else{
                  
                  withProgress(message = 'Saving files', value = 0, {
                    n<-length(rValue$Features_list)
                  for(i in 1:length(rValue$Features_list)){
                    
                    incProgress(1/n, detail = paste(names(rValue$Features_list)[i],
                                                    round(i*100/n,0)," %..."))
                    
                    print(basename(paste0(c(directory_Save_Files,
                                            paste(#as.character(i),
                                              sub(names(rValue$Features_list)[i], pattern = "Massif-list",
                                              replacement = "",
                                              fixed = TRUE),
                                              path_Save_Files[length(path_Save_Files)],
                                              collapse = "", sep = "-")), collapse = "/")))
                    write.table(rValue$Features_list[[i]], file = paste0(c(directory_Save_Files,
                                                                           paste(#as.character(i),
                                                                             sub(names(rValue$Features_list)[i], pattern = "Massif-list",
                                                                                 replacement = "",
                                                                                 fixed = TRUE),
                                                                             path_Save_Files[length(path_Save_Files)],
                                                                             collapse = "", sep = "-")), collapse = "/"),
                                sep = ",", row.names = FALSE)
                    print("Save ok")
                  }
                    
                   })
                  
                }
                
                
              }
            }
            
            
          })
        
        
    observe({
      
      if(is.null(rValue$Features_list)){
        
        output$FateauresList_Plot<-renderPlot({})
      }
      
        
        
    })  

        
    
}

# Run the application 
shinyApp(ui = ui, server = server)
