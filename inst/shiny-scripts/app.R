library(shiny)
library(shinyalert)

# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel(tags$h1(tags$b("mixMVPLN:"),"Mixtures of Matrix Variate
                     Poisson-log Normal With Variational-EM Framework
                     for Clustering Three-way Count Data")),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      tags$p("Description: This is a Shiny App that is part of the mixMVPLN R package
             (Silva et al., 2022). Most of the functions available via the package are made
             available with Shiny App. The mixMVPLN is an R package for performing
             clustering using mixtures of matrix variate Poisson-log normal (mixMVPLN)
             distribution provided a three-way count dataset. The observations of the
             dataset will be clustered into subgroups. The app permits to model select
             using Bayesian information criterion (BIC), Integrated Complete Likelihood (ICL)
             criterion, and Akaike Information Criterion (AIC) values. Results will
             be available as dotted line plots of information criteria value versus
             cluster size, heatmaps of input datset with cluster membership information,
             alluvial plot of observations by cluster, and barplot of posterior
             probabilities."),

      # br() element to introduce extra vertical spacing ----
      br(),
      br(),

      # input
      tags$b("Instructions: Below, upload a three-way count dataset
              in rds format, and enter or select values to perform
              cluster analysis. Default values are shown. Then press
              'Run'. Navigate through the different tabs to the right
              to explore the results."),

      # br() element to introduce extra vertical spacing ----
      br(),
      br(),
      br(),
      # input
      shinyalert::useShinyalert(),  # Set up shinyalert
      uiOutput("tab2"),
      actionButton(inputId = "data1",
                   label = "Dataset 1 Details"),
      fileInput(inputId = "file1",
                label = "Dataset: Select a three-way count dataset to analyze.
                File should be in single R object (.rds) format containing
                three-way data of units (n), defined by occasions/layers (r) and
                the variables/responses (p) save as a list in R. You may
                download an example dataset above and explore first.",
                accept = c(".rds")),
      textInput(inputId = "ngmin",
                label = "ngmin: Enter the minimum value of components or clusters
                for clustering. This should be a positive integer.", "1"),
      textInput(inputId = "ngmax",
                label = "ngmax: Enter the maximum value of components or clusters.
                for clustering. This should be a positive integer, bigger
                than ngmin and less than or equal to number of total units (n)
                in the dataset.", "2"),
      selectInput(inputId = 'typeinitMethod',
                  label = 'initMethod: Select the initialization method.',
                  choices = c("kmeans",
                              "random",
                              "medoids",
                              "clara",
                              "fanny")),
      textInput(inputId = "nInitIterations",
                label = "nInitIterations: Enter the number of initial iterations.
                This should be a positive integer.", "1"),
      selectInput(inputId = 'typenormalize',
                  label = 'Normalization: Select whether to perform normalization
                  or not. Currently, normalization is performed by
                  calculating normalization factors via TMM method of edgeR
                  package (Robinson, et al., 2010). The option Yes is recommended
                  for raw RNA sequencing count data.',
                  choices = c("'Yes' ",
                              "'No' ")),

      # br() element to introduce extra vertical spacing ----
      br(),

      # actionButton
      actionButton(inputId = "button2",
                   label = "Run"),

      # br() element to introduce extra vertical spacing -
      br(),

    ), # End of side pannel


    # Main panel for displaying outputs
    mainPanel(

      # Output: Tabet
      tabsetPanel(type = "tabs",
                  tabPanel("Input Summary",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Summary of Count Dataset:"),
                           br(),
                           verbatimTextOutput('textOut')),
                  tabPanel("Cluster Results",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Summary of Clustering Results:"),
                           br(),
                           verbatimTextOutput('clustering')),
                  tabPanel("Information Criteria Plot",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Model Selection Results:"),
                           br(),
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput('BICvalues'), plotOutput('ICLvalues')),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput('AIC3values'), plotOutput('AICvalues')),
                           )),
                  tabPanel("Alluvial Plot",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Alluvial Plot Showing Observation Memberships by Information Criteria for Input Dataset:"),
                           h5("Note, below the x-axis values are in the order of BIC, ICL, AIC, AIC3.
                              Colors are assigned based on cluster membership of model selected via BIC."),
                           br(),
                           fluidRow(
                             splitLayout(cellWidths = c("100%"), plotOutput("alluvialPlot")),
                             h5("Note, below the x-axis values are in the order of ICL, BIC, AIC, AIC3.
                              Colors are assigned based on cluster membership of model selected via ICL."),
                             splitLayout(cellWidths = c("100%"), plotOutput("alluvialPlot2")),
                             h5("Note, below the x-axis values are in the order of AIC3, ICL, BIC, AIC
                              Colors are assigned based on cluster membership of model selected via AIC3."),
                             splitLayout(cellWidths = c("100%"), plotOutput("alluvialPlot3")),
                             h5("Note, below the x-axis values are in the order of AIC, AIC3, ICL, BIC
                              Colors are assigned based on cluster membership of model selected via AIC."),
                             splitLayout(cellWidths = c("100%"), plotOutput("alluvialPlot4")),
                           )),
                  tabPanel("Barplot",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Barplot of Posterior Probabilities with Cluster Memberships:"),
                           h5("Note, the plots are in the order of models selected by: BIC (top, left), ICL (top, right) and AIC (bottom, left), AIC3 (bottom, right)."),
                           br(),
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("barPlotBIC"), plotOutput('barPlotICL')),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("barPlotAIC3"), plotOutput('barPlotAIC'))
                           ))

      )
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {

  # Reactive expression to generate the requested distribution ----
  # This is called whenever the inputs change. The output functions
  # defined below then use the value computed from this expression


  # Step I: save input csv as a reactive
  matrixInput <- reactive({
    if (! is.null(input$file1))
      as.list(readRDS(input$file1$datapath))
  })


  startclustering <- eventReactive(eventExpr = input$button2, {
    withProgress(message = 'Clustering', value = 0, {
      # Number of times we'll go through the loop

      mixMVPLN::mvplnVGAclus(
        dataset = matrixInput(),
        membership = "none",
        gmin = as.numeric(input$ngmin),
        gmax = as.numeric(input$ngmax),
        initMethod = as.character(input$typeinitMethod),
        nInitIterations = as.numeric(input$nInitIterations),
        normalize = "Yes")
    })
  })

  # Textoutput
  output$textOut <- renderText({
    if (! is.null(startclustering))

      a1 <- paste("Number of occasions/layers (r):", nrow(startclustering()$dataset[[1]]), "\n")

      a2 <- paste("Number of variables/responses (p):", ncol(startclustering()$dataset[[1]]), "\n")

      a3 <- paste("Number of units (n):", length(startclustering()$dataset), "\n")

      paste(a1, a2, a3, sep = "\n")
 })


  # Step II: clustering
  output$clustering <- renderText({
    if (! is.null(startclustering))

    aa <- paste("BIC model selected is:", startclustering()$BICAll$BICmodelselected, "\n")

    bb <- paste("ICL model selected is:", startclustering()$ICLAll$ICLmodelselected, "\n")

    cc <- paste("AIC model selected is:", startclustering()$AICAll$AICmodelselected, "\n")

    dd <- paste("AIC3 model selected is:", startclustering()$AIC3All$AIC3modelselected, "\n")
    paste(aa, bb, cc, dd, sep = "\n")
  })

  # Step III: visualize



  # plot ICL value
  output$ICLvalues <- renderPlot({
    if (! is.null(startclustering))
      if (length(startclustering()$loglikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$ICLAll$allICLvalues), type = "p",
               xlab = "G", ylab = "ICL value",
               main = paste("G vs ICL value"), xaxt="n")
          axis(1, at = seq(as.numeric(0), as.numeric(input$ngmax), by = 1))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$ICLAll$allICLvalues),
               type = "p", xlab = "G", ylab = "ICL value",
               main = paste("G vs ICL value"), xaxt="n")
          axis(1, at = seq(as.numeric(0), as.numeric(input$ngmax), by = 1))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$ICLAll$allICLvalues, type = "l",
             lty = 2, xlab = "G", ylab = "ICL value",
             main = paste("G vs ICL value"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })


  output$BICvalues <- renderPlot({
    if (! is.null(startclustering))
      if (length(startclustering()$loglikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$BICAll$allICLvalues), type = "p",
               xlab = "G", ylab = "BIC value",
               main = paste("G vs BIC value"), xaxt="n")
          axis(1, at = seq(as.numeric(0), as.numeric(input$ngmax), by = 1))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$BICAll$allBICvalues),
               type = "p", xlab = "G", ylab = "BIC value",
               main = paste("G vs BIC value"), xaxt="n")
          axis(1, at = seq(as.numeric(0), as.numeric(input$ngmax), by = 1))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$BICAll$allBICvalues, type = "l",
             lty = 2, xlab = "G", ylab = "BIC value",
             main = paste("G vs BIC value"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })


  output$AICvalues <- renderPlot({
    if (! is.null(startclustering))
      if (length(startclustering()$loglikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$AICAll$allAICvalues), type = "p",
               xlab = "G", ylab = "AIC value",
               main = paste("G vs AIC value"), xaxt="n")
          axis(1, at = seq(as.numeric(0), as.numeric(input$ngmax), by = 1))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$AICAll$allAICvalues),
               type = "p", xlab = "G", ylab = "AIC value",
               main = paste("G vs AIC value"), xaxt="n")
          axis(1, at = seq(as.numeric(0), as.numeric(input$ngmax), by = 1))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$AICAll$allAICvalues, type = "l",
             lty = 2, xlab = "G", ylab = "AIC value",
             main = paste("G vs AIC value"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })


  output$AIC3values <- renderPlot({
    if (! is.null(startclustering))
      if (length(startclustering()$loglikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$AIC3All$allAIC3values), type = "p",
               xlab = "G", ylab = "AIC3 value",
               main = paste("G vs AIC3 value"), xaxt="n")
          axis(1, at = seq(as.numeric(0), as.numeric(input$ngmax), by = 1))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$AIC3All$allAIC3values),
               type = "p", xlab = "G", ylab = "AIC3 value",
               main = paste("G vs AIC3 value"), xaxt="n")
          axis(1, at = seq(as.numeric(0), as.numeric(input$ngmax), by = 1))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$AIC3All$allAIC3values, type = "l",
             lty = 2, xlab = "G", ylab = "AIC3 value",
             main = paste("G vs AIC3 value"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })


  # plot bar - ICL
  barPlottingICL <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      if ((as.numeric(input$ngmax) - as.numeric(input$ngmin) + 1) == 1) {
        MPLNClust::mplnVisualizeBar(
          vectorObservations = 1:length(matrixInput()),
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$ICLAll$ICLmodelselectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$ICLAll$ICLmodelselected)
        MPLNClust::mplnVisualizeBar(
          vectorObservations = 1:length(matrixInput()),
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$ICLAll$ICLmodelselectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - ICL
  output$barPlotICL <- renderPlot({
    barPlottingICL()
  })


  # plot bar - BIC
  barPlottingBIC <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      if ((as.numeric(input$ngmax) - as.numeric(input$ngmin) + 1) == 1) {
        MPLNClust::mplnVisualizeBar(
          vectorObservations = 1:length(matrixInput()),
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$BICAll$BICmodelselectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$BICAll$BICmodelselected)
        MPLNClust::mplnVisualizeBar(
          vectorObservations = 1:length(matrixInput()),
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$BICAll$BICmodelselectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - BIC
  output$barPlotBIC <- renderPlot({
    barPlottingBIC()
  })



  # plot bar - AIC3
  barPlottingAIC3 <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      if ((as.numeric(input$ngmax) - as.numeric(input$ngmin) + 1) == 1) {
        MPLNClust::mplnVisualizeBar(
          vectorObservations = 1:length(matrixInput()),
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AIC3All$AIC3modelselectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$AIC3All$AIC3modelselected)
        MPLNClust::mplnVisualizeBar(
          vectorObservations = 1:length(matrixInput()),
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AIC3All$AIC3modelselectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - AIC3
  output$barPlotAIC3 <- renderPlot({
    barPlottingAIC3()
  })


  # plot bar - AIC
  barPlottingAIC <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      if ((as.numeric(input$ngmax) - as.numeric(input$ngmin) + 1) == 1) {
        MPLNClust::mplnVisualizeBar(
          vectorObservations = 1:length(matrixInput()),
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AICAll$AICmodelselectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$AICAll$AICmodelselected)
        MPLNClust::mplnVisualizeBar(
          vectorObservations = 1:length(matrixInput()),
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AICAll$AICmodelselectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - AIC
  output$barPlotAIC <- renderPlot({
    barPlottingAIC()
  })








  # Alluvial plot
  alluvialPlotting <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      MPLNClust::mplnVisualizeAlluvial(nObservations = length(matrixInput()),
                            firstGrouping =
                              as.numeric(startclustering()$BICAll$BICmodelselectedLabels),
                            secondGrouping =
                              as.numeric(startclustering()$ICLAll$ICLmodelselectedLabels),
                            thirdGrouping =
                              as.numeric(startclustering()$AICAll$AICmodelselectedLabels),
                            fourthGrouping =
                              as.numeric(startclustering()$AIC3All$AIC3modelselectedLabels),
                            fileName = 'alluvial',
                            printPlot = FALSE)
  })

  alluvialPlotting2 <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      MPLNClust::mplnVisualizeAlluvial(nObservations = length(matrixInput()),
                            firstGrouping =
                              as.numeric(startclustering()$ICLAll$ICLmodelselectedLabels),
                            secondGrouping =
                              as.numeric(startclustering()$BICAll$BICmodelselectedLabels),
                            thirdGrouping =
                              as.numeric(startclustering()$AICAll$AICmodelselectedLabels),
                            fourthGrouping =
                              as.numeric(startclustering()$AIC3All$AIC3modelselectedLabels),
                            fileName = 'alluvial',
                            printPlot = FALSE)
  })

  alluvialPlotting3 <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      MPLNClust::mplnVisualizeAlluvial(nObservations = length(matrixInput()),
                            firstGrouping =
                              as.numeric(startclustering()$AIC3All$AIC3modelselectedLabels),
                            secondGrouping =
                              as.numeric(startclustering()$ICLAll$ICLmodelselectedLabels),
                            thirdGrouping =
                              as.numeric(startclustering()$BICAll$BICmodelselectedLabels),
                            fourthGrouping =
                              as.numeric(startclustering()$AICAll$AICmodelselectedLabels),
                            fileName = 'alluvial',
                            printPlot = FALSE)
  })


  alluvialPlotting4 <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      MPLNClust::mplnVisualizeAlluvial(nObservations = length(matrixInput()),
                            firstGrouping =
                              as.numeric(startclustering()$AICAll$AICmodelselectedLabels),
                            secondGrouping =
                              as.numeric(startclustering()$AIC3All$AIC3modelselectedLabels),
                            thirdGrouping =
                              as.numeric(startclustering()$ICLAll$ICLmodelselectedLabels),
                            fourthGrouping =
                              as.numeric(startclustering()$BICAll$BICmodelselectedLabels),
                            fileName = 'alluvial',
                            printPlot = FALSE)
  })


  # Alluvial Plot
  output$alluvialPlot <- renderPlot({
    alluvialPlotting()
  })

  output$alluvialPlot2 <- renderPlot({
    alluvialPlotting2()
  })

  output$alluvialPlot3 <- renderPlot({
    alluvialPlotting3()
  })

  output$alluvialPlot4 <- renderPlot({
    alluvialPlotting4()
  })




  # URLs for downloading data
  url2 <- a("Download Example Dataset 1", href="https://drive.google.com/file/d/1Hbk6hH9QQ-3yR5YkgHCIUtlWkGWS5ipS/view?usp=sharing")
  output$tab2 <- renderUI({
    tagList("Dataset:", url2)
  })

  observeEvent(input$data1, {
    # Show a modal when the button is pressed
    shinyalert(title = "Example Dataset 1",
               text = "This is a simulated dataset generated from mixtures of multivariate Poisson log-normal
               distributions with G = 2 components. It has a size of n = 1000 observations along rows and d = 6
               samples along columns. Data was generated January, 2022. To save the file, click on link, then click 'Download' from the top right side.
               Citation: Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2019). A multivariate Poisson-log normal
               mixture model for clustering transcriptome sequencing data. BMC Bioinformatics. 2019;20(1):394. URL https://pubmed.ncbi.nlm.nih.gov/31311497/",
               type = "info")
  })


}

# Create Shiny app ----
shinyApp(ui, server)
