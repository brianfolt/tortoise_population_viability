# Gopher tortoise PVA app

# This app can be deployed on GitHub pages using instructions from:
# https://medium.com/@rami.krispin/deploy-shiny-app-on-github-pages-b4cbd433bdc

# Run this code through R:
# library(shinylive)
# library(httpuv)
# shinylive::export(appdir = "tort_pva_app", destdir = "docs")
# httpuv::runStaticServer("docs/", port=8008) # to test
# # or visit Rstudio code editor to test manually
# # commit and push the changes
# # i have already set up the app to by live at https://brianfolt.github.io/tortoise_population_viability/

# Call a package
library(popbio)

# Specify two functions from Bruce Kendall's mpmtools package:
#    'subdiag()' function and 'make_Leslie_matrix()'
subdiag <- function(A, sx) {
  if (is.null(dim(A))) {
    A <- matrix(0, A, A)
  }
  stopifnot(is.matrix(A),
            length(sx) == 1 | length(sx) == (nrow(A) - 1),
            nrow(A) == ncol(A))
  n <- nrow(A)
  B <- matrix(0, n, n)
  B[-1, -n] <- diag(sx, n - 1, n - 1)
  return(B + A)
}

make_Leslie_matrix <- function(x, sx = NULL, mx = NULL, model = c("pre", "post")) {
  if(is.data.frame(x)) {
    stopifnot(c("x", "sx", "mx") %in% names(x))
    sx <- x$sx
    mx <- x$mx
    x <- x$x
  }
  if (length(x) == 1) x <- 0:x
  
  # Make prebreeding census model
  n <- length(x) - 1
  A <- subdiag(n, sx[2:n])
  A[n, n] <- sx[n + 1]
  A[1, ] <- sx[1] * mx[-1]
  
  if (model[1] == "post") A <- suppressWarnings(pre_to_post(sx[1], A))
  
  return(A)
}

##### Specify the UI -----
ui <- fluidPage(
  
    tags$head(HTML("<title>Tortoise PVA Tool</title>")),
  
    strong(h3("Tortoise PVA Tool")),
    #includeMarkdown("tort_pva_app/header.Rmd"),
    p(strong("Gopher tortoise"), " populations experience ", strong("varying demographic conditions"), "across the species' range in the southeastern United States. This page comprises a ", strong("flexible tool"), " that allows users to ", strong("simulate tortoise population dynamics"), " under varying conditions of demographic rates. Specifically, populations experience latitudinal variation in maturity age and fecundity, where more southern populations have faster somatic growth rates, reach sexual maturity at young ages, and lay larger clutches of eggs likely due to increased energy assimilation. To accommodate variation in life history, the user can ", strong("adjust mean estimates"), " of maturity age, fecundity (clutch size), and survival rates of different life history stages (nests, hatchlings, juveniles, adults). The juvenile stage includes all 1-year old animals up to the year prior to the maturity age. The software flexibly 'unwinds' the demographic rates to appropriate ages and projects the population using an age-based model. The model is a ",
    strong("female-only model"), " and assumes a ", strong("pre-breeding census.")),
    sidebarLayout(
      sidebarPanel(
        h4(p(strong("Demographic Rates"))),
        numericInput("ma", "Age of maturity:", 18, min = 5, max = 30),
        numericInput("bp", "Probability of laying eggs:", 0.97, 0, 1, step=0.01),
        numericInput("f", "Clutch size:", 6, 0, 25, step=1),
        numericInput("ns", "Probability of nest survival:", 0.35, 0, 1, step=0.01),
        numericInput("ve", "Probability of egg viability:", 0.85, 0, 1, step=0.01),
        numericInput("pf", "Probability of female:", 0.5, 0, 1, step=0.01),
        numericInput("s_h", "Probability of hatchling survival:", 0.13, 0, 1, step=0.01),
        numericInput("s_j", "Probability of juvenile survival:", 0.8, 0, 1, step=0.01),
        numericInput("s_a", "Probability of adult survival:", 0.98, 0, 1, step=0.01),
        numericInput("n", "Initial population size:", 50, 0, 10000, step=1),
        numericInput("pp", "Population persistence threshold (minimum adult females):", 3, 1, 100, step=1),
        br(),
        h4(p(strong("Simulation Inputs"))),
        numericInput("nyears", "Projection interval (years)", 50, 1, 100),
        numericInput("nreps", "Number of simulation replicates", 10, 5, 50),
        br(),
        actionButton("run", "Perform simulation",
                     class = "btn-success",
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
        br(), br(),
        p("It may take 20-30 seconds to display results.")
        
        ),
      
      mainPanel(
        br(),
        plotOutput("nplot"),
        br(),
        verbatimTextOutput("lambda_title"),
        verbatimTextOutput("lambda"),
        br(),
        verbatimTextOutput("gt_title"),
        verbatimTextOutput("gentime"),
        br(),
        verbatimTextOutput("ssd_title"),
        verbatimTextOutput("ssd"),
        br(),
        verbatimTextOutput("rv_title"),
        verbatimTextOutput("rv"),
        br(),
        verbatimTextOutput("elas_title"),
        verbatimTextOutput("elas"),
        br(),
        verbatimTextOutput("persist_title"),
        verbatimTextOutput("persist")
        ),
      
    ), #sidebar
    
    #includeMarkdown("footer.Rmd")
    p("See Folt et al. 2022 ",em("Global Ecology and Conservation"), " for a review of geographic variation of tortoise demographic rates. Note that paper used estimates of 'apparent survival', whereas this model assumes true survival.")
    
  )
  
server <- function(input, output, session) {
  
    dat <- eventReactive(input$run, {
      
      # Use 'req()' to require certain inputs to perform simulations
      req(input$ma, input$bp, input$f, input$ns, input$ve, input$pf,
          input$s_h, input$s_j, input$s_a,
          # Simulation inputs
          input$nyears, input$nreps, input$n, input$pp)
      
      withProgress(message = "Simulating population growth...", {
      
      library(statmod)
    
      ### 1) Specify User Inputs 
      
      # Population parameters
      ma <- input$ma
      bp <- input$bp
      f <- input$f
      ns <- input$ns
      ve <- input$ve
      pf <- input$pf
      s_h <- input$s_h
      s_j <- input$s_j
      s_a <- input$s_a
      
      # Create a demography schedule, with juvenile and mature age classes
      # The model has a pre-breeding census
      ma <- ma + 1
      demog_sched <- data.frame(x = 1:ma,
                                sx = rep(NA, length(1:ma)),
                                mx = rep(NA, length(1:ma)))
      
      # Specify productivity of juvenile and adult females
      demog_sched[1:(ma-1), "mx"] <- 0
      
      # Specify productivity of juvenile and adult females
      demog_sched[ma, "mx"] <- bp * f * ns * ve * pf *s_h

      # Specify juvenile and adult survival rates
      demog_sched[1:(ma-1), "sx"] <- s_j
      demog_sched[ma, "sx"] <- s_a
      
      
      # Construct a Leslie matrix from this demography schedule
      A <- make_Leslie_matrix(demog_sched)
      
      ## Demographic features
      
      # Calculate the asymptotic growth rate of the population governed by this 
      #   demography schedule:
      lam <- popbio::lambda(A)
      ssd <- popbio::stable.stage(A) # Stable stage (age) distribution
      names(ssd) <- paste0("Age_", 1:length(ssd))
      generation.time <- popbio::generation.time(A) # Generation time
      rv <- popbio::reproductive.value(A)  # Reproductive value
      names(rv) <- paste0("Age_", 1:length(rv))
      elas <- popbio::elasticity(A) # Elasticity values
      elas.values <- round(c(diag(elas[-1,]), elas[ma-1,ma-1], elas[1,ma-1]), 3)
      names(elas.values) <- c(paste0("S_", 1:(length(elas.values)-1)), "Fecundity")
      
      # Save the demographic variables together
      demo_vars <- list(lam, ssd, generation.time, rv, elas.values)
      
      
      ### Project and track population structure
      
      # Starting population size
      n <- input$n
      pp <- input$pp
      
      
      # Simulation parameters
      nyears <- input$nyears
      nreps <- input$nreps
      
      nstages <- dim(A)[1]
      
      # Matrices to save population size and structure
      N_tot <- matrix(0, nyears, nreps)
      N_stages <- array(0, c(nyears, nreps, nstages))
      # array to save population structure across years, reps, and ages
      
      ## Initial population size and structure
      # What is the stable-stage distribution (SSD) of the matrix
      ssd <- stable.stage(A)
      
      # How many individuals in the population?
      n <- 50
      
      # Spread individuals across the SSD
      n_ssd <- ssd * n
      
      # Use poisson draws to randomly populate numbers per stage for each simulation
      #   replicate
      n_i <- matrix(NA, nreps, dim(A)[1])
      for (i in 1:nreps){
        n_i[i,]  <- rpois(length(n_ssd), n_ssd)
      }
      # rowSums(n_i) # usually between 40-60 females to start
      
      ## Partition matrix for simulation
      x <- splitA(A)
      x_T <- x$T
      x_F <- x$F
      
      ## Run the simulation!
      for (j in 1:nreps){ # for each replicate
        n_i_j <- n_i[j,] # specify initial population size, randomly drawn above
        for (i in 1:nyears){ # for each year
          # If it's the first year, specify n_i_j; else,
          #   specify N_stage from previous year
          if(i == 1){
            n_sim <- n_i_j
          } else {
              n_sim <- N_stages[i-1,j,]
              }
          N_stages[i,j,] <- multiresultm(n_sim, x_T, x_F) # project pop
          N_tot[i,j] <- sum(N_stages[i,j,]) # save N
        } #year
      } #rep
      
      # Persistence probability in the final year
      N_adult_t <- N_stages[i,,ma-1]
      persistence <- ifelse(N_adult_t > pp-1, 1, 0)
      persist.prob <- sum(persistence)/nreps
      demo_vars <- c(demo_vars, persist.prob)
      
      # Save a vector of demographic rate titles
      demo_titles <- c(c("Population growth rate (lambda):"),
                       c("Generation time (years):"),
                       c("Stable age distribution:"),
                       c("Reproductive value of ages:"),
                       c("Elasticities for demographic rates:"),
                       c("Population persistence probability in the final year:")
                       )
      
      ### Save the data so they can be used
      list(N_tot, demo_vars, demo_titles)
      
  }) # end-withProgress
}) # end-dat
    
## Create an output object where we will save all the graphs and tables
##  from the Basic Tool
graphs <- reactiveValues(nplot = NULL, table1 = NULL)

## Summary graph 1 - abundance ------------
output$nplot <- renderPlot({
  
  # Call the data from the simulation run
  N_tot <- dat()[[1]]
  
  # Load ggplot2
  library(ggplot2)
  
  # Convert the matrix to a data frame
  df <- as.data.frame(N_tot)
  
  # Add a column for time (years)
  df$Time <- 1:50  # Assuming 50 years based on the matrix dimensions
  
  # Reshape the data to long format
  df_long <- tidyr::gather(df, key = "Population", value = "Count", -Time)
  
  # Create a ggplot with a logarithmic y-axis
  graphs$nplot <- ggplot(df_long, aes(x = Time, y = Count, color = Population)) +
    geom_line(linewidth = 1.5) +
    labs(x = "Time (years)", y = "Population size (females)") +
    ylim(0, max(N_tot)) +
    theme_minimal() +
    theme(legend.position = "none",  # Remove the legend
          axis.title = element_text(size = 18),
          axis.line = element_line(size = 1.25),
          axis.text = element_text(size = 14, color = 1))
         
  
  # Print the graph
  graphs$nplot

})

## Summary text - lambdaa
output$lambda_title <- renderText({
  
  # Recall the parameter name
  demo_titles <- dat()[[3]]
  demo_titles[1]

})


## Summary text - lambdaa
output$lambda <- renderText({
  
  # Recall some parameters
  demo_vars <- dat()[[2]]
  
  # Save lambda
  round(demo_vars[[1]], 3)

})

## Summary text - generation time
output$gt_title <- renderText({
  
  # Recall the parameter name
  demo_titles <- dat()[[3]]
  demo_titles[2]
  
})

## Summary text - generation time
output$gentime <- renderText({
  
  # Recall some parameters
  demo_vars <- dat()[[2]]
  
  # Save GT
  round(demo_vars[[3]], 1)
  
})

## Summary text - stable stage distribution
output$ssd_title <- renderText({
  
  # Recall the parameter name
  demo_titles <- dat()[[3]]
  demo_titles[3]
  
})

## Summary text - stable stage distribution
output$ssd <- renderText({
  
  # Recall some parameters
  demo_vars <- dat()[[2]]
  
  # SSD
  ssd <- demo_vars[[2]]
  
  # Save
  t(data.frame(names(ssd), round(ssd, 2)))
  
})


## Summary text - repro value
output$rv_title <- renderText({
  
  # Recall the parameter name
  demo_titles <- dat()[[3]]
  demo_titles[4]
  
})

## Summary text - reproductive value
output$rv <- renderText({
  
  # Recall some parameters
  demo_vars <- dat()[[2]]
  
  # RV
  rv <- demo_vars[[4]]
  
  # Save 
  t(data.frame(names(rv), round(rv, 2)))  
  
})


## Summary text - elasticities
output$elas_title <- renderText({
  
  # Recall the parameter name
  demo_titles <- dat()[[3]]
  demo_titles[5]
  
})

## Summary text - elasticities
output$elas <- renderText({
  
  # Recall some parameters
  demo_vars <- dat()[[2]]
  
  # Elasticities
  elas <- demo_vars[[5]]
  
  # Save text
  t(data.frame(names(elas), round(elas, 2)))

})

## Summary text - elasticities
output$persist_title <- renderText({
  
  # Recall the parameter name
  demo_titles <- dat()[[3]]
  demo_titles[6]
  
})

## Summary text - elasticities
output$persist <- renderText({
  
  # Recall some parameters
  demo_vars <- dat()[[2]]
  
  # Elasticities
  persist <- demo_vars[[6]]
  
  # Save text
  round(persist, 3)
  
})


}


# Create Shiny app ----
shinyApp(ui = ui, server = server)