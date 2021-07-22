#### CenGenAPP 2018-2020
#### Report bugs to: gabrielsantperebaro@gmail.com


require(shiny)
library(shinyjs)
library(shinythemes)
library(shinybusy)

library(DT)
library(dplyr)
library(ggplot2)

library(Matrix)
library(expss)


options(shiny.maxRequestSize = 400 * 1024 ^ 2)
source("Functions.R")

options(scipen = 0)
options(digits = 2)

utr <-
  c(
    "WBGene00023498",
    "WBGene00023497",
    "WBGene00004397",
    "WBGene00006843",
    "WBGene00004010",
    "WBGene00006789",
    "WBGene00001135",
    "WBGene00001079",
    "WBGene00001135",
    "WBGene00006783",
    "WBGene00000501",
    "WBGene00006788",
    "WBGene00001555",
    "WBGene00206533",
    "WBGene00011964",
    "WBGene00018172",
    "WBGene00016259",
    "WBGene00023407"
  )

msg <- filter(gene_list, gene_id %in% utr)$gene_name %>% paste(., collapse = ", ")

## UI ----
ui <- fluidPage(
  # tags$head(includeHTML(("google-analytics-script2.html"))),
  
  theme = "Theme.min.css",
  tags$head(tags$style(
    HTML(".shiny-output-error-validation {color: red;}")
  )),
  tags$style(HTML("
    .tabbable > .nav > li > a[data-value='Enriched Genes by cell type'] {background-color: grey;   color:white}
    .tabbable > .nav > li > a[data-value='Find Differential Expression between Cell Types'] {background-color: grey;  color:white}
    .tabbable > .nav > li > a[data-value='Single cell plot'] {background-color: grey; color:white}
  ")),
  
  # App title ----
  titlePanel(
    "Discovery and analysis of the C. elegans Neuronal Gene Expression Network -- CeNGEN"
  ),
  tags$a(href="http://www.cengen.org/single-cell-rna-seq/", "CeNGENApp Help and Documentation"),
  hr(),
  add_busy_spinner(
    spin = "double-bounce",
    color = "orange",
    timeout = 100,
    position = c("top-right"),
    onstart = TRUE,
    margins = c(10, 10),
    height = "50px",
    width = "50px"
  ),
  
  # Main panel showing plots ----
  tabsetPanel(
    type = "pills",
    id = "tabs",
    ### Cell types Panel ----
    tabPanel(
      "Gene expression by cell type",
      fluidPage(
        hr(),
        h6(
          "Find all genes expressed in a given cell type, or all cell types expressing a given gene (or group of genes)."
        ),
        h6("Select one of four thresholds for expression:"),
        h6("1 (least stringent) to 4 (most stringent) or select unfiltered data"),
        h6("Choose All Cells Unfiltered to query the entire unfiltered dataset, including non-neuronal cells"),
        hr(),
        fluidRow(
          column(1),
          column(
            4,
            selectInput(
              inputId = "Tcell_name",
              label = "Select cell type",
              choices = colnames(L4.all.TPM.raw_th)[1:169],
              selected = "ADA"
            ),
            selectInput(
              inputId = "Tcell_cut",
              label = "Select threshold",
              choices = c(1:4, "Unfiltered", "All Cells Unfiltered"),
              selected = 2
            ),
            actionButton("TCell", "Expressed genes", icon = icon("hand-o-right"))
            
          ),
          #column(width = 1, offset = 0, style='padding:5px;'),
          column(
            4,
            textInput(
              inputId = "Tgene_name",
              label = "Type gene name",
              value = "zig-4"
            ),
            selectInput(
              inputId = "Tgene_cut",
              label = "Select threshold",
              choices = c(1:4, "Unfiltered", "All Cells Unfiltered"),
              selected = 2
            ),
            actionButton("TGene", "Which cell types", icon = icon("hand-o-right"))
          ),
          column(
            2,
            offset = 0,
            style = 'padding:5px;',
            textAreaInput(
              inputId = "Tgene_name_batch",
              label = "Query multiple genes for download",
              value = "flp-1\nflp-2,flp-3,WBGene00001447\nWBGene00001448\nflp-6\nflp-7\nflp-8\nflp-9\nflp-10\nflp-11\nflp-12\nflp-13\nflp-14\nflp-15\nflp-16\nflp-17\nflp-18\nflp-19\nflp-20\nflp-21\nflp-22\nflp-23\nflp-24\nflp-25\nflp-26\nflp-27\nflp-28\nflp-32\nflp-33\nflp-34",
              width = "500px",
              height = "100px"
            ),
            selectInput(
              inputId = "Tgene_cut_batch",
              label = "Select threshold",
              choices = c(1:4, "Unfiltered", "All Cells Unfiltered" ),
              selected = 2
            ),
            downloadButton("TGeneBatch", "Download batch"),
            span(textOutput("textb"), style =
                   "color:red")
            
          )
        ),
        br(),
        h6( paste0("WARNING: Expression values for ",msg," are unreliable as they have been overexpressed to generate transgenic strains."), style="color:orange"),
        fluidRow(
          column(1),
          column(
            4,
            br(),
            span(textOutput("Error1"), style ="color:red"),
            DT::dataTableOutput("Tcell_name_table"),
            br(),
            uiOutput("get_download_gene")
          ),
          #column(width = 1, offset = 0, style='padding:5px;'),
          column(
            4,
            br(),
            DT::dataTableOutput("Tgene_name_table"),
            br(),
            #downloadButton('downloadCell', "Download table"),
            uiOutput("get_download_cell"),
            span(textOutput("text1"), style =
                   "color:red")
          )
        )
      )
    ),
    
    ### Find markers ----
    tabPanel(
      "Find markers based on percentage of expression",
      fluidPage(
        hr(),
        h6(
          "Find genes expressed in one group of cell types and not in another group based on the percentages of cells expressing the gene."
        ),
        h6( paste0("WARNING: Expression values for ",msg," are unreliable as they have been overexpressed to generate transgenic strains."), style="color:orange"),
        br(),
        fluidRow(
          #column(1),
          column(
            4,
            textInput(
              inputId = "String1",
              label = "Group 1",
              value = "AWC_ON,AWC_OFF"
            ),
            textInput(
              inputId = "Expressed",
              label = "Minimum percentage of cells expressing the gene",
              value = 65
            ),
            actionButton("Filter", "Run query", icon = icon("hand-o-right"))
            
          ),
          #column(width = 1, offset = 0, style='padding:5px;'),
          column(
            4,
            textInput(
              inputId = "String2",
              label = "Group 2",
              value = "SMD,SIB"
            ),
            textInput(
              inputId = "NonExpressed",
              label = "Maximum percentage of cells expressing the gene",
              value = 2
            ),
            downloadButton('downloadQuery', "Download table")
          )
        ),
        fluidRow(
          #column(1),
          column(4,
                 br(), DT::dataTableOutput("YesExpressed")),
          column(4,
                 br(),
                 DT::dataTableOutput("NoExpressed")),
          
          column(4,
                 br(),
                 DT::dataTableOutput("Result")
                 #span(textOutput("text3"), style="color:red"))
          )
          
        )
      )
    ),
    ### Enriched types Panel ----
    tabPanel(
      "Enriched Genes by cell type",
      fluidPage(
        hr(),
        h6("This functionality is not available in the Lite App.")
      )
    ),
    ### DEX panel ----
    tabPanel(
      "Find Differential Expression between Cell Types",
      fluidPage(
        hr(),
        h6("This functionality is not available in the Lite App.")
      )
    ),
    ### Single cell panel ----
    tabPanel(
      "Single cell plot",
      fluidPage(
        hr(),
        h6("This functionality is not available in the Lite App.")
      ),
    ),
    
    ### Heatmap ----
    tabPanel(
      "Heatmaps of gene expression",
      fluidPage(
        hr(),
        h6(
          "Display a heatmap showing relative expression and proportion of cells expressing a gene or group of genes across all neurons. This function uses data from threshold 2. Color shows relative scaled expression for each gene across neuron types, and is not comparable between genes."
        ),
        h6( paste0("WARNING: Expression values for ",msg," are unreliable as they have been overexpressed to generate transgenic strains."), style="color:orange"),
        textAreaInput(
          inputId = "genelist",
          label = "Introduce a list of genes",
          value = "flp-1\nflp-2,flp-3,WBGene00001447\nWBGene00001448\nflp-6\nflp-7\nflp-8\nflp-9\nflp-10\nflp-11\nflp-12\nflp-13\nflp-14\nflp-15\nflp-16\nflp-17\nflp-18\nflp-19\nflp-20\nflp-21\nflp-22\nflp-23\nflp-24\nflp-25\nflp-26\nflp-27\nflp-28\nflp-32\nflp-33\nflp-34",
          width = "500px",
          height = "100px"
        ),
        
        
        selectInput(
          inputId = "dataset_heatmap",
          label = "Choose dataset: Neurons (threshold 2), All cells (unfiltered)",
          choices = c("Neurons only", "All cell types")
        ),    
        
        
        actionButton(
          "PlotHeatmap","Plot heatmap from list",
          icon = icon("hand-o-right")
          
          
        ),
        hr(),
        fluidRow(
          column(3,fileInput("file1", NULL,
                             accept = c(
                               "text/csv",
                               "text/comma-separated-values,text/plain",
                               "txt")
          )),
          #column(3,actionButton("resetFile", "Clear uploaded file")),
          column(3,   actionButton(
            "PlotHeatmap2",
            "Plot heatmap from file",
            icon = icon("hand-o-right")
          ))
        ),
        
        hr(),
        h6(
          "You can identify circles by clicking on them."
        ),
        
        div(style="height:30px;width:800px;padding-left:10px;padding-right:10px;background-color:#ffffff;",fluidRow(verbatimTextOutput("vals", placeholder = TRUE))),
        #uiOutput("dynamic"),
        br(),
        br(),
        br(),
        plotOutput("heatmap", width = "100%",hover = "plot_hover"),
        
        hr(),
        downloadLink("downloadheatmap", "Download plot")
        
        
      )
    )
  )
)

