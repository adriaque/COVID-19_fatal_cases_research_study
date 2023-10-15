# Load required libraries
library(shiny)
library(dplyr)
library(readr)


paths<- read.csv('./pathways_names_file.csv')
# Define UI

ui <- fluidPage(
  
  titlePanel("COVID-19 Fatal cases (Lung and Colon) research study"),
  tabsetPanel(
    tabPanel("Part 1: DEA and Plots", 
             h3('Differential Expression Analysis'),
             sidebarLayout(
               sidebarPanel(
                 
                 # Input: Choose a factor for group comparison
                 selectInput("Contrast", "Select a Contrast for Group Comparison",
                             choices = c('COVID-19 vs Normal', 'Lung vs Colon')),
                 
                 conditionalPanel(
                   condition = "input.Contrast == 'COVID-19 vs Normal' ",  # Only show when Lung vs Colon is selected
                   selectInput("Tissue", "Select a Reference Level for the Factor",
                               choices = c('Lung','Colon'))),
                 
                 conditionalPanel(
                   condition = "input.Contrast == 'Lung vs Colon' ",  # Only show when COVID-19 vs Normal factor is selected
                   selectInput("Disease", "Select a Reference Level for the Factor",
                               choices = c('COVID-19','Normal'))),
                 
                 # Input: Adjusted p-value cutoff
                 numericInput("pvalueCutoff", "Adjusted p-value Cutoff", 0.05),
                 
                 # Input: Log2 fold change cutoff
                 numericInput("log2FCcutoff", "Log2 Fold Change Cutoff", 1),
                 
                 # Run button
                 actionButton("runButton_1", "Run Analysis")
               ),
               
               
               mainPanel(
                 # Display the differential expression results
                 tableOutput("diffExprTable"),
                 
                 style = "max-width: 800px; max-height: 400px; overflow: auto;"
               )
             ),
             fluidRow(
               column(6,
                      
                      h3('PCA Plot'),
                      plotOutput('pcaPlot')),
               
               column(6, 
                      h3('Volcano Plot'),
                      plotOutput('volcanoPlot'))),
             
             
             
    ),
    tabPanel("Part 2: Enrichment Analysis", 
             h3("Gene Ontology Enrichment Analysis"),
             
             # Run button
             sidebarLayout(
               sidebarPanel(
                 
                 selectInput("updown", "Choose if you want to see Up or Down regulated Ontologies:", choices = c('Up-regulated', 'Down-Regulated')),
                 selectInput("ontology", "Choose the ontology term you want to plot:", choices = c('Molecular Function (MF)' , 'Cellular Component (CC)','Biological Process (BP)')),
                 actionButton("runButton_2", "Run Analysis")),
               mainPanel(
                 plotOutput('GOGraph')
                 
               )
             ),
             h3("Over-representation analysis in WikiPathways"),
             fluidRow(plotOutput("wikiPath"))
             
             ),
    tabPanel("Part 3: Pathway Analysis in Hipathia", 
             h3("Hipahtia Report"),
             
             sidebarLayout(
               sidebarPanel(
                 
                 selectInput("path", "Choose the KEGG Pathway you want to visualise:", choices = paths$all_pathway_names),
                 actionButton("runButton_3", "Run Analysis")),
               mainPanel(
                 uiOutput("hipathia")
               )),
             
             h3("Top features altered in the selected pathways pathway"),
             
             sidebarLayout(
               sidebarPanel(
                 
                 numericInput("number_genes", "Choose how many features you want to visualise per category:", 10),
                 actionButton("runButton_4", "Run Analysis")),
               mainPanel(
                 ' '
               )
               
               
             ),
             fluidRow(plotOutput("DAtop"))
    )
  )
)




# Define the server logic
server <- function(input, output) {
  
  #Load files and pre-process
  coldata <- read_tsv('./E-ENAD-46-experiment-design.tsv')
  lung_data <- read_delim('./E-ENAD-46-raw-counts.tsv',delim = '\t')
  lung_data_good <- lung_data[,3:40]
  row.names(lung_data_good) <- lung_data$`Gene ID`
  rm(lung_data)
  
  
  # Create a data frame representing the design matrix
  design <- data.frame(row.names = coldata$Run, Condition = coldata$`[disease]`,Tissue= coldata$`[organism part]`)
  
  
  # Create a DESeqDataSet object
  library(DESeq2)
  dds_data <- DESeqDataSetFromMatrix(countData = lung_data_good, colData = design, design = ~ Condition + Tissue)
  
  rm(coldata)
  
  # Perform differential expression analysis
  dds <- DESeq(dds_data)
  rm(dds_data)
  gc()
  
  observeEvent(input$runButton_1, {
    
    # Perform differential expression analysis and display results
    output$diffExprTable <- renderTable({
      req(input$Contrast)
      req(input$pvalueCutoff)
      req(input$log2FCcutoff)
      
      
      ContrastOfInterest <- input$Contrast
      if (ContrastOfInterest=='COVID-19 vs Normal')
      {
        req(input$Tissue)
        tissueLevel <- input$Tissue
      }
      
      if (ContrastOfInterest=='Lung vs Colon')
      {
        req(input$Disease)
        diseaseLevel <- input$Disease
      }
      
      
      
      
      #Create contrasts object depending on the contrast selected
      if (ContrastOfInterest=='COVID-19 vs Normal')
      {
        contrast <- c("Condition", "COVID-19", "normal")
        
        # For the contrast COVID vs Normal; select the tissue
        if (tissueLevel=='Lung'){
          dds <- dds[dds$Tissue == 'lung', ] #covid vs normal in lung
          contrast_hipathia <<-  "Reference_Lung - COVID_Lung"
          
        }
        else{
          dds <- dds[dds$Tissue == 'colon', ] # covid vs normal in colon
          contrast_hipathia <<-  "Reference_Colon - COVID_Colon"
        }
        
      }
      else #Now the situation where the contrast is lung vs colon
      {
        contrast <- c("Tissue", "lung", "colon")
        
        # Modify dds 
        if (diseaseLevel=='COVID-19'){
          dds <- dds[dds$Condition == 'COVID-19', ] #lung vs colon in covid
          contrast_hipathia <<-  "COVID_Lung - COVID_Colon"
          
        } else{
          dds <- dds[dds$Condition == 'normal', ] #lung vs colon in normal
          contrast_hipathia <<-  "Reference_Lung - Reference_Colon"
        }
      }
      
      
      pvalue_cut <-input$pvalueCutoff
      logfold_cut <- input$log2FCcutoff
      
      res <<- data.frame(results(dds,contrast))
      # Filter results based on p-value and fold change cutoffs
      res_filtered_1 <<- res %>% filter(padj<  pvalue_cut)
      res_filtered <- res_filtered_1 %>% filter(log2FoldChange < -logfold_cut | log2FoldChange > logfold_cut)
      
      return(res_filtered)
      
    },include.rownames = TRUE)
    
    output$pcaPlot <- renderPlot({
      library(ggplot2)
      # Perform PCA
      vsd <- vst(dds)
      pca_result <- prcomp(t(assay(vsd)))
      
      # Calculate the percentage of explained variance for each principal component
      explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
      
      # Create a PCA plot
      pca_data <- as.data.frame(pca_result$x)
      rownames(pca_data) <- colnames(lung_data_good)
      pca_data$Condition <- dds$Condition
      pca_data$Tissue <- dds$Tissue
      
      pcaPlot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, shape = Tissue)) +
        geom_point(size = 3) +
        scale_shape_manual(values = c("circle", "triangle")) +
        theme_minimal() +
        geom_text(aes(label = rownames(pca_data)), size = 3, nudge_x = 0.1, nudge_y = -5) +
        xlab(paste("PCA 1 (",round(explained_var[1],2),'%)')) +
        ylab(paste("PCA 2 (",round(explained_var[2],2),'%)'))
      
      return(pcaPlot)
      
      
    })
    
    output$volcanoPlot <- renderPlot({
      
      req(input$pvalueCutoff)
      req(input$log2FCcutoff)
      
      # Set log2 fold change and pvalue cutoffs
      log2fc_cutoff <- input$log2FCcutoff
      significance_cutoff <- input$pvalueCutoff
      
      
      # Prepare the data, including applying the significance cutoff
      res <- res %>%
        mutate(significant = case_when(
          padj < significance_cutoff & log2FoldChange > log2fc_cutoff ~ "Upreg",
          padj < significance_cutoff & log2FoldChange < -log2fc_cutoff ~ "Downreg",
          TRUE ~ "Not Significant"
        ))
      
      # Define custom legend labels
      legend_labels <- c("Upreg" = "Up-Regulated Genes",
                         "Downreg" = "Down-Regulated Genes",
                         "Not Significant" = "Non-DE Genes")
      
      # Customize volcano plot
      ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(aes(color = significant), size = 1) +
        scale_color_manual(values = c("Upreg" = "red", 'Downreg' ='blue',"Not Significant" = "grey"), 
                           labels = legend_labels, 
                           name = "Gene Expression Significance") +  # Specify custom legend title
        geom_hline(yintercept = -log10(significance_cutoff), linetype = "dashed", color = "gray") +
        geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "gray") +  # Add fold change cutoff
        labs(
          x = "Log2 Fold Change",
          y = "-log10(p-value)") +
        theme_minimal() +
        theme(legend.text = element_text(size = 14))
      
      
      
    })
  })
  
  
  library(biomaRt)
  mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  
  observeEvent(input$runButton_2, {
    library(clusterProfiler)
    output$GOGraph <- renderPlot({
      req(input$updown)
      req(input$log2FCcutoff)
      req(input$ontology)
      
      #Get the ontology
      if (input$ontology =='Biological Process (BP)'){
        ont='BP'
      }else if(input$ontology=='Cellular Component (CC)'){
        ont='CC'
      }else {
        ont='MF'
      }
      
      
      if (input$updown=='Up-regulated'){res_filtered_1=res_filtered_1 %>% filter(log2FoldChange>input$log2FCcutoff)}
      else {res_filtered_1=res_filtered_1 %>% filter(log2FoldChange< -(input$log2FCcutoff))}
      
      
      
      #Change the input up or down regulated genes
      test_input_ENSEMBL <- rownames(res_filtered_1)
      rm(res_filtered_1)
      
      # Get gene symbols
      entrez_names <<- getBM(c("entrezgene_id"), "ensembl_gene_id", test_input_ENSEMBL, mart)
      
      test <- enrichGO(entrez_names$entrezgene_id, keyType = 'ENTREZID' ,ont = ont, OrgDb = 'org.Hs.eg.db', pvalueCutoff = 1, qvalueCutoff = 1)
      plot <-dotplot(test, showCategory =10)
      rm(test)
      return(plot)
      
    })
    
    output$wikiPath <- renderPlot({
      
      enrich_result <- enrichWP(gene=entrez_names$entrezgene_id ,organism = "Homo sapiens" )
      wikiPath<-dotplot(enrich_result)
      return(wikiPath)
      
    })
  })
  
  #Pre-processing for hypathia visualisation
  #Translation ENSEMBL to ENTREZ
  library(hipathia)
  trans_lung <- translate_data(as.matrix(lung_data_good), "hsa") 
  
  #Pre-processing
  trans_lung <- normalize_data(trans_lung)
  
  groups_names <- c('Reference_Lung','COVID_Lung','COVID_Lung','COVID_Lung','COVID_Lung','COVID_Lung','COVID_Lung','COVID_Lung','Reference_Lung','Reference_Lung','Reference_Lung','Reference_Lung','Reference_Lung','Reference_Lung','Reference_Lung','Reference_Lung','Reference_Lung','COVID_Lung','COVID_Lung','Reference_Colon','COVID_Colon','COVID_Colon','COVID_Colon','COVID_Colon','COVID_Colon','COVID_Colon','COVID_Colon','COVID_Colon','COVID_Colon','Reference_Colon','Reference_Colon','Reference_Colon','Reference_Colon','Reference_Colon','Reference_Colon','Reference_Colon','Reference_Colon','Reference_Colon')
  groups <- data.frame(groups_names)
  row_names <- colnames(lung_data_good)
  rownames(groups)<- row_names
  
  observeEvent(input$runButton_3,{
    
    output$hipathia <- renderUI({
      req(input$path)
      selected_path <- input$path
      seledcted_hsa <- paths$pathway_ids[paths$all_pathway_names == selected_path]
      
      data <-SummarizedExperiment(assays=SimpleList(raw=trans_lung),
                                  colData=groups)
      
      pathways_only <- load_pathways(species = "hsa", pathways_list = seledcted_hsa)
      
      hidata <- hipathia(data, pathways_only, uni.terms = TRUE, GO.terms = TRUE,decompose = FALSE, verbose=TRUE)
      
      
      DAdata <- DAcomp(hidata, "groups_names", contrast_hipathia, path.method = "limma",fun.method = "limma")
      DAdata2 <<- DAdata
      # Port randomise
      port <- runif(1, min = 2000, max = 4000)
      port <- round(port)
      
      
      #continue the analysis
      HPreport <- DAreport(DAdata, pathways_only)
      visualize_report(HPreport,port = port)
      iframe_code <- paste0('<iframe src="http://127.0.0.1:', port, '" width="800" height="600"></iframe>')
      HTML(iframe_code)
    })
  })
  
  observeEvent(input$runButton_4,{
    
    output$DAtop <- renderPlot({
      req(input$number_genes)
      n<-input$number_genes
      DAtop<-DAtop(DAdata2,n = n)
      return(DAtop)
      
    })
  })
}


# Run the Shiny app
shinyApp(ui = ui, server = server)
