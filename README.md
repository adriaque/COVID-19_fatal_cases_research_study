# COVID-19_fatal_cases_research_study
The COVID-19 Fatal Cases Research Study is a data analysis and visualisation project aimed at understanding the gene expression patterns, enriched ontology terms and differences in metabolic pathways associated with COVID-19 in lung and colon tissues.

This repository contains the code for a Shiny web application that performs a research study on COVID-19 fatal cases, specifically in lung and colon tissues. The project is divided into three parts: 

1. **Differential Expression Analysis and Plots:** Analyze and visualize gene expression data.
2. **Enrichment Analysis:** Perform gene ontology and pathway enrichment analysis.
3. **Pathway Analysis in Hipathia:** Visualize KEGG pathways and altered features.

## Getting Started

To run this Shiny app on your local machine, follow these steps:

1. **Prerequisites:** You need R and RStudio installed.

2. **Clone this Repository:**
   ```sh
   git clone https://github.com/adriaque/COVID-19_fatal_cases_research_study.git
   
3. Set Working Directory: Open the project in RStudio and set your working directory to the project folder.
   
4. Install Required Packages: Install the necessary R packages by running the following commands in RStudio:
   ```sh
   install.packages(c("shiny", "DESeq2", "readr", "dplyr", "readxl", "clusterProfiler", "biomaRt", "org.Hs.eg.db", "hipathia", "ggplot2"))
   
5. Run the Shiny App: Open the app.R file in RStudio and click the "Run App" button to launch the Shiny application.
Usage

## Demo Video

In this section tyuo can see a desmostration on how to run and use all the app features.

https://github.com/adriaque/COVID-19_fatal_cases_research_study/assets/113828670/cbbc0bcb-3098-47db-9db6-f6a49ac0ff68

## Contact and Contributing

Feel free to contribute to this project by opening issues, suggesting improvements, or creating pull requests.
Please contact sby sending an email to alejandroadriaquelozano@gmail.com

## License

This project is licensed under the MIT License - see the LICENSE file for details.



