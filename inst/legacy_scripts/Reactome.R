#### DATABASE CNVEXPRESSION ANNOTATION ####

## 1- REACTOME ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ReactomePA", "clusterProfiler", "org.Hs.eg.db", "dplyr", "tibble"))

library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)

# Filtrar genes diferencialmente expressos

# Ajustar o corte conforme necessário
genes_signif <- CNVExpression %>%
  filter(!is.na(entrez_id), AdjPValue < 0.05) %>%
  distinct(entrez_id, .keep_all = TRUE)

# Extrair apenas os IDs Entrez
#entrez_ids <- genes_signif$entrez_id

# Enriquecimento com ReactomePA
reactome_results <- enrichPathway(gene = genes_signif$entrez_id,
                                  organism = "human",
                                  pvalueCutoff = 0.1,
                                  readable = TRUE)
# Visualização simples
head(reactome_results)

# Gráfico de barras dos caminhos enriquecidos
barplot(reactome_results, showCategory = 20, title = "Reactome Significant Genes")

# Gráfico de bolhas
dotplot(reactome_results, showCategory = 20, title = "Reactome Significant Genes")

##cnvsfinally
reactome_cnvsFinally <- enrichPathway(gene = cnvsFinally$entrez_id,
                                  organism = "human",
                                  pvalueCutoff = 0.1,
                                  readable = TRUE)
barplot(reactome_cnvsFinally, showCategory = 20, title = "Reactome Significant Genes")

# Gráfico de bolhas
dotplot(reactome_cnvsFinally, showCategory = 20, title = "Reactome Significant Genes")

##
library(enrichplot)
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)


# Visualização com enrichMap:
reactome_cnvsFinally <- pairwise_termsim(reactome_cnvsFinally)

emapplot(reactome_cnvsFinally)
dotplot(reactome_cnvsFinally)
cnetplot(reactome_cnvsFinally)

##
cnetplot(reactome_cnvsFinally, 
         foldChange = cnvsFinally$logFC,  # opcional: muda cor dos genes com base em logFC
         showCategory = 10,       # quantas vias mostrar
         circular = FALSE,        # layout circular ou em rede
         colorEdge = TRUE)        # cor das conexões com base em peso

