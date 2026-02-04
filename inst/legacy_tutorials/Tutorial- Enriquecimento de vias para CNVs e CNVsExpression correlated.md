---
title: "Enriquecimento de Vias: integração de CNVs e Expressão na LMA"
author: "Bruno Rodrigo Assunção"
date: "2024-02-14"
output: github_document
---

## Corpathway: Enriquecimento de vias de Genes correlacionados com CNVs e relacionados a LMA.

Este é um fluxo de análise para investigar o enriquecimento de vias moleculares para os genes envolvidos em vias biológicas da Leucêmia Mielóide Aguda. Genes os quais nesse conjunto experimental possuem alterações no número de cópias gênicas apresentam correlação estatísticas com a expressão gênicas.

## Fundo

O conjunto de dados utilizados para esta análise parte da experimentação prévia realizada por ASSUNCAO B.R (2023) [CNeExp](https://rpubs.com/Brunor/1149286 "CNeExp"). Deste modo, temos dois data frames para filtragem e anotação de vias.

Estes são:

1.  anotação de cnvs e os seus genes de um experimento tcga.

2.  correlação entre estes cnvs e expressão gênica.

#### **00 - Input**

```{r include=FALSE}
all_cnvs_exp <- read.delim("C:/Users/bruno/Desktop/cnvs_all_annotated.csv")
chr_correlations <- read.csv("C:/Users/bruno/Desktop/Apêndice A - correlação de spearman entre segmentação (CNVs) expressão de rna (RNASeq).csv")
positive_correlations <- read.csv("C:/Users/bruno/Desktop/positive_correlations.csv")
```

#### 1 - Seleção dos genes de interesse

##### 1.1 Lista dos genes de interesse

```{r include=FALSE}
genesAML <- c("KIT", "KITLG", "GRB2", " SOS","RAS","STAT","PI3","RAF", "JAK","PIP","PKB","MTOR","MEK","MAPK","IKK","BAD","ERK","p7026k","PIM","EBP", "CEBP","PU", "RUNX1", "ETO", "PML","TRIM", "MYL", "PML", "RAR","CSF","GMCSF","PER2","CSF","CD14","CD64","CYC","MYC","DUSP","PPAR","TCF","JUP")

genesAML <- unique(genesAML)
```

Como esses genes possuem diferentes alelos e assim diferentes nomenclaturas, vamos utilizar um grep para buscar todos os genes na tabela de cnvs, utilizando a lista acima como referência. Deste modo, é possível obter todos os símbolos gênicos.

##### 1.2 Grep da lista gênica de interesse.

Utilizamos um loop for() e a função grep() em R

```{r include=FALSE}
# Loop para pegar as diferentes nomenclaturas para os simbolos dos genes no df
genesAML_new <- NULL
for (i in genesAML) {
  selected<- grep(pattern = i,x = all_cnvs_exp$SYMBOL ,value = T)
  genesAML_new <- c(genesAML_new, selected)
}
```

Do modo acima, a nossa lista foi atualizada para todas as variações presentes a partir da nossa lista de entrada. Então, podemos finalmente selecionar os cnvs envolvendo os nossos genes de interesse presentes na tabela de cnvs.

##### 1.3 Selecionar cnvs do data frame

```{r include=FALSE}
selectedGenesCNVs <- all_cnvs_exp[all_cnvs_exp$SYMBOL %in% genesAML_new,]
```

##### 1.4 Selecionar correlações

```{r}
selectedGenesCorrelation <- chr_correlations[chr_correlations$gene %in% genesAML,]
```

##### 1.5 Filtrar apenas regiões definidas como deleção ou duplicação.

```{r include=FALSE}
selectedGenesCNVs <- selectedGenesCNVs[selectedGenesCNVs$type != "neutral",]
```

##### 1.6 Tabela com a contagem de cnvs por genes

```{r include=FALSE}
selectedGenesCNVs <- selectedGenesCNVs[selectedGenesCNVs$type != "neutral",]
```

#### **2- Análise de enriquecimento de vias biológicas**

##### 2.1 Converter simbolos genicos para entrez ids

Primeiro devemos modificar os simbólos de genes para entrez id. Visto que esse é o input esperado pelo pacote reactomePA. Para isso vamos utilizar a função criada localmente symbols_to_entrez_id(), a qual esta depositada no pacote [CNeEXP](https://github.com/brunorsci/CNeEXP/tree/bruno "CNeEXP").

**Função criada:**

```{r}
symbols_to_entrez_id <-function(input){
  library(org.Hs.eg.db) # mudar a forma de carregar essa library
  entrez_ids=NULL
  for (i in input) {
    x <-AnnotationDbi::select(x = org.Hs.eg.db, 
                              keys=i, 
                              keytype = "SYMBOL", 
                              columns = "ENTREZID")
    entrez_ids=rbind(entrez_ids, x)}
    }
```

**Correlações positivas**

```{r include=FALSE}
inputPosCor <- symbols_to_entrez_id(positive_correlations$gene)
```

```{r echo=TRUE}
inputPosCor[1:5,]
```

**Correlações Negativas**

``` r
inputNegCor <- symbols_to_entrez_id(negative_correlations$gene)
```

##### 2.2 Realizando a análise de enriquecimento

```{r}
library(ReactomePA)
reactomePositiveCorrelations <- enrichPathway(
  sort(na.omit(inputPosCor$ENTREZID)), ## espera entrez_id
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  #universe,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE
)
```

``` r
library(ReactomePA)
reactomeNegativeCorrelations <- enrichPathway(
  sort(na.omit(inputNegCor$ENTREZID)), ## espera entrez_id
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  #universe,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE
)
```

**Armazenar os resultados em um data frame**

```{r}
reactomePositiveCorrelations <-as.data.frame(reactomePositiveCorrelations)
```

``` r
reactomeNegativeCorrelations <- as.data.frame(reactomeNegativeCorrelations)
```

```{r echo=TRUE}
reactomePositiveCorrelations[1:5,-1]
```
