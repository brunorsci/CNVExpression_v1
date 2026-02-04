---
title: "Correlação entre Copy Number Variations e Gene expression"
author: "Bruno Rodrigo Assunção"
date: "2024-02-10"
output: github_documento
---

#### WORKFLOW- análise de variantes do número de cópias e integração com Sequencimento de RNA em dados de Leucemia Mielóide Aguda

Analisamos amostras de Leucemia Mielóide do projeto ICGC Data Portal TCGA-LAML obtidas do banco de dados ICGC Data Potal <https://dcc.icgc.org/>. O objetivo foi processar e integrar a análise de "Copy Number Variations" e Sequênciamento de RNA (expressão de genes). Este workflow computacional criar novas funções, bem como utilizar diferentes ferramentas para trabalhar com dados de projetos TCGA e concatena os algoritimos computacionais em um pacote denominado CNVcExp.

[Bibliotecas necessárias:]{.underline}

-dplyr

-biomaRt

-org.Hs.eg.db

-pbapply

#### 1- Input dos arquivos de cnvs e expressão de RNA (RNASeq).

```{r}
cnvs <- read.delim("C:/Users/bruno/Desktop/working/LMA-CNV/Data/icgc/tcga lma us/copy_number_somatic_mutation.LAML-US.tsv") 
exp_seq <- read.delim("C:/Users/bruno/Desktop/working/LMA-CNV/Data/icgc/tcga lma us/exp_seq.LAML-US.tsv")
```

#### 2- Modificar identificadores conforme o TCGA Barcode

``` r
library(stringr)

# Substituir hífens por underscores
samples <- str_replace_all(cnvs[["submitted_sample_id"]], "-", "_")

# Função para modificar o identificador TCGA
mudar_identificador_tcga <- function(x, sep) {
  sapply(x, function(i) {
    partes <- unlist(strsplit(i, sep))
    paste(partes[1:4], collapse = sep)
  })
}

# Aplicando a função
samples <- mudar_identificador_tcga(samples, "_")

# A numeração ao final indica se é uma amostra de tumor ou de tecido normal. 
samples <- gsub("3A", "3", x = samples) 
samples <- gsub("11A", "11", x = samples)
# Adicionar samples ids ao objeto cnv
cnvs <- cbind(samples, cnvs)
#Remover amostras "Normal"
keep = grep("-03A-", cnvs$submitted_sample_id) 
keep = cnvs[keep,] cnvs=keep
```

Acima, mantemos apenas as amostras de tumor através da interpretação do padrão do TCGA Barcode, no qual nesse estudo as amostras tumorais são identificadas pelo sufixo "-03A-". Verifique se a tabela "cnv" contém algum identificador com o sufixo "11A" ou "11B", o qual nesse projeto é utilizado para identificar as amostras de tecido normal.

**Selecionar colunas e renomea-las**

``` r

cnvs <- cnvs %>%
  select(samples, donor_id = icgc_donor_id, chr = chromosome, start = chromosome_start, 
    end = chromosome_end, segment_mean)
```

#### 4- Remover os cromossomos X e Y e mudar os valores 1 para chr1.

Utilizamos a função para modificar a nomenclatura da nomeação dos cromossomos e também para remover os cromossomos sexuais, as quais não desejamos em nossa análise.

``` r
cnvs <- rm_sexual(df = cnvs, col = "chr", remove_X = TRUE, remove_Y = TRUE)
```

Isso é util a depender do genomas de referência que irá utilizar mais a frente, adepender da referência a nomenclatura dos cromossomos estará de uma forma ou de outra.

**4.1 Podemos utilizar gsub() para remover o caracter "chr" da nomeação dos cromossomos.**

``` r
cnvs$seqnames <- gsub(x = cnvs$seqnames, "chr", "")
```

#### 5 Anotar CNVS: duplicações ou deleções

#### 5.1- Transformar segment_mean para números inteiros de cópias

```{r}
chr_cnv <- convert_segment_mean_to_cn(chr_cnv, col = "segment_mean", 
                                      cn_col = "copy_number", 
                                      annotation_col = "annotation")

```

Esta função converte valores de Segment Mean para número de cópias usando a fórmula:

`CN = 2 * (2 ^ SM) e arredonda o resultado para o inteiro mais próximo.`

Também adiciona uma anotação com base no CN calculado:

\- CN = 0 -\> "deleção total"

\- CN = 1 -\> "deleção parcial"

\- CN = 2 -\> "neutral"

\- CN \> 2 -\> "duplicação"

#### 6- Anotar genes que se sobrepõem as regiões de cnvs

Utilizamos um genoma de referência para anotar genes que se encontram nas regiões genomicas as quais identificamos cnvs (ganhos ou perdas de segmentos de DNA).

O genoma de referêmcia "GRCh37" foi utilizado, de acordo ao genoma de referência utilizado para mapear as amostras deste estudo. Os genomas de referência podem ser obtidos no banco de dados UCSC genome.

``` r
cnvs_ann_result <- cnv_annotate_genes(input = chr_cnv, ref_genome = "hg19", organism = "hsapiens_gene_ensembl")
```

Foi adicionada a coluna gene_ids ao dataframe com os cnvs. Essa coluna contêm os ENSEMBL ids dos genes sobrepostos a cada região em formato de lista para cada grange. Algumas regiões não possuem genes, o que gera uma coluna vazia.

**6.1 Remover linhas não identificadas (que não sobrepõem genes)**

``` r
cnvs_ann_genes=NULL 
for(i in 1:10){ 
  if (S4Vectors::isEmpty(cnvs$gene_ids[i])==FALSE){ 
    y=cnv[i,] cnvs_ann_genes=rbind(cnvs_ann_genes,y)} 
  }
```

**6.2 Expandir dataframe: separar genes em linhas individuais**

**Use a função unnest para expandir a coluna gene_ids**

``` r
cnvs_ann_genes_new <- cnvs_ann_genes %>% tidyr::unnest(cols = gene_ids) %>% filter(gene_ids != "") # Filtras as colunas vazias
```

#### 7- Anotar simbolos de genes

**7.1 Conectar ao database biomaRt**

``` r
library(biomaRt) ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
```

Objeto com os identificadores de genes

``` r
gene_ids <- cnvs_ann_genes_new$gene_ids
```

**7.2 Obter os simbolos genicos, ensemble e entrezid**

``` r
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"), filters = "ensembl_gene_id", values = gene_ids, mart = ensembl, filter)
```

**7.3 Remover duplicatas**

``` r
colnames(gene_symbols) <- c("gene_ids", "symbols", "entrezgene_id") gene_symbols <- gene_symbols %>% distinct(gene_ids, .keep_all = TRUE) #gene_symbols <- gene_symbols[!duplicated(gene_symbols$gene_ids), ]
```

**7.4 Juntar ao dataframe original com left join**

Esta função manterá as colunas que não tiveram simbolos gênicos anotados.

``` r
cnvs_ann_genes_new <- cnvs_ann_genes_new %>% 
  left_join(gene_symbols, by = c("gene_ids" = "gene_ids"))
```

Os dados com as informações que desejamos estão no objeto cnvs_ann_genes_new.

*Observação: se o conjunto de dados for muito grande/pesado, é possível realizar os processos de anotação de genes dividindo o dataframe. Além de outras maneiras computacionalmente mais eficazes através das utilização de funções como apply, lapply, pbapply, dentre outras. Basta modificar as funções acima para que elas aceitem apenas um valor como entrada. Temos essas outras maneiras de fazer isso, mas não está disponibilizados nessa demonstração.*

**Alternativa 1:** Dividir os data frames e realizar a anotação em cada parte individualmente.

Esta função divide os dataframes de acordo a quantidade de partes que se definir. Como exemplo estaremos dividindo-o em 100 partes iguais. Cada parte ficará disponível no ambiente global da sessão a ser utilizada.

[*Função 1:*]{.underline} `percentiles_df_divide()`.

**Divide o data frame.**

``` r
percentiles_df_divide <- function(df, n_parts){
  n_rows <- nrow(df) 
  perc <- seq(0, 100, length.out = n_parts + 1) 
  parts <- vector("list", n_parts)
  for (i in 1:n_parts){ 
    start_row <- floor(perc[i] * n_rows / 100) + 1 
    end_row <- floor(perc[i + 1] * n_rows / 100) 
    parts[[i]] <- df[start_row:end_row, ]
    } 
  return(parts)
  }
```

[Função 2]{.underline}: `cnvr_to_genes_divided()`.

Realiza a anotação de genes que se sobrepõem as regiões de cnvs.

``` r
cnv_genes <- cnvr_to_genes_divided(df = cnv, n=100)
```

Para esse passo acima, também se aplica o que está disposto na explicação abaixo.

[Função 3]{.underline}: `cnv_individual_cromossomes()`. Optamos por utilizar essa função e realizar a anotação de genes em cada cromossomo individualmente.

Esta função separa cada cromossomo para realizar a anotação individualmente. Falta a implementação para expandir o data frame de acordo a colunas gene_ids, adicionada durante o processo de sobreposição de genes, como realizado em etapas anteriores. A falta dessa etapa pode gerar problemas durante a anotação dos simbolos de genes e entrezid. Após aplicar essas função, check os resultados, verifique a coluna gene_ids e realize a etapa de expansão do data frame.

``` r
cnv_individual_chromossomes <- function(n=22,
                                        dataf=cnv,
                                        chr_colx="seqnames", 
                                        genes_colz="gene_ids"){ 
  pb <- txtProgressBar(min = 0, max = n,style = 3, width = 50, char = "=")
  
  for (i in 1:n) { 
    assign(paste0("chr",i), cnv_extract_chromossomes(chr_df=dataf, chr_col=chr_colx)) 
    assign(paste0("chr",i), cnvr_to_genes(df = paste0("chr",i), col = genes_colz)) 
    assign(paste0("chr",i), paste0("chr",i)[,-10])} 
  }
```

#### 8- Preparar a tabela de expressão

Mais a frente iremos realizar o cálculo da correlação da expressão de genes e as mudanças no número de cópias. Para isso, precisamos que a tabela com as expressões esteja em um formato especifíco, como veremos a seguir.

**8.1 Criar a coluna samples, coluna samples+genes**

Criar esse novo identificador torna possível individualizar/ relacionar cada gene corretamente a sua amostra.

-   mudar o identificador do tcga para a tabela de expressão como feito para a tabela de cnvs em etapas anteriores\*

``` r
samples <- chr_exp[["submitted_sample_id"]]

samples=gsub(x = samples, "-", "_")

samples <- pbapply(as.data.frame(samples), 1, FUN = mudar_identificador_tcga)

chr_exp <- cbind(samples, chr_exp)
```

Vamos criar os identificadores samples+ genes: podemos fazer para ambas as tabelas, cnvs e expressão, utilizando uma das maneiras a seguir.

``` r
chr_cnv <- sample_genes_ids(x1=cnvs, col_sample="samples", col_gene="SYMBOL") 
chr_exp <- sample_genes_ids(x1=chr_exp, col_sample="samples", col_gene="gene_id")
```

[Ou assim, utilizado pbaply: mais rápida.]{.underline}

As funções internas são disponibilizadas mais abaixo do código de exemplificação.

``` r
library(pbapply) 
chr1_cnv= dplyr::select(chr1_cnv, "samples", "SYMBOL", "segment_mean", "type" ) 
chr1_cnv_new <- pbapply::pbapply(chr1_cnv, 1, FUN = sample_genes_ids_new) 
chr1_exp_new <- pbapply::pbapply(chr1_exp, 1, FUN = sample_genes_ids_new)
```

Observar a coluna que corresponde aos simbolos dos genes em cada objeto e imputar corretamente na função.

**Funções**:`sample_genes_ids()`

Essa função junta amostras e simbolos de genes/nomes em única coluna.

``` r
sample_genes_ids <- function(x1=chr1_cnv, col_sample="samples", col_gene="SYMBOL"){
  
  pb <- txtProgressBar(min = 0, max = nrow(x1), style = 3, width = 50, char = "=")
  
  sample_genes <- NULL 
    
  for (rows in 1:nrow(x1)){ 
      ids<-paste0(x1[[col_sample]][rows], "_", x1[[col_gene]][rows])
      sample_genes <-rbind(sample_genes, ids)
      close(pb)
      setTxtProgressBar(pb, rows)}
  x2 <-cbind(sample_genes, x1) 
  rownames(x2)<-NULL return(x2) 
}
```

[Alternativa a função acima]{.underline}. Essa é mais rápida. sample_genes_ids_new().

Para ser utilizada com pbaply ou outras funções semelhantes. É a função acima simplificada.

``` r
sample_genes_ids_new <- function(x){ 
  # OBS.: The function expect sample names in vector (column) 1 and genes names in      vector (column) 2 
   
   # Loop to create the labels 
   sample_genes <- paste0(chr_cnv[1], "_", chr_cnv[2]) 
   
   # return a data frame with the column containing samples+genes names 
   rownames(sample_genes)<-NULL 
   return(sample_genes)
  }
```

#### 9- Correlação

**9.1 Criar a tabela para cálculo da correlação**

``` r
cordata <- merge(x = chr_cnv, y = chr_exp, by = "sample_genes")
```

Observe que utilizamos o identificador sample + genes que criamos.

**9.2 Calcular a correlação entre a variação na expressão de acordo a variação do número de cópias.**

Aqui aplicamos o método de "spearman", visto que esse não exije que os dados sigam uma distribuição normal.

``` r
correlation_results_chr1 <- cnv_exp_correlation(df = cordata, methodx = "spearman")
```

**Função**:`cnv_exp_correlation()`

Essa é uma função para testar a correlação entre as regiões de cnvs para cada gene individualmente anotado e a expressão de cada gene. Os métodos de correlação podem ser "pearson", "kendall", "spearman".

-   df = dataframe com uma coluna com os genes e as contagens de segment e read count

-   gene_ids_col= coluna que contem os genes

``` r
 cnv_exp_correlation <- function(df=z, 
                                 x_col="segment_mean",
                                 y_col="normalized_read_count", 
                                 gene_ids_col="gene_id", 
                                 methodx="spearman"){ 
   
  pb <- txtProgressBar(min = 0, 
                        max = length(df[[gene_ids_col]]),
                        style = 3, 
                        width = 50, 
                        char = "=") 
  correlation_results<-NULL 
 
   for (gene in unique(df[[gene_ids_col]])){ 
     
     setTxtProgressBar(pb, gene) 
     x1 <- df[df[[gene_ids_col]]==gene,] 
     
     cor <- cor.test(x=x1[[x_col]], 
                     y=x1[[y_col]], 
                     method = methodx) 
     
     statistic <- cor$statistic 
     estimate <- cor$estimate 
     p_value <- cor$p.value 
     method <- cor$method 
     res <- cbind(gene, statistic, estimate , p_value, method) 
     rownames(res) <- NULL 
     correlation_results <- rbind(correlation_results, res)
     } 
   close(pb) 
 return(correlation_results)
   } 
```

[Quanto aos resultados obtidos pela função acima:]{.underline}

-   S é o valor do teste estatístico (S = 10.871)

-   p-value é o nível de significância do teste estatístico (p-value = 0.4397)

-   Hipotese alternativa é um caracter descrevendo a hipotese alternativa (VERDADEIRO rho não é igual a 0)

-   Estimativa da amostra é o coeficiente de correlação. Para correlação de Spearman o coeficiente de correlação é denominado rho (Cor.coeff=0.4564)

#### 10- Visualização gráfica

**Avaliar a correlação da tabela e plotar apenas genes únicos.**

Em z teremos a tabela concatenada com apenas um gene, x é o dado normalizado da contagens das reads no experimento de expressão de rna e y= a media dos segmentos do experimento de cnvs.

**Plot 1: Scatter Plots**

``` r
library("ggpubr") 

ggscatter(z, 
          x = "normalized_read_count", 
          y = "segment_mean", 
          add = "reg.line", 
          conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.method = "pearson", 
          xlab = "CNVs segment mean", 
          ylab = "Normalized read counts")
```
