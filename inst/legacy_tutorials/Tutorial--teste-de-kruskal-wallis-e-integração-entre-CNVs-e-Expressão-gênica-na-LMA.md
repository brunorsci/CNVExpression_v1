CNKExp: Teste de kruskal-wallis na integração entre CNVs e Expressão
Gênica na LMA
================
Bruno Rodrigo Assunção
2024-02-16

### CNKExp: Teste de kruskal-wallis na integração entre CNVs e Expressão Gênica na LMA

Neste tutorial exemplificamos a aplicação do teste de hipótese
kruskal-wallis na integração de “Copy Number Variation” e Expressão de
Genes. Essa integração foi demonstrada no tutorial
[CNeEXP](https://rpubs.com/Brunor/1149286 "CNeEXP"). Os testes de
avaliação da distribuição normal foi previamente realizada no referido
tutorial acima, o que demonstrou que os dados não seguem uma
distribuição normal, deste modo, sendo indicado o teste aqui
exemplificado. A doença em estudo é a Leucemia Mielóide Aguda (LMA).

#### 00 - Input dos dados

#### 01 - Filtragem dos genes de interesse

Aqui iremos filtrar apenas os genes de interesse, isto é, os
relacionados a Leucemia Mielóide Aguda.

##### **1.1 Genes de interesse:**

``` r
genesAML <- c("KIT", "KITLG", "GRB2", " SOS","RAS","STAT","PI3","RAF","JAK","PIP","PKB","MTOR","MEK","MAPK","IKK","BAD","ERK","p7026k","PIM","EBP", "CEBP","PU", "RUNX1", "ETO", "PML","TRIM", "MYL", "PML", "RAR","CSF","GMCSF","PER2","CSF","CD14","CD64","CYC","MYC","DUSP","PPAR","TCF","JUP")

genesAML <- unique(genesAML)
```

##### **1.2 loop for para “pegar” as variações de nomenclatura dos simbolos dos genes de interesse:**

##### 1.3 Apenas unicos valores de genes:

**Vamos visualizar esses genes:**

``` r
genesAML_new[1:10]
```

    ##  [1] NA NA NA NA NA NA NA NA NA NA

#### 02- Filtrando o data frame para selecionar apenas os genes de interesse:

##### 2.1 Loop *for()* para o teste kruskal-wallis:

``` r
for (i in genesAML_new) {
  inputKruskal <- all_cnvs_analysis[all_cnvs_analysis$SYMBOL %in% i,]
  if (length(unique(all_cnvs_analysis[all_cnvs_analysis$SYMBOL %in% i,]$type))>1) {
    testeK       <- kruskal.test(x = inputKruskal$normalized_read_count,
                                 inputKruskal$type)
    assign(paste0("res_Kruskal_",i), data.frame(gene=i,
                                                statistic=testeK$statistic,
                                                p_value = testeK$p.value,
                                                method = testeK$method))
    }else{
      print(paste("The observation in some 'type' values only have one group:",i))
    }
}
resultsKruskal <- NULL
for (variable in grep(pattern = "res_K", value = T, ls())) {
  resultsKruskal <- rbind(resultsKruskal, get(variable))
  row.names(resultsKruskal) <-NULL
}
```

O loop for acima retorna os resultados do teste de kruskal-wallis para
todos os genes fornecidos como entrada (genesAML_new). Ele verifica se
tem algum gene no qual há apenas um grupo, ou seja, regiões neutras sem
variação de cópias. Os grupos fornecidos são três: deleções, duplicações
e regiões neutras.

##### 2.2 Gerando a tabela com os resultados acima:

Este trecho foi executado no código acima. Não é necessário executa-lo
novamente.

``` r
resultsKruskal <- NULL
for (variable in grep(pattern = "res_K", value = T, ls())) {
  resultsKruskal <- rbind(resultsKruskal, get(variable))
  row.names(resultsKruskal) <-NULL
}
```

##### 2.3 Podemos filtrar apenas os resultados com valor de P significativos

##### 2.4 Visualizar os resultados significativos

``` r
resultsKruskal[resultsKruskal$p_value <=0.05,][1:5,]
```

    ## NULL

Então, na tabela denominada resultsKruskal_Significativos temos os
resultados desta análise.
