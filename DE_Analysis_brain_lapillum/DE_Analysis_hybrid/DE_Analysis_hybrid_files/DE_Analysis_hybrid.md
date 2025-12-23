DE_analysis_L.Rmd
================
2025-12-23

Once gene counts have been obtained and summarised within the
RSEM_gene_counts file, differential expression analyses can be performed
on the data. For this the R Studio software is utilised. ***Please
note***, this code is repeated for each individual tissue and species
combination, used to compare expression between temperature treatments.
**This example is for the brain tissue of the common triplefin**.

# Loading necessary libraries

``` r
library(tidyr)
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomicRanges

    ## Loading required package: Seqinfo

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(pheatmap)
library(stringr)
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following object is masked from 'package:matrixStats':
    ## 
    ##     count

    ## The following objects are masked from 'package:GenomicRanges':
    ## 
    ##     intersect, setdiff, union

    ## The following object is masked from 'package:Seqinfo':
    ## 
    ##     intersect

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:generics':
    ## 
    ##     explain

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(matrixStats)
library(edgeR)
```

    ## Loading required package: limma

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:DESeq2':
    ## 
    ##     plotMA

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

``` r
library(grid)
library("RColorBrewer")
library(ggrepel)
```

# Loading and subsetting the count data

The raw RSEM gene count data can be accessed here
[RSEM_gene_counts.txt](https://github.com/breanariordan/triplefinRNA/blob/main/RSEM_gene_counts.txt)

``` r
# changing column names
original=read.csv("RSEM_gene_counts_hybrid.txt", sep="", head=T)
head(original)
```

    ##                   X10L1brain X10L1heart X10L2brain X10L2heart X10L3brain
    ## TRINITY_DN0_c0_g1      82.94     343.34      84.17     287.24      82.33
    ## TRINITY_DN0_c0_g2     126.62     502.20     129.65     337.51      97.98
    ## TRINITY_DN0_c0_g3      49.93     145.42      68.57     183.32     133.14
    ## TRINITY_DN0_c1_g1     102.89    1042.53     113.97     795.60     239.41
    ## TRINITY_DN0_c1_g2       3.78       5.46      10.54       4.81       2.43
    ## TRINITY_DN0_c2_g1     403.89    1714.76     366.53     816.07     220.29
    ##                   X10L3heart X10L4brain X10L4heart X10L5brain X10L5heart
    ## TRINITY_DN0_c0_g1     378.11      41.66     246.92      68.60     371.04
    ## TRINITY_DN0_c0_g2     345.32     103.16     365.21     136.99     343.94
    ## TRINITY_DN0_c0_g3     265.31     115.94     281.90     132.17     256.06
    ## TRINITY_DN0_c1_g1    1206.30     181.10    1095.26     266.42    1113.32
    ## TRINITY_DN0_c1_g2       8.14      10.57       9.59       2.82      13.71
    ## TRINITY_DN0_c2_g1    1449.88     340.80    1063.21     248.17    1118.59
    ##                   X10N1brain X10N1heart X10N2brain X10N2heart X10N3brain
    ## TRINITY_DN0_c0_g1      45.66      77.43      15.04      43.84      46.87
    ## TRINITY_DN0_c0_g2     803.41    1086.03     336.18     621.80     483.55
    ## TRINITY_DN0_c0_g3     231.96     170.20     149.65     217.25     193.63
    ## TRINITY_DN0_c1_g1     371.08    1230.13     159.73     988.80     156.01
    ## TRINITY_DN0_c1_g2      70.12      58.61      39.93      64.65      47.41
    ## TRINITY_DN0_c2_g1     330.78    1248.87     383.84    1474.00     283.67
    ##                   X10N3heart X10N4brain X10N4heart X10N5brain X10N5heart
    ## TRINITY_DN0_c0_g1     108.02      27.98     213.03      92.51     812.92
    ## TRINITY_DN0_c0_g2     778.83     506.65     753.76     543.31     898.59
    ## TRINITY_DN0_c0_g3     198.13     153.04     158.86     122.32     179.50
    ## TRINITY_DN0_c1_g1     902.53     160.39    1252.95     198.95    1417.81
    ## TRINITY_DN0_c1_g2      43.04      18.65      47.23      22.79      58.12
    ## TRINITY_DN0_c2_g1    1114.65     331.00    1721.08     378.00    2302.30
    ##                   X22L1brain X22L1heart X22L2brain X22L2heart X22L3brain
    ## TRINITY_DN0_c0_g1      70.94     519.92      86.85     413.27      50.44
    ## TRINITY_DN0_c0_g2     134.56     260.71      94.64     239.66     119.20
    ## TRINITY_DN0_c0_g3     305.84     257.71     152.10      82.98     144.10
    ## TRINITY_DN0_c1_g1     206.86     932.93     352.81     436.80     229.86
    ## TRINITY_DN0_c1_g2      19.91      24.18       5.93       4.49      10.77
    ## TRINITY_DN0_c2_g1     323.42    1389.82     358.73    1266.90     354.23
    ##                   X22L3heart X22L4brain X22L4heart X22L5brain X22L5heart
    ## TRINITY_DN0_c0_g1     414.02     108.86     725.57     127.21     470.38
    ## TRINITY_DN0_c0_g2     270.07     111.10     166.78      80.11     219.69
    ## TRINITY_DN0_c0_g3     233.35     120.25     267.30     276.37     187.00
    ## TRINITY_DN0_c1_g1     851.85     243.88     630.48     140.90     712.77
    ## TRINITY_DN0_c1_g2      15.40      10.87      18.70      10.13       6.26
    ## TRINITY_DN0_c2_g1    1734.54     357.77    1324.54     321.64    1413.69
    ##                   X22N1brain X22N1heart X22N2brain X22N2heart X22N3brain
    ## TRINITY_DN0_c0_g1     149.10     585.91      36.54     205.47      85.67
    ## TRINITY_DN0_c0_g2     323.82     517.87     263.39     385.22     252.29
    ## TRINITY_DN0_c0_g3     222.39     224.35     129.88     107.58     176.78
    ## TRINITY_DN0_c1_g1     224.97    1036.06     290.10    1448.48     196.73
    ## TRINITY_DN0_c1_g2      71.05      81.54      50.08      50.52      59.12
    ## TRINITY_DN0_c2_g1     381.05    2025.81     311.91    1496.57     442.51
    ##                   X22N3heart X22N4brain X22N4heart X22N5brain X22N5heart
    ## TRINITY_DN0_c0_g1     275.16      79.58     208.93     131.50     425.49
    ## TRINITY_DN0_c0_g2     489.11     322.50     510.96     268.99     464.07
    ## TRINITY_DN0_c0_g3     143.30     117.02     168.23     152.96     211.98
    ## TRINITY_DN0_c1_g1     782.52     441.63    1358.12     453.26    1330.39
    ## TRINITY_DN0_c1_g2      60.68      30.58      52.32      57.65      91.20
    ## TRINITY_DN0_c2_g1    1319.18     291.43    1222.18     328.01    1717.70

``` r
#rownames(original)<-original[,1]
#original<-original[,-1]

colnames(original) <- c ("10L1brain", "10L1heart", "10L2brain", "10L2heart", "10L3brain", "10L3heart", "10L4brain", "10L4heart", "10L5brain", "10L5heart", "10N1brain", "10N1heart", "10N2brain", "10N2heart", "10N3brain","10N3heart", "10N4brain", "10N4heart", "10N5brain", "10N5heart", "22L1brain", "22L1heart", "22L2brain", "22L2heart", "22L3brain", "22L3heart", "22L4brain", "22L4heart", "22L5brain", "22L5heart", "22N1brain", "22N1heart", "22N2brain", "22N2heart", "22N3brain", "22N3heart", "22N4brain", "22N4heart", "22N5brain", "22N5heart")

# subsetting data to target specific tissues and species, this will change depending on what you are looking at
subset_datalapillum <- original[, grepl("L", colnames(original))]
subset_datalapillumbrain <- subset_datalapillum[, grepl("brain",colnames(subset_datalapillum))]

# creating temperature variable for comparison
sampleslapbrain <- cbind(colnames(subset_datalapillumbrain), rep(c("10", "22"), each = 5))
rownames(sampleslapbrain) <- sampleslapbrain[,1]
sampleslapbrain <- as.data.frame(sampleslapbrain[, -1])
colnames(sampleslapbrain) <- c("treatment")
sampleslapbrain$treatment <- factor(rep(c("10", "22"), each = 5))
```

# Basic filtering

``` r
des <- model.matrix(~sampleslapbrain$treatment)
dge <- DGEList(counts = subset_datalapillumbrain)
keep <- filterByExpr(subset_datalapillumbrain, des)
dge <- dge[keep, keep.lib.sizes=TRUE]

dge <- calcNormFactors(dge)
dim(dge)
```

    ## [1] 62074    10

# Exploratory analyses and visualisation after filtering

``` r
# checking number of reads per gene
quantile(apply(dge$counts,1,sum)/10,seq(0,1,0.05))
```

    ##           0%           5%          10%          15%          20%          25% 
    ##      6.40000     11.50000     13.80000     16.34495     19.51860     23.45450 
    ##          30%          35%          40%          45%          50%          55% 
    ##     28.66030     35.34365     43.70920     55.06585     70.25650     90.70965 
    ##          60%          65%          70%          75%          80%          85% 
    ##    118.97000    158.52860    214.36820    290.17000    403.51720    581.31905 
    ##          90%          95%         100% 
    ##    898.49860   1644.91080 441798.87600

``` r
hist(apply(dge$counts,1,sum)/10, xlim = c(0,2000), breaks = 100000, main = "number of reads per gene per sample", xlab = "number of reads per sample", ylab = "frequency")
```

![](DE_Analysis_hybrid_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# look at the most expressed genes
which(apply(dge$counts,1,sum)/10>1e5)
```

    ## TRINITY_DN1100_c0_g1 TRINITY_DN2141_c0_g2 TRINITY_DN224_c14_g1 
    ##                 2626                21845                23199 
    ##  TRINITY_DN285_c1_g1 TRINITY_DN3350_c1_g1  TRINITY_DN340_c0_g1 
    ##                30187                34919                35360 
    ##  TRINITY_DN340_c1_g2 
    ##                35368

# Performing the DE analysis

``` r
logCPM <- cpm(dge, log=TRUE)
design <- model.matrix(~sampleslapbrain$treatment)

v <- voom(dge, design, plot=TRUE, normalize="quantile")
```

![](DE_Analysis_hybrid_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# visualising clustering of samples
plotMDS(v, pch = 19, col = rep(c( "#377EB8","#E41A1C"), each = 5), cex = 2, frame = F, main = "", cex.axis = 1)
```

![](DE_Analysis_hybrid_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
#legend(x=-2.3, y=2.4, legend = c("10", "22"), col = c("#377EB8","#E41A1C"), pch = 19, cex = 0.7, title = "Temperature treatment (°C), ")
```

10 is blue, 22 is red

``` r
# over/ underexpressed gene counts

fit <- lmFit(v, design)
fit <- eBayes(fit, trend=TRUE)
results <- topTable(fit, coef=ncol(design),n=Inf)

sum(results$adj.P.Val<0.05)
```

    ## [1] 12079

``` r
DE_results <- results[results$adj.P.Val<0.05,]
DE_counts<-dge[which(rownames(dge)%in%rownames(results)),]$counts

#Also take account logFC

DEdown<-DE_results[which(DE_results$logFC<(-1)),]
DEup<-DE_results[which(DE_results$logFC>(1)),]
DE_results<-rbind(DEdown,DEup)
```

# Creating the figures

``` r
# volcano plot

results <- results %>%  mutate(
  Expression = case_when(logFC >= 1 & adj.P.Val <= 0.05 ~ "Up-regulated",
  logFC <= -1 & adj.P.Val <= 0.05 ~ "Down-regulated",
  TRUE ~ "Unchanged"))

ggvolcano <- ggplot(results, aes(logFC,-log10(adj.P.Val))) +
  geom_point(aes(color = Expression), size = 3/5) +
  scale_color_manual(values = c(  "#E41A1C","gray50","#377EB8")) +
  theme_minimal() + theme(legend.text=element_text(size=6)) +
  xlab(expression("log"[2]*"FC")) + 
  ylim(0,6.5) +
  xlim(-8,8) +
  ylab(expression("-log"[10]*"FDR")) +
  guides(colour = guide_legend(override.aes = list(size=6),legend.key.size=6,title="")) 
```

    ## Warning in guide_legend(override.aes = list(size = 6), legend.key.size = 6, : Arguments in `...` must be used.
    ## ✖ Problematic argument:
    ## • legend.key.size = 6
    ## ℹ Did you misspell an argument name?

``` r
ggvolcano
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](DE_Analysis_hybrid_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

# Creating tables for Gene Ontology analyses

``` r
# Printing DE results

final_results_DE <- merge(DE_results, DE_counts, by=0)
colnames(final_results_DE)[1] <- "GeneID"
write.csv(final_results_DE, "lapbrain_DE_results.csv", row.names = T, col.names = T)
```

    ## Warning in write.csv(final_results_DE, "lapbrain_DE_results.csv", row.names =
    ## T, : attempt to set 'col.names' ignored

``` r
# GO writing

write(rownames(DE_results), "GOlapbrainallDE.txt", sep="n")

write.table(rownames(DE_results) [which(DE_results$logFC>0)], "GOlapbrain_upregulated.txt", row.names=F, col.names=F, quote=F)
write.table(rownames(DE_results) [which(DE_results$logFC<0)], "GOlapbrain_downregulated.txt", row.names=F, col.names=F, quote=F)

write(rownames(dge$counts), "GOlapbrain_background.txt", sep="\n")
write(rownames(dge$counts)[-which(rownames(dge$counts)%in%rownames(DE_results))], "GOla-brain_background_no_overlap.txt", sep="n")

final_results<-merge(results, dge$counts, by=0)
colnames(final_results)[1]<-"GeneID"
write.table(final_results, "lapbrain_ALLgenes_results.txt", row.names=T, col.names=T, sep="\t")
```
