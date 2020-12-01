# MeRIPtools
Tool sets to analyze high throughput data for RNA modifications

### Install the R package from Github

Depends: GenomicFeatures, Rsamtools, ggplot2, doParallel, foreach,grid,rtracklayer,GenomicAlignments,reshape2,Rcpp,RcppArmadillo,
Guitar, stringr,vcfR,gamlss, broom, DESeq2

	install.packages("devtools")
	library(devtools)
	install_github("scottzijiezhang/MeRIPtools")
	library("MeRIPtools")

## Manual page

Please refer to [manual page](https://scottzijiezhang.github.io/MeRIPtoolsManual/) for detailed instructions.  

## Citation 
MeRIPtools is a tool sets that implemented functions for peak calling, QTL calling, differential methylation analysis, visualization. 

If you used MeRIPtools in your publication, please cite:
*Something to be published...*

**Note** MeRIPtools also have wrapper functions to call functions from other R packages to do specific analysis.  
If you used the `plotMetaGene` or `MetaGene` function, please cite the original R package [`Guitar`](https://bioconductor.org/packages/release/bioc/html/Guitar.html)  
Cui X, Wei Z, Zhang L, Liu H, Sun L, Zhang s, Huang Y, Meng J (2016). “Guitar: an R/Bioconductor package for gene annotation guided transcriptomic analysis of RNA related genomic features.” BioMed Research International. 


****************

# 代码解读

## 1. reads counting
该软件首先将基因区域分成 bins,然后统计每个 bin 上的 counts。

可以通过 `MeRIPtools::countReads` 函数获得 `MeRIP` 对象，对象信息如下

* MeRIP@geneBins    
包含两列的数据框，行名为`基因名.bin`    
  - 基因ID    
  - 根据事先设定的 bin 大小得到的 bin 终止长度    
 
  e.g.
```
                                     gene bin
ENSG00000000003.10,20  ENSG00000000003.10  20
ENSG00000000003.10,60  ENSG00000000003.10  60
ENSG00000000003.10,99  ENSG00000000003.10  99
ENSG00000000003.10,139 ENSG00000000003.10 139
ENSG00000000003.10,178 ENSG00000000003.10 178
ENSG00000000003.10,218 ENSG00000000003.10 218
```

## 2.peak calling
在 `method.R` 文件中，函数 `callPeakFisher` 对每一个 bin 都进行了 fisher's exact test。

```
fisher_result <- t( mapply(.fisher_exact_test, batch_ip[,j], batch_input[,j], overall_ip[j], overall_input[j]) )
```
其中，
* batch_ip[,j] 表示某个基因的第 j 个 bin 在 IP 中的 reads count
* overall_ip[j] 表示某个基因所有 bin 在 IP 中 reads count 的中值
* input 同理

### 改进

由于该软件并没有输出检验之后的关键信息，例如 p-value/fdr 等，因此，这里我添加了更多信息 
  —— 查看 `callPeakFisher` 函数，另外相应地为 `MeRIP.peak` 类增加了一个 slot `peakCallResultMore`

