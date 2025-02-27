# HMCC
**HMCC-package**: A hybrid model for classifying cancer phenotype. HMCC used for identifying cancerous samples based on the combination of statistical tests and ML algorithms.
**HMCC-r-package** provides a user-friendly R-package interface. 

```r
#Installation process

# development version from GitHub:
require("devtools")
devtools::install_github("BandhanSarker/HMCC")
```

## Usage
```r
#Identify key features/genes
GS(Data,n_groups, DE_ANOVA = TRUE, DE_KW = TRUE,DE_Ttest = TRUE, DE_Sam = TRUE, DE_limma = TRUE,
  P = 0.05, Adjp = "none", N_top = 15)
```
## Arguments
 
* **Data**   Microarray data matrix, whose row contain features/genes and column contain samples
* **n_groups**  Number of groups eg. n1=n2=3; n_groups <- c(3,3)
* **DE_ANOVA**  Analysis of variance (ANOVA) by default TRUE
* **DE_KW**    Kruskal–Wallis (KW) test) by default TRUE
* **DE_Ttest** T-test) by default TRUE
* **DE_Sam** Significance analysis of microarrays (SAM) by default TRUE
* **DE_limma** Linear models for microarray data (LIMMA) by default TRUE
* **P** pvalue by default=0.05
* **Adjp** adjusted correction method see also: p.adjust
* **N_top** Number of significant features/genes selection by default 15
 
## Value 
Top significant features/genes matrix

## Examples
```r
top.df<-HMCC::GS(bladder,n_groups =c(3,9))
col <- colorRampPalette(c("green", "black", "red"))(100)
heatmap(as.matrix(t(top.df[,-16])),scale = "none",col = col,Rowv = NA,Colv=NA)

#Classification
attach(top.df)

#SVM
library(e1071)
model.svm <- svm(DataLabels ~ ., data = top.df)
summary(model.svm)
x <- subset(top.df, select = -DataLabels)
pred.svm <- predict(model.svm, x,type = "class")
table(predict = pred.svm, actual =DataLabels)

#random forest
library(randomForest)
model.rf <- randomForest(DataLabels~.,data=top.df)
pred.rf <- predict(model.rf, x,type = "class")
table(predict = pred.rf, actual =DataLabels)

```
 
## Author(s): 
Md. Matiur Rahaman, Bandhan Sarker, Muhammad Habibulla Alamin, and Farzana Ferdousi

## References
