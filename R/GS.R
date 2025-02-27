#library()
library(samr)
library(limma)
#######
#' @title Bladder cancer microarray data
#' @description Bladder cancer microarray data
#' @format microarray \code{data_matrix},whose row contains 54675 features and column contains 12 samples, where 3 normal groups and 9 tumor groups
#' @references
#'L. Mengual, M. Burset, E. Ars, J.J. Lozano, H. Villavicencio, M.J. Ribal, A. Alcaraz, DNA Microarray Expression Profiling of Bladder Cancer Allows Identification of Noninvasive Diagnostic Markers, Journal of Urology 182 (2009) 741–748. https://doi.org/10.1016/j.juro.2009.03.084.
#' @examples
#' data(bladder)
"bladder"

#' Identify Key features
#' @description This function is used to identify key features/genes.
#' @param Data  Microarray data matrix, whose row contain features/genes and column contain samples
#' @param n_groups Number of groups eg. n1=n2=3; n_groups <- c(3,3)
#' @param DE_ANOVA Analysis of variance (ANOVA) by default TRUE
#' @param DE_KW Kruskal–Wallis (KW) test) by default TRUE
#' @param DE_Ttest T-test) by default TRUE
#' @param DE_Sam Significance analysis of microarrays (SAM) by default TRUE
#' @param DE_limma Linear models for microarray data (LIMMA) by default TRUE
#' @param P  pvalue by default=0.05
#' @param Adjp adjusted correction method see also: \code{\link[stats]{p.adjust}}
#' @param N_top Number of significant features/genes selection by default 15
#' @return Top significant features/genes matrix
#' @import samr
#' @import limma
#' @import stats
#' @import RColorBrewer
#' @import e1071
#' @import randomForest
#' @seealso \code{\link[e1071]{svm}}, \code{\link[randomForest]{randomForest}}
#' @examples
#' top.df<-HMCC::GS(bladder,n_groups =c(3,9))
#' col <- colorRampPalette(c("green", "black", "red"))(100)
#' heatmap(as.matrix(t(top.df[,-16])),scale = "none",col = col,Rowv = NA,Colv=NA)
#'
#' #Classification
#' attach(top.df)
#'
#' #SVM
#' library(e1071)
#' model.svm <- svm(DataLabels ~ ., data = top.df)
#' summary(model.svm)
#' x <- subset(top.df, select = -DataLabels)
#' pred.svm <- predict(model.svm, x,type = "class")
#' table(predict = pred.svm, actual =DataLabels)
#'
#' #random forest
#' library(randomForest)
#' model.rf <- randomForest(DataLabels~.,data=top.df)
#' pred.rf <- predict(model.rf, x,type = "class")
#' table(predict = pred.rf, actual =DataLabels)
#'
#' @export
GS<-function(Data,n_groups,DE_ANOVA=TRUE,DE_KW=TRUE,DE_Ttest=TRUE,
             DE_Sam=TRUE,DE_limma=TRUE,P=0.05,Adjp="none",N_top=15){
  message("HMCC::INFO:Step:1-Rank Transformation... ...")
  #-----------------Rank Transformation-------------
  #-------------------------------------------------
  rtData<-function (DtMat)
  {
    y<-matrix(0,nr=nrow(DtMat),nc=ncol(DtMat))
    for (ii in 1:nrow(DtMat))
    {
      xx<-DtMat[ii,]
      r=rank(xx)
      y[ii,]<-r
    }
    return(y)
  }
  tdata<-rtData(Data)
  ####
  colnames(tdata)<-colnames(Data)
  rownames(tdata)<-rownames(Data)
  dataT<-as.matrix(tdata)
  #####
  # Function to generate DataLabels dynamically
  generate_data_labels <- function(n_groups) {
    DataLabels <- c()
    for (i in seq_along(n_groups)) {
      DataLabels <- c(DataLabels, rep(i, n_groups[i]))
    }
    return(DataLabels)
  }

  # Call the functions
  DataLabels <- generate_data_labels(n_groups)
  ###########
  message("HMCC::INFO:Step:2-statistical test optimization... ...")

  pValue.ANOVA<-NULL ; pValue.KW<-NULL; pTtest<-NULL;

  for (j1 in 1:dim(Data)[1])
  {
    DataYY<- data.frame(YY =dataT[j1,], FactorLevels = factor(DataLabels))
    if(DE_ANOVA){
      pValue.ANOVA[j1] <-  anova(lm(YY~FactorLevels,data=DataYY))[[5]][1]
    }
    if(DE_KW){
      pValue.KW[j1] <- kruskal.test(YY~FactorLevels, data=DataYY)[[3]]
    }
    if(DE_Ttest){
      pTtest[j1]<-t.test(YY~FactorLevels, data=DataYY)[[3]]
    }
  }
  #---------------------------------
  DE_results <- list()
  #
  if(DE_ANOVA){
    message("HMCC::INFO:ANOVA Done!")
    pValue_ANOVA.adj <- p.adjust(pValue.ANOVA, method=Adjp)#BH
    ANOVA<-which(pValue_ANOVA.adj<P)
    ANOVA.DE<-dataT[ANOVA,]
    DE_results$ANOVA <-rownames(ANOVA.DE)
  }
  #
  if(DE_KW){
    message("HMCC::INFO:KW-test Done!")
    pValue_KW.adj <- p.adjust(pValue.KW, method=Adjp)#BH
    KW<-which(pValue_KW.adj<P)
    KW.DE<-dataT[KW,]
    DE_results$KWtest <-rownames(KW.DE)
  }
  #
  if(DE_Ttest){
    message("HMCC::INFO:t-test Done!")
    pTtest.adj <- p.adjust(pTtest, method=Adjp)#BH
    Ttest<-which(pTtest.adj<P)
    T.DE<-dataT[Ttest,]
    DE_results$ttest <-rownames(T.DE)
  }
  #
  #=================================
  ######## samr
  if(DE_Sam){
    message("HMCC::INFO:SAM running... ...")
    y<-DataLabels# For Unpaired
    data.list <- list(x=dataT, y=y, geneid=rownames(dataT), genenames=rownames(dataT), logged2=F)
    if(length(n_groups)==2){
      SAMP <- samr(data.list, resp.type="Two class unpaired", s0=NULL, s0.perc=NULL, nperms=100, center.arrays=FALSE, testStatistic="wilcoxon")
    }else{
      SAMP <- samr(data.list, resp.type="Multiclass", s0=NULL, s0.perc=NULL, nperms=100, center.arrays=FALSE, testStatistic="standard")
    }
    Pvalue_Sam=samr.pvalues.from.perms(SAMP$tt, SAMP$ttstar)
    Pvalue_Sam.adj <- p.adjust(Pvalue_Sam, method=Adjp)#BH
    sam<-which(Pvalue_Sam.adj<P)
    sam.DE<-dataT[sam,]
    DE_results$sam <-rownames(sam.DE)
  }
  ########################
  #limma
  if(DE_limma){
    message("HMCC::INFO:LIMMA running... ...")
    groupid<-DataLabels
    design <- model.matrix(~ -1+factor(groupid))
    if(length(n_groups)==2){
      colnames(design) <- c("group1", "group2")
    }else{
      dd<-c()
      for(i in seq_along(n_groups)){
        d<-paste0("group",i)
        dd<-c(dd,d)
      }
      colnames(design) <-dd
    }
    limma.fit <- lmFit(dataT, design)
    contrast.matrix <- makeContrasts(group2-group1, levels=design)
    limma.fit2 <- contrasts.fit(limma.fit, contrast.matrix)
    limma.fit2 <- eBayes(limma.fit2)
    Pvalue_limma<- p.adjust(limma.fit2$p.value, method=Adjp)
    limma<-which(Pvalue_limma<P)
    L.DE<-dataT[limma,]
    DE_results$limma <-rownames(L.DE)
  }
  ######################
  ###### Overlap gene selection
  if (length(DE_results) > 1) {
    message("HMCC::INFO:Step3:Overlap gene selection...")
    common_genes <- Reduce(intersect, DE_results)
  } else {
    stop("Required at least two DE method")
  }

  ##############
  DataTT<-Data
  DataTT<-data.frame(DataTT,Pvalue=rep(0,nrow(DataTT)))
  # prepare reData
  if (DE_limma && length(common_genes) > 0) {
    DataTT$Pvalue<-Pvalue_limma
  } else if(DE_Sam && length(common_genes) > 0) {
    DataTT$Pvalue<-Pvalue_Sam.adj
  }else if (DE_KW && length(common_genes) > 0) {
    DataTT$Pvalue<-pValue_KW.adj
  } else if(DE_ANOVA && length(common_genes) > 0) {
    DataTT$Pvalue<-pValue_ANOVA.adj
  } else if (DE_Ttest && length(common_genes) > 0) {
    DataTT$Pvalue<-pTtest.adj
  } else {
    stop("No valid method selected for p-values!")
  }
  #
  df.p<-DataTT[common_genes,]
  df.p<-df.p[order(df.p$Pvalue),][1:200,]
  #For binary class
  message("HMCC::INFO:Step4:filtering with FC ...")
  if(length(n_groups)==2){
    N<-apply(df.p[,(1:n_groups[1])],1,mean)
    C<-apply(df.p[,((n_groups[1]+1):sum(n_groups))],1,mean)
    FC<-abs(N-C)
  }else{
    ####### multiple class
    # Function to calculate pairwise differences and max difference
    calculate_max_diff <- function(data1, n_groups) {
      # Calculate the cumulative sum of group sizes
      cumsum_groups <- cumsum(n_groups)

      # Initialize a list to store group means
      group_means <- list()

      # Calculate means for each group
      for (i in seq_along(n_groups)) {
        if (i == 1) {
          start_col <- 1
        } else {
          start_col <- cumsum_groups[i - 1] + 1
        }
        end_col <- cumsum_groups[i]

        group_means[[i]] <- apply(data1[, start_col:end_col], 1, mean)
      }

      # Calculate pairwise differences
      pairwise_diffs <- list()
      k <- 1
      for (i in seq_along(group_means)) {
        for (j in seq_along(group_means)) {
          if (i < j) {
            pairwise_diffs[[k]] <- abs(group_means[[i]] - group_means[[j]])
            k <- k + 1
          }
        }
      }

      # Combine pairwise differences into a matrix
      midd <- do.call(rbind, pairwise_diffs)

      # Calculate the maximum difference for each row
      mdiff <- apply(midd, 2, max)

      return(mdiff)
    }

    # Call the function
    FC <- calculate_max_diff(df.p, n_groups)
  }
  #####
  mdiff.df<- data.frame(df.p,FC)
  f.df<-mdiff.df[order(mdiff.df$FC, decreasing = TRUE),]
  top<-f.df[1:N_top,][,(1:sum(n_groups))]
  Top<-data.frame(t(top),DataLabels)
  Top$DataLabels<-as.factor(Top$DataLabels)
  message("HMCC::INFO:Computation done!")
  return(Top)
}
#########################################################################################
###########################     End        ##############################################
#########################################################################################
##
