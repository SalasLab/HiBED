#' @title
#' HiBED_deconvolution
#'
#' @name
#' HiBED_deconvolution
#'
#' @description
#' The function estimates proportions up to 7 cell types in brain tissues.
#'
#' @examples
#' #Step 1: Load required libraries
#' library(FlowSorted.Blood.EPIC)
#' library(FlowSorted.DLPFC.450k)
#' library(minfi)
#' #Step 2: Load the library, example data and preprocess
#' data("HiBED_Libraries")
#' Mset<-preprocessRaw(FlowSorted.DLPFC.450k)
#' Examples_Betas<-getBeta(Mset)
#' #Step 3: Run HiBED and show results
#' HiBED_result<-HiBED_deconvolution(Examples_Betas, h=2)
#' head(HiBED_result)
#'
#' @param
#' Beta Methylation beta matrix from brain samples.
#'
#' @param
#' h Specify the layer of deconvolution in the hierarchical model. Default is 2.
#'
#' @return
#' A matrix with predicted cell proportions in brain tissues.
#'
#' @import  FlowSorted.Blood.EPIC
#'
#' @import  dplyr
#'
#' @import  tibble
#'
#' @import  FlowSorted.DLPFC.450k
#'
#' @import  minfi
#'
#' @export

HiBED_deconvolution <- function(Beta, h=2){

  Library_Layer1<-as.data.frame(HiBED_Libraries$Library_Layer1)
  Library_Layer2A<-as.data.frame(HiBED_Libraries$Library_Layer2A)
  Library_Layer2B<-as.data.frame(HiBED_Libraries$Library_Layer2B)
  Library_Layer2C<-as.data.frame(HiBED_Libraries$Library_Layer2C)

  proj1<-as.data.frame(projectCellType_CP(Beta[rownames(Library_Layer1[which(rownames(Library_Layer1) %in% rownames(Beta)),]),],
                                          as.matrix(Library_Layer1[which(rownames(Library_Layer1) %in% rownames(Beta)),]), lessThanOne = TRUE))
  proj1[proj1<1e-05]<-0
  message(nrow(Library_Layer1) - nrow(Library_Layer1[which(rownames(Library_Layer1) %in% rownames(Beta)),]), " out of ",nrow(Library_Layer1),"probes were missing in the Beta matrix for L1")

  proj2A<-as.data.frame(projectCellType_CP(Beta[rownames(Library_Layer2A[which(rownames(Library_Layer2A) %in% rownames(Beta)),]),],
                                           as.matrix(Library_Layer2A[which(rownames(Library_Layer2A) %in% rownames(Beta)),]), lessThanOne = TRUE))
  proj2A[proj2A<1e-05]<-0

  for (i in seq_len(nrow(proj2A))) {
    z<-1/sum(proj2A[i,])
    proj2A[i,]<-z*proj2A[i,]
  }
  message(nrow(Library_Layer2A) - nrow(Library_Layer2A[which(rownames(Library_Layer2A) %in% rownames(Beta)),]), " out of ",nrow(Library_Layer2A),"probes were missing in the Beta matrix for L2A")


  proj2B<-as.data.frame(projectCellType_CP(Beta[rownames(Library_Layer2B[which(rownames(Library_Layer2B) %in% rownames(Beta)),]),],
                                           as.matrix(Library_Layer2B[which(rownames(Library_Layer2B) %in% rownames(Beta)),]), lessThanOne = TRUE))
  proj2B[proj2B<1e-05]<-0

  for (i in seq_len(nrow(proj2B))) {
    z<-1/sum(proj2B[i,])
    proj2B[i,]<-z*proj2B[i,]
  }
  message(nrow(Library_Layer2B) - nrow(Library_Layer2B[which(rownames(Library_Layer2B) %in% rownames(Beta)),]), " out of ",nrow(Library_Layer2B),"probes were missing in the Beta matrix for L2B")

  proj2C<-as.data.frame(projectCellType_CP(Beta[rownames(Library_Layer2C[which(rownames(Library_Layer2C) %in% rownames(Beta)),]),],
                                           as.matrix(Library_Layer2C[which(rownames(Library_Layer2C) %in% rownames(Beta)),]), lessThanOne = TRUE))
  proj2C[proj2C<1e-05]<-0

  for (i in seq_len(nrow(proj2C))) {
    z<-1/sum(proj2C[i,])
    proj2C[i,]<-z*proj2C[i,]
  }
  message(nrow(Library_Layer2C) - nrow(Library_Layer2C[which(rownames(Library_Layer2C) %in% rownames(Beta)),]), " out of ",nrow(Library_Layer2C),"probes were missing in the Beta matrix for L2C")

  proj<-proj1
  h1_proj<-proj
  proj2A<-proj2A[,-which(colnames(proj2A) %in% c("Glial","Neuronal")),
                 drop=FALSE]
  proj2A<-proj2A/rowSums(proj2A)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Endothelial and Stromal"))],
              proj[, c("Endothelial and Stromal")] * proj2A)
  h2A_proj<-proj

  proj2B<-proj2B[,-which(colnames(proj2B) %in% c("Endothelial and Stromal","Neuronal")),
                 drop=FALSE]
  proj2B<-proj2B/rowSums(proj2B)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Glial"))],
              proj[, c("Glial")] * proj2B)
  h2B_proj<-proj

  proj2C<-proj2C[,-which(colnames(proj2C) %in% c("Endothelial and Stromal","Glial")),
                 drop=FALSE]
  proj2C<-proj2C/rowSums(proj2C)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Neuronal"))],
              proj[, c("Neuronal")] * proj2C)
  h2C_proj<-proj

  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))

  proj[is.nan.data.frame(proj)]<-0


  proj$Sum<-round(rowSums(proj),2)

  proj_low<-proj %>% filter(Sum<1)
  ID_low<-rownames(proj_low)
  # empty_list <- vector(mode = "list", length = length(ID_low))
  # names(empty_list)<-ID_low

  proj<-proj[,!colnames(proj)=="Sum"]

  proj[ID_low,]<- h2C_proj[ID_low,]

  proj$Neuronal<-h1_proj$Neuronal
  proj$Glial<-h1_proj$Glial
  proj$`Endothelial and Stromal`<-h1_proj$`Endothelial and Stromal`

  proja<-proj[!rownames(proj)%in%ID_low,]
  proja[is.nan.data.frame(proja)]<-0
  proj[rownames(proja),]<-proja

  if(h=="1"){
    output<-proj[,c("Endothelial and Stromal","Glial","Neuronal")]
  }else{output<-proj[,c("Endothelial","Stromal","Astrocyte","Microglial","Oligodendrocyte","GABA","GLU")]
  }

  output$Sum<-rowSums(output)

  output_low<-output %>% filter(Sum == "NaN")
  ID_low<-rownames(output_low)
  if(length(ID_low)!=0){
    message("NaN indicates noisy deconvolution signal, thus removed")}
  return(output[,!colnames(output)=="Sum"]*100)
}
