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
#' #Step 2: Load example data and preprocess
#' Mset<-minfi::preprocessRaw(FlowSorted.DLPFC.450k)
#' Examples_Betas<-minfi::getBeta(Mset)
#' #Step 3: Run HiBED and show results
#' HiBED_result<-HiBED_deconvolution(Examples_Betas, h=2)
#' head(HiBED_result)
#'
#' @param
#' Beta Methylation beta in the format of matrix or data frame or Mset or
#' SummarizedExperiment from brain samples.
#'
#' @param
#' h Numeric variable.
#' Specify the layer of deconvolution in the hierarchical model. Default is 2.
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
#' @import  AnnotationHub
#'
#' @import  SummarizedExperiment
#'
#' @importFrom minfi preprocessRaw
#'
#' @importFrom minfi getBeta
#'
#' @importFrom utils data
#'
#' @export
#'

HiBED_deconvolution <- function(Beta, h=2){
  data_env <- new.env(parent = emptyenv())
  data("HiBED_Libraries", envir = data_env, package = "HiBED")
  HiBED_Libraries<-data_env[["HiBED_Libraries"]]
  Library_Layer1<-as.data.frame(
    assay(HiBED_Libraries$Library_Layer1,"counts"))
  Library_Layer2A<-as.data.frame(
    assay(HiBED_Libraries$Library_Layer2A,"counts"))
  Library_Layer2B<-as.data.frame(
    assay(HiBED_Libraries$Library_Layer2B,"counts"))
  Library_Layer2C<-as.data.frame(
    assay(HiBED_Libraries$Library_Layer2C,"counts"))

  if ((!is(h, "numeric"))) {
    stop(strwrap(sprintf(
      "object is of class '%s', but needs to be of
                                class 'numeric' to use this function",
      class(h)
    ),
    width = 80, prefix = " ",
    initial = ""
    ))
  }
  if (is(Beta, "MethylSet")) {
    Beta <- getBeta(Beta)
  }
  if ((!is(Beta, "matrix")) && (!is(Beta, "data.frame")) &&
      (!is(Beta, "SummarizedExperiment"))) {
    stop(strwrap(sprintf(
      "object is of class '%s', but needs to be of
                                class 'matrix' 'data.frame' or
                                'SummarizedExperiment' to use this function",
      class(Beta)
    ),
    width = 80, prefix = " ",
    initial = ""
    ))
  }

  if (is(Beta,'SummarizedExperiment')){
    Beta <- assay(Beta)
  }


  proj1<-as.data.frame(projectCellType_CP(Beta[rownames(Library_Layer1[
    which(rownames(Library_Layer1) %in% rownames(Beta)),]),],
                                          as.matrix(Library_Layer1[
                                            which(rownames(Library_Layer1)
                                                  %in% rownames(Beta)),]),
    lessThanOne = TRUE))
  proj1[proj1<1e-05]<-0


  proj2A<-as.data.frame(projectCellType_CP(Beta[
    rownames(Library_Layer2A[which(rownames(Library_Layer2A)
                                   %in% rownames(Beta)),]),],
                                           as.matrix(Library_Layer2A[
                                             which(rownames(Library_Layer2A)
                                                   %in% rownames(Beta)),]),
    lessThanOne = TRUE))
  proj2A[proj2A<1e-05]<-0

  for (i in seq_len(nrow(proj2A))) {
    z<-1/sum(proj2A[i,])
    proj2A[i,]<-z*proj2A[i,]
  }

  proj2B<-as.data.frame(projectCellType_CP(Beta[rownames(Library_Layer2B[
    which(rownames(Library_Layer2B) %in% rownames(Beta)),]),],
                                           as.matrix(Library_Layer2B[
                                             which(rownames(Library_Layer2B)
                                                   %in% rownames(Beta)),]),
    lessThanOne = TRUE))
  proj2B[proj2B<1e-05]<-0

  for (i in seq_len(nrow(proj2B))) {
    z<-1/sum(proj2B[i,])
    proj2B[i,]<-z*proj2B[i,]
  }


  proj2C<-as.data.frame(projectCellType_CP(Beta[rownames(Library_Layer2C[
    which(rownames(Library_Layer2C) %in% rownames(Beta)),]),],
                                           as.matrix(Library_Layer2C[
                                             which(rownames(Library_Layer2C)
                                                   %in% rownames(Beta)),]),
    lessThanOne = TRUE))
  proj2C[proj2C<1e-05]<-0

  for (i in seq_len(nrow(proj2C))) {
    z<-1/sum(proj2C[i,])
    proj2C[i,]<-z*proj2C[i,]
  }


  proj<-proj1
  h1_proj<-proj
  proj2A<-proj2A[,-which(colnames(proj2A) %in% c("Glial","Neuronal")),
                 drop=FALSE]
  proj2A<-proj2A/rowSums(proj2A)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Endothelial and Stromal"))],
              proj[, c("Endothelial and Stromal")] * proj2A)
  h2A_proj<-proj

  proj2B<-proj2B[,-which(colnames(proj2B) %in% c("Endothelial and Stromal",
                                                 "Neuronal")),
                 drop=FALSE]
  proj2B<-proj2B/rowSums(proj2B)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Glial"))],
              proj[, c("Glial")] * proj2B)
  h2B_proj<-proj

  proj2C<-proj2C[,-which(colnames(proj2C) %in% c("Endothelial and Stromal",
                                                 "Glial")),
                 drop=FALSE]
  proj2C<-proj2C/rowSums(proj2C)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Neuronal"))],
              proj[, c("Neuronal")] * proj2C)
  h2C_proj<-proj

  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))

  proj[is.nan.data.frame(proj)]<-0


  proj$SumValue<-round(rowSums(proj),2)

  proj_low<-proj %>% filter(SumValue<1)
  ID_low<-rownames(proj_low)

  proj<-proj[,!colnames(proj)=="SumValue"]

  proj[ID_low,]<- h2C_proj[ID_low,]

  proj$Neuronal<-h1_proj$Neuronal
  proj$Glial<-h1_proj$Glial
  proj$`Endothelial and Stromal`<-h1_proj$`Endothelial and Stromal`

  proja<-proj[!rownames(proj)%in%ID_low,]
  proja[is.nan.data.frame(proja)]<-0
  proj[rownames(proja),]<-proja

  if(h==1){
    output<-proj[,c("Endothelial and Stromal","Glial","Neuronal")]
  }else{output<-proj[,c("Endothelial","Stromal","Astrocyte","Microglial",
                        "Oligodendrocyte","GABA","GLU")]
  }

  output$SumValue<-rowSums(output)

  output_low<-output %>% filter(SumValue == "NaN")
  ID_low<-rownames(output_low)
  if(length(ID_low)!=0){
    message("NaN indicates noisy deconvolution signal, thus removed")}
  return(output[,!colnames(output)=="SumValue"]*100)
}
