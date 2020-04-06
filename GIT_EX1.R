#### SIMPLIFIED MODEL THAT OUTPUTS RESULTS OF DISTANCE MATRICES BASED ON CENTROID AND SYMBOLIC COMPARISONS
#### THIS MODEL WAS VERIFIED AGAINST THE MANUAL CALCULATONS PRODUCED FROM THE COMPANION EXCEL FILE########
##########################################################################################################
set.seed(250)
library(readxl)
library(sqldf)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(ggpubr)
library(RSDA)
library(symbolicDA)
library(MVA)
library(grid)
library(reshape2)
library(lattice)
#######################################################
#########Clear#########################################
rm(list = ls())
cat("\014")
#######################################################
#########SOURCE AND FUNCTIONS #########################
source('C:/Users/Rick/Desktop/R_Auburn/Feb2019/Histo_Fun_Sturges.R')
source('C:/Users/Rick/Desktop/R_Auburn/Feb2019/Histo_Fun_Scott.R')
source('C:/Users/Rick/Desktop/R_Auburn/Feb2019/cprobs_fun.R')
source('C:/Users/Rick/Desktop/R_Auburn/Feb2019/newSobject_fun.R')
source('C:/Users/Rick/Desktop/R_Auburn/Feb2019/Set_Fun.R')
source('C:/Users/Rick/Desktop/R_Auburn/Feb2019/Modal_Fun.R')
source('C:/Users/Rick/Desktop/R_Auburn/Feb2019/Cont_Fun.R')

calculate.probs <- cprobs_fun
newSobject <- newSobject_fun

process.histogram.variable <- Histo_Fun_Sturges
#process.histogram.variable <- Histo_Fun_Scott
process.set.variable <- Set_Fun
process.modal.variable <- Modal_Fun
process.continuum.variable <- Cont_Fun

#########Create DATASET################################
datafile <- "Simul_D1_Code"
#####################################################
excel <- ".xlsx"
path1 <- "C:/Users/Rick/Desktop/R_Auburn/Nov_2018/Dist/"
pathfinal <- paste(path1,datafile, excel, sep = "")
#######################################################
Vendor_Data_Set <- read_excel(pathfinal)
Vendor_Data_Set$Vendor <- as.factor(Vendor_Data_Set$Group)
Vendor_Data_Set = as.data.frame(Vendor_Data_Set)
####################################################################################################################
dist <- stats::dist
C1 <- ncol(Vendor_Data_Set)
# this removes the label (as it is added below)
red_vendor_Set <- Vendor_Data_Set[,-c(4)]
df <- red_vendor_Set
df <- melt(df, id='Vendor')
#### SCALING FUNCTION ##############################################################################################
scaling <- function(x) {
  num <- x
  denom <- max(x)
  return (num/denom)
}

dataTable <- Vendor_Data_Set
concept=c('Vendor')
####################################################################################################################
###########################################CONSTRUCT DISTANCE MATRIX################################################
####################################################################################################################
#CHOOSE ANALYSIS TYPE FROM THE FOUR OPTIONS BELOW. UNCOMMENT THE "variables" and "variable.types" lines to select###
####################################################################################################################

#CENTROIDAL#########################################################################################################
variables=c('V1', 'V2', 'V3')
variables.types=c('$C', '$C', '$C')

#SYMBOLIC###########################################################################################################
#variables=c('V1', 'V2', 'V3')
#variables.types=c('$H', '$H', '$H')

#CENTROIDAL WITH CATEGORICAL########################################################################################
# variables=c('V1', 'V2', 'V3', 'Orange', 'Blue')
# variables.types=c('$C', '$C', '$C', '$C', '$C')

#SYMBOLIC WITH CATEGORICAL##########################################################################################
#variables=c('V1', 'V2', 'V3', 'CAT1' )
#variables.types=c('$H', '$H', '$H', '$M')
####################################################################################################################
if (length(variables) != length(variables.types)) {
  stop("variables and variables.types must have the same length")
}
dataTable <- dataTable[, which(colnames(dataTable) %in% 
                                 c(variables, concept))]
concept <- paste0("[", concept, "]")
variables <- paste0("[", variables, "]")
conceptColumns <- paste(concept, collapse = ", ")

sqldf()


sqldf(paste0("CREATE INDEX main.concept_index ON dataTable (", 
             conceptColumns, ")"))
sym.obj <- sqldf(paste0("SELECT DISTINCT ", conceptColumns, 
                        " FROM main.dataTable ORDER BY ", conceptColumns))
sym.obj.names <- do.call("paste", args = c(sym.obj, sep = "."))
symObjTable <- data.frame(SymObjNames = sym.obj.names)
sqldf("SELECT SymObjNames FROM symObjTable")
meta.data <- list()

for (i in 1:length(variables)) {
  n <- as.numeric(stringr::str_extract(variables.types[[i]], 
                                       "[[:digit:]]"))
  variables.types[[i]] <- stringr::str_replace(variables.types[[i]], 
                                               "[[:digit:]]", "")
  switch(variables.types[[i]], `$C` = {
    meta.data[[i]] <- process.continuum.variable(variables[[i]], 
                                                 conceptColumns)
  }, `$I` = {
    meta.data[[i]] <- process.interval.variable(variables[[i]], 
                                                conceptColumns)
  }, `$H` = {
    meta.data[[i]] <- process.histogram.variable(variables[[i]], 
                                                 concept, dataTable, n)
  }, `$M` = {
    meta.data[[i]] <- process.modal.variable(variables[[i]], 
                                             concept, sym.obj.names)
  }, `$S` = {
    meta.data[[i]] <- process.set.variable(variables[[i]], 
                                           concept, sym.obj.names)
  }, stop("Invalid variable type"))
}

sqldf()

meta.data <- data.frame(meta.data, check.names = F)
rownames(meta.data) <- sym.obj.names

colnames(meta.data)[colnames(meta.data) == "'$C'"] <- "$C"
colnames(meta.data)[colnames(meta.data) == "'$I'"] <- "$I"
colnames(meta.data)[colnames(meta.data) == "'$S'"] <- "$S"
colnames(meta.data)[colnames(meta.data) == "'$M'"] <- "$M"
colnames(meta.data)[colnames(meta.data) == "'$H'"] <- "$H"

new.sym <- newSobject(meta.data)
class(new.sym) <- c("list", "sym.data.table")
new.sym
result.histogram.raw <- new.sym
###############################################################################################
###############################################################################################
hr.dist <- stats::dist(result.histogram.raw[["data"]], method = "euclidean")
hr.dist.scaled <- hr.dist
hr.dist.scaled <- as.numeric(hr.dist.scaled)
hr.dist.scaled <- as.data.frame(hr.dist.scaled)
hr.dist.scaled <- lapply(hr.dist.scaled, scaling)
hr.dist.scaled <- as.data.frame(hr.dist.scaled)
hr.dist.scaled <- as.numeric(hr.dist.scaled$hr.dist.scaled)

r <- nrow(meta.data)
b2 <- matrix(0,r,r)
b2[lower.tri(b2, diag = FALSE)] <- hr.dist.scaled

hr.dist.scaled <- as.dist(b2)

#print(hr.dist, diag = TRUE, digits = 3)
#print(hr.dist.scaled, diag = TRUE, digits = 3)

m2 <- as.matrix(hr.dist.scaled)

m3 <- m2

m3[upper.tri(m3, diag = FALSE)] <- NA
m10<- as.data.frame(m3)
colnames(m10) <- levels(df$Vendor)
rownames(m10) <- levels(df$Vendor)
cf <- format(round(m10, 2), nsmall=2)
cf[is.na(m10)] <- ""
cf


#m2
diag(m2) <- NA

#####################################################################################
c2 <- as.data.frame(colMeans(m2, na.rm = TRUE))
d3 <- median(c2[,1])
median_symbolic <- d3
q1_symbolic <- unname(quantile(c2[,1], .25))
q3_symbolic <- unname(quantile(c2[,1], .75))
iqr_symbolic <- IQR(c2[,1])
upperlimit_symbolic <- q3_symbolic + 1.5*iqr_symbolic
format(round(c2, 2), nsmall=2)

print(median_symbolic)
print(q1_symbolic)
print(q3_symbolic)
print(iqr_symbolic)
print(upperlimit_symbolic)


