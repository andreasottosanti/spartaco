## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = F----------------------------------------------------------------
library(spartaco)
library(ggplot2)
library(scry)

## ---- eval=F------------------------------------------------------------------
#  library(spartaco)
#  library(spatialLIBD)
#  library(scry)
#  library(ggplot2)
#  
#  ehub <- ExperimentHub::ExperimentHub()

## ---- eval=F------------------------------------------------------------------
#  if (!exists("sce")) sce <- fetch_data(type = "sce", eh = ehub) # SingleCellExperiment
#  if (!exists("spe")) spe <- fetch_data(type = "spe", eh = ehub) # SpatialExperiment

## ---- eval = F----------------------------------------------------------------
#  spe <- spe[,colData(spe)$sample_id == "151673"]

## ---- eval=F------------------------------------------------------------------
#  spe <- spe[which(rowSums(counts(spe)) >= 300),]

## ---- eval = F----------------------------------------------------------------
#  spe <- devianceFeatureSelection(spe, assay="counts", sorted=TRUE)
#  DeviancePlot <- data.frame(x = 1:nrow(spe),
#                             y = rowData(spe)$binomial_deviance)
#  ggplot(DeviancePlot, aes(x, y))+geom_line()+
#    labs(x = "ranked genes", y = "deviance", color = "")+theme_classic()+
#    geom_vline(aes(xintercept = 200, color = "Ideal cut-off"), lty = 2)+
#    geom_vline(aes(xintercept = 500, color = "Actual cut-off"), lty = 2)+
#    scale_color_manual( values = c('Ideal cut-off' = 4, 'Actual cut-off' = 2))+
#    theme(legend.position = "bottom")

## ---- echo = F, fig.width = 6, fig.height=4, fig.align="center"---------------
load("DeviancePlot.Rdata")
ggplot(DeviancePlot, aes(x, y))+geom_line()+
  labs(x = "ranked genes", y = "deviance", color = "")+theme_classic()+
  geom_vline(aes(xintercept = 200, color = "Ideal cut-off"), lty = 2)+
  geom_vline(aes(xintercept = 500, color = "Actual cut-off"), lty = 2)+
  scale_color_manual( values = c('Ideal cut-off' = 4, 'Actual cut-off' = 2))+
  theme(legend.position = "bottom")

## ---- eval=F------------------------------------------------------------------
#  spe <- nullResiduals(spe, assay="counts", fam = "binomial", type = "deviance")
#  spe <- spe[1:500,]

## ---- eval = F----------------------------------------------------------------
#  results <- spartaco_multirun(spe, assay = "binomial_deviance_residuals", K = 2, R = 9, nstart = 5, max.iter = 3000, mc.cores = 5)

## ---- eval = F----------------------------------------------------------------
#  multiresults <- parallel::mclapply(7:12, function(r)
#    spartaco_multirun(data = spe, assay = "binomial_deviance_residuals", K = 2, R = r, nstart = 5, max.iter = 3000, mc.cores = 5), mc.cores = length(7:12))

