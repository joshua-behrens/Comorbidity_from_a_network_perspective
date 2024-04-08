library(graphicalVAR)
library(Amelia)
library(abind)
library(dgof)
library(dplyr)
load("momentary.Rda")

#Specify variables
ids <- unique(momentary$Name)
dep_vars <- c("sad", "guilty", "hopeless",
              "lonely", "look_forward", "satisf_self", "enjoy_activity_now")
anx_vars <- c("anxious", "overwhelmed", "paranoid",
              "trembling", "racing_heart", "control_loss", "concentration",
              "worry")
vars <- c(dep_vars, anx_vars)
dayvar <- "nday"
beepvar <- "beep"
filtervars <- c(vars, dayvar, beepvar, "Name")
filtered_momentary <- momentary[filtervars]
ids_to_remove <- c()#add ids to remove for compliance here

#Removing participants with <50% compliance
ids <- setdiff(ids, ids_to_remove)

#Imputing missing values
amelia_results <- amelia(x=filtered_momentary, idvars="Name", ts=beepvar, ords=vars)
imp1 <- amelia_results$imputations[[1]]
imp2 <- amelia_results$imputations[[2]]
imp3 <- amelia_results$imputations[[3]]
imp4 <- amelia_results$imputations[[4]]
imp5 <- amelia_results$imputations[[5]]
temp_array <- abind(imp1, imp2, imp3, imp4, imp5, along=3)
imputed_momentary <- apply(temp_array, 1:2, mean)
imputed_momentary <- as.data.frame(imputed_momentary)

#Correlation matrix
correlation_momentary <- momentary[vars]
correlation_momentary <- na.omit(correlation_momentary)
correlation_matrix <- cor(correlation_momentary)
correlation_matrix <- as.data.frame(correlation_matrix)

#Testing multivariate normality assumption
normality_assumption <- list()
for (v in 1:length(vars)) {
ksres <- ks.test(filtered_momentary[[v]], ecdf(1:3), 3, 2)
normality_assumption <- append(normality_assumption, ksres$p.value)
}

#Testing stationarity assumption
imputed_momentary <- imputed_momentary %>% group_by(Name) %>% 
  mutate(time = dplyr::row_number())
detrended_momentary <- as.data.frame(imputed_momentary)
detrends <- as.data.frame(
  matrix(0, nrow = length(unique(detrended_momentary$Name)), ncol = length(vars)))
for (i in 1:length(ids)) {
  pp <- detrended_momentary[detrended_momentary$Name == ids[i],]
  for (v in 1:length(vars)) {
    ff <- as.formula(paste0(vars[[v]], " ~ time"))
    fit <- lm(ff, data = pp)
    if(anova(fit)$P[1] < 0.05) {
      detrends[i,v] <- 1
      pp[[vars[v]]][!is.na(pp[[vars[[v]]]])] <- residuals(fit)
      detrended_momentary[detrended_momentary$Name == 
                      ids[i],vars[v]] <- pp[,vars[v]]
    } else {
      detrends[i,v] <- 0
    }
  }
}

#Create objects to store effects
array_pdc <- array(0,dim=c(15,15,length(ids)))
array_pcc <- array(0,dim=c(15,15,length(ids)))

#Calculate Network per participant
for(i in 1:length(ids)){
  print(i)
  data <- detrended_momentary[detrended_momentary$Name == ids[i],]
  print(unique(data$Name))
  model <- graphicalVAR(data, vars=vars, dayvar=dayvar, beepvar=beepvar)
  
  #Saving image of network
  filename <- paste("networks/",ids[i],".pdf", sep="")
  pdf(filename)
  plot(model)
  dev.off()

  #Storing effects
  array_pcc[,,i] <- model$PCC #contemporaneous
  array_pdc[,,i] <- model$PDC #temporal
}

#Saving to RDS files
saveRDS(array_pcc, "contemporaneous_individual_array")
saveRDS(array_pdc, "temporal_individual_array")
array_pcc <- readRDS("contemporaneous_individual_array")
array_pdc <- readRDS("temporal_individual_array")

#Array analysis
#Splitting arrays into depression and anxiety arrays
#PCC
dep_array_pcc <- array_pcc[-c(8:15),-c(8:15),]
anx_array_pcc <- array_pcc[-c(1:7),-c(1:7),]
inter_array_pcc <- array_pcc[-c(1:7),-c(8:15),]

#PDC
dep_array_pdc <- array_pdc[-c(8:15),-c(8:15),]
anx_array_pdc <- array_pdc[-c(1:7),-c(1:7),]
inter_array_pdc <- array_pdc[-c(1:7),-c(8:15),]

#Averaging triangles of array
#PCC
num_slices <- dim(dep_array_pcc)[3]
for(s in length(num_slices)) {
  dep_array_pcc[,,s] <- (dep_array_pcc[,,s] + t(dep_array_pcc[,,s]))/2
}
num_slices <- dim(anx_array_pcc)[3]
for(s in length(num_slices)) {
  anx_array_pcc[,,s] <- (anx_array_pcc[,,s] + t(anx_array_pcc[,,s]))/2
}

#PDC
num_slices <- dim(dep_array_pdc)[3]
for(s in length(num_slices)) {
  dep_array_pdc[,,s] <- (dep_array_pdc[,,s] + t(dep_array_pdc[,,s]))/2
}
num_slices <- dim(anx_array_pdc)[3]
for(s in length(num_slices)) {
  anx_array_pdc[,,s] <- (anx_array_pdc[,,s] + t(anx_array_pdc[,,s]))/2
}

#Removing upper triangles and diagonals
#PCC
num_slices <- dim(dep_array_pcc)[3]
for(s in length(num_slices)) {
  dep_array_pcc[,,s][upper.tri(dep_array_pcc[,,s][,], diag=T)] <- 0
}
num_slices <- dim(anx_array_pcc)[3]
for(s in length(num_slices)) {
  anx_array_pcc[,,s][upper.tri(anx_array_pcc[,,s][,], diag=T)] <- 0
}

#PDC
num_slices <- dim(dep_array_pdc)[3]
for(s in length(num_slices)) {
  dep_array_pdc[,,s][upper.tri(dep_array_pdc[,,s][,], diag=F)] <- 0
}
num_slices <- dim(anx_array_pdc)[3]
for(s in length(num_slices)) {
  anx_array_pdc[,,s][upper.tri(anx_array_pdc[,,s][,], diag=F)] <- 0
}

#Counting edges
#PCC
edges_dep_pcc <- (length(dep_array_pcc[dep_array_pcc < -0.123]) + 
                    length(dep_array_pcc[dep_array_pcc > 0.123]))
edges_anx_pcc <- (length(anx_array_pcc[anx_array_pcc < -0.123]) + 
                    length(anx_array_pcc[anx_array_pcc > 0.123]))
edges_inter_pcc <- (length(inter_array_pcc[inter_array_pcc < -0.123]) + 
                    length(inter_array_pcc[inter_array_pcc > 0.123]))

#PDC
edges_dep_pdc <- (length(dep_array_pdc[dep_array_pdc < -0.123]) + 
                    length(dep_array_pdc[dep_array_pdc > 0.123]))
edges_anx_pdc <- (length(anx_array_pdc[anx_array_pdc < -0.123]) + 
                    length(anx_array_pdc[anx_array_pdc > 0.123]))
edges_inter_pdc <- (length(inter_array_pdc[inter_array_pdc < -0.123]) + 
                      length(inter_array_pdc[inter_array_pdc > 0.123]))

#Goodness-of-fit test
#PCC
pcc_testresult_1 <- prop.test(c(edges_dep_pcc, edges_anx_pcc),
                              c((21*length(ids)), (28*length(ids))),
                              p = NULL, alternative = "two.sided",correct = TRUE)
pcc_testresult_2 <- prop.test(c(edges_anx_pcc, edges_inter_pcc),
                              c((28*length(ids)), (56*length(ids))),
                              p = NULL, alternative = "two.sided",correct = TRUE)
pcc_testresult_3 <- prop.test(c(edges_dep_pcc, edges_inter_pcc),
                              c((21*length(ids)), (56*length(ids))),
                              p = NULL, alternative = "two.sided",correct = TRUE)
#PDC
pdc_testresult_1 <- prop.test(c(edges_dep_pdc, edges_anx_pdc),
                              c((21*length(ids)), (28*length(ids))),
                              p = NULL, alternative = "two.sided",correct = TRUE)
pdc_testresult_2 <- prop.test(c(edges_anx_pdc, edges_inter_pdc),
                              c((28*length(ids)), (56*length(ids))),
                              p = NULL, alternative = "two.sided",correct = TRUE)
pdc_testresult_3 <- prop.test(c(edges_dep_pdc, edges_inter_pdc),
                              c((21*length(ids)), (56*length(ids))),
                              p = NULL, alternative = "two.sided",correct = TRUE)

#Calculating percentages
#PCC
perc_dep_pcc <- (edges_dep_pcc / (21*length(ids))) * 100
perc_anx_pcc <- (edges_anx_pcc / (28*length(ids))) * 100
perc_inter_pcc <- (edges_inter_pcc / (56*length(ids))) * 100

#PDC
perc_dep_pdc <- (edges_dep_pdc / (28*length(ids))) * 100
perc_anx_pdc <- (edges_anx_pdc / (36*length(ids))) * 100
perc_inter_pdc <- (edges_inter_pdc / (56*length(ids))) * 100

#Calculating results and placing in dataframe
results <- data.frame(first_column  = c("",
                                        "% of edges withing MDD",
                                        "% of edges within GAD",
                                        "% of edges between", "", 
                                        "p-value: Edges within MDD ~ Edges within GAD", 
                                        "p-value: Edges within GAD ~ Edges between", 
                                        "p-value: Edges between ~ Edges within MDD"),
                      second_column = c("Contemporaneous", perc_dep_pcc, 
                                        perc_anx_pcc, perc_inter_pcc, "", 
                                        pcc_testresult_1$p.value,
                                        pcc_testresult_2$p.value, 
                                        pcc_testresult_3$p.value),
                      third_column = c("Temporal", perc_dep_pdc, perc_anx_pdc, 
                                       perc_inter_pdc, "", 
                                       pdc_testresult_1$p.value, 
                                       pdc_testresult_2$p.value, 
                                       pdc_testresult_3$p.value))
View(results)
print("View any of the following additional results by calling that object: correlation_matrix, normality_assumption.")
