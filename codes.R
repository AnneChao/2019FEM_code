#------------------------------------------------------
#Before running the code, you should install and library the following packages:
install.packages("FD")
install.packages("Rcpp")
install.packages("ggplot2")
install.packages("iNEXT")
library(FD)
library(Rcpp)
library(ggplot2)
library(readxl)
library(dplyr)
library(iNEXT)
#------------------------------------------------------
#Read data.
# Data_established = read.csv("matrix_abundance_established_site.csv")
# Data_regenerating = read.csv("matrix_abundance_regenerating_site.csv")
Data_established = read.csv("matrix_abundance_established_site_NO_EXOTICS_26oct2018.csv")
Data_regenerating = read.csv("matrix_abundance_regenerating_site_NO_EXOTICS_26oct2018.csv")

#Read trait matrix.
# trait_established <- read_excel("matrix_traits_FRAGSeRESTS_established_15ago2018.xlsx",sheet = 1,
#                                 col_types = c("text", "text", "text", "numeric", "numeric"))
# trait_regenerating <- read_excel("matrix_traits_FRAGSeRESTS_regenerants_15ago2018.xlsx",sheet = 1,
#                                 col_types = c("text", "text", "text", "numeric", "numeric"))
trait_established <- read_excel("matrix_traits_FRAGSeRESTS_established_NO_EXOTICS_26oct2018.xlsx",sheet = 1,
                                col_types = c("text", "text", "text", "numeric", "numeric"))
trait_regenerating <- read_excel("matrix_traits_FRAGSeRESTS_regenerants_NO_EXOTICS_26oct2018.xlsx",sheet = 1,
                                 col_types = c("text", "text", "text", "numeric", "numeric"))
trait_established <- data.frame(trait_established)
trait_regenerating <- data.frame(trait_regenerating)

#Turn the NA character into NA value; trait_regenerating has no such case.
trait_established$seed_disp[which(trait_established$seed_disp=="NA")]= NA

#For trait, remove the species and trait which is "NA" in the trait: seed_disp (only 1 species).
Data_established = Data_established[-which(is.na(trait_established$seed_disp)),]
trait_established = trait_established[-which(is.na(trait_established$seed_disp)),]

#For traits, If there is "NA" in seed.size or wood.density, I replace it with mean value of other species.
trait_established[which(is.na(trait_established[,4])),4] = mean(trait_established[,4],na.rm = T)
trait_established[which(is.na(trait_established[,5])),5] = mean(trait_established[,5],na.rm = T)

trait_regenerating[which(is.na(trait_regenerating[,4])),4] = mean(trait_regenerating[,4],na.rm = T)
trait_regenerating[which(is.na(trait_regenerating[,5])),5] = mean(trait_regenerating[,5],na.rm = T)

#Regard "cerrado" and "exotic" as the same category,so I replace "cerrado" by "exotic".
# trait_established[trait_established[,2]=="cerrado", 2] = "exotic"
# trait_regenerating[trait_regenerating[,2]=="cerrado", 2] = "exotic"

#Use gower distance to change traits matrix into species-pairwise distance matrix.
D_established = as.matrix(gowdis(trait_established[,-1]))
D_regenerating = as.matrix(gowdis(trait_regenerating[,-1]))

#Remove useless column in data 20181029-yhc
Data_established = Data_established[,-3]
Data_regenerating = Data_regenerating[,-3]
#------------------------------------------------------
#The following two functions are used in the main function FD_MLE.
cppFunction(
  "NumericMatrix empirical(NumericVector ai,NumericMatrix dij,float q,float Q){
  const int S = ai.size();
  NumericMatrix temp(S,S);
  for(int i=0;i<S;i++){
  for(int j = 0;j<S;j++){
  temp(i,j) = dij(i,j)*pow((ai[i]*ai[j]/Q),q);
  }
  }
  return(temp);
  }")
cppFunction(
  "NumericMatrix empirical_q1(NumericVector ai,NumericMatrix dij,float Q){
  const int S = ai.size();
  NumericMatrix temp(S,S);
  for(int i=0;i<S;i++){
  for(int j = 0;j<S;j++){
  temp(i,j) = dij(i,j)*(ai[i]*ai[j]/Q)*log(ai[i]*ai[j]/Q);
  }
  }
  return(temp);
  }")

#------------------------------------------------------------------------------
# Functional diversity (MLE:empirical) See Chiu and Chao (2014, Plos One) paper
#-------------------------------------------------------------------------------
#' FD_MLE(q, data, Dij) is a function which computes FD of order q based on abundance data.
#' @param q a numeric or a vector of diversity order. The suggested range for q is [0, 3].
#' @param data a vector of species sample frequencies.
#' @param Dij a matrix of distance matrix.
#' @return a numerical vector of FD.
FD_MLE = function(q, data ,Dij){
  Xi <- data[data!=0]
  distance <- Dij[data!=0, data!=0]
  a <- Xi/sum(Xi)

  Q = sum(distance*(a %*% t(a)))

  Emp <- function(q){
    if(q==1){
      Empirical = exp(-sum(empirical_q1(a, as.matrix(distance), Q)))
    }else{
      Empirical = sum(empirical(a, as.matrix(distance), q, Q))^(1/(1-q))
    }
    Empirical
  }
  sapply(q, Emp)
}

#-----------------------------------------------------------------------------------------------
# Functional dissimilarity (empirical) See Chiu and Chao (2014, Plos One) paper with modification
#------------------------------------------------------------------------------------------------
#' FD_Beta(q, data, dij, CU, method) is a function which computes functional dissimilarity measure of order q based on abundance data.
#' @param q a numeric or a vector of diversity order; The suggested range for q is [0, 3].
#' @param data a S*N matrix of species sample frequencies with S species (rows), N communities (columns).
#' @param dij a matrix of distance matrix.
#' @param CU a character to choose method, "C" for 1-CqN (Sorensen) ; "U" for 1-UqN (Jaccard).
#' @param method a character to choose method, "relative" or "absolute".
#' @return a numerical vector of functional dissimilarity.
FD_Beta = function(q, data, dij, CU, method){
  if(method == "relative") data <- sapply(1:ncol(data), FUN = function(i) data[,i]/sum(data[,i]))
  N = ncol(data)
  S = nrow(dij)
  m = dij

  M=matrix(rep(t(m), N), ncol = ncol(m), byrow = TRUE)
  M=matrix(rep(M, N), nrow = nrow(M), byrow = F)

  alpha = (1/(N)^2)*FD_MLE(q, unlist(data), M)
  gamma = FD_MLE(q, rowSums(data), dij)
  beta = gamma/alpha
  beta[beta<1] = 1

  MAX.dis = matrix(1, N*S, N*S)
  for(i in 1:N){
    MAX.dis[(1:S)+((i-1)*S), (1:S)+((i-1)*S)] = dij
  }

  MAX.gamma =  FD_MLE(q, unlist(data), MAX.dis)
  MAX.beta = MAX.gamma/alpha
  #MAX.alpha =  FD_MLE(q, unlist(data), MAX.dis)
  #MAX.beta = gamma/MAX.alpha
  if(CU=="C"){
    out = (1-beta^(1-q))/(1-MAX.beta^(1-q))
  }else{
    out = (1-(beta)^(q-1))/(1-MAX.beta^(q-1))
  }
  out[q==1] = log(beta[q==1])/log(MAX.beta[q==1])
  out
}

#--------------------------------------------------------------------
# Plot the output of iNEXT function (See Fig. 3 of 2019 FEM) paper
#--------------------------------------------------------------------
#' p_inext(x, type, se, facet.var, color.var, grey) is a function to plot the output of iNEXT
#' @param x is the output of iNEXT function
p_inext = function (x, type = 1, se = TRUE, facet.var = "order", color.var = "site", grey = FALSE) {
  TYPE <- c(1, 2, 3)
  SPLIT <- c("none", "order", "site", "both")
  if (is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  if (is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var,
                                               SPLIT) == -1)
    stop("invalid facet variable")
  if (is.na(pmatch(color.var, SPLIT)) | pmatch(color.var,
                                               SPLIT) == -1)
    stop("invalid color variable")
  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  if (facet.var == "order")
    color.var <- "site"
  if (facet.var == "site")
    color.var <- "order"
  options(warn = -1)
  z <- fortify(x, type = type)
  z$order=factor(z$order,levels = c("0", "1","2"))
  levels(z$order)<- c("italic(q)==0", "italic(q)==1","italic(q)==2")
  options(warn = 0)
  if (ncol(z) == 7) {
    se <- FALSE
  }
  datatype <- unique(z$datatype)
  if (color.var == "none") {
    if (levels(factor(z$order)) > 1 & "site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object consists multiple sites and orders, change setting as both")
      color.var <- "both"
      z$col <- z$shape <- paste(z$site, z$order, sep = "-")
    }
    else if ("site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as order")
      color.var <- "site"
      z$col <- z$shape <- z$site
    }
    else if (levels(factor(z$order)) > 1) {
      warning("invalid color.var setting, the iNEXT object consists multiple sites, change setting as site")
      color.var <- "order"
      z$col <- z$shape <- factor(z$order)
    }
    else {
      z$col <- z$shape <- rep(1, nrow(z))
    }
  }
  else if (color.var == "order") {
    z$col <- z$shape <- factor(z$order)
  }
  else if (color.var == "site") {
    if (!"site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- z$site
  }
  else if (color.var == "both") {
    if (!"site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- paste(z$site, z$order, sep = "-")
  }
  #my adjustment
  data.sub <- z[which(z$method == "observed"), ]
  z = z[which(z$method != "observed"),]
  z$lty <- factor(z$method, c("interpolated", "extrapolated"),
                  c("interpolation", "extrapolation"))
  z$col <- factor(z$col)


  g <- ggplot(z, aes_string(x = "x", y = "y", colour = "col")) +  theme_bw()+
    geom_point(aes_string(shape = "shape"), size = 5, data = data.sub)+
    scale_color_manual(values = c("#CC0099","darkgreen","#56B4E9"))
  g <- g + geom_line(aes_string(linetype = "lty"), lwd = 1.5) +
    guides(linetype = guide_legend(title = "Method"), colour = guide_legend(title = "Guides"),
           fill = guide_legend(title = "Guides"), shape = guide_legend(title = "Guides")) +
    theme(panel.grid.major = element_line(colour = "gray80"),
          panel.grid.minor = element_line(colour = "gray80"),
          legend.position = "bottom", legend.title = element_blank(),
          legend.key.width = unit(2,"cm"),legend.key.height = unit(1,"cm"),
          text = element_text(size = 18))
  if (type == 2L) {
    g <- g + labs(x = "Number of sampling units", y = "Sample coverage")
    if (datatype == "abundance")
      g <- g + labs(x = "Number of individuals", y = "Sample coverage")
  }
  else if (type == 3L) {
    g <- g + labs(x = "Sample coverage", y = "Species diversity")
  }
  else {
    g <- g + labs(x = "Number of sampling units", y = "Species diversity")
    if (datatype == "abundance")
      g <- g + labs(x = "Number of individuals", y = "Species diversity")
  }
  if (se)
    g <- g + geom_ribbon(aes_string(ymin = "y.lwr", ymax = "y.upr",
                                    fill = "factor(col)", colour = "NULL"), alpha = 0.2)+
    scale_fill_manual(values = c("#CC0099","darkgreen","#56B4E9"))
  if (facet.var == "order") {
    if (length(levels(factor(z$order))) == 1 & type != 2) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple orders.")
    }
    else {
      g <- g + facet_wrap(~order,labeller=label_parsed, nrow = 1)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Guides",
                                              ncol = length(levels(factor(z$order))), byrow = TRUE),
                        fill = guide_legend(title = "Guides"))
      }
    }
  }
  if (facet.var == "site") {
    if (!"site" %in% names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites.")
    }
    else {
      g <- g + facet_wrap(~site, nrow = 1)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Guides",
                                              nrow = length(levels(factor(z$order)))), fill = guide_legend(title = "Guides"))
      }
    }
  }
  if (facet.var == "both") {
    if (length(levels(factor(z$order))) == 1 | !"site" %in%
        names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites or orders.")
    }
    else {
      g <- g + facet_wrap(site ~ order)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Guides",
                                              nrow = length(levels(factor(z$site))), byrow = TRUE),
                        fill = guide_legend(title = "Guides"))
      }
    }
  }
  if (grey) {
    g <- g + theme_bw(base_size = 18) + scale_fill_grey(start = 0,
                                                        end = 0.4) + scale_colour_grey(start = 0.2, end = 0.2) +
      guides(linetype = guide_legend(title = "Method"),
             colour = guide_legend(title = "Guides"), fill = guide_legend(title = "Guides"),
             shape = guide_legend(title = "Guides")) + theme(legend.position = "bottom",
                                                             legend.title = element_blank())
  }
  return(g)
}

#-----------------------------------------------
# Plot the dissimilarity curve
#-----------------------------------------------
#' Get_plot(out1, out2, out3, CU, ER) is a function to plot the output of FD_Beta or beta_diversity_MLE.
#' @param out1 a 3*N matrix, 3 means q=0,1,2 ; N means 2 to N subgroups for "ALL".
#' @param out2 a 3*N matrix, 3 means q=0,1,2 ; N means 2 to N subgroups for "Restoration".
#' @param out3 a 3*N matrix, 3 means q=0,1,2 ; N means 2 to N subgroups for "Fragments".
#' @param CU a character to specify method, "C" for Sorensen type non-overlap measure (1-CqN);
#'                                         "U" for Jaccard type non-overlap measure (1-UqN).
#' @param ER a character to specify data, "E" for Established ; "R" for Regenerating.
#' @param dtype a character to specify the dissimilarity type, "Tax" for taxonmic ; "Fun" for functional
#' @return a plot.
Get_plot = function(out1, out2, out3, CU, ER, dtype){

  result.all = data.frame(N=c(2:32,2:14,2:18),
                          Dissimilarity=c(out1[1,],out2[1,],out3[1,],
                                          out1[2,],out2[2,],out3[2,],
                                          out1[3,],out2[3,],out3[3,]),
                          Type=rep(c("ALL","Restoration","Fragments"),c(31,13,17)),
                          q=rep(c("q=0","q=1","q=2"),each=61))
  #result.all = result.all[(result.all$N %% 2 == 0),]
  result.all$Type=factor(result.all$Type,levels = c("Fragments","Restoration","ALL"))
  result.all$q=factor(result.all$q,levels = c("q=0", "q=1","q=2"))
  levels(result.all$q)<- c("italic(q)==0", "italic(q)==1","italic(q)==2")

  p = ggplot(result.all, aes(x=N,y=Dissimilarity, color=Type, linetype=Type))+ theme_bw()+
    theme(text=element_text(size=25),legend.position="bottom",legend.title = element_blank(),
          legend.key.width = unit(2,"cm"),legend.key.height = unit(1,"cm"),
          axis.text.x = element_text(color = "black"),legend.text=element_text(size=25),
          axis.text.y = element_text(color = "black"),
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour = "gray"))+
    facet_wrap(~q,labeller=label_parsed)+
    scale_linetype_manual(values=c("1111","solid", "3232"),breaks=c("ALL","Fragments","Restoration"))+
    scale_color_manual(values = c("darkgreen","#56B4E9","#CC0099"),breaks=c("ALL","Fragments","Restoration"))+
    geom_smooth(formula = y~log(x) ,se = F,fullrange = F, size=2)
  # geom_line(aes(linetype=Type),size=2)
  # geom_line(data = subset(result.all, Type == "Fragments"), linetype = "1111" ,size = 2)+
  #   geom_line(data = subset(result.all, Type == "Restoration"), size = 2)+
  #   geom_line(data = subset(result.all, Type == "ALL"), size = 2)+
  #panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
  #panel.border = element_rect( colour = "black", size=1)

  if(dtype == "Tax"){
    p = p + ylab("Taxonomic dissimilarity")
  }else{
    p = p + ylab("Functional dissimilarity")
  }

  if(CU=="C"){
    if(ER=="E"){
      # p + labs(title= bquote(paste("Established (", 1-C[qN],")" )))
      p + labs(title= "Established (Sorensen-type)")
    }else{
      # p + labs(title= bquote(paste("Regenerating (",1-C[qN],")" )))
      p + labs(title= "Regenerating (Sorensen-type)")
    }
  }else{
    if(ER=="E"){
      # p + labs(title= bquote(paste("Established (", 1-U[qN],")" )))
      p + labs(title= "Established (Jaccard-type)")
    }else{
      # p + labs(title= bquote(paste("Regenerating (",1-U[qN],")" )))
      p + labs(title= "Regenerating (Jaccard-type)")
    }
  }
}

####Plot Fig 3. in 2019 FEM paper (based on iNEXT output)
data_est = data.frame(ALL = c(32,rowSums(Data_established[,3:34]>0)),
                      Rsetoration = c(14,rowSums(Data_established[,21:34]>0)),
                      Fragments = c(18,rowSums(Data_established[,3:20]>0)))
data_reg = data.frame(ALL = c(32,rowSums(Data_regenerating[,3:34]>0)),
                      Rsetoration = c(14,rowSums(Data_regenerating[,21:34]>0)),
                      Fragments = c(18,rowSums(Data_regenerating[,3:20]>0)))

outpout_est = iNEXT(x = data_est,q = c(0,1,2),datatype = "incidence_freq")
outpout_reg = iNEXT(x = data_reg,q = c(0,1,2),datatype = "incidence_freq")
p_inext(x = outpout_est)
p_inext(x = outpout_reg)

####Plot Appendix S6: Figure S2
#B is the number of combinations of sites for each chosen number of point on x-axis. 
#If the number of combinations exceeds 100, then cut off at B = 100
B = 100                                          
#((a) Sorensen-type functional dissimilarity measure) Compare all, Restoration Fragments base on established, 1-CqN, relative.
result.ALL_established1 = sapply(2:32, function(x) {
  if(choose(n = 32,k = x)>B){
    chosen = replicate(B,sort(sample(1:32,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:32,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:32,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    FD_Beta(c(0,1,2),Data_established[,(chosen[i,]+2)], D_established, "C", "relative")
  }) %>% rowMeans; })
result.R_established1 = sapply(2:14, function(x) {
  if(choose(n = 14,k = x)>B){
    chosen = replicate(B,sort(sample(1:14,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:14,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:14,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    FD_Beta(c(0,1,2),Data_established[,(chosen[i,]+20)], D_established, "C", "relative")
  }) %>% rowMeans; })
result.F_established1 = sapply(2:18, function(x) {
  if(choose(n = 18,k = x)>B){
    chosen = replicate(B,sort(sample(1:18,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:18,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:18,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    FD_Beta(c(0,1,2),Data_established[,(chosen[i,]+2)], D_established, "C", "relative")
  }) %>% rowMeans; })
Get_plot(result.ALL_established1, result.R_established1, result.F_established1, "C", "E", "Fun")

#((a) Sorensen-type functional dissimilarity measure) Compare All, Restoration Fragments base on regenerating, 1-CqN, relative.
result.ALL_regenerating1 <- sapply(2:32, function(x) {
  if(choose(n = 32,k = x)>B){
    chosen = replicate(B,sort(sample(1:32,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:32,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:32,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    FD_Beta(c(0,1,2),Data_regenerating[,(chosen[i,]+2)], D_regenerating, "C", "relative")
  }) %>% rowMeans; })
result.R_regenerating1 = sapply(2:14, function(x) {
  if(choose(n = 14,k = x)>B){
    chosen = replicate(B,sort(sample(1:14,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:14,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:14,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    FD_Beta(c(0,1,2),Data_regenerating[,(chosen[i,]+20)], D_regenerating, "C", "relative")
  }) %>% rowMeans; })
result.F_regenerating1 = sapply(2:18, function(x) {
  if(choose(n = 18,k = x)>B){
    chosen = replicate(B,sort(sample(1:18,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:18,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:18,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    FD_Beta(c(0,1,2),Data_regenerating[,(chosen[i,]+2)], D_regenerating, "C", "relative")
  }) %>% rowMeans; })
Get_plot(result.ALL_regenerating1, result.R_regenerating1, result.F_regenerating1, "C", "R", "Fun")

#### Plot Fig 6. in the main text 
#((a) Jaccard-type functional dissimilarity measure) Compare All, Restoration Fragments base on established, 1-UqN, relative.
result.ALL_established2 = sapply(2:32, function(x) {
  if(choose(n = 32,k = x)>B){
    chosen = replicate(B,sort(sample(1:32,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:32,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:32,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    FD_Beta(c(0,1,2),Data_established[,(chosen[i,]+2)], D_established, "U", "relative")
  }) %>% rowMeans; })
result.R_established2 = sapply(2:14, function(x) {
  if(choose(n = 14,k = x)>B){
    chosen = replicate(B,sort(sample(1:14,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:14,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:14,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    FD_Beta(c(0,1,2),Data_established[,(chosen[i,]+20)], D_established, "U", "relative")
  }) %>% rowMeans; })
result.F_established2 = sapply(2:18, function(x) {
  if(choose(n = 18,k = x)>B){
    chosen = replicate(B,sort(sample(1:18,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:18,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:18,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    FD_Beta(c(0,1,2),Data_established[,(chosen[i,]+2)], D_established, "U", "relative")
  }) %>% rowMeans; })
Get_plot(result.ALL_established2, result.R_established2, result.F_established2, "U", "E", "Fun")

#((b) Jaccard-type functional dissimilarity measure) Compare All, Restoration Fragments base on regenerating, 1-UqN, relative.
result.ALL_regenerating2 <- sapply(2:32, function(x) {
  if(choose(n = 32,k = x)>B){
    chosen = replicate(B,sort(sample(1:32,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:32,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:32,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    FD_Beta(c(0,1,2),Data_regenerating[,(chosen[i,]+2)], D_regenerating, "U", "relative")
  }) %>% rowMeans; })
result.R_regenerating2 = sapply(2:14, function(x) {
  if(choose(n = 14,k = x)>B){
    chosen = replicate(B,sort(sample(1:14,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:14,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:14,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    FD_Beta(c(0,1,2),Data_regenerating[,(chosen[i,]+20)], D_regenerating, "U", "relative")
  }) %>% rowMeans; })
result.F_regenerating2 = sapply(2:18, function(x) {
  if(choose(n = 18,k = x)>B){
    chosen = replicate(B,sort(sample(1:18,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:18,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:18,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    FD_Beta(c(0,1,2),Data_regenerating[,(chosen[i,]+2)], D_regenerating, "U", "relative")
  }) %>% rowMeans; })
Get_plot(result.ALL_regenerating2, result.R_regenerating2, result.F_regenerating2, "U", "R", "Fun")

#-----------------------------------------------
# Taxonmic diversity (empirical Hill numbers)
#-----------------------------------------------
#' Diversity_profile_MLE(data, q) is a function which computes taxonmic diversity of order q based on abundance data.
#' @param data a vector of species sample frequencies.
#' @param q a numeric or a vector of diversity order; The suggested range for q is [0, 3].
#' @return a numerical vector of taxonmic diversity.
Diversity_profile_MLE <- function(data, q){
  pi <- if (sum(data) != 1) data/sum(data) else data
  pi <- pi[pi>0]
  Sub <- function(q){
    if (q == 1) exp(-sum(pi*log(pi))) else exp(1/(1-q)*log(sum(pi^q)))
  }
  sapply(q, Sub)

}

#-------------------------------------------------------------------
# Taxonmic dissimilarity (empirical) See Chao and Chiu (2016, MEE)
#------------------------------------------------------------------
#' beta_diversity_MLE(q, data, CU, method) is a function which computes taxonmic dissimilarity measure of order q based on abundance data.
#' @param q a numeric or a vector of diversity order; The suggested range for q is [0, 3].
#' @param data a S*N matrix of species sample frequencies with S species, N communities.
#' @param CU a character to choose method, "C" for 1-CqN (Sorensen); "U" for 1-UqN (Jaccard).
#' @param method a character to choose method, "relative" or "absolute".
#' @return a numerical vector of taxonmic dissimilarity.
beta_diversity_MLE <- function(q, data, CU, method){
  if(method == "relative") data <- sapply(1:ncol(data), FUN = function(i) data[,i]/sum(data[,i]))
  N <- ncol(data)
  gamma_MLE = Diversity_profile_MLE(rowSums(data), q)
  alpha_MLE = Diversity_profile_MLE(unlist(data), q)/N
  beta_MLE = gamma_MLE/alpha_MLE
  beta_MLE[beta_MLE<1] = 1

  if(CU=="C"){
    out = (1-beta_MLE^(1-q))/(1-N^(1-q))
  }else{
    out = (1-(beta_MLE)^(q-1))/(1-N^(q-1))
  }
  out[q==1] = log(beta_MLE[q==1])/log(N)
  out
}

#### Plot Appendix S5: Figure S1
#B is the number of combinations of sites for each chosen number of point on x-axis. 
#If the number of combinations exceeds 100, then cut off at B = 500
B <- 500
#((a) Sorensen dissimilarity measure) Compare all, Restoration Fragments base on established, 1-CqN, relative
result.ALL_established5 = sapply(2:32, function(x) {
  if(choose(n = 32,k = x)>B){
    chosen = replicate(B,sort(sample(1:32,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:32,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:32,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    beta_diversity_MLE(c(0,1,2),Data_established[,(chosen[i,]+2)], "C", "relative")
  }) %>% rowMeans; })
result.R_established5 = sapply(2:14, function(x) {
  if(choose(n = 14,k = x)>B){
    chosen = replicate(B,sort(sample(1:14,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:14,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:14,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    beta_diversity_MLE(c(0,1,2),Data_established[,(chosen[i,]+20)], "C", "relative")
  }) %>% rowMeans; })
result.F_established5 = sapply(2:18, function(x) {
  if(choose(n = 18,k = x)>B){
    chosen = replicate(B,sort(sample(1:18,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:18,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:18,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    beta_diversity_MLE(c(0,1,2),Data_established[,(chosen[i,]+2)], "C", "relative")
  }) %>% rowMeans; })
Get_plot(result.ALL_established5, result.R_established5, result.F_established5, "C", "E", "Tax")

#((a) Sorensen dissimilarity measure) Compare all, Restoration Fragments base on regenerating, 1-CqN, relative
result.ALL_regenerating5 = sapply(2:32, function(x) {
  if(choose(n = 32,k = x)>B){
    chosen = replicate(B,sort(sample(1:32,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:32,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:32,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    beta_diversity_MLE(c(0,1,2),Data_regenerating[,(chosen[i,]+2)], "C", "relative")
  }) %>% rowMeans; })
result.R_regenerating5 = sapply(2:14, function(x) {
  if(choose(n = 14,k = x)>B){
    chosen = replicate(B,sort(sample(1:14,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:14,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:14,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    beta_diversity_MLE(c(0,1,2),Data_regenerating[,(chosen[i,]+20)], "C", "relative")
  }) %>% rowMeans; })
result.F_regenerating5 = sapply(2:18, function(x) {
  if(choose(n = 18,k = x)>B){
    chosen = replicate(B,sort(sample(1:18,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:18,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:18,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    beta_diversity_MLE(c(0,1,2),Data_regenerating[,(chosen[i,]+2)], "C", "relative")
  }) %>% rowMeans; })
Get_plot(result.ALL_regenerating5, result.R_regenerating5, result.F_regenerating5, "C", "R", "Tax")

#### Plot Fig 5. in the main text
#((a) Jaccard dissimilarity measure) Compare all, Restoration Fragments base on established, 1-UqN, relative
result.ALL_established6 = sapply(2:32, function(x) {
  if(choose(n = 32,k = x)>B){
    chosen = replicate(B,sort(sample(1:32,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:32,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:32,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    beta_diversity_MLE(c(0,1,2),Data_established[,(chosen[i,]+2)], "U", "relative")
  }) %>% rowMeans; })
result.R_established6 = sapply(2:14, function(x) {
  if(choose(n = 14,k = x)>B){
    chosen = replicate(B,sort(sample(1:14,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:14,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:14,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    beta_diversity_MLE(c(0,1,2),Data_established[,(chosen[i,]+20)], "U", "relative")
  }) %>% rowMeans; })
result.F_established6 = sapply(2:18, function(x) {
  if(choose(n = 18,k = x)>B){
    chosen = replicate(B,sort(sample(1:18,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:18,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:18,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    beta_diversity_MLE(c(0,1,2),Data_established[,(chosen[i,]+2)], "U", "relative")
  }) %>% rowMeans; })
Get_plot(result.ALL_established6, result.R_established6, result.F_established6, "U", "E", "Tax")

#((b) Jaccard dissimilarity measure) Compare all, Restoration Fragments base on regenerating, 1-UqN, relative
result.ALL_regenerating6 = sapply(2:32, function(x) {
  if(choose(n = 32,k = x)>B){
    chosen = replicate(B,sort(sample(1:32,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:32,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:32,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    beta_diversity_MLE(c(0,1,2),Data_regenerating[,(chosen[i,]+2)], "U", "relative")
  }) %>% rowMeans; })
result.R_regenerating6 = sapply(2:14, function(x) {
  if(choose(n = 14,k = x)>B){
    chosen = replicate(B,sort(sample(1:14,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:14,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:14,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    beta_diversity_MLE(c(0,1,2),Data_regenerating[,(chosen[i,]+20)], "U", "relative")
  }) %>% rowMeans; })
result.F_regenerating6 = sapply(2:18, function(x) {
  if(choose(n = 18,k = x)>B){
    chosen = replicate(B,sort(sample(1:18,size = x,replace = F))) %>% t()
    chosen = unique(chosen)
    while(nrow(chosen)<B){
      new = replicate(B-nrow(chosen),sort(sample(1:18,size = x,replace = F))) %>% t()
      chosen = rbind(chosen,new)
      chosen = unique(chosen)
    }
  }else{
    chosen = t(combn(x = 1:18,m = x))
  }
  sapply(1:nrow(chosen), function(i){
    beta_diversity_MLE(c(0,1,2),Data_regenerating[,(chosen[i,]+2)], "U", "relative")
  }) %>% rowMeans; })
Get_plot(result.ALL_regenerating6, result.R_regenerating6, result.F_regenerating6, "U", "R", "Tax")
