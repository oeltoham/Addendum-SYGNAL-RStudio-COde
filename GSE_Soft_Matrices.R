library(GEOquery)
setwd("~/GSE_Soft_Matrices") #sets directory
getwd()
getGEO(filename = paste("C:", "Users", "oelto", "Documents", "GSE_Soft_Matrices", paste("GSE", gse1, "_series_matrix.txt.gz",sep=""), sep = "/", collapse = NULL), GSEMatrix = T, getGPL = T) #extracts GSE matrix file from directory
gseNums <- c(15460, 14860, 10141, 15641, 19417, 19422, 19750, 19987, 22138, 25097, 26253, 26566, 26939, 27155, 29174, 29354, 29695, 31448, 32062, 32225, 32984, 33630, 35158, 39366, 39582, 40435, 45725, 4573, 49278, 56303, 65858, 71118, 71729, 72094, 8607, 9891, 10846, 9843, 46517, 21034, 19949)
#logtrans = rep(FALSE,36)
#names(logtrans) = gseNums

logtrans = rep(FALSE,41)
names(logtrans) = gseNums
for (gse1 in c(15460, 4573, 8607, 15641, 19750, 19987, 25097, 26253, 26566, 29695, 31448, 32225, 32984, 39582, 40435, 45725, 65858, 71118, 71729, 72094)) #sets specific GSEs as true for log transformation 
  {
  logtrans[as.character(gse1)] = TRUE
  
}

for(gse1 in gseNums){
  print(gse1)
  gset1 = getGEO(filename = paste("C:", "Users", "oelto", "Documents", "GSE_Soft_Matrices", paste("GSE", gse1, "_series_matrix.txt.gz",sep=""),  sep = "/"), GSEMatrix = T, getGPL = T)
  png(filename = paste("C:", "Users", "oelto", "Documents", paste("GSE", gse1, "_LogTrans3.png", sep = ""), sep = "/"), width = 1024, height = 576) #saves png file of each boxplot
  if(as.logical(logtrans[as.character(gse1)])==T) {
    boxplot(log2(exprs(gset1)), main = paste('GSE',gse1,sep=''), outline = F)
  }  else {
    boxplot(exprs(gset1), main = paste('GSE',gse1,sep=''), outline = F)
  }
  dev.off()
}

