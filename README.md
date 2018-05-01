# Addendum-SYGNAL-RStudio-COde
Code for RStudio and Everything done for SYGNAL Project
Boxplot and normalization code:
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

Pheno data & expression set code:

# Load libraries
library(GEOquery)
library(genefilter)

# Set working directory
setwd('C:/Users/cplaisie/Dropbox (ASU)/ASU/Students/Omar')

# Load up data
gset <- getGEO(filename="data/GSE26253_series_matrix.txt.gz", GSEMatrix=TRUE, AnnotGPL=T)
pData(gset)
write.csv(pData(gset), file = "GSE26253PHENO.csv")

# Use genefilter to filter out gene names with redundant Entrez IDs
library(genefilter)
library("illuminaHumanWGDASLv3.db")
library(MASS)
gN = rownames(exprs(gset))
testStat = rowVars(exprs(gset))/rowMeans(exprs(gset))
filteredNames = findLargest(gN, testStat, data = "illuminaHumanWGDASLv3")

# Apply to expression set
exprs1 = exprs(gset)
exprs1 = exprs1[filteredNames,]
rownames(exprs1) = names(filteredNames)
write.csv(exprs1, 'GSE26253_expression_nonredundant.csv')

Replication analysis code:
#!/usr/bin/env Rscript

suppressMessages(library(WGCNA))
#suppressMessages(library(multicore))
suppressMessages(library(getopt))
suppressMessages(library(parallel))
suppressMessages(library(impute))
suppressMessages(library(survival))

# read command line arguments
spec = matrix(c(
  'dataset','d',1,'character',
  'outdir', 'o', 1, 'character',
  'help', 'h', 0, 'logical'
), byrow=TRUE, ncol=4)

opt <- getopt(spec)

if (is.null(opt$dataset) || is.null(opt$outdir) || !is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

#allowWGCNAThreads(opt$cores)

getEigengene <- function (expr, colors, impute = TRUE, nPC = 1, align = "along average", 
                          excludeGrey = FALSE, grey = ifelse(is.numeric(colors), 0, 
                                                             "grey"), subHubs = FALSE, trapErrors = FALSE, returnValidOnly = trapErrors, 
                          softPower = 1, scale = TRUE, verbose = 0, indent = 0) 
{
  spaces = indentSpaces(indent)
  if (verbose == 1) 
    printFlush(paste(spaces, "moduleEigengenes: Calculating", 
                     nlevels(as.factor(colors)), "module eigengenes in given set."))
  if (is.null(expr)) {
    stop("moduleEigengenes: Error: expr is NULL. ")
  }
  if (is.null(colors)) {
    print("moduleEigengenes: Error: colors is NULL. ")
    stop()
  }
  if (is.null(dim(expr)) || length(dim(expr)) != 2) 
    stop("moduleEigengenes: Error: expr must be two-dimensional.")
  #if (dim(expr)[2] != length(colors)) 
  #    stop("moduleEigengenes: Error: ncol(expr) and length(colors) must be equal (one color per gene).")
  #if (is.factor(colors)) {
  #    nl = nlevels(colors)
  #    nlDrop = nlevels(colors[, drop = TRUE])
  #    if (nl > nlDrop) 
  #        stop(paste("Argument 'colors' contains unused levels (empty modules). ", 
  #            "Use colors[, drop=TRUE] to get rid of them."))
  #}
  if (softPower < 0) 
    stop("softPower must be non-negative")
  alignRecognizedValues = c("", "along average")
  if (!is.element(align, alignRecognizedValues)) {
    printFlush(paste("ModulePrincipalComponents: Error:", 
                     "parameter align has an unrecognised value:", align, 
                     "; Recognized values are ", alignRecognizedValues))
    stop()
  }
  maxVarExplained = 10
  if (nPC > maxVarExplained) 
    warning(paste("Given nPC is too large. Will use value", 
                  maxVarExplained))
  nVarExplained = min(nPC, maxVarExplained)
  modlevels = 1:length(colors)
  PrinComps = data.frame(matrix(NA, nrow = dim(expr)[[1]], 
                                ncol = length(modlevels)))
  averExpr = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modlevels)))
  varExpl = data.frame(matrix(NA, nrow = nVarExplained, ncol = length(modlevels)))
  validMEs = rep(TRUE, length(modlevels))
  validAEs = rep(FALSE, length(modlevels))
  isPC = rep(TRUE, length(modlevels))
  isHub = rep(FALSE, length(modlevels))
  validColors = colors
  names(PrinComps) = paste(moduleColor.getMEprefix(), modlevels, 
                           sep = "")
  names(averExpr) = paste("AE", modlevels, sep = "")
  for (i in c(1:length(modlevels))) {
    if (length(colors[[i]])>0) {
      if (verbose > 1) 
        printFlush(paste(spaces, "moduleEigengenes : Working on ME for module", 
                         modlevels[i]))
      modulename = modlevels[i]
      restrict1 = colors[[modulename]]
      #if (verbose > 2) 
      #    printFlush(paste(spaces, " ...", sum(restrict1), 
      #        "genes"))
      datModule = as.matrix(t(expr[ ,restrict1]))
      n = dim(datModule)[1]
      p = dim(datModule)[2]
      pc = try({
        if (nrow(datModule) > 1 && impute) {
          seedSaved = FALSE
          if (exists(".Random.seed")) {
            saved.seed = .Random.seed
            seedSaved = TRUE
          }
          if (verbose > 5) 
            printFlush(paste(spaces, " ...imputing missing data"))
          datModule = impute.knn(as.matrix(datModule), 
                                 k = min(10, nrow(datModule) - 1))
          try({
            if (!is.null(datModule$data)) 
              datModule = datModule$data
          }, silent = TRUE)
          if (seedSaved) 
            .Random.seed <<- saved.seed
        }
        if (verbose > 5) 
          printFlush(paste(spaces, " ...scaling"))
        if (scale) 
          datModule = t(scale(t(datModule)))
        if (verbose > 5) 
          printFlush(paste(spaces, " ...calculating SVD"))
        svd1 = svd(datModule, nu = min(n, p, nPC), nv = min(n, 
                                                            p, nPC))
        if (verbose > 5) 
          printFlush(paste(spaces, " ...calculating PVE"))
        veMat = cor(svd1$v[, c(1:min(n, p, nVarExplained))], 
                    t(datModule), use = "p")
        varExpl[c(1:min(n, p, nVarExplained)), i] = apply(veMat^2, 
                                                          1, mean, na.rm = TRUE)
        svd1$v[, 1]
      }, silent = TRUE)
      if (class(pc) == "try-error") {
        if ((!subHubs) && (!trapErrors)) 
          stop(pc)
        if (subHubs) {
          if (verbose > 0) {
            printFlush(paste(spaces, " ..principal component calculation for module", 
                             modulename, "failed with the following error:"))
            printFlush(paste(spaces, "     ", pc, spaces, 
                             " ..hub genes will be used instead of principal components."))
          }
          isPC[i] = FALSE
          pc = try({
            scaledExpr = scale(t(datModule))
            covEx = cov(scaledExpr, use = "p")
            modAdj = abs(covEx)^softPower
            kIM = (apply(modAdj, 1, sum, na.rm = TRUE))^3
            if (max(kIM, na.rm = TRUE) > 1) 
              kIM = kIM - 1
            kIM[is.na(kIM)] = 0
            hub = which.max(kIM)
            alignSign = sign(covEx[, hub])
            alignSign[is.na(alignSign)] = 0
            isHub[i] = TRUE
            pcxMat = scaledExpr * matrix(kIM * alignSign, 
                                         nrow = nrow(scaledExpr), ncol = ncol(scaledExpr), 
                                         byrow = TRUE)/sum(kIM)
            pcx = apply(pcxMat, 1, sum, na.rm = TRUE)
            varExpl[1, i] = mean(cor(pcx, t(datModule), 
                                     use = "p")^2, na.rm = TRUE)
            pcx
          }, silent = TRUE)
        }
      }
      if (class(pc) == "try-error") {
        if (!trapErrors) 
          stop(pc)
        if (verbose > 0) {
          printFlush(paste(spaces, " ..ME calculation of module", 
                           modulename, "failed with the following error:"))
          printFlush(paste(spaces, "     ", pc, spaces, 
                           " ..the offending module has been removed."))
        }
        warning(paste("Eigengene calculation of module", 
                      modulename, "failed with the following error \n     ", 
                      pc, "The offending module has been removed.\n"))
        validMEs[i] = FALSE
        isPC[i] = FALSE
        isHub[i] = FALSE
        validColors[restrict1] = grey
      }
      else {
        PrinComps[, i] = pc
        ae = try({
          if (isPC[i]) 
            scaledExpr = scale(t(datModule))
          averExpr[, i] = apply(scaledExpr, 1, mean, na.rm = TRUE)
          if (align == "along average") {
            if (verbose > 4) 
              printFlush(paste(spaces, " .. aligning module eigengene with average expression."))
            if (cor(averExpr[, i], PrinComps[, i], use = "p") < 
                0) 
              PrinComps[, i] = -PrinComps[, i]
          }
          0
        }, silent = TRUE)
        if (class(ae) == "try-error") {
          if (!trapErrors) 
            stop(ae)
          if (verbose > 0) {
            printFlush(paste(spaces, " ..Average expression calculation of module", 
                             modulename, "failed with the following error:"))
            printFlush(paste(spaces, "     ", ae, spaces, 
                             " ..the returned average expression vector will be invalid."))
          }
          warning(paste("Average expression calculation of module", 
                        modulename, "failed with the following error \n     ", 
                        ae, "The returned average expression vector will be invalid.\n"))
        }
        validAEs[i] = !(class(ae) == "try-error")
      }
    }
  }
  allOK = (sum(!validMEs) == 0)
  if (returnValidOnly && sum(!validMEs) > 0) {
    PrinComps = PrinComps[, validMEs]
    averExpr = averExpr[, validMEs]
    varExpl = varExpl[, validMEs]
    validMEs = rep(TRUE, times = ncol(PrinComps))
    isPC = isPC[validMEs]
    isHub = isHub[validMEs]
    validAEs = validAEs[validMEs]
  }
  allPC = (sum(!isPC) == 0)
  allAEOK = (sum(!validAEs) == 0)
  list(eigengenes = PrinComps, averageExpr = averExpr, varExplained = varExpl, 
       nPC = nPC, validMEs = validMEs, validColors = validColors, 
       allOK = allOK, allPC = allPC, isPC = isPC, isHub = isHub, 
       validAEs = validAEs, allAEOK = allAEOK)
}

# Read in genes for each cluster
g1 = read.csv('gene_members/tfbs_db/cluster.members.genes.txt',header=F)
biclustMembership = list()
for(j in 1:length(g1[,1])) {
  biclustMembership[[paste('TFBS_DB_',j,sep='')]] = strsplit(as.character(g1[j,]),split=' ')[[1]][-1]
}
g2 = read.csv('gene_members/pita/cluster.members.genes.txt',header=F)
for(j in 1:length(g2[,1])) {
  biclustMembership[[paste('PITA_',j,sep='')]] = strsplit(as.character(g2[j,]),split=' ')[[1]][-1]
}
g3 = read.csv('gene_members/targetscan/cluster.members.genes.txt',header=F)
for(j in 1:length(g3[,1])) {
  biclustMembership[[paste('TargetScan_',j,sep='')]] = strsplit(as.character(g3[j,]),split=' ')[[1]][-1]
}

# Read in a second dataset
ratSec <- read.delim( file=paste('replication_',opt$dataset,'/',opt$dataset,'_ratios.csv',sep=''), sep=",", as.is=T, header=T,row.names=1 )
ratSec = ratSec[which(apply(ratSec,1,sum)!=0),]
rownames( ratSec ) <- toupper( rownames(ratSec) )
biclustMembership.sec = list()
for(j in names(biclustMembership)) {
  biclustMembership.sec[[j]] = intersect(unlist(biclustMembership[[j]]), rownames(ratSec))
}
ratSec <- as.matrix(ratSec)
p1 = read.csv(paste('replication_',opt$dataset,'/',opt$dataset,'_phenotypes.csv',sep=''),header=T,row.names=1)

# Calculate the residuals for all clusters in the second dataset
ks = length(biclustMembership)
outNames = c('n.rows','overlap.rows','pc1.var.exp','avg.pc1.var.exp','pc1.perm.p','survival','survival.p','survival.age','survival.age.p','survival.age.sex','survival.age.sex.p')
m1 = matrix(ncol=length(outNames),nrow=ks,dimnames=list(names(biclustMembership),outNames))
permutations = 100
for(k in names(biclustMembership)) {
  # Get and add number of rows and columns
  k.rows.sec = biclustMembership.sec[[k]]
  print(k)
  if(length(k.rows.sec)>1) {
    m1[k,1] = length(biclustMembership[[k]])
    m1[k,2] = length(k.rows.sec)
    # Use eigengenes
    testEm.rows = list()
    testEm.rows[[1]] = k.rows.sec
    for( i in 2:(permutations+1)) {
      testEm.rows[[i]] = sample(rownames(ratSec),m1[k,2])
    }
    eg1 = getEigengene(t(ratSec),testEm.rows) #,verbose=10)
    var.exp = t(eg1$varExplained)[,1]
    m1[k,3] = var.exp[1]
    m1[k,4] = mean(var.exp[2:length(var.exp)],na.rm=TRUE)
    m1[k,5] = length(which(na.omit(var.exp[2:length(var.exp)]) >= m1[k,3]))/length(na.omit(var.exp[2:length(var.exp)]))
    pc.1 = t(eg1$eigengenes)[1,]
    # Survival analysis
    d2 = data.frame(p1[colnames(ratSec),],pc.1)
    scph1 = summary(coxph(Surv(OS.time,OS==1) ~ pc.1, data=d2))
    m1[k,6] = scph1$coef[1,4]
    m1[k,7] = scph1$coef[1,5]
    if('age' %in% colnames(d2)) {
      scph2 = summary(coxph(Surv(OS.time,OS==1) ~ pc.1 + age, data=d2))
      m1[k,8] = scph2$coef[1,4]
      m1[k,9] = scph2$coef[1,5]
      if('gender' %in% colnames(d2)) {
        scph3 = summary(coxph(Surv(OS.time,OS==1) ~ pc.1 + age + gender, data=d2))
        m1[k,10] = scph3$coef[1,4]
        m1[k,11] = scph3$coef[1,5]
      }
    }
  }
  print(m1[k,])
}
write.csv(m1,file=paste(opt$outdir,'/replicationPvalues_STAD_',opt$dataset,'.csv',sep=''))
