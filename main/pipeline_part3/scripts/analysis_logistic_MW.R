###----- ANALYSIS WITH R LOGISTIC.R -----#
library(gplots, lib.loc = "/proj/sens2018106/softwares/R_packages/")
library(ROCR, lib.loc = "/proj/sens2018106/softwares/R_packages/")

source("scripts/logistic.R")

args       <- commandArgs(trailingOnly=TRUE)
cypfile    <- args[1]
pcfile     <- args[2]
samplefile <- args[3]
assocfile  <- args[4]
phenofile  <- args[5]
outfile    <- args[6]
outfolder  <- args[7]

# read in cyp data
cyp           <- read.table(cypfile, header=F, sep=",", stringsAsFactors = FALSE)
colnames(cyp) <- c("sample","Hap1","Hap2","ActScore","PredPheno","SV","CNV","Cand1","Cand2","ChangePoint")

# read in PC data
pcs           <- read.table(pcfile, stringsAsFactors = FALSE)
colnames(pcs) <- c("sample","sample2","T","U","V","W","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

# read in samplefile
sample              <- read.table(samplefile,fill = T, header = T, encoding = "UTF-8", stringsAsFactors = FALSE)
colnames(sample)[1] <- "sample"


# read in file with which associations to test
assoc             <- read.table(assocfile, encoding = "UTF-8", stringsAsFactors = FALSE, col.names = c("tests"))
assoc$tests.fixed <- chartr("~", ".", assoc$tests)
assoc$tests.fixed <- chartr("-", ".", assoc$tests.fixed)
assoc$tests.fixed <- chartr("/", ".", assoc$tests.fixed)


# read in phenofile
pheno              <- read.table(phenofile, header = T, encoding = "UTF-8", stringsAsFactors = FALSE, na.strings = c("NA", -9))
colnames(pheno)[1] <- "sample"


# merge cyp and PC data with sample and pheno data
int1   <- merge(cyp,pcs,by="sample")
int2   <- merge(int1,sample,by="sample")
data.1 <- int2[,c("sample","Hap1","Hap2","ActScore","PredPheno","SV","CNV","Cand1","Cand2","sex","cohort","batch","suspectsubstance","adrtype","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
data.2 <- merge(data.1,pheno, by="sample")

# make PredPheno categorical
data.2$PredPhenoBin <- NA
data.2$PredPhenoBin[data.2$PredPheno=="PM"] <- 0
data.2$PredPhenoBin[data.2$PredPheno=="IM"] <- 1
data.2$PredPhenoBin[data.2$PredPheno=="NM"] <- 2
data.2$PredPhenoBin[data.2$PredPheno=="UM"] <- 3

# run glm on all tests specified
sink(outfile, append=TRUE, split=TRUE)
for (test in assoc$tests.fixed) {
  for (i in test) {
    print(i)
    model <- glm(recode(data.2[[i]])~data.2$PredPhenoBin+data.2$PC1+data.2$PC2+data.2$PC3+data.2$PC4+data.2$PC5+data.2$PC6+data.2$PC7+data.2$PC8+data.2$PC9+data.2$PC10, family=binomial, na.action = na.exclude)
    print(summary(model))
    
    p     <- predict(model, type="response", na.action = na.exclude)
    pr    <- prediction(p, recode(data.2[[i]]))
    prf   <- performance(pr, measure = "tpr", x.measure = "fpr")
    auc   <- performance(pr, measure = "auc")
    auc   <- auc@y.values[[1]]
    
    png(paste0(outfolder, i, '.png'))
    plot(prf, lwd=2, col='blue', main = i)
    abline(0,1, col='red', lwd=2, lty=2)
    text(0.9,0.1, labels=paste('AUC:',round(auc,digits=2)))
    dev.off()
  }
}
sink()



