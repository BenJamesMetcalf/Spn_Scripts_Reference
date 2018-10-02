

## Function PBP_AA_TO_MIC2

PBP_AA_TO_MIC2<- function(cwd){
  dbdir="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/share/PBP_AA_to_MIC/newDB/"
  libpath="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/share/PBP_AA_to_MIC/Rlib"
  #x1=.libPaths()
  x2=c(libpath, "/usr/lib64/R/library", "/usr/share/R/library")
  .libPaths(x2)
  #.libPaths(c(.libPaths(), libpath))
  library("methods")
  library("randomForest")
  library('iterators')
  library('foreach')
  library('glmnet')
  
  ## Data processing
  setwd(cwd)

  colAA=5:918
  m1=read.csv(paste(dbdir, "train_PT_AA_table.csv", sep=""), colClasses="character") 
  n1=dim(m1)[1] 
  x1=matrix(unlist(m1[,colAA]), ncol=length(colAA), nrow=n1, byrow=F)
 
  m2=read.csv("Sample_PBP_AA_table.csv", colClasses="character")
  #m2=m2[, c(1, 9:925)]
  AF=paste(m2[, 2], m2[, 3], m2[, 4])=="0 0 0"
  m2=m2[AF, ]
    # Only keep isolates with all PBPs found
  n2=dim(m2)[1]
  m2[,559]="N"

  ## Median MIC Prediction
  APT=rep(NA, n2)
  DIS1A=rep(NA, n2)
  DIS2B=rep(NA, n2)
  DIS2X=rep(NA, n2)
	
  for (j1 in 1:n2)
  {
    print (c("processing sample", m2[j1,1]))
    x2=matrix(rep(m2[j1, colAA], n1), ncol=length(colAA), nrow=n1, byrow=T)
    x3=apply((x2 != x1), 1, sum, na.rm=T)
    x4=which(x3==min(x3, na.rm=T))[1]
    APT[j1]=m1[x4, 1]

    x5 = (m2[j1, colAA] != m1[x4, colAA])        
    DIS1A[j1]=sum(x5[grep("PBP1A_", colnames(x5))], na.rm=T)
    DIS2B[j1]=sum(x5[grep("PBP2B_", colnames(x5))], na.rm=T)
    DIS2X[j1]=sum(x5[grep("PBP2X_", colnames(x5))], na.rm=T)
  }

  m3=data.frame(cbind(sampleID=m2[, 1], APT, DIS1A, DIS2B, DIS2X))

  m3$PT=as.character(m3$APT)
  m3$NDIST=DIS1A+DIS2B+DIS2X
  m3$PT[m3$NDIST>0]="new"
  m4=read.csv(paste(dbdir, "Ref_PBPtype_MIC.csv", sep=""), colClasses="character") 
  m5=merge(m3, m4, by="APT", all.x=T)
  m6=m5[, c(2,6,1,3:5, 7:13)]
  x1=colnames(m6)[8:13]
  colnames(m6)[8:13]=paste(x1, "_MIC_MM", sep="")
  #write.csv(m6, file="Sample_PBPtype_MIC2_Prediction.csv", row.names=F)
  
  BLAclass=c("PEN", "AMO", "MER", "TAX", "CFT", "CFX")
  BKclass1=c(0.06, 2, 0.25, 0.5, 0.5, 0.5)
  BKclass2=c(4,    4, 0.5,  1,   1,   1)
  
  # Median MIC fro BK1
  m3=m6[, c(1, 8:13)]
  colnames(m3)[1]="sampleID"
  colnames(m3)[2:7]=paste(BLAclass, "_BK1_MM", sep="")
  m7=merge(m6, m3, by="sampleID", all.x=T)
  m6=m7
 

  ## RandomForest and Elastic-Net
  ## Re-level predictors
  m1=read.csv(paste(dbdir, "NR_Varsites_Randomforest.csv", sep=""), colClasses="character") 
  colAA=colnames(m1)
  
  m2=read.csv("Sample_PBP_AA_table.csv", colClasses="character")
  #m2=m2[, c(1, 9:925)]
  colnames(m2)[1]="sampleID"

  AF=paste(m2[, 2], m2[, 3], m2[, 4])=="0 0 0"
  m2=m2[AF, ]
    # Only keep isolates with all PBPs found
  m2.1=m2[colAA]  
  n2=dim(m2.1)[2]


  m3=read.csv("/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/share/PBP_AA_to_MIC/scripts/BLOSUM62.csv", colClasses="character", header=T) 
  colnames(m3)[2:25]=m3$src
  m2.1[m2.1=="*"]="-"
  m2.2=m2.1
  for (j1 in 1:n2)
  {
    x1=levels(factor(m1[, j1]))
    x1.1=x1[nchar(x1)==1]
    x1.2=x1
    if (length(x1.1)==0) {m2.2[, j1]=factor(m1[1, j1], levels=x1)} else {
    x1=x1.1
    x2=levels(factor(m2.1[, j1]))
    x3=x2[ !(x2 %in% x1)]
    if (length(x3)>0) 
    {
      for (j2 in x3)
      {
        AA=substr(j2, 1, 1)
        x4=(x1[order(m3[m3$src==AA, x1], decreasing =T)])[1]
        x5=m2.1[, j1]
        m2.1[j2==x5, j1]=x4
        print(c(j1, x4))
      }       
    }
    m2.2[, j1]=factor(m2.1[, j1], levels=x1.2) }
  }

  #Model fitting Random Forest MIC
  setwd(dbdir)

  RFdbMIC=c("PEN_RandomForestMIC", "AMO_RandomForestMIC", "MER_RandomForestMIC",
	  "TAX_RandomForestMIC", "CFT_RandomForestMIC", "CFX_RandomForestMIC")  
  m3=m2$sampleID
  for (j1 in 1:6)
  {
    rm(fit1)
    load(RFdbMIC[j1])
    m3=cbind(m3, round(2^predict(fit1, newdata = m2.2), 2))
  }
  colnames(m3)[1]="sampleID"
  colnames(m3)[2:7]=paste(BLAclass, "_MIC_RF", sep="")
  m7=merge(m6, m3, by="sampleID", all.x=T)
  m6=m7
  
  #Model fitting Random Forest BK1

  RFdbBK1=c("PEN_RandomForestBK1", "AMO_RandomForestBK1", "MER_RandomForestBK1",
	  "TAX_RandomForestBK1", "CFT_RandomForestBK1", "CFX_RandomForestBK1")  
  m3=m2$sampleID
  for (j1 in 1:6)
  {
    rm(fit1)
    load(RFdbBK1[j1])
    m3=cbind(m3, round(predict(fit1, newdata = m2.2, type="prob"), 3)[, 2])
  }

  colnames(m3)[1]="sampleID"
  colnames(m3)[2:7]=paste(BLAclass, "_BK1_RF", sep="")
  m7=merge(m6, m3, by="sampleID", all.x=T)
  m6=m7  

  #Model fitting Elastic Net MIC

  ENdbMIC=c("PEN_ElasticNetMIC", "AMO_ElasticNetMIC", "MER_ElasticNetMIC",
	  "TAX_ElasticNetMIC", "CFT_ElasticNetMIC", "CFX_ElasticNetMIC")  

  xfactors <- model.matrix(~ ., data=m2.2)[, -1]
  x1 <- as.matrix(data.frame(xfactors))
  if (dim(x1)[2]==1) {x1=t(x1)}  
  m3=m2$sampleID
  for (j1 in 1:6)
  {
    rm(fit1)
    load(ENdbMIC[j1])
    m3=cbind(m3, round(2^predict(fit1, x1, type="response"), 2))
  }

  colnames(m3)[1]="sampleID"
  colnames(m3)[2:7]=paste(BLAclass, "_MIC_EN", sep="")
  m7=merge(m6, m3, by="sampleID", all.x=T)
  m6=m7

  #Model fitting Elastic Net BK1
  
  ENdbBK1=c("PEN_ElasticNetBK1", "AMO_ElasticNetBK1", "MER_ElasticNetBK1",
	  "TAX_ElasticNetBK1", "CFT_ElasticNetBK1", "CFX_ElasticNetBK1")  

  xfactors <- model.matrix(~ ., data=m2.2)[, -1]
  x1 <- as.matrix(data.frame(xfactors)) 
  if (dim(x1)[2]==1) {x1=t(x1)}  
  m3=m2$sampleID
  for (j1 in 1:6)
  {
    rm(fit1)
    load(ENdbBK1[j1])
    m3=cbind(m3, round(predict(fit1, x1, type="response"), 3))
  }

  colnames(m3)[1]="sampleID"
  colnames(m3)[2:7]=paste(BLAclass, "_BK1_EN", sep="")
  m7=merge(m6, m3, by="sampleID", all.x=T)
  m6=m7

  setwd(cwd)
  write.csv(m6, file="Sample_PBPtype_MIC2_Prediction.csv", row.names=F) 
}

## End of Function PBP_AA_TO_MIC2

args <- commandArgs(TRUE)
cwd = args[1]

print (cwd)
PBP_AA_TO_MIC2(cwd)



