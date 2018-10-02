
MIC_format <- function(PBPtype_MIC2_Prediction)
{

  
#PBPtype_MIC2_Prediction="/scicomp/home/yqh8/download/GPS/Spn_Typing_Analysis/NICD_859/SPN_Typing_Output/PBPtoMIC/Sample_PBPtype_MIC2_Prediction.csv"  
  fout=paste(PBPtype_MIC2_Prediction, "MIC_formatted_with_SIR.csv", sep="_")
  m3=read.csv(PBPtype_MIC2_Prediction, colClasses="character")
  m3a=m3[, c(1, 20:25)]  
  n1=dim(m3a)[1]
  SampleID=m3a[, 1]

  #Format PEN
  i=2
  x1=as.integer(round(log2(as.numeric(m3a[, i]))))
  sx1=rep("=", n1)
  sx1[x1<=-5]="<="
  x1[x1<=-5]=-5
  sx1[x1>4]=">"
  x1[x1>4]=4
  x1a=rep(NA, n1)
  x1a[x1==-5]="0.03"
  x1a[x1==-4]="0.06"
  x1a[x1==-3]="0.12"
  x1a[x1==-2]="0.25"
  x1a[x1==-1]="0.5"
  x1a[x1==0]="1"
  x1a[x1==1]="2"
  x1a[x1==2]="4"
  x1a[x1==3]="8"
  x1a[x1==4]="16"
 
  PEN_SIGN=sx1
  PEN=x1a

  SIR_M=rep("U", n1)
  SIR_M[x1<=-4]="S"
  SIR_M[x1>-4]="R"  
  WGS_PEN_SIR_Meningitis=SIR_M  

  SIR_NM=rep("U", n1)
  SIR_NM[x1<=1]="S"
  SIR_NM[x1==2]="I"
  SIR_NM[x1>=3]="R"
  WGS_PEN_SIR_Nonmeningitis=SIR_NM

  #Format AMO
  i=3
  x1=as.integer(round(log2(as.numeric(m3a[, i]))))
  sx1=rep("=", n1)
  sx1[x1<=-5]="<="
  x1[x1<=-5]=-5
  sx1[x1>3]=">"
  x1[x1>3]=3
  x1a=rep(NA, n1)
  x1a[x1==-5]="0.03"
  x1a[x1==-4]="0.06"
  x1a[x1==-3]="0.12"
  x1a[x1==-2]="0.25"
  x1a[x1==-1]="0.5"
  x1a[x1==0]="1"
  x1a[x1==1]="2"
  x1a[x1==2]="4"
  x1a[x1==3]="8"
  x1a[x1==4]="16"
  AMO_SIGN=sx1
  AMO=x1a

  SIR_NM=rep("U", n1)
  SIR_NM[x1<=1]="S"
  SIR_NM[x1==2]="I"
  SIR_NM[x1>=3]="R"
  WGS_AMO_SIR=SIR_NM

  #Format MER
  i=4
  x1=as.integer(round(log2(as.numeric(m3a[, i]))))
  sx1=rep("=", n1)
  sx1[x1<=-4]="<="
  x1[x1<=-4]=-4
  sx1[x1>0]=">"
  x1[x1>0]=0
  x1a=rep(NA, n1)
  x1a[x1==-5]="0.03"
  x1a[x1==-4]="0.06"
  x1a[x1==-3]="0.12"
  x1a[x1==-2]="0.25"
  x1a[x1==-1]="0.5"
  x1a[x1==0]="1"
  x1a[x1==1]="2"
  x1a[x1==2]="4"
  x1a[x1==3]="8"
  x1a[x1==4]="16"
  MER_SIGN=sx1
  MER=x1a

  SIR_NM=rep("U", n1)
  SIR_NM[x1<=-2]="S"
  SIR_NM[x1==-1]="I"
  SIR_NM[x1>=0]="R"
  WGS_MER_SIR=SIR_NM


  #Format TAX
  i=5
  x1=as.integer(round(log2(as.numeric(m3a[, i]))))
  sx1=rep("=", n1)
  sx1[x1<=-4]="<="
  x1[x1<=-4]=-4
  sx1[x1>4]=">"
  x1[x1>4]=4
  x1a=rep(NA, n1)
  x1a[x1==-5]="0.03"
  x1a[x1==-4]="0.06"
  x1a[x1==-3]="0.12"
  x1a[x1==-2]="0.25"
  x1a[x1==-1]="0.5"
  x1a[x1==0]="1"
  x1a[x1==1]="2"
  x1a[x1==2]="4"
  x1a[x1==3]="8"
  x1a[x1==4]="16"
  TAX_SIGN=sx1
  TAX=x1a

  SIR_M=rep("U", n1)
  SIR_M[x1<=-1]="S"
  SIR_M[x1==0]="I"
  SIR_M[x1>=1]="R"
  WGS_TAX_SIR_Meningitis=SIR_M


  SIR_NM=rep("U", n1)
  SIR_NM[x1<=0]="S"
  SIR_NM[x1==1]="I"
  SIR_NM[x1>=2]="R"
  WGS_TAX_SIR_Nonmeningitis=SIR_NM


  #Format CFT
  i=6
  x1=as.integer(round(log2(as.numeric(m3a[, i]))))
  sx1=rep("=", n1)
  sx1[x1<=-1]="<="
  x1[x1<=-1]=-1
  sx1[x1>3]=">"
  x1[x1>3]=3
  x1a=rep(NA, n1)
  x1a[x1==-5]="0.03"
  x1a[x1==-4]="0.06"
  x1a[x1==-3]="0.12"
  x1a[x1==-2]="0.25"
  x1a[x1==-1]="0.5"
  x1a[x1==0]="1"
  x1a[x1==1]="2"
  x1a[x1==2]="4"
  x1a[x1==3]="8"
  x1a[x1==4]="16"
  CFT_SIGN=sx1
  CFT=x1a
 
  SIR_M=rep("U", n1)
  SIR_M[x1<=-1]="S"
  SIR_M[x1==0]="I"
  SIR_M[x1>=1]="R"
  WGS_CFT_SIR_Meningitis=SIR_M

  SIR_NM=rep("U", n1)
  SIR_NM[x1<=0]="S"
  SIR_NM[x1==1]="I"
  SIR_NM[x1>=2]="R"
  WGS_CFT_SIR_Nonmeningitis=SIR_NM

  #Format CFX
  i=7
  x1=as.integer(round(log2(as.numeric(m3a[, i]))))
  sx1=rep("=", n1)
  sx1[x1<=-1]="<="
  x1[x1<=-1]=-1
  sx1[x1>1]=">"
  x1[x1>1]=1
  x1a=rep(NA, n1)
  x1a[x1==-5]="0.03"
  x1a[x1==-4]="0.06"
  x1a[x1==-3]="0.12"
  x1a[x1==-2]="0.25"
  x1a[x1==-1]="0.5"
  x1a[x1==0]="1"
  x1a[x1==1]="2"
  x1a[x1==2]="4"
  x1a[x1==3]="8"
  x1a[x1==4]="16"
  CFX_SIGN=sx1
  CFX=x1a

  SIR_NM=rep("U", n1)
  SIR_NM[x1<=-1]="S"
  SIR_NM[x1==0]="I"
  SIR_NM[x1>=1]="R"
  WGS_CFX_SIR=SIR_NM

  m4=cbind(SampleID, 
          PEN_SIGN, PEN, WGS_PEN_SIR_Meningitis, WGS_PEN_SIR_Nonmeningitis, 
          AMO_SIGN, AMO, WGS_AMO_SIR, 
          MER_SIGN, MER, WGS_MER_SIR,
          TAX_SIGN, TAX, WGS_TAX_SIR_Meningitis, WGS_TAX_SIR_Nonmeningitis,
          CFT_SIGN, CFT, WGS_CFT_SIR_Meningitis, WGS_CFT_SIR_Nonmeningitis,
          CFX_SIGN, CFX, WGS_CFX_SIR)

  write.csv(m4, fout, row.names=F)

}

args <- commandArgs(TRUE)
fin = args[1]

#print (c(train_file, cwd))
MIC_format(fin)

