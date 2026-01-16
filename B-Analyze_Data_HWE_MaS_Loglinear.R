rm(list=ls())

source("../run_analysis_functions_V3.R")


## Read in the datasets
Hybrid_Asia=read.csv("Data/Hybrid_Asia_counts_loglin_NoE.csv")
Hybrid_Europe=read.csv("Data/Hybrid_Europe_counts_loglin_NoE.csv")
Hybrid_Both=read.csv("Data/Hybrid_BothContinent_counts_loglin_NoE.csv")

## SNPs
snpnames=colnames(Hybrid_Asia)[7:ncol(Hybrid_Both)]

rs=unlist(strsplit(snpnames,split="_"))[2*(1:(length(snpnames)/2))-1]
gene=unlist(strsplit(snpnames,split="_"))[2*(1:(length(snpnames)/2))]

## Population genetic tests

# HWE 
HWE_res_Asia=matrix(nrow=length(snpnames),ncol=6)
HWE_res_Europe=HWE_res_Asia
HWE_res_Both=HWE_res_Asia

for (i in 7:ncol(Hybrid_Both))
{
  
  # GENEVA first
  # Asia
  dat=subset(Hybrid_Asia[,c(3,4,i)],Hybrid_Asia$D==1)
  colnames(dat)[3]="count"
  HWE_res_Asia[i-6,1:3]=getHWE(dat)
  
  # Europe
  dat=subset(Hybrid_Europe[,c(3,4,i)],Hybrid_Europe$D==1)
  colnames(dat)[3]="count"
  HWE_res_Europe[i-6,1:3]=getHWE(dat)
  
  # Both
  dat=subset(Hybrid_Both[,c(3,4,i)],Hybrid_Both$D==1)
  colnames(dat)[3]="count"
  HWE_res_Both[i-6,1:3]=getHWE(dat)
  
  # 1000G second
  # Asia
  dat=subset(Hybrid_Asia[,c(3,4,i)],Hybrid_Asia$D==0)
  colnames(dat)[3]="count"
  HWE_res_Asia[i-6,4:6]=getHWE(dat)
  
  # Europe
  dat=subset(Hybrid_Europe[,c(3,4,i)],Hybrid_Europe$D==0)
  colnames(dat)[3]="count"
  HWE_res_Europe[i-6,4:6]=getHWE(dat)
  
  # Both
  dat=subset(Hybrid_Both[,c(3,4,i)],Hybrid_Both$D==0)
  colnames(dat)[3]="count"
  HWE_res_Both[i-6,4:6]=getHWE(dat)
  
}

HWE_res=data.frame(gene,rs,HWE_res_Asia,HWE_res_Europe,HWE_res_Both)
colnames(HWE_res)=c("GENE","SNP","GENEVA.AS.M","GENEVA.AS.F","GENEVA.AS.Both",
                    "GENEVA.EU.M","GENEVA.EU.F","GENEVA.EU.Both",
                    "GENEVA.All.M","GENEVA.All.F","GENEVA.All.Both",
                    "1000G.AS.M","1000G.AS.F","1000G.AS.Both",
                    "1000G.EU.M","1000G.EU.F","1000G.EU.Both",
                    "1000G.All.M","1000G.All.F","1000G.All.Both")
write.csv(HWE_res,"Results/HWE_AllSNPs.csv",row.names=F)




# MaS
MaS_res_Asia=matrix(nrow=length(snpnames),ncol=2)
MaS_res_Europe=MaS_res_Asia
MaS_res_Both=MaS_res_Asia

for (i in 7:ncol(Hybrid_Both))
{
  # GENEVA first
  # Asia
  dat=subset(Hybrid_Asia[,c(3,4,i)],Hybrid_Asia$D==1)
  colnames(dat)[3]="count"
  MaS_res_Asia[i-6,1]=MaSlrt(dat)
  
  # Europe
  dat=subset(Hybrid_Europe[,c(3,4,i)],Hybrid_Europe$D==1)
  colnames(dat)[3]="count"
  MaS_res_Europe[i-6,1]=MaSlrt(dat)
  
  # Both
  dat=subset(Hybrid_Both[,c(3,4,i)],Hybrid_Both$D==1)
  colnames(dat)[3]="count"
  MaS_res_Both[i-6,1]=MaSlrt(dat)
  
  # 1000G first
  # Asia
  dat=subset(Hybrid_Asia[,c(3,4,i)],Hybrid_Asia$D==0)
  colnames(dat)[3]="count"
  MaS_res_Asia[i-6,2]=MaSlrt(dat)
  
  # Europe
  dat=subset(Hybrid_Europe[,c(3,4,i)],Hybrid_Europe$D==0)
  colnames(dat)[3]="count"
  MaS_res_Europe[i-6,2]=MaSlrt(dat)
  
  # Both
  dat=subset(Hybrid_Both[,c(3,4,i)],Hybrid_Both$D==0)
  colnames(dat)[3]="count"
  MaS_res_Both[i-6,2]=MaSlrt(dat)
  
}

MaS_res=data.frame(gene,rs,MaS_res_Asia,MaS_res_Europe,MaS_res_Both)
colnames(MaS_res)=c("GENE","SNP","GENEVA.AS.MaS","1000G.AS.MaS",
                    "GENEVA.EU.MaS","1000G.EU.MaS",
                    "GENEVA.All.MaS","1000G.All.MaS")
write.csv(MaS_res,"Results/MaS_AllSNPs.csv",row.names=F)


# Loglinear tests

# Look at test of C and M in two groups separate and compare to 
# combined - essentially created PopStrat
Loglin_res_Asia=matrix(nrow=length(snpnames),ncol=6)
colnames(Loglin_res_Asia)=c("HW.M","HW.C","MS.M","MS.C","MaS.M","MaS.C")
Loglin_res_Europe=Loglin_res_Asia
Loglin_res_Both=Loglin_res_Asia

for (i in 7:ncol(Hybrid_Both))
{
  
  # Asia
  
  dat=subset(Hybrid_Asia[,c(1:6,i)],Hybrid_Asia$D==1)
  colnames(dat)[7]="count"
  res=runLoglin("HW", c("M","C"), dat)
  Loglin_res_Asia[i-6,1:2]=res$pvals
  
  res=runLoglin("MS", c("M","C"), dat)
  Loglin_res_Asia[i-6,3:4]=c(res$pvals)
  
  dat=Hybrid_Asia[,c(1:6,i)]
  colnames(dat)[7]="count"
  res=runLoglin("MaS", c("M","C"), dat)
  Loglin_res_Asia[i-6,5:6]=c(res$pvals)
  
  # Europe
  dat=subset(Hybrid_Europe[,c(1:6,i)],Hybrid_Europe$D==1)
  colnames(dat)[7]="count"
  res=runLoglin("HW", c("M","C"), dat)
  Loglin_res_Europe[i-6,1:2]=res$pvals
  
  res=runLoglin("MS", c("M","C"), dat)
  Loglin_res_Europe[i-6,3:4]=c(res$pvals)
  
  dat=Hybrid_Europe[,c(1:6,i)]
  colnames(dat)[7]="count"
  res=runLoglin("MaS", c("M","C"), dat)
  Loglin_res_Europe[i-6,5:6]=c(res$pvals)
  
  
  # Both
  dat=subset(Hybrid_Both[,c(1:6,i)],Hybrid_Both$D==1)
  colnames(dat)[7]="count"
  res=runLoglin("HW", c("M","C"), dat)
  Loglin_res_Both[i-6,1:2]=res$pvals
  
  res=runLoglin("MS", c("M","C"), dat)
  Loglin_res_Both[i-6,3:4]=c(res$pvals)
  
  dat=Hybrid_Both[,c(1:6,i)]
  colnames(dat)[7]="count"
  res=runLoglin("MaS", c("M","C"), dat)
  Loglin_res_Both[i-6,5:6]=c(res$pvals)
  
}




Loglin_res=data.frame(gene,rs,Loglin_res_Asia,Loglin_res_Europe,Loglin_res_Both)
colnames(Loglin_res)=c("GENE","SNP",paste0("AS.",colnames(Loglin_res_Asia)),
                       paste0("EU.",colnames(Loglin_res_Europe)),
                       paste0("Both.",colnames(Loglin_res_Both)))
write.csv(Loglin_res[,c(1,2,2*(2:10)-1)],"Results/Loglin_AllSNPs_M.csv", row.names=F)
write.csv(Loglin_res[,c(1,2,2*(2:10))],"Results/Loglin_AllSNPs_C.csv",row.names=F)



# MAF for s4646703
i=85
dat=subset(Hybrid_Asia[,c(3,4,i)],Hybrid_Asia$D==0)
colnames(dat)[3]="count"
mom=tapply(dat$count, dat$M, sum)
dad=tapply(dat$count, dat$F, sum)
both=dad+mom
(2*both[3]+both[2])/(2*sum(both))

dat=subset(Hybrid_Asia[,c(3,4,i)],Hybrid_Asia$D==1)
colnames(dat)[3]="count"
mom=tapply(dat$count, dat$M, sum)
dad=tapply(dat$count, dat$F, sum)
both=dad+mom
(2*both[3]+both[2])/(2*sum(both))


dat=subset(Hybrid_Europe[,c(3,4,i)],Hybrid_Europe$D==0)
colnames(dat)[3]="count"
mom=tapply(dat$count, dat$M, sum)
dad=tapply(dat$count, dat$F, sum)
both=dad+mom
(2*both[3]+both[2])/(2*sum(both))

dat=subset(Hybrid_Europe[,c(3,4,i)],Hybrid_Europe$D==1)
colnames(dat)[3]="count"
mom=tapply(dat$count, dat$M, sum)
dad=tapply(dat$count, dat$F, sum)
both=dad+mom
(2*both[3]+both[2])/(2*sum(both))


