

rm(list=ls())

source("../simulate_data_functions_V3.R")
source("../run_analysis_functions_V3.R")


# Read in the 1000G and replace the disease column with 2
Asia_dat=read.table("Data/1000G_Asian.ped",colClasses="character")
Asia_dat[,6]=rep(2,nrow(Asia_dat))
write.table(Asia_dat,"Data/1000G_Asian_Disease2.ped",quote=F,col=F,row=F)
system("cp Data/1000G_Asian.map Data/1000G_Asian_Disease2.map")

Europe_dat=read.table("Data/1000G_Europe.ped",colClasses="character")
Europe_dat[,6]=rep(2,nrow(Europe_dat))
write.table(Europe_dat,"Data/1000G_Europe_Disease2.ped",quote=F,col=F,row=F)
system("cp Data/1000G_Europe.map Data/1000G_Europe_Disease2.map")


## Set up the datasets with the controls 

## Use PREMIM to get overall counts
system("~/StatGenTools/emim-v3.22-code/premim Data/GENEVA_Asian.ped")
system("mv caseparenttrios.dat Data/GENEVA_Asian_trio_counts.dat")

system("~/StatGenTools/emim-v3.22-code/premim Data/1000g_Asian_Disease2.ped")
system("mv caseparenttrios.dat Data/1000G_Asian_trio_counts.dat")

system("~/StatGenTools/emim-v3.22-code/premim Data/GENEVA_Europe.ped")
system("mv caseparenttrios.dat Data/GENEVA_Europe_trio_counts.dat")

system("~/StatGenTools/emim-v3.22-code/premim Data/1000G_Europe_Disease2.ped")
system("mv caseparenttrios.dat Data/1000G_Europe_trio_counts.dat")

system("rm case*.dat con*.dat emimparams.dat premim.log")


# EMIM counts compared to ours:
convertmat=data.frame(EMIM=1:15,LOGLIN=c(15:11,7,6,10:8,5:1))


# Create GENEVA datasets for log-linear
precols=mtmat()[,1:5]

#Snps in both
snpnames_EUR=read.table("Data/GENEVA_Europe.map")[,2]
snpnames_AS=read.table("Data/GENEVA_Asian.map")[,2]
snpnames_both=intersect(snpnames_AS,snpnames_EUR)

# Euro 
GENEVA_Europe=read.table("Data/GENEVA_Europe_trio_counts.dat",skip=1)
GENEVA_Europe=GENEVA_Europe[,-1]
GENEVA_Europe=GENEVA_Europe[,convertmat$LOGLIN]
GENEVA_Europe=data.frame(t(GENEVA_Europe))
colnames(GENEVA_Europe)=snpnames_EUR
GENEVA_Europe=GENEVA_Europe[,snpnames_both]

# Asian 
GENEVA_Asia=read.table("Data/GENEVA_Asian_trio_counts.dat",skip=1)
GENEVA_Asia=GENEVA_Asia[,-1]
GENEVA_Asia=GENEVA_Asia[,convertmat$LOGLIN]
GENEVA_Asia=data.frame(t(GENEVA_Asia))
colnames(GENEVA_Asia)=snpnames_AS
GENEVA_Asia=GENEVA_Asia[,snpnames_both]


# Combine them and write out the files
GENEVA_Both=GENEVA_Asia+GENEVA_Europe
write.csv(data.frame(precols,GENEVA_Asia),
          "Data/GENEVA_Asia_counts_loglin_NoE.csv",row.names=F)
write.csv(data.frame(precols,GENEVA_Europe),
          "Data/GENEVA_Europe_counts_loglin_NoE.csv",row.names=F)
write.csv(data.frame(precols,GENEVA_Both),
          "Data/GENEVA_Both_counts_loglin_NoE.csv",row.names=F)



## Create datasets with the controls

# Verify the SNPs are the same based on location
snps_EUR_1000G=read.table("Data/1000G_Europe.map")
snps_EUR_GENEVA=read.table("Data/GENEVA_Europe.map")
snps_EUR_1000G[,4]-snps_EUR_GENEVA[,4]

snps_AS_1000G=read.table("Data/1000G_Asian.map")
snps_AS_GENEVA=read.table("Data/GENEVA_Asian.map")
snps_AS_1000G[,4]-snps_AS_GENEVA[,4]


# They are the same. We can use the GENEVA maps since they have
# SNP names and rs numbers
snpnames_EUR=read.table("Data/GENEVA_Europe.map")[,2]
snpnames_AS=read.table("Data/GENEVA_Asian.map")[,2]
snpnames_both=intersect(snpnames_AS,snpnames_EUR)

# Euro 
Europe_1000G=read.table("Data/1000G_Europe_trio_counts.dat",skip=1)
Europe_1000G=Europe_1000G[,-1]
Europe_1000G=Europe_1000G[,convertmat$LOGLIN]
Europe_1000G=data.frame(t(Europe_1000G))
colnames(Europe_1000G)=snpnames_EUR
Europe_1000G=Europe_1000G[,snpnames_both]

# Asian 
Asia_1000G=read.table("Data/1000G_Asian_trio_counts.dat",skip=1)
Asia_1000G=Asia_1000G[,-1]
Asia_1000G=Asia_1000G[,convertmat$LOGLIN]
Asia_1000G=data.frame(t(Asia_1000G))
colnames(Asia_1000G)=snpnames_AS
Asia_1000G=Asia_1000G[,snpnames_both]


# Add 1000G to GENEVA in one dataset
D=c(rep(1,15), rep(0,15))
precols=rbind(precols,precols)
write.csv(data.frame(precols,D,rbind(GENEVA_Asia,Asia_1000G)),
          "Data/Hybrid_Asia_counts_loglin_NoE.csv",row.names=F)

write.csv(data.frame(precols,D,rbind(GENEVA_Europe,Europe_1000G)),
          "Data/Hybrid_Europe_counts_loglin_NoE.csv",row.names=F)

# Create a stratified dataset
All_1000G=Asia_1000G+Europe_1000G
write.csv(data.frame(precols,D,rbind(GENEVA_Both,All_1000G)),
          "Data/Hybrid_BothContinent_counts_loglin_NoE.csv",row.names=F)






