## This script is useqd to run all the scenarios with mating asymmetry
## and no population substructure.

## Kelly Burkett
## April 18, 2023
## Modification Notes
## - Added HW/PS/MaS tests; Now write out summary dataset with mean risk and 
##   type 1/power across all simulations (June 1, 2023)
## - Re-run with lower main and interaction effects and changed to 
##   run_analysis_functions_V4 which returns "M" effect from unstratified 
##   analyses (July 6, 2023)

options(warn=1)

## Source the analysis functions
source("simulate_data_functions_V3.R")
source("run_analysis_functions_V3.R")
library(Haplin)

## Set directories
prefix="NoPS_MaS_C1C2only_LowerEffects/"
setwd(paste0("../Run_directories/",prefix))
outdir=paste0("../../Results/",prefix)

## Save the seed or use an existing seed
useseed=FALSE
seedname="seed_1.txt"
date=paste(strsplit(date(),split=" ")[[1]],collapse="_")


if (!useseed){
  if(!exists(".Random.seed")) set.seed(NULL)
  s1 <- .Random.seed
  write(s1,file=paste0(outdir,"seed_",date,".txt"),ncol=length(s1))
  write(runif(1),file=paste0(outdir,"rng_check",date,".txt"))

} else {
  seeds=scan(paste0(seedname,what=integer()))
  .Random.seed=seeds
  write(seeds,file=paste0(outdir,"seed_",date,"_2.txt"),ncol=length(seeds))
  write(runif(1),file=paste0(outdir,"rng_check",date,"_2.txt"))
}


## Setup parameters for simulations
alpha.type1=0.05
numsim=10000 # Number of replicates per set of parameters
NumCase=1000 # Number of case trios
NumControl=c(0,500,1000) # Number of control trios
maf=0.3 # Minor allele frequency
ef=0.3  # Environmental variable frequency

# Each element of the list is the risk ratio, of the C, M and ME
# effect per risk allele. Assuming a multiplicative model. 
effectmat=list(c(1,1,1),c(1.2,1,1),c(1.2,1.2,1),c(1.2,1,1.4))
masmat=list(c(0.85,0.85,1),c(0.95,0.95,1),c(1.05,1.05,1),c(1.15,1.15,1))

namesres=c("Controls","Ceffect","Meffect","MEeffect",
                  "C1","C2","C4",
                   "EMIM.HW.C","EMIM.HW.M","EMIM.HW.ME",
                   "EMIM.MS.C","EMIM.MS.M","EMIM.MS.ME",
                   "EMIM.MaS.C","EMIM.MaS.M","EMIM.MaS.ME",
                   "Haplin.C","Haplin.M","Haplin.ME",
                   "Loglin.HW.C","Loglin.HW.M","Loglin.HW.ME",
                   "Loglin.MS.C","Loglin.MS.M","Loglin.MS.ME",
                   "Loglin.MaS.C","Loglin.MaS.M","Loglin.MaS.ME")
namestests=c("HW.Mothers.Case","HW.Fathers.Case","HW.Both.Case",
                   "HW.Mothers.Control","HW.Fathers.Control","HW.Both.Control",
                   "PS.HW","PS.MS","PS.MaS","MS")

resvec=vector(length=length(namesres))
names(resvec)=namesres
pvalvec=vector(length=length(namesres))
names(pvalvec)=namesres
testvec=vector(length=length(namestests))
names(testvec)=namestests


write(namesres,paste0(outdir,"Mean_RiskEstimates.dat"),ncol=length(namesres),append=F)
write(namesres,paste0(outdir,"Type1_Power_Risk.dat"),ncol=length(namesres),append=F)
write(namestests,paste0(outdir,"Type1_Power_HW_PS_MaS.dat"),ncol=length(namesres),append=F)



# Loop to run all scenarios
Sys.time()

for (i in 1:length(NumControl)){
  
  for (j in 1:length(effectmat)){
    
    for (k in 1:length(masmat)){

    cat(paste("\n\nStarting Controls:", NumControl[i],"; Effects:", effectmat[[j
]][1]," ",effectmat[[j]][2]," ", effectmat[[j]][3],"; MaS: ", masmat[[k]][1]," ", masmat[[k]][2],"\n\n"))

      
      mas=masmat[[k]]
  
      C=effectmat[[j]][1]
      M=effectmat[[j]][2]
      ME=effectmat[[j]][3]

      # Summarized result and pval matrix
      resvec[1]=NumControl[i]
      resvec[2:4]=effectmat[[j]]
      resvec[5:7]=mas

      pvalvec[1]=NumControl[i]
      pvalvec[2:4]=effectmat[[j]]
      pvalvec[5:7]=mas
       
      testvec[1]=NumControl[i]
      testvec[2:4]=effectmat[[j]]
      testvec[5:7]=mas


      # Run simulations
      if (NumControl[i]==0){ # No control trios
        res=runTrioReps(nreps=numsim, ntrios=NumCase,
                        R=c(1,C,C^2), S=c(1,M,M^2), V=c(1,ME,ME^2),mtCoef=mas,
                        effects=c("M","C","E:M"),
                        includeE=TRUE, prE=ef, envint="Mother", 
                        includePopStrat=FALSE,
                        includeControl=FALSE,
                        Loglin=TRUE, loglinMT=c("HW","MS"), 
                        Haplin=TRUE, haplinControls=FALSE,
                        EMIM=TRUE, EMIMMT=c("HW","MS"), 
                        PStest=FALSE, HWtest=TRUE, MStest=TRUE)
      } else { # Control trios
        NumTrios=NumCase+NumControl[i]
        propControl=NumControl[i]/NumTrios
        res=runTrioReps(nreps=numsim, ntrios=NumTrios,
                        R=c(1,C,C^2), S=c(1,M,M^2), V=c(1,ME,ME^2),mtCoef=mas,
                        effects=c("M","C","E:M"),
                        includeE=TRUE, prE=ef, envint="Mother", 
                        includePopStrat=FALSE,
                        includeControl=TRUE, prControl=propControl,
                        Loglin=TRUE, loglinMT=c("HW","MS","MaS"), 
                        Haplin=TRUE, haplinControls=TRUE,
                        EMIM=TRUE, EMIMMT=c("HW","MS","MaS"), 
                        PStest=FALSE, HWtest=TRUE, MStest=TRUE)
      }
    
# Write out summarized type 1 error/power
      
      # EMIM
      pvalvec["EMIM.HW.C"]=sum(res$pvalEMIM[,"HW.C"]<alpha.type1)/numsim
      pvalvec["EMIM.HW.M"]=sum(res$pvalEMIM[,"HW.M"]<alpha.type1)/numsim
      pvalvec["EMIM.HW.ME"]=sum(res$pvalEMIM[,"HW.E:M"]<alpha.type1)/numsim
      
      pvalvec["EMIM.MS.C"]=sum(res$pvalEMIM[,"MS.C"]<alpha.type1)/numsim
      pvalvec["EMIM.MS.M"]=sum(res$pvalEMIM[,"MS.M"]<alpha.type1)/numsim
      pvalvec["EMIM.MS.ME"]=sum(res$pvalEMIM[,"MS.E:M"]<alpha.type1)/numsim
      
      if (NumControl[i]>0){
        
        pvalvec["EMIM.MaS.C"]=sum(res$pvalEMIM[,"MaS.C"]<alpha.type1)/numsim
        pvalvec["EMIM.MaS.M"]=sum(res$pvalEMIM[,"MaS.M"]<alpha.type1)/numsim
        pvalvec["EMIM.MaS.ME"]=sum(res$pvalEMIM[,"MaS.E:M"]<alpha.type1)/numsim
        
      }
      
      # Haplin
      pvalvec["Haplin.C"]=sum(res$pvalHaplin[,"C"]<alpha.type1)/numsim
      pvalvec["Haplin.M"]=sum(res$pvalHaplin[,"M"]<alpha.type1)/numsim
      pvalvec["Haplin.ME"]=sum(res$pvalHaplin[,"E:M"]<alpha.type1)/numsim
      
      
      # Loglinear
      pvalvec["Loglin.HW.C"]=sum(res$pvalLoglin[,"HW.C"]<alpha.type1)/numsim
      pvalvec["Loglin.HW.M"]=sum(res$pvalLoglin[,"HW.M"]<alpha.type1)/numsim
      pvalvec["Loglin.HW.ME"]=sum(res$pvalLoglin[,"HW.E:M"]<alpha.type1)/numsim
      
      pvalvec["Loglin.MS.C"]=sum(res$pvalLoglin[,"MS.C"]<alpha.type1)/numsim
      pvalvec["Loglin.MS.M"]=sum(res$pvalLoglin[,"MS.M"]<alpha.type1)/numsim
      pvalvec["Loglin.MS.ME"]=sum(res$pvalLoglin[,"MS.E:M"]<alpha.type1)/numsim
      
      if (NumControl[i]>0){
        
        pvalvec["Loglin.MaS.C"]=sum(res$pvalLoglin[,"MaS.C"]<alpha.type1)/numsim
        pvalvec["Loglin.MaS.M"]=sum(res$pvalLoglin[,"MaS.M"]<alpha.type1)/numsim
        pvalvec["Loglin.MaS.ME"]=sum(res$pvalLoglin[,"MaS.E:M"]<alpha.type1)/numsim
        
      }
      
      
      # HW/MS/PS tests
      
      testvec["MS"]=sum(res$pvalMaS<alpha.type1)/numsim
      testvec["PS.HW"]=NA
      testvec["PS.MS"]=NA
      testvec["PS.MaS"]=NA
      
      if (NumControl[i]>0){
	# Don't provide when there are 10,000 simulations
        #testvec["PS.HW"]=sum(res$pvalPS[,"HW"]<alpha.type1)/numsim
        #testvec["PS.MS"]=sum(res$pvalPS[,"MS"]<alpha.type1)/numsim
        #testvec["PS.MaS"]=sum(res$pvalPS[,"MaS"]<alpha.type1)/numsim
        testvec["HW.Mothers.Case"]=sum(res$pvalHW[,"Mothers.Case"]<alpha.type1)/numsim
        testvec["HW.Fathers.Case"]=sum(res$pvalHW[,"Fathers.Case"]<alpha.type1)/numsim
        testvec["HW.Both.Case"]=sum(res$pvalHW[,"Both.Case"]<alpha.type1)/numsim
        testvec["HW.Mothers.Control"]=sum(res$pvalHW[,"Mothers.Control"]<alpha.type1)/numsim
        testvec["HW.Fathers.Control"]=sum(res$pvalHW[,"Fathers.Control"]<alpha.type1)/numsim
        testvec["HW.Both.Control"]=sum(res$pvalHW[,"Both.Control"]<alpha.type1)/numsim
      } else {
        testvec["HW.Mothers.Case"]=sum(res$pvalHW[,"Mothers"]<alpha.type1)/numsim
        testvec["HW.Fathers.Case"]=sum(res$pvalHW[,"Fathers"]<alpha.type1)/numsim
        testvec["HW.Both.Case"]=sum(res$pvalHW[,"Both"]<alpha.type1)/numsim
      }
      
      # Write out summarized coefficient estimates
      
      # EMIM
      resvec["EMIM.HW.C"]=mean(res$resEMIM[,"HW.C"])
      resvec["EMIM.HW.M"]=mean(apply(res$resEMIM[,c("HW.M[E=0]","HW.M[E=1]")],1,mean))
      resvec["EMIM.HW.ME"]=mean(res$resEMIM[,"HW.M[E=1]/M[E=0]"])
      
      resvec["EMIM.MS.C"]=mean(res$resEMIM[,"MS.C"])
      resvec["EMIM.MS.M"]=mean(apply(res$resEMIM[,c("MS.M[E=0]","MS.M[E=1]")],1,mean))
      resvec["EMIM.MS.ME"]=mean(res$resEMIM[,"MS.M[E=1]/M[E=0]"])
      
      if (NumControl[i]>0){
        
        resvec["EMIM.MaS.C"]=mean(res$resEMIM[,"MaS.C"])
        resvec["EMIM.MaS.M"]=mean(apply(res$resEMIM[,c("MaS.M[E=0]","MaS.M[E=1]")],1,mean))
        resvec["EMIM.MaS.ME"]=mean(res$resEMIM[,"MaS.M[E=1]/M[E=0]"])
        
      }
      
      # Haplin
      resvec["Haplin.C"]=mean(res$resHaplin[,"C"])
      resvec["Haplin.M"]=mean(apply(res$resHaplin[,c("M[E=0]","M[E=1]")],1,mean))
      resvec["Haplin.ME"]=mean(res$resHaplin[,"M[E=1]/M[E=0]"])
      
      
      # Loglinear
      resvec["Loglin.HW.C"]=mean(res$resLoglin[,"HW.C"])
      resvec["Loglin.HW.M"]=mean(res$resLoglin[,"HW.M"])
      resvec["Loglin.HW.ME"]=mean(res$resLoglin[,"HW.E:M"])
      
      resvec["Loglin.MS.C"]=mean(res$resLoglin[,"MS.C"])
      resvec["Loglin.MS.M"]=mean(res$resLoglin[,"MS.M"])
      resvec["Loglin.MS.ME"]=mean(res$resLoglin[,"MS.E:M"])
      
      if (NumControl[i]>0){
        
        resvec["Loglin.MaS.C"]=mean(res$resLoglin[,"MaS.C"])
        resvec["Loglin.MaS.M"]=mean(res$resLoglin[,"MaS.M"])
        resvec["Loglin.MaS.ME"]=mean(res$resLoglin[,"MaS.E:M"])
        
      }

      write(round(resvec,5),paste0(outdir,"Mean_RiskEstimates.dat"),ncol=length(resvec), append=T)
      write(round(pvalvec,5),paste0(outdir,"Type1_Power_Risk.dat"),ncol=length(resvec), append=T)
      write(round(testvec,5),paste0(outdir,"Type1_Power_HW_PS_MaS.dat"),ncol=length(resvec), append=T)         

    }
  }
}

Sys.time()
