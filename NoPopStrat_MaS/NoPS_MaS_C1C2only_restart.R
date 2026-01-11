## This script is used to fill in scenarios with mating asymmetry
## and no population substructure.

## Kelly Burkett
## April 18, 2023
## Modification Notes
## - Added HW/PS/MaS tests; Now write out summary dataset with mean risk and 
##   type 1/power across all simulations (June 1, 2023)
## - Modified the script to runs subsets at a time (July 11, 2025)

options(warn=1)

## Source the analysis functions
source("simulate_data_functions_V3.R")
source("run_analysis_functions_V3.R")
library(Haplin)

## Set directories
prefix="NoPS_MaS_C1C2only_rerun/"
setwd(paste0("../Run_directories/",prefix))
outdir=paste0("../../Results/",prefix)

## Setup parameters for simulations
alpha=0.05
numsim=2000 # Number of replicates per set of parameters
NumCase=1000 # Number of case trios
NumControl=500 # Number of control trios
maf=0.3 # Minor allele frequency
ef=0.3  # Environmental variable frequency

# Each element of the list is the risk ratio, of the C, M and ME
# effect per risk allele. Assuming a multiplicative model. 
effectmat=list(c(1,1,1))
masmat=list(c(1.05,1.05,1))



numsettings=length(NumControl)*length(effectmat)*length(masmat)
resmat=matrix(nrow=numsettings,ncol=28)
colnames(resmat)=c("Controls","Ceffect","Meffect","MEeffect",
                  "C1","C2","C4",
                   "EMIM.HW.C","EMIM.HW.M","EMIM.HW.ME",
                   "EMIM.MS.C","EMIM.MS.M","EMIM.MS.ME",
                   "EMIM.MaS.C","EMIM.MaS.M","EMIM.MaS.ME",
                   "Haplin.C","Haplin.M","Haplin.ME",
                   "Loglin.HW.C","Loglin.HW.M","Loglin.HW.ME",
                   "Loglin.MS.C","Loglin.MS.M","Loglin.MS.ME",
                   "Loglin.MaS.C","Loglin.MaS.M","Loglin.MaS.ME")
pvalmat=matrix(nrow=numsettings,ncol=28)
colnames(pvalmat)=c("Controls","Ceffect","Meffect","MEeffect",
                    "C1","C2","C4",
                    "EMIM.HW.C","EMIM.HW.M","EMIM.HW.ME",
                    "EMIM.MS.C","EMIM.MS.M","EMIM.MS.ME",
                    "EMIM.MaS.C","EMIM.MaS.M","EMIM.MaS.ME",
                    "Haplin.C","Haplin.M","Haplin.ME",
                    "Loglin.HW.C","Loglin.HW.M","Loglin.HW.ME",
                    "Loglin.MS.C","Loglin.MS.M","Loglin.MS.ME",
                   "Loglin.MaS.C","Loglin.MaS.M","Loglin.MaS.ME")
testmat=matrix(nrow=numsettings,ncol=17)
colnames(testmat)=c("Controls","Ceffect","Meffect","MEeffect",
                    "C1","C2","C4",
                    "HW.Mothers.Case","HW.Fathers.Case","HW.Both.Case",
                    "HW.Mothers.Control","HW.Fathers.Control","HW.Both.Control",
                    "PS.HW","PS.MS","PS.MaS","MS")

counter=1



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
      resmat[counter,1]=NumControl[i]
      resmat[counter,2]=C
      resmat[counter,3]=M
      resmat[counter,4]=ME
      resmat[counter,5:7]=mas
      
      pvalmat[counter,1]=NumControl[i]
      pvalmat[counter,2]=C
      pvalmat[counter,3]=M
      pvalmat[counter,4]=ME
      pvalmat[counter,5:7]=mas
      
      testmat[counter,1]=NumControl[i]
      testmat[counter,2]=C
      testmat[counter,3]=M
      testmat[counter,4]=ME
      testmat[counter,5:7]=mas
      
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
                        PStest=TRUE, HWtest=TRUE, MStest=TRUE)
      }
    
      # Write out separate datasets with the effect estimates and p-values 
      # across all replicates for the 3 different methods
      write.table(res$resLoglin, paste0(outdir,"LogLin_Estimates_Controls",
                                        NumControl[i],"_C",
                                        C,"_M",M,"_ME",ME,"_Mas",mas[1]),
                  quote=F, row=F)
      write.table(res$resHaplin, paste0(outdir,"Haplin_Estimates_Controls",
                                        NumControl[i],"_C",
                                        C,"_M",M,"_ME",ME,"_Mas",mas[1]),
                  quote=F, row=F)
      write.table(res$resEMIM, paste0(outdir,"EMIM_Estimates_Controls",
                                      NumControl[i],"_C",
                                      C,"_M",M,"_ME",ME,"_Mas",mas[1]),
                  quote=F, row=F)
      write.table(res$pvalLoglin, paste0(outdir,"LogLin_pval_Controls",
                                         NumControl[i],"_C",
                                         C,"_M",M,"_ME",ME,"_Mas",mas[1]),
                  quote=F, row=F)
      write.table(res$pvalHaplin, paste0(outdir,"Haplin_pval_Controls",
                                        NumControl[i],"_C",
                                        C,"_M",M,"_ME",ME,"_Mas",mas[1]),
                  quote=F, row=F)
      write.table(res$pvalEMIM, paste0(outdir,"EMIM_pval_Controls",
                                      NumControl[i],"_C",
                                      C,"_M",M,"_ME",ME,"_Mas",mas[1]),
                  quote=F, row=F)
      
      write.table(round(res$pvalHW,6), paste0(outdir,"HW_pval_Controls",
                                              NumControl[i],"_C",
                                              C,"_M",M,"_ME",ME,"_Mas",mas[1]),
                  quote=F, row=F)
      write.table(round(res$pvalMaS,6), paste0(outdir,"MaS_pval_Controls",
                                         NumControl[i],"_C",
                                         C,"_M",M,"_ME",ME,"_Mas",mas[1]),
            quote=F,row=F)
      
      if (NumControl[i]>0){
        write.table(round(res$pvalPS,6), paste0(outdir,"PS_pval_Controls",
                                                NumControl[i],"_C",
                                                C,"_M",M,"_ME",ME,"_Mas",mas[1]),
                    quote=F, row=F)
      }
      
      # Write out summarized type 1 error/power
      
      # EMIM
      pvalmat[counter,"EMIM.HW.C"]=sum(res$pvalEMIM[,"HW.C"]<alpha)/numsim
      pvalmat[counter,"EMIM.HW.M"]=sum(res$pvalEMIM[,"HW.M"]<alpha)/numsim
      pvalmat[counter,"EMIM.HW.ME"]=sum(res$pvalEMIM[,"HW.E:M"]<alpha)/numsim
      
      pvalmat[counter,"EMIM.MS.C"]=sum(res$pvalEMIM[,"MS.C"]<alpha)/numsim
      pvalmat[counter,"EMIM.MS.M"]=sum(res$pvalEMIM[,"MS.M"]<alpha)/numsim
      pvalmat[counter,"EMIM.MS.ME"]=sum(res$pvalEMIM[,"MS.E:M"]<alpha)/numsim
      
      if (NumControl[i]>0){
        
        pvalmat[counter,"EMIM.MaS.C"]=sum(res$pvalEMIM[,"MaS.C"]<alpha)/numsim
        pvalmat[counter,"EMIM.MaS.M"]=sum(res$pvalEMIM[,"MaS.M"]<alpha)/numsim
        pvalmat[counter,"EMIM.MaS.ME"]=sum(res$pvalEMIM[,"MaS.E:M"]<alpha,na.rm=T)/numsim
        
      }
      
      # Haplin
      pvalmat[counter,"Haplin.C"]=sum(res$pvalHaplin[,"C"]<alpha,na.rm=T)/numsim
      pvalmat[counter,"Haplin.M"]=sum(res$pvalHaplin[,"M"]<alpha,na.rm=T)/numsim
      pvalmat[counter,"Haplin.ME"]=sum(res$pvalHaplin[,"E:M"]<alpha,na.rm=T)/numsim
      
      
      # Loglinear
      pvalmat[counter,"Loglin.HW.C"]=sum(res$pvalLoglin[,"HW.C"]<alpha,na.rm=T)/numsim
      pvalmat[counter,"Loglin.HW.M"]=sum(res$pvalLoglin[,"HW.M"]<alpha,na.rm=T)/numsim
      pvalmat[counter,"Loglin.HW.ME"]=sum(res$pvalLoglin[,"HW.E:M"]<alpha,na.rm=T)/numsim
      
      pvalmat[counter,"Loglin.MS.C"]=sum(res$pvalLoglin[,"MS.C"]<alpha,na.rm=T)/numsim
      pvalmat[counter,"Loglin.MS.M"]=sum(res$pvalLoglin[,"MS.M"]<alpha,na.rm=T)/numsim
      pvalmat[counter,"Loglin.MS.ME"]=sum(res$pvalLoglin[,"MS.E:M"]<alpha,na.rm=T)/numsim
      
      if (NumControl[i]>0){
        
        pvalmat[counter,"Loglin.MaS.C"]=sum(res$pvalLoglin[,"MaS.C"]<alpha)/numsim
        pvalmat[counter,"Loglin.MaS.M"]=sum(res$pvalLoglin[,"MaS.M"]<alpha,na.rm=T)/numsim
        pvalmat[counter,"Loglin.MaS.ME"]=sum(res$pvalLoglin[,"MaS.E:M"]<alpha,na.rm=T)/numsim
        
      }
      
      
      # HW/MS/PS tests
      
      testmat[counter,"MS"]=sum(res$pvalMaS<alpha)/numsim
      
      if (NumControl[i]>0){
        testmat[counter,"PS.HW"]=sum(res$pvalPS[,"HW"]<alpha, na.rm=T)/numsim
        testmat[counter,"PS.MS"]=sum(res$pvalPS[,"MS"]<alpha,na.rm=T)/numsim
        testmat[counter,"PS.MaS"]=sum(res$pvalPS[,"MaS"]<alpha,na.rm=T)/numsim
        testmat[counter,"HW.Mothers.Case"]=sum(res$pvalHW[,"Mothers.Case"]<alpha,na.rm=T)/numsim
        testmat[counter,"HW.Fathers.Case"]=sum(res$pvalHW[,"Fathers.Case"]<alpha, na.rm=T)/numsim
        testmat[counter,"HW.Both.Case"]=sum(res$pvalHW[,"Both.Case"]<alpha, na.rm=T)/numsim
        testmat[counter,"HW.Mothers.Control"]=sum(res$pvalHW[,"Mothers.Control"]<alpha,na.rm=T)/numsim
        testmat[counter,"HW.Fathers.Control"]=sum(res$pvalHW[,"Fathers.Control"]<alpha,na.rm=T)/numsim
        testmat[counter,"HW.Both.Control"]=sum(res$pvalHW[,"Both.Control"]<alpha,na.rm=T)/numsim
      } else {
        testmat[counter,"HW.Mothers.Case"]=sum(res$pvalHW[,"Mothers"]<alpha,na.rm=T)/numsim
        testmat[counter,"HW.Fathers.Case"]=sum(res$pvalHW[,"Fathers"]<alpha,na.rm=T)/numsim
        testmat[counter,"HW.Both.Case"]=sum(res$pvalHW[,"Both"]<alpha,na.rm=T)/numsim
      }
      
      # Write out summarized coefficient estimates
      
      # EMIM
      resmat[counter,"EMIM.HW.C"]=mean(res$resEMIM[,"HW.C"],na.rm=T)
      resmat[counter,"EMIM.HW.M"]=mean(apply(res$resEMIM[,c("HW.M[E=0]","HW.M[E=1]")],1,mean,na.rm=T))
      resmat[counter,"EMIM.HW.ME"]=mean(res$resEMIM[,"HW.M[E=1]/M[E=0]"],na.rm=T)
      
      resmat[counter,"EMIM.MS.C"]=mean(res$resEMIM[,"MS.C"],na.rm=T)
      resmat[counter,"EMIM.MS.M"]=mean(apply(res$resEMIM[,c("MS.M[E=0]","MS.M[E=1]")],1,mean,na.rm=T))
      resmat[counter,"EMIM.MS.ME"]=mean(res$resEMIM[,"MS.M[E=1]/M[E=0]"],na.rm=T)
      
      if (NumControl[i]>0){
        
        resmat[counter,"EMIM.MaS.C"]=mean(res$resEMIM[,"MaS.C"],na.rm=T)
        resmat[counter,"EMIM.MaS.M"]=mean(apply(res$resEMIM[,c("MaS.M[E=0]","MaS.M[E=1]")],1,mean,na.rm=T))
        resmat[counter,"EMIM.MaS.ME"]=mean(res$resEMIM[,"MaS.M[E=1]/M[E=0]"],na.rm=T)
        
      }
      
      # Haplin
      resmat[counter,"Haplin.C"]=mean(res$resHaplin[,"C"],na.rm=T)
      resmat[counter,"Haplin.M"]=mean(apply(res$resHaplin[,c("M[E=0]","M[E=1]")],1,mean,na.rm=T))
      resmat[counter,"Haplin.ME"]=mean(res$resHaplin[,"M[E=1]/M[E=0]"],na.rm=T)
      
      
      # Loglinear
      resmat[counter,"Loglin.HW.C"]=mean(res$resLoglin[,"HW.C"],na.rm=T)
      resmat[counter,"Loglin.HW.M"]=mean(res$resLoglin[,"HW.M"],na.rm=T)
      resmat[counter,"Loglin.HW.ME"]=mean(res$resLoglin[,"HW.E:M"],na.rm=T)
      
      resmat[counter,"Loglin.MS.C"]=mean(res$resLoglin[,"MS.C"],na.rm=T)
      resmat[counter,"Loglin.MS.M"]=mean(res$resLoglin[,"MS.M"],na.rm=T)
      resmat[counter,"Loglin.MS.ME"]=mean(res$resLoglin[,"MS.E:M"],na.rm=T)
      
      if (NumControl[i]>0){
        
        resmat[counter,"Loglin.MaS.C"]=mean(res$resLoglin[,"MaS.C"],na.rm=T)
        resmat[counter,"Loglin.MaS.M"]=mean(res$resLoglin[,"MaS.M"],na.rm=T)
        resmat[counter,"Loglin.MaS.ME"]=mean(res$resLoglin[,"MaS.E:M"],na.rm=T)
        
      }
      
      counter=counter+1
    
    }
  }
  write.table(resmat,paste0(outdir,"Mean_RiskEstimates.dat"),quote=F,row=F)
  write.table(pvalmat,paste0(outdir,"Type1_Power_Risk.dat"),quote=F,row=F)
  write.table(testmat,paste0(outdir,"Type1_Power_HW_PS_MaS.dat"),quote=F,row=F)
}

Sys.time()
