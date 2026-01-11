## This script is used to run all the scenarios with population substructure
## but no MaS

## Kelly Burkett
## May 4, 2023
## Modification Notes
## - June 14 : Added code to write out summarized results on each simulation

options(warn=1)

## Source the analysis functions
source("simulate_data_functions_V3.R")
source("run_analysis_functions_V3.R")
library(Haplin)

## Set directories
prefix="PS_MSHWE/"
setwd(paste0("../Run_directories/",prefix))
outdir=paste0("../../Results/",prefix)

## Setup parameters for simulations
alpha.type1=0.05
numsim=2000 # Number of replicates per set of parameters
NumCase=1000 # Number of case trios
effectmat=list(c(1,1,1))
masmat=list(c(1,1,1))

# Parameters that vary in this simulation
NumControl=c(0,500,1000) # Number of control trios
maf=list(c(0.3,0.3),c(0.25,0.35),c(0.2,0.4)) # Minor allele frequency
ef=maf  # Environmental variable frequency
alpha=list(c(0.1,0.1),c(0.15,0.05))
omega=list(c(0.5,0.5),c(0.4,0.6),c(0.6,0.4))

## Set up output files
namesres=c("Controls","maf.1","maf.2","ef.1","ef.2","alpha.1","alpha.2","omega.1","omega.2",
                   "EMIM.HW.C","EMIM.HW.M","EMIM.HW.ME",
                   "EMIM.MS.C","EMIM.MS.M","EMIM.MS.ME",
                   "EMIM.MaS.C","EMIM.MaS.M","EMIM.MaS.ME",
                   "Haplin.C","Haplin.M","Haplin.ME",
                   "Loglin.HW.C","Loglin.HW.M","Loglin.HW.ME",
                   "Loglin.MS.C","Loglin.MS.M","Loglin.MS.ME",
                   "Loglin.MaS.C","Loglin.MaS.M","Loglin.MaS.ME")
namestests=c("Controls","maf.1","maf.2","ef.1","ef.2","alpha.1","alpha.2","omega.1","omega.2",
                    "HW.Mothers.Case","HW.Fathers.Case","HW.Both.Case",
                    "HW.Mothers.Control","HW.Fathers.Control","HW.Both.Control",
                    "PS.HW","PS.MS","PS.MaS","MS")

resvec=vector(length=length(namesres))
pvalvec=vector(length=length(namesres))
testvec=vector(length=length(namestests))

resvec=vector(length=length(namesres))
names(resvec)=namesres
pvalvec=vector(length=length(namesres))
names(pvalvec)=namesres
testvec=vector(length=length(namestests))
names(testvec)=namestests



           
# Loop to run all scenarios. Note that some combinations of parameters will not
# be run in order to reduce run time and uninteresting results. The boolean variable
# runscenario controls this behaviour

Sys.time()
i=3 # Num Control
j=2 # MAF
k=1 # EF
l=2 # alpha
m=3 # omega

          runscenario=TRUE

          if ( (l==1)&&(m!=1) ) {
            runscenario=FALSE
          }


          if ( ((j==1)&&(l==1)&&(k>1)) || ((j==1)&&(l>1)&&(k==1)) || ((l==1)&&(k==1)&&(j>1)) ) {
            runscenario=FALSE
          }

         if ( (l==1)&&(j==1)&&(k==1)&&(m>1) ) {
            runscenario=FALSE
          }

         if (runscenario==TRUE){

	   cat(paste0("\n\n","Starting Controls",NumControl[i],
                                                     "_maf",paste(maf[[j]],collapse="-"),
                                                     "_ef",paste(ef[[k]],collapse="-"),
                                                     "_alpha",paste(alpha[[l]],collapse="-"),
                                                     "_omega",paste(omega[[m]],collapse="-"),"\n"))

           # Summarized result and pval matrix
           resvec[1]=NumControl[i]
           resvec[2:3]=maf[[j]]
           resvec[4:5]=ef[[k]]
           resvec[6:7]=alpha[[l]]
           resvec[8:9]=omega[[m]]

           pvalvec[1]=NumControl[i]
           pvalvec[2:3]=maf[[j]]
           pvalvec[4:5]=ef[[k]]
           pvalvec[6:7]=alpha[[l]]
           pvalvec[8:9]=omega[[m]]

           testvec[1]=NumControl[i]
           testvec[2:3]=maf[[j]]
           testvec[4:5]=ef[[k]]
           testvec[6:7]=alpha[[l]]
           testvec[8:9]=omega[[m]]

          ratioByPop=getRatio(alpha=alpha[[l]],maf=maf[[j]],omega=omega[[m]])
          
          if (NumControl[i]==0){ # No control trios
            res=runTrioReps(nreps=numsim, ntrios=NumCase,
                            effects=c("M","C","E:M"),
                            includeE=TRUE, prE=ef[[k]], envint="Mother", 
                            includePopStrat=TRUE, maf=maf[[j]],
                            includeControl=FALSE, numPop=2,
                            prCase.byPop=ratioByPop$prPopGivenCase,
                            Loglin=TRUE, loglinMT=c("HW","MS"), 
                            Haplin=TRUE, haplinControls=FALSE,
                            EMIM=TRUE, EMIMMT=c("HW","MS"), 
                            PStest=FALSE, HWtest=TRUE, MStest=TRUE)
          } else { # Control trios
            NumTrios=NumCase+NumControl[i]
            propControl=NumControl[i]/NumTrios
            res=runTrioReps(nreps=numsim, ntrios=NumTrios,
                            effects=c("M","C","E:M"),
                            includeE=TRUE, prE=ef[[k]], envint="Mother", 
                            includePopStrat=TRUE, maf=maf[[j]],
                            includeControl=TRUE, numPop=2, prControl=propControl,
                            prCase.byPop=ratioByPop$prPopGivenCase,
                            prControl.byPop=ratioByPop$prPopGivenControl,
                            Loglin=TRUE, loglinMT=c("HW","MS","MaS"), 
                            Haplin=TRUE, haplinControls=TRUE,
                            EMIM=TRUE, EMIMMT=c("HW","MS","MaS"), 
                            PStest=TRUE, HWtest=TRUE, MStest=TRUE)
          }
    
    
          # Write out separate datasets with the effect estimates and p-values 
          # across all replicates for the 3 different methods
          write.table(round(res$resLoglin,4), paste0(outdir,"LogLin_Estimates_Controls",
                                                     NumControl[i],
                                                     "_maf",paste(maf[[j]],collapse="-"),
                                                     "_ef",paste(ef[[k]],collapse="-"),
                                                     "_alpha",paste(alpha[[l]],collapse="-"),
                                                     "_omega",paste(omega[[m]],collapse="-"),".dat"),
                        quote=F, row=F)
          write.table(round(res$resHaplin,4), paste0(outdir,"Haplin_Estimates_Controls",
                                                     NumControl[i],
                                                     "_maf",paste(maf[[j]],collapse="-"),
                                                     "_ef",paste(ef[[k]],collapse="-"),
                                                     "_alpha",paste(alpha[[l]],collapse="-"),
                                                     "_omega",paste(omega[[m]],collapse="-"),".dat"),
                      quote=F, row=F)
          write.table(round(res$resEMIM,4), paste0(outdir,"EMIM_Estimates_Controls",
                                             NumControl[i],
                                             "_maf",paste(maf[[j]],collapse="-"),
                                             "_ef",paste(ef[[k]],collapse="-"),
                                             "_alpha",paste(alpha[[l]],collapse="-"),
                                             "_omega",paste(omega[[m]],collapse="-"),".dat"),
                quote=F, row=F)
          write.table(round(res$pvalLoglin,6), paste0(outdir,"LogLin_pval_Controls",
                                                NumControl[i],
                                                "_maf",paste(maf[[j]],collapse="-"),
                                                "_ef",paste(ef[[k]],collapse="-"),
                                                "_alpha",paste(alpha[[l]],collapse="-"),
                                                "_omega",paste(omega[[m]],collapse="-"),".dat"),
                quote=F, row=F)
          write.table(round(res$pvalHaplin,6), paste0(outdir,"Haplin_pval_Controls",
                                                NumControl[i],
                                                "_maf",paste(maf[[j]],collapse="-"),
                                                "_ef",paste(ef[[k]],collapse="-"),
                                                "_alpha",paste(alpha[[l]],collapse="-"),
                                                "_omega",paste(omega[[m]],collapse="-"),".dat"),
                quote=F, row=F)
          write.table(round(res$pvalEMIM,6), paste0(outdir,"EMIM_pval_Controls",
                                              NumControl[i],
                                              "_maf",paste(maf[[j]],collapse="-"),
                                              "_ef",paste(ef[[k]],collapse="-"),
                                              "_alpha",paste(alpha[[l]],collapse="-"),
                                              "_omega",paste(omega[[m]],collapse="-"),".dat"),
                quote=F, row=F)
          write.table(round(res$pvalHW,6), paste0(outdir,"HW_pval_Controls",
                                            NumControl[i],
                                            "_maf",paste(maf[[j]],collapse="-"),
                                            "_ef",paste(ef[[k]],collapse="-"),
                                            "_alpha",paste(alpha[[l]],collapse="-"),
                                            "_omega",paste(omega[[m]],collapse="-"),".dat"),
                quote=F, row=F)
          write(round(res$pvalMaS,6), paste0(outdir,"MaS_pval_Controls",
                                       NumControl[i],
                                       "_maf",paste(maf[[j]],collapse="-"),
                                       "_ef",paste(ef[[k]],collapse="-"),
                                       "_alpha",paste(alpha[[l]],collapse="-"),
                                       "_omega",paste(omega[[m]],collapse="-"),".dat"),
          ncol=1)
          
          if (NumControl[i]>0){
            write.table(round(res$pvalPS,6), paste0(outdir,"PS_pval_Controls",
                                              NumControl[i],
                                              "_maf",paste(maf[[j]],collapse="-"),
                                              "_ef",paste(ef[[k]],collapse="-"),
                                              "_alpha",paste(alpha[[l]],collapse="-"),
                                              "_omega",paste(omega[[m]],collapse="-"),".dat"),
                  quote=F, row=F)
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
      
      if (NumControl[i]>0){
        testvec["PS.HW"]=sum(res$pvalPS[,"HW"]<alpha.type1)/numsim
        testvec["PS.MS"]=sum(res$pvalPS[,"MS"]<alpha.type1)/numsim
        testvec["PS.MaS"]=sum(res$pvalPS[,"MaS"]<alpha.type1)/numsim
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





    
Sys.time()
