# R script
# Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>
# Mon Apr 9 2018 

# based on  https://informatics.nescent.org/wiki/R_Hackathon_1/Ancestral_State_Reconstruction#Reconstructing_Ancestral_States_for_Discrete_Variables
# script runs two likelihood ratio tests, one with the ctenophore sister tree 
#   and one with the sponge sister tree.

library(ape)

mydir         <- "[input directory name]"

cucumbertree  <- "cucumber.nex"
charmat       <- "charmat_nomissindata.dat"

# if there is missing data we have to remove it because missing data causes a problem with ace
unknowns <- c(insert "names" of missing data)

setwd(mydir)

#lrt: likelihood ratio test
mylrt <- function(mytree) {
    anitree <- read.nexus(mytree)
    anidata <- read.table(charmat)

    for (unk in unknowns) {
        anitree <- tip(anitree, unk)
    }

    ERreconstruction <- ace(anidata$V2, anitree, type="discrete", model="ER")
    SYMreconstruction <- ace(anidata$V2, anitree, type="discrete", model="SYM")
    ARDreconstruction <- ace(anidata$V2, anitree, type="discrete", model="ARD")

    cat("RESULTS FOR ",mytree,"\n")

    erlnl <- ERreconstruction$loglik
    cat("ER log likelihood: ",erlnl,"\n")
    
    symlnl <- SYMreconstruction$loglik
    cat("SYM log likelihood: ",symlnl,"\n")
    
    ardlnl <- ARDreconstruction$loglik
    cat("ARD log likelihood: ",ardlnl,"\n")
    
    #  For a three-state character, 
    #     ER is a one parameter model, 
    #     SYM a three parameter model,
    #     ARD a six parameter model.

    # df = 5 (i.e. ARD(6) - ER(1))
    erard_stat <- 2*abs(ERreconstruction$loglik - ARDreconstruction$loglik)
    erard <- 1-pchisq(2*abs(ERreconstruction$loglik - ARDreconstruction$loglik), 5)
    cat("ER vs. ARD: ", erard, ", stat = ", erard_stat, "\n")
    cat("    if <= 0.05 than ARD is significantly better than ER\n")

    # df = 3 (i.e. ARD(6) - SYM(3))
    symard <- 1-pchisq(2*abs(SYMreconstruction$loglik - ARDreconstruction$loglik), 3)
    symard_stat <- 2*abs(SYMreconstruction$loglik - ARDreconstruction$loglik)
    cat("SYM vs. ARD: ",symard, ", stat = ", symard_stat, "\n")
    cat("    if <= 0.05 than ARD is significantly better than SYM\n")

    # df = 2 (i.e. SYM(3) - ER(1))
    ersym <- 1-pchisq(2*abs(ERreconstruction$loglik - SYMreconstruction$loglik), 2)
    ersym_stat <- 2*abs(ERreconstruction$loglik - SYMreconstruction$loglik)
    cat("ER vs. SYM: ",ersym, ", stat = ", ersym_stat, "\n")
    cat("    if <= 0.05 than SYM is significantly better than ER\n")
    cat("\n---------------------------------------------------------------\n\n")
}    
    
mylrt(cucumbertree)
sessionInfo()
