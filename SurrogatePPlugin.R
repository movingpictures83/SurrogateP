###################################################
### code chunk number 1: sva.Rnw:5-6
###################################################
options(width=65)


###################################################
### code chunk number 3: input
###################################################
library(sva)
library(bladderbatch)
#data(bladderdata)
library(pamr)
library(limma)


dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {
pheno = read.csv(paste(pfix, parameters["pheno", 2], sep="/"))
edata = as.matrix(read.csv(paste(pfix, parameters["edata", 2], sep="/")))
mod = model.matrix(~as.factor(cancer), data=pheno)
mod0 = model.matrix(~1,data=pheno)
n.sv = as.integer(parameters["numsv", 2])#num.sv(edata,mod,method="leek")

	###################################################
### code chunk number 9: input
###################################################
svobj = sva(edata,mod,mod0,n.sv=n.sv)

###################################################
### code chunk number 6: input
###################################################
mod = model.matrix(~as.factor(cancer), data=pheno)


###################################################
### code chunk number 7: input
###################################################
mod0 = model.matrix(~1,data=pheno)


###################################################
### code chunk number 8: input
###################################################
#n.sv = num.sv(edata,mod,method="leek")
#n.sv


###################################################
### code chunk number 11: input
###################################################
modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)

pValuesSv = f.pvalue(edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

write.csv(qValuesSv, outputfile)

}
