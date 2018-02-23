#!/usr/bin/Rscript --vanilla

# Load required libraries
# library(fExtremes,quietly=TRUE,verbose=FALSE)
# Workaround to load library without messages.
msg.trap <- capture.output( suppressMessages( library( fExtremes ) ) )
library(boot,quietly=TRUE)

#source("/u1/qudo046/ADG/adglib/uq/bootstrap/mcvar.R")
#source("~/EXPORT/BIN/mcvar.R")
#source("../mcvar.R")
MCVaR<-
function (x, alpha = 0.05, type = "sample", tail = c("lower", 
    "upper")) 
{
    x = as.matrix(x)
    tail = match.arg(tail)
    VaR = VaR(x, alpha, type, tail)
#   if (tail == "upper") 
#       alpha = 1 - alpha
    if (alpha==1)
        return(VaR)
    if (type == "sample") {
        MCVaR = NULL
        for (i in 1:ncol(x)) {
            X = as.vector(x[, i])
            MCVaR = c(MCVaR, VaR[i] + 0.5 * mean(((-VaR[i] + X) + 
                abs(-VaR[i] + X)))/(1-alpha))
        }
    }
    MCVaR
}


## Collect arguments
args <- commandArgs(TRUE)
#print(args)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      The R Script
 
      Arguments:
      --Nboot=someValue   - numeric, Bootsrap iterations (if equal to 0, no bootstrap analysis is performed)
      --DB=someValue      - character, Data table
      --Vname=someValue   - character, Name of the variable to be extracted and elaborated from Data table
      --alpha=[0,1]       - numeric, Percentile value
      --subset=<n>        - Resample without replacement the original data frame randomly picking n elements
      --plot=[TRUE,FALSE] - if TRUE pdf plots are produced in output (default: TRUE)
      --help              - print this text
 
      Example:
      ./compute_boot --Nboot=10000 --Vname=response_fn_1 --alpha=0.99 --DB=uq_tabular.dat --subset=1000\n\n")
 
  q(save="no")
}
 
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

#print(argsDF)
#print(argsL)
#print( parseArgs(args))

## plot
if(is.null(argsL$plot)) {
  plot=TRUE   # Default
} else {
  plot=as.logical(argsL$plot)
}
cat("plot = ",plot,"\n")


## Nboot default
if(is.null(argsL$Nboot)) {
  Nboot=100
} else {
  Nboot=as.numeric(argsL$Nboot)
}
cat("Nboot = ",Nboot,"\n")
 
## DB default
if(is.null(argsL$DB)) {
  DB="DATABASE.des"
} else {
  DB=argsL$DB
}
cat("DB    = ",DB,"\n")
 
## Vname default
if(is.null(argsL$Vname)) {
  Vname="OB_00000"
} else {
  Vname=argsL$Vname
}
cat("Vname = ",Vname,"\n")

## alpha default
if(is.null(argsL$alpha)) {
  alpha=0.95
} else {
  alpha=as.numeric(argsL$alpha)
}
cat("alpha = ",alpha,"\n")

DATABASE <- read.table(DB, header=T, quote="\"")
#attach(DATABASE)

if(is.null(argsL$subset)) {
  f=DATABASE[,Vname]
} else {
  nsubset=as.numeric(argsL$subset)
  nr<-nrow(DATABASE)
  f<-DATABASE[sample(1:nr,nsubset),Vname]
}
cat("DB size = ",length(f),"\n")
# To access a column, use:
# my_data[ , cond]
# or
# my_data[[cond]]
# The ith row can be accessed with:
# my_data[i, ]
# Combine both to obtain the desired value:
# my_data[i, cond]
# or
# my_data[[cond]][i]


#print(f)
cat(" VaR = ",VaR(f,alpha),"\n")
cat("CVaR = ",MCVaR(f,alpha),"\n")
if(Nboot>0) {
  varboot<- function(f,inds,alpha){VaR(f[inds],alpha)}
  b1<-boot(f,varboot,Nboot,alpha=alpha)
  cat("\n\n VaR_Bootstrap:\n")
  print(b1)
  #boot.ci(b1, type = "all")
  b1ci<-boot.ci(b1, type =  "perc")
  print(b1ci)
  if(is.null(b1ci)) {
    cat(" Cannot calculate confidence intervals\n")
    cat("VaR.ci ----- ",VaR(f,alpha),"  ",VaR(f,alpha),"\n")
  } else {
    cat("VaR.ci ",b1ci$percent[1],"  ",b1ci$percent[4],"  ",b1ci$percent[5],"\n")
  }
  cvarboot<- function(f,inds,alpha){MCVaR(f[inds],alpha)}
  b2<-boot(f,cvarboot,Nboot,alpha=alpha)
  cat("\n\nCVaR_Bootstrap:\n")
  print(b2)
  #boot.ci(b2, type = "all")
  #boot.ci(b2, type =  c("norm","basic", "perc"))
  b2ci<-boot.ci(b2, type =  "perc")
  print(b2ci)
  if(is.null(b1ci)) {
    cat(" Cannot calculate confidence intervals\n")
    cat("CVaR.ci ----- ",CVaR(f,alpha),"  ",CVaR(f,alpha),"\n")
  } else {
    cat("CVaR.ci ",b2ci$percent[1],"  ",b2ci$percent[4],"  ",b2ci$percent[5],"\n\n")
  }
}
if(plot) {
  plot(ecdf(f))
  if(Nboot>0) {
    plot(b1)
    plot(b2)
  }
}
q(save="no")
