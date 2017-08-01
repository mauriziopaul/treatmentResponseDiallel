##---------------------------------------------------------------------------------------------------------------------
## Title: Functions for the treatmentResponseDiallel
## Author: Paul L. Maurizio
## Institution: University of North Carolina at Chapel Hill
## Program: Bioinformatics and Computational Biology Ph.D. Curriculum
## Advisors: Will Valdar, Mark Heise
## Date Created: 2015-06-18
## Date Updated: 2017-03-22
##---------------------------------------------------------------------------------------------------------------------

#' treatmentResponseDiallel: A package for analysis of infection/treatment response in a diallel cross of inbred lines.
#' The treatment response diallel (tRD) package provides the important function:
#' run.tr.diallel
#' 
#' @docType package

#' @name treatmentResponseDiallel
NULL

#' @import BayesDiallel
NULL

#' @importFrom MESS auc
NULL

#' @import cmdline
NULL

#' @importFrom coda nchain as.mcmc varnames
NULL

#' @import configfile
NULL

#' @import lattice
NULL

#' @import methods
NULL

#' @import tools
NULL

#' @importFrom grDevices dev.off pdf png rgb
NULL

#' @importFrom graphics abline axis barplot frame lines par plot points segments text title
NULL

#' @importFrom stats complete.cases density rnorm sd
NULL

#' @importFrom utils data head read.csv str write.csv write.csv2
NULL

#' @section require namespaces:

#requireNamespace("BayesDiallel", quietly=TRUE)
requireNamespace("cmdline", quietly=TRUE)
requireNamespace("coda", quietly=TRUE)
requireNamespace("configfile", quietly=TRUE)
requireNamespace("lattice", quietly=TRUE)
requireNamespace("MESS", quietly=TRUE) ## required for AUC
requireNamespace("methods", quietly=TRUE)
requireNamespace("tools", quietly=TRUE)
#requireNamespace("MCMCglmm", quietly=TRUE)
#requireNamespace("WVmisc", quietly=TRUE)
#requireNamespace("BayesSpike")

#' @section treatmentResponseDiallel functions:

#' @title removeDots: remove dots, underscores from a string
#' @description This removes dots and underscores from a string, replacing them with a space.
#' 
#' @param string the dot-containing string or a vector of character strings
#' @param ... additional arguments
#' @return Returns string without dot/underscore.
#' @examples
#' removeDots("This_is_a.test.string.")
#' @export
removeDots <- function(string, ...){
  if(!(is.character(string))){
  	cat("ERROR:", as.character(string), "is not a string! \n")
		return(cat(NULL))
	}
	new.string <- gsub(".", " ", string, fixed=TRUE)
	new.string <- gsub("_", " ", new.string, fixed=TRUE)
	return(new.string)
}

#' @title col.switch: Convert CC names to colors 
#' @description Convert CC founder abbreviated strain names to standard CC colors
#' 
#' @param arg the argument, in the form of "AA", "BB", etc.
#' @param ... additional arguments
#' @return Defines standard CC founder colors in hex-code
#' @examples
#' col.switch("A", "B", "C", "DD", "EE", "FF")
#' @export
col.switch <- function(arg, ...){
 	founder.colors <- c("#F0F000", "#808080", "#F08080", "#1010F0", 
 		"#00A0F0", "#00A000", "#F00000", "#9000E0")
 	ret <- NULL
 	for(i in 1:length(arg)){
 		ret[i] <- switch(as.character(arg[i]),
 				AA=founder.colors[1],
 				BB=founder.colors[2],
 				CC=founder.colors[3],
 				DD=founder.colors[4],
 				EE=founder.colors[5],
 				FF=founder.colors[6],
 				GG=founder.colors[7],
 				HH=founder.colors[8],
 				A=founder.colors[1],
 				B=founder.colors[2],
 				C=founder.colors[3],
 				D=founder.colors[4],
 				'E'=founder.colors[5],
 				F=founder.colors[6],
 				G=founder.colors[7],
 				H=founder.colors[8],
 				AJ=founder.colors[1],
 				B6=founder.colors[2],
 				'129'=founder.colors[3],
 				NOD=founder.colors[4],
 				NZO=founder.colors[5],
 				CAST=founder.colors[6],
 				PWK=founder.colors[7],
 				WSB=founder.colors[8],
 				'white'="white",
 				'black'="black",
 				"black")
 	}
 	return(ret)
}

#' @title lettersToNumbers
#' @description Convert CC founder strain name letters to numbers
#' 
#' @param arg the argument, in the form of "A", "B", ... "H", or vector of these.
#' @param ... additional arguments
#' @return Defines standard CC founder numbers from single letters.
#' @examples
#' lettersToNumbers("A", "B", "C", "D", "E", "F")
#' @export
lettersToNumbers <- function(arg, ...){
  ret <- NULL
  for(i in 1:length(arg)){
    ret[i] <- switch(as.character(arg[i]),
                     A=1, B=2, C=3, D=4,
                     E=5, F=6, G=7, H=8)
  }
  return(ret)
}

#' @title removeDots
#' @description Remove dots and underscores from strings
#' 
#' @param x A string or vector of strings 
#' @param ... additional arguments
#' @return Returns string (or vector of strings) without dots or underscores
#' @examples
#' removeDots("This.is.a.test_string")
#' @export
removeDots <- function(x){
	if(!(is.character(x))){
		cat("ERROR:", as.character(x), "is not a string! \n")
		return(cat(NULL))
	}
	new.x <- gsub(".", " ", x, fixed=TRUE)
	new.x <- gsub("_", " ", new.x, fixed=TRUE)
	return(new.x)
}

mcmc.stack <- function (coda.object, ...){
	## This function is from Will; also part of BayesDiallel
    if (inherits(coda.object, "mcmc")) {
        return(coda.object)
    }
    if (!inherits(coda.object, "mcmc.list")) {
        stop("Non-mcmc object passed to function\n")
    }
    chain <- coda.object[[1]]
    for (i in 2:nchain(coda.object)) {
        chain <- rbind(chain, coda.object[[i]])
    }
    as.mcmc(chain)
}

#' @title mcmc.stack.and.burn: Stack MCMC chains, after burning each
#' @description This is a modification of mcmc.stack that also removes specified number of burn-ins from chains.
#' 
#' @param coda.object this is the mcmc object
#' @param burnin this is the common amount to burn off each chain
#' @param ... additional arguments
#' @return returns single chain, stacked and burned MCMC object
#' @examples
#' NULL
#' @export
mcmc.stack.and.burn <- function (coda.object, burnin, ...){
	burnin <- burnin+1
    if (inherits(coda.object, "mcmc")) {
        return(coda.object)
    }
    if (!inherits(coda.object, "mcmc.list")) {
        stop("Non-mcmc object passed to function\n")
    }
    len <- dim(coda.object[[1]])[1]
    chain <- coda.object[[1]][c(burnin:len),]
    for (i in 2:nchain(coda.object)) {
        chain <- rbind(chain, coda.object[[i]][c(burnin:len),])
    }
    as.mcmc(chain)
}

#' @title var.translate: translate variable names
#' @description Translate variable names to simpler/more readable standards.
#' 
#' @param variable.name this is the name or names of variables to be translated
#' @param ... additional arguments
#' @return returns translation result
#' @examples
#' strings <- c("tau:Gender:aj", "SymCrossjk", ":j:4")
#' var.translate(strings)
#' @export
var.translate <- function(variable.name, ...){
	var1 <- variable.name
	var1 <- sub(" - mean\\(.+?\\)", "", var1)
	var1 <- sub("^Sigma:1", "Sigma", var1)
	#var1 <- sub("^tau:RandomEffect:1", "tau:RandomEffect", var1)
	var1 <- sub(":j:1", ":AJ", var1)
	var1 <- sub("j:1", "j:AJ", var1)
	var1 <- sub(":j:2", ":B6", var1)
	var1 <- sub(":j:3", ":129", var1)
	var1 <- sub(":j:4", ":NOD", var1)
	var1 <- sub(":j:5", ":NZO", var1)
	var1 <- sub(":j:6", ":CAST", var1)
	var1 <- sub(":j:7", ":PWK", var1)
	var1 <- sub(":j:8", ":WSB", var1)
	var1 <- sub(":j8", ":WSB", var1)
	var1 <- sub("j:2", "j:B6", var1)
	var1 <- sub("j:3", "j:129", var1)
	var1 <- sub("j:4", "j:NOD", var1)
	var1 <- sub("j:5", "j:NZO", var1)
	var1 <- sub("j:6", "j:CAST", var1)
	var1 <- sub("j:7", "j:PWK", var1)
	var1 <- sub("j:8", "j:WSB", var1)
	var1 <- sub("j8", "j:WSB", var1)
	var1 <- sub(";k:1", ";AJ", var1)
	var1 <- sub(";k:2", ";B6", var1)
	var1 <- sub(";k:3", ";129", var1)
	var1 <- sub(";k:4", ";NOD", var1)
	var1 <- sub(";k:5", ";NZO", var1)
	var1 <- sub(";k:6", ";CAST", var1)
	var1 <- sub(";k:7", ";PWK", var1)
	var1 <- sub(";k:8", ";WSB", var1)
	var1 <- sub("tau:Gender:aj", "tau:female:additive", var1)
	var1 <- sub("tau:Gender:motherj", "tau:female:maternal", var1)
	var1 <- sub("tau:Gender:dominancej", "tau:female:inbreeding", var1)
	var1 <- sub("tau:Gender:ASymCrossjkDkj", "tau:female:w", var1)
	var1 <- sub("tau:Gender:SymCrossjk", "tau:female:v", var1)
	var1 <- sub("Gender:aj", "f:additive", var1)
	var1 <- sub("Gender:motherj", "f:maternal", var1)
	var1 <- sub("Gender:dominancej", "f:inbreeding", var1)
	var1 <- sub("Gender:BetaHybrid:Av", "f:inbreed.overall", var1)
	var1 <- sub("Gender:ASymCrossjkDkj", "f:w", var1)
	var1 <- sub("Gender:SymCrossjk", "f:v", var1)
	var1 <- sub("aj", "additive", var1)
	var1 <- sub("motherj", "maternal", var1)
	var1 <- sub("dominancej", "inbreeding", var1)
	var1 <- sub("BetaInbred:Av", "inbreed.overall", var1)
	var1 <- sub("ASymCrossjkDkj", "w", var1)
	var1 <- sub("SymCrossjk", "v", var1)
	var1 <- sub("BetaInbred:Gender:Av", "female.inbred", var1)
	var1 <- sub("Gender:Av", "female.overall", var1)
	var1 <- sub("RandomEffect", "Batch", var1)
	return(var1)
}

#' @title var.translate.glmm: translate variable names from BayesDiallel MCMCglmm model
#' @description Translate BayesDiallel MCMCglmm variable names to standard formats.
#' 
#' @param variable.vec a vector of the name(s) of variables to be translated
#' @param ... additional arguments
#' @return returns translation result
#' @examples
#' strings <- c("(Intercept)", "is.female", "units")
#' var.translate.glmm(strings)
#' @export

var.translate.glmm <- function(variable.vec, ...){
	var1 <- variable.vec
	var1 <- sub("^t.", "", var1)
	var1 <- sub("^units", "Sigma", var1)
	var1 <- sub("^(Intercept)", "mu", var1)
	var1 <- sub(".NA.1$", "", var1)
	#var1 <- sub("^tau:RandomEffect:1", "tau:RandomEffect", var1)
	var1 <- sub("^add.mat", "additive:", var1)
	var1 <- sub("^mat.mat", "maternal:", var1)
	var1 <- sub("^inbred.mat", "inbred:", var1)
	var1 <- sub("^jk.mat", "v:", var1)
	var1 <- sub("^asymm.mat", "w:", var1)
	var1 <- sub("^f.add.mat", "f:additive:", var1)
	var1 <- sub("^f.mat.mat", "f:maternal:", var1)
	var1 <- sub("^f.inbred.mat", "f:inbred:", var1)
	var1 <- sub("^f.jk.mat", "f:v:", var1)
	var1 <- sub("^f.asymm.mat", "f:w:", var1)
	var1 <- sub(":1", ":AJ", var1)
	var1 <- sub(":2", ":B6", var1)
	var1 <- sub(":3", ":129", var1)
	var1 <- sub(":4", ":NOD", var1)
	var1 <- sub(":5", ":NZO", var1)
	var1 <- sub(":6", ":CAST", var1)
	var1 <- sub(":7", ":PWK", var1)
	var1 <- sub(":8", ":WSB", var1)
	# var1 <- sub("A[A-H]", "j:AJ", var1)
	# var1 <- sub("B[A-H]", "j:B6", var1)
	# var1 <- sub("C[A-H]", "j:129", var1)
	# var1 <- sub("D[A-H]", "j:NOD", var1)
	# var1 <- sub("E[A-H]", "j:NZO", var1)
	# var1 <- sub("F[A-H]", "j:CAST", var1)
	# var1 <- sub("G[A-H]", "j:PWK", var1)
	# var1 <- sub("H[A-H]", "j:WSB", var1)
	# var1 <- sub("[A-H]A", ";AJ", var1)
	# var1 <- sub("[A-H]B", ";B6", var1)
	# var1 <- sub("[A-H]C", ";129", var1)
	# var1 <- sub("[A-H]D", ";NOD", var1)
	# var1 <- sub("[A-H]E", ";NZO", var1)
	# var1 <- sub("[A-H]F", ";CAST", var1)
	# var1 <- sub("[A-H]G", ";PWK", var1)
	# var1 <- sub("[A-H]H", ";WSB", var1)
	var1 <- sub("^f:additive:.$", "tau:female:additive", var1)
	var1 <- sub("^f:maternal:.$", "tau:female:maternal", var1)
	var1 <- sub("^f:inbred:.$", "tau:female:inbreeding", var1)
	var1 <- sub("^f:v:.$", "tau:female:v", var1)
	var1 <- sub("^f:w:.$", "tau:female:w", var1)
	var1 <- sub("^additive:.$", "tau:additive", var1)
	var1 <- sub("^maternal:.$", "tau:maternal", var1)
	var1 <- sub("^inbred:.$", "tau:inbreeding", var1)
	var1 <- sub("^v:.$", "tau:v", var1)
	var1 <- sub("^w:.$", "tau:w", var1)
	var1 <- sub("Gender:aj", "f:additive", var1) # fix
	var1 <- sub("Gender:motherj", "f:maternal", var1) # fix 
	var1 <- sub("Gender:dominancej", "f:inbreeding", var1) # fix
	var1 <- sub("Gender:BetaHybrid:Av", "f:inbreed.overall", var1) # fix
	var1 <- sub("Gender:ASymCrossjkDkj", "f:w", var1) # fix
	var1 <- sub("Gender:SymCrossjk", "f:v", var1) # fix
	var1 <- sub("batch.mat.$", "tau:Batch", var1)
	var1 <- sub("aj", "additive", var1)
	var1 <- sub("motherj", "maternal", var1)
	var1 <- sub("dominancej", "inbreeding", var1)
	var1 <- sub("BetaInbred:Av", "inbreed.overall", var1)
	var1 <- sub("ASymCrossjkDkj", "w", var1)
	var1 <- sub("SymCrossjk", "v", var1)
	var1 <- sub("inbred:is.female", "female.inbred", var1)
	var1 <- sub("is.female", "female.overall", var1)
	return(var1)
}

#' @title dot.plotter: Make dot plots indicating treated and control per cross
#' @description Create plots, with dots, indicating treated and control groups in each jk cross.
#' 
#' @param data.frame the data frame containing elements
#' @param title the title of the plot
#' @param batch.num the number of the experimental batch
#' @param trt.string a string indicating the treatment group
#' @param ctrl.string a string indicating the control group
#' @param ... additional arguments
#' @return make pdf plots of treated and control, by batch
#' @examples
#' ## not run
#' @export
dot.plotter <- function(data.frame, title, trt.string, ctrl.string, batch.num=TRUE, ...){
	## Note: modify this to define a plot directory
	bottom.margin 	<- 0.5
	name.margin 	<- 1.1
	title.margin 	<- 3.1
	right.margin 	<- 1.0
	pdf(paste0(title,".pdf"), width=12, height=12)
	mar <- c(bottom.margin, name.margin, title.margin, right.margin)
	par(mar=mar)
	par(mfrow=c(8,8))

	for(j in c(1:8)){
		for(k in c(1:8)){
			this.cross <- paste0(LETTERS[j], LETTERS[k])
			plot(1, type="n", main=NULL, xlab="", ylab="",
				ylim=c(-5,5), xlim=c(-0.5,12), xaxt='n', yaxt='n', cex=0.75)
			title(this.cross, adj=1)
			if(1==j && 1==k){title(title, adj=0)}
			trt.dat <- subset(data.frame, Strain==this.cross & Trt==trt.string)
			mock.dat <- subset(data.frame, Strain==this.cross & Trt==ctrl.string)
			if(nrow(trt.dat)>0){
				try(points(c(rep(-2,nrow(trt.dat)))~ c(1:nrow(trt.dat)),
					col="red", pch=16, cex=1.5))
				if(TRUE==batch.num){
					try(text(c(rep(-3.75,nrow(trt.dat)))~ c(1:nrow(trt.dat)), 
						pch=16, cex=0.6, labels=sort(trt.dat$Block)))
				}
			}
			if(nrow(mock.dat)>0){
				try(points(c(rep(2,nrow(mock.dat)))~ c(1:nrow(mock.dat)),
					col="black", pch=16, cex=1.5))
				if(TRUE==batch.num){
					try(text(c(rep(3.75,nrow(mock.dat)))~ c(1:nrow(mock.dat)),
						pch=16, cex=0.6, labels=sort(mock.dat$Block)))
				}
			}
		}
	}
	dev.off()
}

#' @title batched.plotter: Make batched dot plots
#' @description Create batched dot plots indicating treated and control for each jk cross.
#' 
#' @param data.frame the data frame containing elements
#' @param title the title of the plot
#' @param sexed define whether or not to indicate sex-specificity
#' @param trt.string a string indicating the treatment group
#' @param ctrl.string a string indicating the control group
#' @param ... additional arguments
#' @return make pdf plots of treated and control, by batch, and by sex
#' @examples
#' ## not run
#' @export
batched.plotter <- function(data.frame, title, trt.string, ctrl.string, sexed=c(FALSE, TRUE)[1], ...){
	bottom.margin 	<- 0.5
	name.margin 	<- 1.1
	title.margin 	<- 3.1
	right.margin 	<- 1.0
	pdf(paste0(title,".pdf"), width=12, height=12)
	mar <- c(bottom.margin, name.margin, title.margin, right.margin)
	par(mar=mar)

	batch.list <- sort(unique(data$Block)) ## 52 batches numbered between 1 and 92
	batch.length <- length(batch.list)
	for(b in batch.list){
		par(mfrow=c(8,8))
		for(j in c(1:8)){
			for(k in c(1:8)){
				this.cross <- paste0(LETTERS[j], LETTERS[k])
				plot(1, type="n", main=NULL, xlab="", ylab="",
					ylim=c(-5,5), xlim=c(-0.5,16), xaxt='n', yaxt='n')
				title(this.cross, adj=1)
				if(1==j && 1==k){title(paste0("Batch:", b), adj=0)}
				trt.dat <- subset(data.frame, Strain==this.cross & Trt==trt.string & Block==b)
				mock.dat <- subset(data.frame, Strain==this.cross & Trt==ctrl.string & Block==b)
				if(nrow(trt.dat)>0){
					pch <- 16
					try(pch <- if(TRUE==sexed){pch=c(1,16)[c(as.integer(trt.dat$Sex))]})
					if(1==nrow(trt.dat)){
							x <- 8
						}else{
							x <- 16*(c(1:nrow(trt.dat)-1)/(nrow(trt.dat)-1))
						}
					try(points(c(rep(-2,nrow(trt.dat)))~ x,
						col="red", pch=pch, cex=1))
					try(text(c(rep(-3.75,nrow(trt.dat)))~ x, 
						cex=0.5, labels=sort(trt.dat$Day)))
				}
				if(nrow(mock.dat)>0){
					pch <- 16
					try(pch <- if(TRUE==sexed){pch=c(1,16)[c(as.integer(mock.dat$Sex))]})
					if(1==nrow(mock.dat)){
							x <- 8
						}else{
							x <- 16*(c(1:nrow(mock.dat)-1)/(nrow(mock.dat)-1))
						}
					try(points(c(rep(2,nrow(mock.dat)))~ x,
						col="black", pch=pch, cex=1))
					try(text(c(rep(3.75,nrow(mock.dat)))~ x, 
						cex=0.5, labels=sort(mock.dat$Day)))				
				}

			}
		}
	}
	dev.off()
}

#' @title batch.plotter.wrapper: Make dot plots by batch
#' @description Generate dot plots, separated by batch.
#' 
#' @param data the data frame being used, in this case, from FluDiData
#' @param trt.string a string indicating the treatment group
#' @param ctrl.string a string indicating the control group
#' @param ... additional arguments
#' @return a wrapper for making pdf plots of treated and control, by batch
#' @examples
#' ## not run
#' @export
batch.plotter.wrapper <- function(data, trt.string, ctrl.string, ...){
	data.D2 <- subset(data, Day=="D2")
	data.D4 <- subset(data, Day=="D4")
	data.m <- subset(data, Sex=="M")
	data.f <- subset(data, Sex=="F")
	data.D2.m <- subset(data.D2, Sex=="M")
	data.D2.f <- subset(data.D2, Sex=="F")
	data.D4.m <- subset(data.D4, Sex=="M")
	data.D4.f <- subset(data.D4, Sex=="F")

	data.frame <- data.D2.m
	title <- "Batches_D2m"
	dot.plotter(data.frame=data.frame, title=title, trt.string=trt.string, ctrl.string=ctrl.string)

	data.frame <- data.D2.f
	title <- "Batches_D2f"
	dot.plotter(data.frame=data.frame, title=title, trt.string=trt.string, ctrl.string=ctrl.string)

	data.frame <- data.D4.m
	title <- "Batches_D4m"
	dot.plotter(data.frame=data.frame, title=title, trt.string=trt.string, ctrl.string=ctrl.string)

	data.frame <- data.D4.f
	title <- "Batches_D4f"
	dot.plotter(data.frame=data.frame, title=title, trt.string=trt.string, ctrl.string=ctrl.string)

	data.frame <- data.m
	title <- "Batches_m"
	batched.plotter(data.frame=data.frame, sexed=FALSE, title=title, trt.string=trt.string, ctrl.string=ctrl.string)

	data.frame <- data.f
	title <- "Batches_f"
	batched.plotter(data.frame=data.frame, sexed=FALSE, title=title, trt.string=trt.string, ctrl.string=ctrl.string)

	data.frame <- data
	title <- "Batches"
	batched.plotter(data.frame=data.frame, sexed=TRUE, title=title, trt.string=trt.string, ctrl.string=ctrl.string)
}

#' @title batch.plotter.wrapper.new: Make dot plots by batch
#' @description Make dot plots, separated by batch (variant).
#' 
#' @param data the data frame being used, in this case, from FluDiData
#' @param phenotype the phenotype for which you are making plots
#' @param trt.string a string indicating the treatment group
#' @param ctrl.string a string indicating the control group
#' @param ... additional arguments
#' @return a wrapper for making pdf plots of treated and control, by batch
#' @examples
#' ## not run
#' @export
batch.plotter.wrapper.new <- function(data, phenotype, trt.string, ctrl.string, ...){
	data <- data[!is.na(data[,phenotype]),]
	data.m <- droplevels(subset(data, Sex=="M"))
	data.f <- droplevels(subset(data, Sex=="F"))

	title <- paste(phenotype, "m", sep="_")
	dot.plotter(data.frame=data.m, title=title, trt.string=trt.string, ctrl.string=ctrl.string)

	title <- paste(phenotype, "f", sep="_")
	dot.plotter(data.frame=data.f, title=title, trt.string=trt.string, ctrl.string=ctrl.string)

	title <- paste("Batches", phenotype, "m", sep="_")
	batched.plotter(data.frame=data.m, sexed=FALSE, title=title, trt.string=trt.string, ctrl.string=ctrl.string)

	title <- paste("Batches", phenotype, "f", sep="_")
	batched.plotter(data.frame=data.f, sexed=FALSE, title=title, trt.string=trt.string, ctrl.string=ctrl.string)

	title <- paste("Batches", phenotype, sep="_")
	batched.plotter(data.frame=data, sexed=TRUE, title=title, trt.string=trt.string, ctrl.string=ctrl.string)
}


##-----------------------------------------------------------------------------------------------
## TITLE: MAKE MATCHES
##-----------------------------------------------------------------------------------------------

#' @title make.matches: make matches
#' @description Makes one or more matches between treated and control mice in the flu diallel 
#' based on specified classes and match strategy. In the manuscript, matches do not differ
#' between imputations (same quartets always selected), but imputed phenotype values differ.
#' 
#' @param data data set
#' @param reps number of repetitions or matches
#' @param trt.colname treatment column name
#' @param trt.string treated
#' @param ctrl.string control
#' @param fdir file directory
#' @param force logical indicating whether to force matching despite warnings/errors
#' @param strategy this gives the treatment-to-control ratio (one-to-one, two-to-one), or mean match
#' @param select.phenotype logical indicating whether to select a specific phenotype
#' @param ignore.batch logical indicating whether to match independently of batch (TRUE), or within batch (FALSE)
#' @param ignore.day logical indicating whether to match independently of day assignment (TRUE), or within day assignment (FALSE)
#' @param batch character string giving batch column name
#' @param matchname name of match file
#' @param ... additional arguments
#' @return returns matched treatment response files
#' @examples
#' ## This example may take a couple of minutes
#' data(FluDiData)
#' write.csv(FluDiData, file="FluDiData.csv")
#' filename <- "FluDiData.csv"
#' data <- read.csv("FluDiData.csv")
#' reps <- 1
#' trt.string <- "FLU"
#' ctrl.string <- "MOCK"
#' fdir <- "." 
#' strategy <- "mean"
#' make.matches(data=data, reps=reps, trt.string=trt.string, 
#'              ctrl.string=ctrl.string, fdir=fdir, strategy=strategy)
#' @export
make.matches <- function(data, reps, trt.colname="Trt", trt.string, ctrl.string, 
	fdir, force=TRUE,
	strategy=c("one-to-one", "two-to-one", "three-to-one", "mean")[1], 
	select.phenotype = FALSE, ignore.batch=FALSE, ignore.day=TRUE, batch="Block",
	matchname="",
	...){
	if("mean"==strategy){
		reps <- 1
		print("With mean matching, there is only one imputation. (rep=1)")
	}

	raw.names <- c("D0", "D1", "D2", "D3", "D4", "D5")
	alt.names <- c("d0", "d1", "d2", "d3", "d4", "d5")
	pct.names <- c("pct_D1", "pct_D2", "pct_D3", "pct_D4", "pct_D5")
	auc.names <- c("auc0", "auc1", "auc2")

	## Generate categories of animals
	phenotypes <- c(raw.names, alt.names, pct.names, auc.names)
	phenotypes <- phenotypes[phenotypes %in% names(data)]
	if(strategy %in% c("two-to-one", "three-to-one")){
		if(FALSE==select.phenotype){
			print("For the two/three-to-one strategy, you must select a phenotype\n")
			break()
		}
		data <- data[complete.cases(data[,as.character(select.phenotype)]),]
	}
	if(FALSE==ignore.batch){
		categ <- apply(data, MARGIN=1, FUN=function(x){
			paste(x[["Strain"]], x[["Sex"]], as.character(as.integer(x[[batch]])), sep="_")})
	}else{
		categ <- apply(data, MARGIN=1, FUN=function(x){
			paste(x[["Strain"]], x[["Sex"]], sep="_")})
	}
	data <- cbind(data, categ)
	unique.cat <- sort(unique(categ))
	L <- length(unique.cat)
	
	factors <- c("dam", "sire", "Strain", "Sex", "Mx1Diplo", "Mx1Diplo6", 
		"Mx1Dam", "Mx1Sire", "categ", batch)
	factors <- factors[factors %in% colnames(data)]
	filenames <- rep(NA, reps)

	for(r in 1:reps){
		numb <- sprintf("%04d", r)

		if("mean"==strategy){
			filenames <- file.path(fdir, "TreatDi_Mean_MP.csv")
			if(file.exists(filenames)){
				cat("Mean match file already exists! \n")
				if(!(TRUE==force)){next}
			}
		}else{
			if(1==reps){
				filenames=matchname
			}else{
				filenames[r]	<- file.path(fdir, paste("TreatDi_MP",numb,".csv", sep=""))
				if(file.exists(filenames[r])){
					cat("Match file for MIMP number ", r, " already exists! \n")
					if(!(TRUE==force)){next}
				}			
			}
		}

		matched.data 	<- NULL
		total.ctrl.discarded <- 0
		total.trt.discarded <- 0		
		total.ctrl.matched 	<- 0
		total.trt.matched 	<- 0
		MAT 			<- NULL

		for(i in unique.cat){
			dat 		<- subset(data, categ==i)
			dat.trt 	<- subset(dat, dat[,trt.colname]==trt.string, droplevels=FALSE)
			dat.ctrl 	<- subset(dat, dat[,trt.colname]==ctrl.string, droplevels=FALSE)
			num.ctrl	<- nrow(dat.ctrl)
			num.trt 	<- nrow(dat.trt)
			mat     	<- NULL

			## Specify the minimum number of flu-infected animals necessary for matching
			if(strategy %in% c("one-to-one", "mean")){trt.min=1}else{
				if("two-to-one"==strategy){
					trt.min=2
				}else{
					if("three-to-one"==strategy){
							trt.min=3
						}else{
							cat("Error in strategy specification."); break()
					}
				}
			}
			## Check that you can make at least one match
			if("one-to-one"==strategy & !(num.ctrl < 1) & !(num.trt < 1)){
				n.matches  		<- min(c(num.trt, num.ctrl))
				dat.ctrl.temp 	<- dat.ctrl
				dat.trt.temp 	<- dat.trt					
				for(k in c(1:n.matches)){
					while((nrow(dat.ctrl.temp) > 0) & (nrow(dat.trt.temp) > 0)){
						sample.ctrl 	<- dat.ctrl.temp[sample(nrow(dat.ctrl.temp), 1, replace=FALSE),]
						sample.trt 		<- dat.trt.temp[sample(nrow(dat.trt.temp), 1, replace=FALSE), ]
						phen.ctrl 		<- sample.ctrl[, phenotypes, drop=FALSE]
						phen.trt 		<- sample.trt[, phenotypes, drop=FALSE]
						ctrl.id 		<- sample.ctrl[,"ID"]
						trt.id			<- sample.trt[,"ID"]
						dat.ctrl.temp 	<- dat.ctrl.temp[!(dat.ctrl.temp$ID %in% sample.ctrl[,"ID"]),]
						dat.trt.temp 	<- dat.trt.temp[!(dat.trt.temp$ID %in% sample.trt[,"ID"]),]
						matched.diff 	<- phen.trt - phen.ctrl
						matched.data 	<- rbind(matched.data, matched.diff)
						mat 			<- cbind(dat[1,factors], ctrl.id=ctrl.id, trt.id=trt.id, match.num=r)
						MAT 			<- rbind(MAT, mat)
						total.ctrl.matched <- total.ctrl.matched + 1
						total.trt.matched <- total.trt.matched + 1
					}	
				}
				total.trt.discarded <- total.trt.discarded + nrow(dat.trt.temp)
				total.ctrl.discarded <- total.ctrl.discarded + nrow(dat.ctrl.temp)
			}
			if("mean"==strategy & !(num.ctrl < 1) & !(num.trt < 1)){
				phen.ctrl 		<- apply(X=dat.ctrl[,phenotypes], MARGIN=2, FUN=mean, na.rm=TRUE)
				phen.trt 		<- apply(X=dat.trt[,phenotypes], MARGIN=2, FUN=mean, na.rm=TRUE)
				ctrl.id 		<- apply(as.matrix(dat.ctrl[,"ID"]),MARGIN=2,paste,collapse=";")
				trt.id 			<- apply(as.matrix(dat.trt[,"ID"]),MARGIN=2,paste,collapse=";")
				n.discarded 	<- 0
				matched.diff 	<- phen.trt - phen.ctrl
				matched.data 	<- rbind(matched.data, matched.diff)
				mat 			<- cbind(dat[1,factors], ctrl.id=ctrl.id, trt.id=trt.id, match.num=r)
				MAT <- rbind(MAT, mat)
				total.ctrl.discarded <- 0 # does not account for categories skipped
				total.trt.discarded <- 0
				total.ctrl.matched <- total.ctrl.matched + num.ctrl
				total.trt.matched <- total.trt.matched + num.trt
			}
			if(("two-to-one"==strategy) & !(num.ctrl < 1) & !(num.trt < 2)){ 
				n.matches 		<- floor(num.trt/(num.ctrl*2))
				dat.ctrl.temp 	<- dat.ctrl
				dat.trt.temp 	<- dat.trt					
				for(k in c(1:n.matches)){
					while((nrow(dat.ctrl.temp) > 0) & (nrow(dat.trt.temp) > 1)){
						sample.ctrl 	<- dat.ctrl.temp[sample(nrow(dat.ctrl.temp), 1, replace=FALSE),]
						sample.trt 		<- dat.trt.temp[sample(nrow(dat.trt.temp), 2, replace=FALSE), ]
						phen.ctrl 		<- sample.ctrl[, phenotypes, drop=FALSE]
						phen.trt 		<- sample.trt[, phenotypes, drop=FALSE]
						ctrl.id 		<- sample.ctrl[,"ID"]
						trt.id.1 		<- sample.trt[1,"ID"]
						trt.id.2 		<- sample.trt[2,"ID"]
						dat.ctrl.temp 	<- dat.ctrl.temp[!(dat.ctrl.temp$ID %in% sample.ctrl[,"ID"]),]
						dat.trt.temp 	<- dat.trt.temp[!(dat.trt.temp$ID %in% sample.trt[,"ID"]),]
						matched.diff 	<- (phen.trt[1,]+phen.trt[2,])/2 - phen.ctrl
						matched.data 	<- rbind(matched.data, matched.diff)
						mat 			<- cbind(dat[1,factors], ctrl.id=ctrl.id, 
											trt.id=paste(trt.id.1, trt.id.2, sep=";"), match.num=r)
						MAT 			<- rbind(MAT, mat)
						total.ctrl.matched <- total.ctrl.matched + 1
						total.trt.matched <- total.trt.matched + 2
					}	
				}
				total.trt.discarded <- total.trt.discarded + nrow(dat.trt.temp)
				total.ctrl.discarded <- total.ctrl.discarded + nrow(dat.ctrl.temp)
			}
			if(("three-to-one"==strategy) & !(num.ctrl < 1) & !(num.trt < 3)){ 
				n.matches 		<- floor(num.trt/(num.ctrl*3))
				dat.ctrl.temp 	<- dat.ctrl
				dat.trt.temp 	<- dat.trt					
				for(k in c(1:n.matches)){
					while((nrow(dat.ctrl.temp) > 0) & (nrow(dat.trt.temp) > 1)){
						sample.ctrl 	<- dat.ctrl.temp[sample(nrow(dat.ctrl.temp), 1, replace=FALSE),]
						sample.trt 		<- dat.trt.temp[sample(nrow(dat.trt.temp), 3, replace=FALSE), ]
						phen.ctrl 		<- sample.ctrl[, phenotypes, drop=FALSE]
						phen.trt 		<- sample.trt[, phenotypes, drop=FALSE]
						ctrl.id 		<- sample.ctrl[,"ID"]
						trt.id.1 		<- sample.trt[1,"ID"]
						trt.id.2 		<- sample.trt[2,"ID"]
						trt.id.3 		<- sample.trt[3,"ID"]
						dat.ctrl.temp 	<- dat.ctrl.temp[!(dat.ctrl.temp$ID %in% sample.ctrl[,"ID"]),]
						dat.trt.temp 	<- dat.trt.temp[!(dat.trt.temp$ID %in% sample.trt[,"ID"]),]
						matched.diff 	<- (phen.trt[1,]+phen.trt[2,]+phen.trt[3,])/3 - phen.ctrl
						matched.data 	<- rbind(matched.data, matched.diff)
						mat 			<- cbind(dat[1,factors], ctrl.id=ctrl.id, 
											trt.id=paste(trt.id.1, trt.id.2, trt.id.3, sep=";"), match.num=r)
						MAT 			<- rbind(MAT, mat)
						total.ctrl.matched <- total.ctrl.matched + 1
						total.trt.matched <- total.trt.matched + 3
					}	
				}
				total.trt.discarded <- total.trt.discarded + nrow(dat.trt.temp)
				total.ctrl.discarded <- total.ctrl.discarded + nrow(dat.ctrl.temp)
			}
		}
		cat("Rep: ", numb, ", Matched(CTRL): ", total.ctrl.matched, ", Matched(TRT): ", total.trt.matched,
				", Discarded(CTRL): ", total.ctrl.discarded, ", Discarded(TRT): ", total.trt.discarded, "\n", sep="")
		post.dat 		<- cbind(MAT, matched.data)
		names(post.dat)	<- c(names(dat[1,factors]), "ctrl.id", "trt.id", "match.num", phenotypes)
		if("mean"==strategy){
			write.csv(post.dat, filenames, row.names=FALSE)
		}else{
			if(1==reps){
				write.csv(post.dat, matchname, row.names=FALSE)
				return(matchname)
			}else{
				write.csv(post.dat, filenames[r], row.names=FALSE)
			}

		}
	}
	return(filenames)
}

#' @title csv.maker: Save analysis to csv
#' @description Save BayesDiallel analysis to csv files.
#' 
#' @param fname the name of the file path
#' @param TreatDi the BayesDiallel object
#' @param burnin specify the burnin
#' @param ... additional arguments
#' @return Save csv files
#' @examples
#' ## not run
#' @export
csv.maker <- function(fname, TreatDi, burnin, ...){
	# Make CSV Files
	burnin <- as.numeric(as.character(burnin))
	HPDTable <- TreatDi$AllDiallelObs[[1]]$HPDTable
	write.csv(HPDTable, file.path(fname, "HPDTable.csv"), row.names=TRUE)
	MyPostSum <- PosteriorPredSummary(TreatDi$AllDiallelObs[[1]], TreatDi, 
		burnin = 500, AFD=TreatDi, keep=TRUE);
	write.csv(MyPostSum, file.path(fname, "MyPostSum.csv"), row.names=TRUE)
	MySum <- summary(TreatDi$AllDiallelObs[[1]]$cent.chains, burnin = 500, AFD=TreatDi)
	write.csv(MySum$statistics, file.path(fname, "MySumStats.csv"), row.names=TRUE)
	write.csv(MySum$quantiles, file.path(fname, "MySumQuants.csv"), row.names=TRUE)
	psq <- TreatDi$AllDiallelObs[[1]]$PosteriorPSqTable
	psq.sub <- TreatDi$AllDiallelObs[[1]]$PosteriorPSqSubTable
	dsq <- TreatDi$AllDiallelObs[[1]]$PosteriorDSqTable
	dsq.sub <- TreatDi$AllDiallelObs[[1]]$PosteriorDSqSubTable

	try(colnames(psq[[2]])[c(3,1,5)] <- c("mean", "lower.bound", "upper.bound")) ## mean, in this case median
	try(colnames(psq.sub[[2]])[c(3,1,5)] <- c("mean", "lower.bound", "upper.bound"))
	try(colnames(dsq[[2]])[c(3,1,5)] <- c("mean", "lower.bound", "upper.bound"))
	try(colnames(dsq.sub[[2]])[c(3,1,5)] <- c("mean", "lower.bound", "upper.bound"))
	write.csv(psq[[2]], file=file.path(fname, "PSqTable.csv"), row.names=TRUE)
	write.csv(psq[[1]], file=file.path(fname, "PSqTable01.csv"), row.names=TRUE)
	write.csv(psq.sub[[2]], file=file.path(fname, "PSqSubTable.csv"), row.names=TRUE)
	write.csv(dsq[[2]], file=file.path(fname, "DSqTable.csv"), row.names=TRUE)
	write.csv(dsq.sub[[2]], file=file.path(fname, "DSqSubTable.csv"), row.names=TRUE)
	try(BSAFDMIP <- TreatDi$BSAFD$MIP)
	try(write.csv(BSAFDMIP, file.path(fname, "MIP.csv"), row.names=TRUE))
}

#' @title stackMatchReps: Combine results across replicates
#' @description Stack multiple imputation matched replicates.
#' 
#' @param filenames stack mp imputations based on filenames
#' @param fdir file directory
#' @param ... additional arguments
#' @return saves csv and plot
#' @examples
#' data(FluDiData)
#' write.csv(FluDiData, file="FluDiData.csv")
#' filename <- "FluDiData.csv"
#' data <- read.csv("FluDiData.csv")
#' reps <- 2
#' trt.string <- "FLU"
#' ctrl.string <- "MOCK" 
#' strategy <- "one-to-one"
#' fdir <- getwd()
#' filenames	<- make.matches(data=data, reps=reps, trt.string=trt.string, 
#'                              ctrl.string=ctrl.string, strategy=strategy, fdir=".")
#' MIMPstack	<- stackMatchReps(filenames, fdir=fdir)
#' @export
stackMatchReps <- function(filenames, fdir, ... ){
	raw.names 	<- c("D0", "D1", "D2", "D3", "D4")
	pct.names 	<- c("pct_D1", "pct_D2", "pct_D3", "pct_D4")
	auc.names 	<- c("auc0", "auc1", "auc2")
	phens 		<- c(raw.names, pct.names, auc.names)
	stack 		<- NULL
	layers 		<- NULL

	for(i in 1:length(filenames)){
		layers[[i]] 	<- read.csv(filenames[i])
		stack 			<- rbind(stack, layers[[i]])
	}

	write.csv(stack, file.path(fdir, "stack_match_reps.csv"), row.names=FALSE)

	pdf(file.path(fdir, "stack_match_reps_plot.pdf"), width=18, height=18)
	par(mfrow=c(4,3))
   	col=rgb(0,0,0,maxColorValue=1, alpha=0.05)
   	for(j in 1:length(phens)){
   		if(phens[j] %in% c("D0", "D1", "D2")){
   		   	ylim=c(0, 0.27)
		   	xlim=c(-15,15)
   		}
   		if(phens[j] %in% c("pct_D1", "pct_D2")){
   			ylim=c(0, 0.22)
		   	xlim=c(-22,18)
   		}
   		if(phens[j] %in% c("pct_D3", "pct_D4")){
	   	   	ylim=c(0, 0.08)
		   	xlim=c(-35,20)
		}
   		if(phens[j] %in% auc.names){
   			ylim=c(0, 0.043)
   			xlim=c(-45, 60)
   		}
   		
  		plot(density(layers[[1]][,phens[j]], na.rm=TRUE), col=col, main=phens[[j]],
  			ylim=ylim, xlim=xlim, xlab="", ylab="")
  		abline(v=0, lty=2, lwd=4)
		for(s in 2:length(layers)){
			lines(density(layers[[s]][,phens[j]], na.rm=TRUE), col=col)
		}
		# lines(density(stack[,phens[j]], na.rm=TRUE), col="red", cex=4)
   	}
	dev.off()

	return(stack)
}

#' @title run.tr.diallel: run tRD analysis.
#' @description Run the treatment response diallel analysis.
#' 
#' @param filename name of data file
#' @param savedir directory to save files in
#' @param treatment indicate subset
#' @param trt.colname indicate treatment columname
#' @param phenotype indicate parameter of interest
#' @param random name of random effect, or NULL
#' @param fixed name of fixed effect, or NULL
#' @param type specify matched pairs or other type of analysis
#' @param BS whether or not to use BayesSpike
#' @param thin thinning factor for chains
#' @param lengthChains MCMC chain length
#' @param numChains number of independent MCMC chains
#' @param burnin amount of burn-in (discarded) samples from warm-up
#' @param force logical indicating whether to proceed despite warnings/errors
#' @param ... additional arguments
#' @return runs the analysis, generating data files.
#' @examples
#' ## This example may take a couple of minutes
#' data(FluDiData)
#' data(PiximusData, envir = environment())
#' write.csv(FluDiData, file="FluDiData.csv")
#' filename <- "FluDiData.csv"
#' args <- c("All", "D0", "NULL", "NULL", "pre", FALSE, "Trt")
#' treatment 		<- args[1]
#' phenotype  		<- args[2]
#' random 			<- args[3]
#' fixed 			<- args[4]
#' type				<- args[5]
#' BS 				<- as.logical(args[6])
#' trt.colname      <- args[7]
#' savedir			<- "."
#' run.tr.diallel(filename=filename, savedir=savedir, treatment=treatment, 
#'	trt.colname=trt.colname, phenotype=phenotype, random=random,
#'	fixed=fixed, type=type, BS=BS, thin=1, lengthChains=3000, numChains=5, 
#'	burnin=500, force = TRUE)
#' @export
run.tr.diallel <- function(	filename, savedir, treatment, trt.colname, 
							phenotype, random, fixed, type, BS=TRUE, 
							thin, lengthChains=3000, numChains, 
							burnin, force=TRUE, ...){
	f <- filename

	data(PiximusData, envir = environment())

	plot.dir <- getwd()
	data <- read.csv(f)

	## Generate female data column
	data$is_female <- sapply(data$Sex, function(x){ifelse("F"==x, 1, 0)})

	## Subset data
	data$starting_weight <- data$D0

	if(!("All"==treatment)){
		data <- subset(data, data[,trt.colname]==treatment)
	}

	if(!dir.exists(savedir)){
		dir.create(savedir, recursive=TRUE)
		print(paste("Creating directories and file:", savedir, sep=""))
	}

	ZeroOutBeforeSample <- 0

	father.strain="sire"
	mother.strain="dam"
	is.female="is_female"

	RandomEffects=if("NULL" %in% random){NULL}else{as.list(random)}
	nameRandomEffects=if("NULL" %in% random){NULL}else{random}
	FixedEffects=if("NULL" %in% fixed){NULL}else{fixed}
	nameFixedEffects=if("NULL" %in% fixed){NULL}else{fixed}

	for(r in 1:length(random)){
		try(if("NULL" %in% random[[r]]){random[[r]] <- NULL})
		try(if(!(is.null(random[[r]]))){
			random.list <- unique(data[,random[[r]]])
			random.list.id <- unique(as.integer(as.factor(data[,random[[r]]])))
			print(paste(random[[r]], ":", random.list, "=", random.list.id, sep=""))
			data[,random[[r]]] <- as.integer(as.factor(data[,random[[r]]]))
		})
	}

	print(dim(data))
	print(str(data))
	print(head(data))
	print(paste("RandomEffects", as.character(random), sep=": "))
	print(paste("FixedEffects", as.character(fixed), sep=": "))

	Models=c('fulls', 'fullu', 'fulls, df_6', 'BSabm', 'BSasbsms', 'B,v,w,w_s,b,m,ms,mu', 
	 'Mu, strain, gender, gender-mother, gender-inbred, symmetric-cross, inbred', 
	 'Mu, inbred, inbred-strain, strain, gender-inbred')[1]

	if(TRUE==as.logical(BS)){DoBayesSpike=TRUE}else{DoBayesSpike=FALSE}
	# if(file.exists(file.path(savedir, "TreatDi_MIMP.RData"))){
	# 	cat("Treat Di MIMP RData file already exists! \n")
	# 	if(!(force==TRUE)){return(NULL)}
	# }

	TreatDi = DiallelAnalyzer(data = data,
	    father.strain=father.strain,
	    mother.strain=mother.strain,
	    phenotype=phenotype,
	    is.female=is.female,
	    FixedEffects=FixedEffects,
	    nameFixedEffects=nameFixedEffects,
	    RandomEffects=RandomEffects,
	    nameRandomEffects=nameRandomEffects,
	    Models=Models,
	    sep="", na.strings="NA",
	    sigmasq.start=1, numChains=numChains, lengthChains=lengthChains,
	    Verbose=TRUE,
	    burnin=burnin, DIC.Only=FALSE,
	    DoBayesSpike=DoBayesSpike,
	    BSSaveDir=savedir, thin=thin, 
	    SaveAFDFile=file.path(savedir, "TreatDi_MIMP.RData"),
	    LogTransform=FALSE,
	    ZeroOutBeforeSample=0, DoGelman=TRUE, DoFirstCenter=TRUE)
	##----------------------
	## STACKED MCMC TABLES
	##----------------------
	# if(file.exists(file.path(savedir, "StackedTreatDi.csv"))){
	# 	cat("StackedTreatDi file already exists! \n")
	# 	if(!(force==TRUE)){return(NULL)}
	# }
	var.labels.all <- varnames(TreatDi$AllDiallelObs[[1]]$centered.chains)
	new.labels <- sapply(X=var.labels.all, FUN=var.translate, USE.NAMES=FALSE)

	stacked.unburned <- mcmc.stack(TreatDi$AllDiallelObs[[1]]$cent.chains)
	varnames(stacked.unburned) <- new.labels
	stacked.unburned <- as.matrix(stacked.unburned)
	write.csv(stacked.unburned, 
		file=gzfile(file.path(savedir, "StackedTreatDiUnburned.csv.gz")),
		row.names=FALSE)

	stacked.burned <- mcmc.stack.and.burn(TreatDi$AllDiallelObs[[1]]$cent.chains, burnin=burnin)
	varnames(stacked.burned) <- new.labels
	stacked.burned <- as.matrix(stacked.burned)
	write.csv(stacked.burned, 
		file=gzfile(file.path(savedir, "StackedTreatDi.csv.gz")),
		row.names=FALSE)

	PPSq <- TreatDi$AllDiallelObs[[1]]$PosteriorPSq()
	stacked.posterior <- PPSq$LPPSq

	stacked.matches.post <- mcmc.stack(stacked.posterior)
 	stacked.matches.post <- as.matrix(stacked.matches.post)
	write.csv2(stacked.matches.post, 
		file=gzfile(file.path(savedir, "StackedTreatDiPostUnburned.csv.gz")),
		row.names=FALSE)

	stacked.matches.post <- mcmc.stack.and.burn(stacked.posterior, burnin=burnin)
 	stacked.matches.post <- as.matrix(stacked.matches.post)
	write.csv2(stacked.matches.post, 
		file=gzfile(file.path(savedir, "StackedTreatDiPost.csv.gz")),
		row.names=FALSE)

	ADO = TreatDi$AllDiallelObs[[1]]
	MyPostSum = PosteriorPredSummary(ADO, TreatDi, burnin = burnin, AFD=TreatDi, keep = TRUE);
	#write.csv(MyPostSum, file.path(fdir, "MyPostSum.csv"))
	#dat <- read.csv("MyPostSum.csv", row.names=1)

	stacked.matches.postpred <- mcmc.stack(ADO$PostKeeper$FakeCoda)
	stacked.matches.postpred <- as.matrix(stacked.matches.postpred)
	write.csv(stacked.matches.postpred,
		file=gzfile(file.path(savedir, "StackedTreatDiPostPred.csv.gz")),
		row.names=FALSE)

	##-----------------------------------------------
	## SAVE CSV FILES
	##-----------------------------------------------
	for(i in c(1:numChains)){
		write.csv(TreatDi$AllDiallelObs[[1]]$cent.chains[[i]], 
			file=file.path(savedir, paste0("TreatDi", i, ".csv")),
			row.names=FALSE)
	}

	##----------------------
	## Make CSV files
	##----------------------
	csv.maker(fname=savedir, TreatDi=TreatDi, burnin=burnin)

	returned <- list()
	returned[[1]] <- savedir
	returned[[2]] <- TreatDi

	return(returned)
}

#' @title inner.plotter: Generate plots for each rep
#' @description Generate plots from AFD object for each rep.
#' 
#' @param AFD Diallel object to make plots out of.
#' @param plotdir refers to the compound path.
#' @param xlim limits for x-axis
#' @param ... additional arguments
#' @return Makes standard plots in the specified plot directory.
#' @examples
#' ## This example may take a couple of minutes
#' data(FluDiData, envir = environment())
#' data(PiximusData, envir = environment())
#' write.csv(FluDiData, file="FluDiData.csv")
#' filename <- "FluDiData.csv"
#' args <- c("All", "D0", "NULL", "NULL", "pre", FALSE, "Trt")
#' treatment 		<- args[1]
#' phenotype  		<- args[2]
#' random 			<- args[3]
#' fixed 			<- args[4]
#' type				<- args[5]
#' BS 				<- args[6]
#' trt.colname      <- args[7]
#' savedir			<- "."
#' returned <- run.tr.diallel(filename=filename, savedir=savedir, treatment=treatment, 
#'	trt.colname=trt.colname, phenotype=phenotype, random=random,
#'	fixed=fixed, type=type, BS=BS, thin=1, lengthChains=3000, numChains=5, 
#'	burnin=500, force = TRUE)
#' plotdir <- returned[[1]]
#' AFD <- returned[[2]]
#' inner.plotter(AFD=AFD, plotdir=plotdir)
#' @export
inner.plotter <- function(AFD, plotdir, xlim=c(-12,12), ...){

	pdf(file.path(plotdir, "HPD_plots.pdf"), width=11, height=8.5)
	par(mfrow=c(1,4))

	var.labels.all <- varnames(AFD$AllDiallelObs[[1]]$centered.chains)
	new.labels <- sapply(X=var.labels.all, FUN=var.translate, USE.NAMES=FALSE)
	var.labels <- grep(pattern="^aj:", x=var.labels.all, value=TRUE)
	var.labels <- c(var.labels, grep(pattern="^motherj:", x=var.labels.all, value=TRUE))
	var.labels <- c(var.labels, grep(pattern="^BetaInbred:Av", x=var.labels.all, value=TRUE))
	var.labels <- c(var.labels, grep(pattern="^dominancej:", x=var.labels.all, value=TRUE))
	new.labels <- sapply(X=var.labels, FUN=var.translate, USE.NAMES=FALSE)
	plot.hpd(AFD$AllDiallelObs[[1]]$centered.chains[,var.labels], names=new.labels, 
		xlim=xlim, main="General effects")

	abline(v=0, col="darkgray")
	segments(x0=xlim[1], x1=xlim[2], y0=9,y1=9, col="darkgray")
	segments(x0=xlim[1], x1=xlim[2], y0=17,y1=17, col="darkgray")

	var.labels <- grep(pattern="^SymCrossjk:j", x=var.labels.all, value=TRUE)
	var.labels <- c(var.labels, grep(pattern="^ASymCrossjkDkj:j", x=var.labels.all, value=TRUE))
	new.labels <- sapply(X=var.labels, FUN=var.translate, USE.NAMES=FALSE)
	plot.hpd(AFD$AllDiallelObs[[1]]$centered.chains[,var.labels], names=new.labels, 
		xlim=xlim, main="Strainpair-specific")
	abline(v=0, col="darkgray");
	segments(x0=-20, x1=xlim[2], y0=28,y1=28, col="darkgray")

	var.labels <- grep(pattern="^Gender:Av", x=var.labels.all, value=TRUE)
	var.labels <- c(var.labels, grep(pattern="Gender:aj:", x=var.labels.all, value=TRUE))
	var.labels <- c(var.labels, grep(pattern="Gender:motherj:", x=var.labels.all, value=TRUE))
	var.labels <- c(var.labels, grep(pattern="^BetaInbred:Gender:Av", x=var.labels.all, value=TRUE))
	var.labels <- c(var.labels, grep(pattern="Gender:dominancej:", x=var.labels.all, value=TRUE))
	new.labels <- sapply(X=var.labels, FUN=var.translate, USE.NAMES=FALSE)
	plot.hpd(AFD$AllDiallelObs[[1]]$centered.chains[,var.labels], names=new.labels, 
		xlim=xlim, main="Sex-specific")
	abline(v=0, col="darkgray")
	segments(x0=xlim[1], x1=xlim[2], y0=9,y1=9, col="darkgray")
	segments(x0=xlim[1], x1=xlim[2], y0=17,y1=17, col="darkgray")

	var.labels <- grep(pattern="^Gender:SymCrossjk:j", x=var.labels.all, value=TRUE)
	var.labels <- c(var.labels, grep(pattern="^Gender:ASymCrossjkDkj:j", x=var.labels.all, value=TRUE))
	new.labels <- sapply(X=var.labels, FUN=var.translate, USE.NAMES=FALSE)
	plot.hpd(AFD$AllDiallelObs[[1]]$centered.chains[,var.labels], names=new.labels, 
		xlim=xlim, main="Sex/strainpair-specific")
	abline(v=0, col="darkgray")
	segments(x0=xlim[1], x1=xlim[2], y0=28,y1=28, col="darkgray")

	dev.off()

	pdf(file.path(plotdir, "HPD_plots_additional.pdf"), width=11, height=8.5)
		par(mfrow=c(1,2))
		var.labels <- grep(pattern="^Mu", x=var.labels.all, value=TRUE)
		var.labels <- c(var.labels, grep(pattern="^Sigma:", x=var.labels.all, value=TRUE))
		var.labels <- c(var.labels, grep(pattern="^FixedEffect:", x=var.labels.all, value=TRUE))
		new.labels <- sapply(X=var.labels, FUN=var.translate, USE.NAMES=FALSE)
		plot.hpd(AFD$AllDiallelObs[[1]]$centered.chains[,var.labels], 
			names=new.labels, cex.label=0.7)
		abline(v=0, col="darkgray")

		var.labels <- grep(pattern="^tau", x=var.labels.all, value=TRUE)
		var.labels <- c(var.labels, grep(pattern="^Batch:", x=var.labels.all, value=TRUE))
		var.labels <- c(var.labels, grep(pattern="^Random:", x=var.labels.all, value=TRUE))
		new.labels <- sapply(X=var.labels, FUN=var.translate, USE.NAMES=FALSE)
		plot.hpd(AFD$AllDiallelObs[[1]]$centered.chains[,var.labels], 
			names=new.labels, cex.label=0.7)
		abline(v=0, col="darkgray")
	dev.off()

	pdf(file.path(plotdir,"Heatmap_plots.pdf"), width=11, height=8.5)
		par(mfrow=c(1,1))
		## Plot females, observed versus fitted
		AFD$AllDiallelObs[[1]]$TwoDiallelPlot(FemalePlot = TRUE,
		MaleAgainstFemale = FALSE);
		## Plot males, observed versus fitted
		AFD$AllDiallelObs[[1]]$TwoDiallelPlot(FemalePlot = FALSE,
		MaleAgainstFemale = FALSE);
		## Plot male predicted versus female predicted
		AFD$AllDiallelObs[[1]]$TwoDiallelPlot(FemalePlot = TRUE,
		MaleAgainstFemale = TRUE);
		## Plot male observed vresus female observed
		AFD$AllDiallelObs[[1]]$TwoDiallelPlot(PlotOnlyObserved = TRUE)
	dev.off()

	new.strain.names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
	print(new.strain.names)

	pdf(file.path(plotdir, "TwoPlot.pdf"), width=11, height=8.5)
	par(mfrow=c(1,1))
	AFD$AllDiallelObs[[1]]$TwoDiallelPlot(PlotObservedVersusExpectedAll=TRUE,
		show.strain.names=TRUE, show.strain.coords=TRUE, 
		new.strain.names=new.strain.names)
	dev.off()

	pdf(file.path(plotdir, "TwoPlot_predicted.pdf"), width=11, height=8.5)
	par(mfrow=c(1,1))
	AFD$AllDiallelObs[[1]]$TwoDiallelPlot(MaleAgainstFemale=TRUE,
		PlotOnlyObserved=FALSE,
		show.strain.names=TRUE, show.strain.coords=TRUE, 
		new.strain.names=new.strain.names)
	dev.off()

	pdf(file.path(plotdir, "TwoPlot_observed.pdf"), width=11, height=8.5)
	par(mfrow=c(1,1))
	AFD$AllDiallelObs[[1]]$TwoDiallelPlot(MaleAgainstFemale=TRUE,
		PlotOnlyObserved=TRUE,
		show.strain.names=TRUE, show.strain.coords=TRUE, 
		new.strain.names=new.strain.names)
	dev.off()

	pdf(file.path(plotdir, "TwoPlot_males.pdf"), width=11, height=8.5)
	par(mfrow=c(1,1))
	AFD$AllDiallelObs[[1]]$TwoDiallelPlot(FemalePlot=FALSE,
		show.strain.names=TRUE, show.strain.coords=TRUE, 
		new.strain.names=new.strain.names)
	dev.off()

	pdf(file.path(plotdir, "TwoPlot_females.pdf"), width=11, height=8.5)
	par(mfrow=c(1,1))
	AFD$AllDiallelObs[[1]]$TwoDiallelPlot(FemalePlot=TRUE,
		show.strain.names=TRUE, show.strain.coords=TRUE, 
		new.strain.names=new.strain.names)
	dev.off()

	print(AFD$strain.map)

	pdf(file.path(plotdir, "Straw.pdf"), width=11, height=8.5)
	AFD$AllDiallelObs[[1]]$PlotStrawPlot()
	dev.off()

	png(file.path(plotdir, "Straw.png"), width=550, height=425)
	AFD$AllDiallelObs[[1]]$PlotStrawPlot()
	dev.off()

	try(plotVarps(fdir=plotdir, fname="PSqTable.csv"))

	try(plotVarps(fdir=plotdir, fname="PSqSubTable.csv"))

	try(plotMips(fdir=plotdir, fname="MIP.csv"))
}

#' @title poe.analyzer: Generate plots of POE
#' @description Generate plots of estimated parent-of-origin effects.
#' 
#' @param stacked stacked chains
#' @param plotdir the directory where the POE plots should be placed
#' @param ... additional arguments
#' @return saves plots of asymmetric cross-pair-specific effects
#' @examples
#' #not.run
#' @export
poe.analyzer <- function(stacked, plotdir, ...){
	its <- nrow(stacked)
	array.f <- array(data=rep(NA, 8*8*its), dim=c(8,8,its))
	array.m <- array(data=rep(NA, 8*8*its), dim=c(8,8,its))

	strain <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")

	for(tt in 1:its){
		for(j in 1:8){
			for(k in 1:j){
				this.layer	 	<- 	stacked[tt,]
				unsexed 		<- 	2*this.layer[paste0("maternal:", strain[j])][[1]] -
									2*this.layer[paste0("maternal:", strain[k])][[1]] +
								  	this.layer[paste0("w:", strain[j], ";", strain[k])][[1]]
				female.comp		<- 	0.5*( 2*this.layer[paste0("f:maternal:", strain[j])][[1]] -
							  		2*this.layer[paste0("f:maternal:", strain[k])][[1]] +
							  		this.layer[paste0("f:w:", strain[j], ";", strain[k])][[1]])
				array.f[j,k,tt] <- 	unsexed + female.comp
				array.m[j,k,tt] <- 	unsexed - female.comp
				array.f[k,j,tt] <- 	-1*(unsexed + female.comp)
				array.m[k,j,tt] <- 	-1*(unsexed - female.comp)
			}
			array.f[j,j,tt] <- 0
			array.m[j,j,tt] <- 0
		}
	}

	xlim=c(-10,10)
	ylim=c(0,0.2)

	pdf(file.path(plotdir, "asymm.plot.females.pdf"), width=20, height=20)
	par(mfrow=c(8,8))
	for(j in 1:8){
		for(k in 1:8){
	#		hist(array.f[j,k,], main=paste0(strain[j],"x",strain[k]), 
	#			xlab="asymm.effect", ylab="freq", xlim=xlim)
			plot(density(array.f[j,k,], na.rm=TRUE), xlim=xlim, xlab="", ylab="", ylim=ylim, 
				main=paste0(strain[j],"x",strain[k], " (f)"))
			abline(v=0, lty=2)
		}
	}
	dev.off()

	pdf(file.path(plotdir, "asymm.plot.males.pdf"), width=20, height=20)
	par(mfrow=c(8,8))
	for(j in 1:8){
		for(k in 1:8){
	#		hist(array.f[j,k,], main=paste0(strain[j],"x",strain[k]), 
	#			xlab="asymm.effect", ylab="freq", xlim=xlim)
			plot(density(array.m[j,k,]), xlim=xlim, xlab="", ylab="", ylim=ylim, 
				main=paste0(strain[j],"x",strain[k]," (m)"))
			abline(v=0, lty=2)
		}
	}
	dev.off()

	pdf(file.path(plotdir, "asymm.ci.plot.pdf"), width=15, height=15)
	par(mfrow=c(8,8))
	for(j in 1:8){
		for(k in 1:8){
			if(j==k){
				frame()
	#			plot(1, type="n", axes=F, xlab="", ylab="")
			}else{
				plot.hpd(as.mcmc(cbind(F=array.f[j,k,], M=array.m[j,k,])), xlim=c(-8,8),
					xlab=paste0(strain[j],"x",strain[k]),
					axes=sides(top=FALSE, bottom=TRUE, left=TRUE, right=FALSE),
				    name.margin=4.1,
				    title.margin=3.1,
	    			bottom.margin=4.1,
	    			right.margin=1.1)
				abline(v=0, col="red")
			}
		}
	}
	dev.off()
}

#' @title hpd.plotter: Make HPD plots
#' @description Generate highest posterior density plots.
#' 
#' @param plotdir the directory to save the HPD plots
#' @param chain.object the MCMC chain object
#' @param batched logical that defines whether to plot batch effect estimates
#' @param fixed defines whether to plot fixed effect estimates 
#' @param xlim gives limits of x-axis
#' @param ... additional arguments
#' @return make HPD plots and save them to plotdir
#' @examples
#' ## not run
#' @export
hpd.plotter <- function(plotdir, chain.object, batched=c(FALSE, TRUE)[1], 
	fixed=c(FALSE, TRUE)[1], xlim=c(-10,10), ...){
	pdf(file.path(plotdir, "HPD_plots.pdf"), width=11, height=8.5)
	par(mfrow=c(1,4))
	var.labels.all <- varnames(chain.object)
	var.labels <- grep(pattern="^additive:", x=var.labels.all, value=TRUE)
	var.labels <- c(var.labels, grep(pattern="^maternal:", x=var.labels.all, value=TRUE))
	var.labels <- c(var.labels, grep(pattern="^inbreed.overall", x=var.labels.all, value=TRUE))
	var.labels <- c(var.labels, grep(pattern="^inbreeding:", x=var.labels.all, value=TRUE))
	plot.hpd(chain.object[,var.labels], names=var.labels, 
		xlim=xlim, main="General effects")
	abline(v=0, col="gray")
	
	var.labels <- grep(pattern="^v:", x=var.labels.all, value=TRUE)
	var.labels <- c(var.labels, grep(pattern="^w:", x=var.labels.all, value=TRUE))
	plot.hpd(chain.object[,var.labels], names=var.labels, 
		xlim=xlim, main="Strainpair-specific")
	abline(v=0, col="gray");

	var.labels <- grep(pattern="^female.overall", x=var.labels.all, value=TRUE)
	var.labels <- c(var.labels, grep(pattern="^f:additive:", x=var.labels.all, value=TRUE))
	var.labels <- c(var.labels, grep(pattern="^f:maternal:", x=var.labels.all, value=TRUE))
	var.labels <- c(var.labels, grep(pattern="^female.inbred", x=var.labels.all, value=TRUE))
	var.labels <- c(var.labels, grep(pattern="^f:inbreeding", x=var.labels.all, value=TRUE))
	plot.hpd(chain.object[,c(var.labels)],
		names=var.labels, xlim=xlim, main="Sex-specific"); 
	abline(v=0, col="gray")

	var.labels <- grep(pattern="^f:v:", x=var.labels.all, value=TRUE)
	var.labels <- c(var.labels, grep(pattern="^f:w:", x=var.labels.all, value=TRUE))
	plot.hpd(chain.object[,c(var.labels)],
		names=var.labels, xlim=xlim, main="Sex/strainpair-specific"); 
	abline(v=0, col="gray")
	dev.off()

	pdf(file.path(plotdir, "HPD_plots_additional.pdf"), width=11, height=8.5)
	par(mfrow=c(1,2))
	if(TRUE==as.logical(batched)){
		par(mfrow=c(1,3))
	}
	var.labels <- grep(pattern="^Mu", x=var.labels.all, value=TRUE)
	var.labels <- c(var.labels, grep(pattern="^Sigma", x=var.labels.all, value=TRUE))
	var.labels <- c(var.labels, grep(pattern="^FixedEffect:", x=var.labels.all, value=TRUE))
	new.labels <- sapply(X=var.labels, FUN=var.translate, USE.NAMES=FALSE)
	plot.hpd(chain.object[,var.labels], 
	         names=new.labels, cex.label=0.7, main="Mu, Sigma^2 and Fixed")
	abline(v=0, col="gray")
	
	var.labels <- grep(pattern="^tau", x=var.labels.all, value=TRUE)
	new.labels <- sapply(X=var.labels, FUN=var.translate, USE.NAMES=FALSE)
	plot.hpd(chain.object[,var.labels], 
	         names=new.labels, cex.label=0.7, main="Tau^2")
	abline(v=0, col="gray")

	if(TRUE==as.logical(batched)){
		var.labels <- grep(pattern="^Batch:", x=var.labels.all, value=TRUE)
		var.labels <- c(var.labels, grep(pattern="^Random", x=var.labels.all, value=TRUE))
		plot.hpd(chain.object[,var.labels],
		        names=var.labels, cex.label=0.7, main="Batch and/or Random")
		abline(v=0, col="gray")
	}
	dev.off()
	## To add: heatmap and straw plots
}

#' @title simplify.labels: Convert labels
#' @description Convert/simplify labels such that they are readable in a table.
#' 
#' @param vector.of.strings defines the arguments to be substituted, if applicable, with new strings
#' @param ... additional arguments
#' @return convert effects/parameters to more readable formats
#' @examples
#' simplify.labels("Gender:motherj")
#' @export
simplify.labels <- function(vector.of.strings, ...){
	var1 <- vector.of.strings
	var1 <- sub("ProbFixed:", "", var1)
	var1 <- sub("Prob:tau:", "", var1)
	var1 <- sub("Gender:aj", "sex-by-additive", var1)
	var1 <- sub("Gender:motherj", "sex-by-maternal", var1)
	var1 <- sub("Gender:dominancej", "sex-by-inbreeding", var1)
	var1 <- sub("Gender:BetaHybrid:Av", "sex-by-inbreed.overall", var1)
	var1 <- sub("Gender:ASymCrossjkDkj", "sex-by-asymmetric.epistatic", var1)
	var1 <- sub("Gender:SymCrossjk", "sex-by-symmetric.epistatic", var1)
	var1 <- sub("aj", "additive", var1)
	var1 <- sub("motherj", "maternal", var1)
	var1 <- sub("dominancej", "inbreeding", var1)
	var1 <- sub("BetaInbred:Av", "inbreed.overall", var1)
	var1 <- sub("ASymCrossjkDkj", "asymmetric.epistatic", var1)
	var1 <- sub("SymCrossjk", "symmetric.epistatic", var1)
	var1 <- sub("BetaInbred:Gender:Av", "female.inbred", var1)
	var1 <- sub("Gender:Av", "female.overall", var1)
	var1 <- sub("RandomEffect", "Batch", var1)
	return(var1)
}

#' @title sub.labels: Replace labels with just the strain information
#' @description Replace labels with just the strain information, removing the prefix portion of the string.
#' 
#' @param vector.of.strings defines the arguments to be substituted, if applicable, with new strings
#' @param ... additional arguments
#' @return returns vector of strings that contains just the strain information
#' @examples
#' simplify.labels("additive:AJ")
#' @export
sub.labels <- function (vector.of.strings, ...){
    var1 <- vector.of.strings
    var1 <- sub("f:additive:", "", var1)
    var1 <- sub("f:maternal:", "", var1)
    var1 <- sub("f:inbreeding:", "", var1)
    var1 <- sub("f:v:", "", var1)
    var1 <- sub("f:w:", "", var1)
    var1 <- sub("additive:", "", var1)
    var1 <- sub("maternal:", "", var1)
    var1 <- sub("inbreeding:", "", var1)
    var1 <- sub("v:", "", var1)
    var1 <- sub("w:", "", var1)
    return(var1)
}

#' @title plotVarps: Plot the variance projections
#' @description Plot variance projection (VarP) confidence intervals.
#' 
#' @param fdir file directory
#' @param fname file name
#' @param reorder optional reordering of arguments
#' @param ... additional arguments
#' @return make variance projection plots and save them to file directory
#' @examples
#' ## not run
#' @export
plotVarps <- function(fdir, fname, reorder=FALSE, ...){
	psq <- read.csv(file.path(fdir, fname))
	psq <- psq[!grepl("FixedEffect",psq$X),]
	psq <- psq[!grepl("RandomEffect",psq$X),]
	if(TRUE==reorder){
		rownames(psq) <- PSq$X
		psq <- psq[c('aj','motherj','BetaInbred:Av','dominancej','SymCrossjk','ASymCrossjkDkj',
					 'Gender:Av',
					 'Gender:aj','Gender:motherj',
					 'BetaInbred:Gender:Av','Gender:dominancej','Gender:SymCrossjk',
					 'Gender:ASymCrossjkDkj','total.explained','Noise'),]
	}
	pdf(file.path(fdir, paste(file_path_sans_ext(fname), "_fig.pdf", sep="")), width=10, height=6)
	par(mar=c(5.1,13.5,2.1,12.1))
	rows <- nrow(psq)
	plot(y=c(0.5, rows+0.5), x=c(-0.01, 1.01), yaxs="i", xaxs="i",
		col="white", ylab="", xlab="VarP C.I.", yaxt="n")
	abline(v=c(0, 0.2, 0.4, 0.6, 0.8, 1), lty=2, col="grey")
	abline(h=2.5, lwd=1.5)
	segments(x0=psq$lower.bound, x1=psq$upper.bound, 
		y0=c(rows:1), y1=c(rows:1), lwd=2)
	axis(side=2, at=c(rows:1), labels=as.vector(simplify.labels(psq$X)), las=1)
	rightlabels <- paste(as.character(sprintf("%.2f", round(100*psq$mean, digits=2))), 
		"% (", as.character(sprintf("%.2f", round(100*psq$lower.bound, digits=2))), 
		"%, ", as.character(sprintf("%.2f", round(100*psq$upper.bound, digits=2))), 
		"%)", sep="")
	axis(side=4, at=c(rows:1), labels=rightlabels, tick=FALSE, las=1)
	points(x=psq$mean, y=c(rows:1), pch=16, col="black", cex=1.5)
	dev.off()
}

#' @title plotMips: Plot the MIPs
#' @description Plot model inclusion probabilities, saved as a csv.
#' 
#' @param fdir file directory
#' @param fname file name
#' @param reorder optional reordering of arguments
#' @param ... additional arguments
#' @return make MIP plots and save them to file directory
#' @examples
#' ## not run
#' @export
plotMips <- function(fdir, fname, reorder=FALSE, ...){
	mip <- read.csv(file.path(fdir, fname))
	pdf(file.path(fdir, paste(file_path_sans_ext(fname), "_fig.pdf", sep="")), width=10, height=6)
	par(mar=c(5.1,13.5,2.1,0.1))
	if(TRUE==reorder){
		reorder <- c('ProbFixed:Mu',
					'ProbFixed:FixedEffect:1',
					'Prob:tau:RandomEffect:1',
					'Prob:tau:aj',
					'Prob:tau:motherj',
					'ProbFixed:BetaInbred:Av','Prob:tau:dominancej',
					'Prob:tau:SymCrossjk',
					'Prob:tau:ASymCrossjkDkj',
					'ProbFixed:Gender:Av',
					'Prob:tau:Gender:aj',
					'Prob:tau:Gender:motherj',
					'ProbFixed:BetaInbred:Gender:Av','Prob:tau:Gender:dominancej',
					'Prob:tau:Gender:SymCrossjk',
					'Prob:tau:Gender:ASymCrossjkDkj')
		rownames(mip) <- mip$X
		mip <- mip[reorder,]
	}
	mip$X <- simplify.labels(mip$X)
	rows <- nrow(mip)
	xx <- barplot(height=rev(mip$x), horiz=TRUE, names.arg=rev(mip$X), las=1,
		xlab="model inclusion probability (MIP)", xlim=c(-0.01,1.19))
	abline(v=c(0, 0.01, 0.02, 0.25, 0.5, 0.75, 0.98, 0.99, 1), lty=2, col="black")
	abline(v=0.5, lty=1, col="black")
#	rightlabels <- as.character(round(mip$x, digits=4))
#	axis(side=4, at=c(rows:1), labels=rightlabels, tick=FALSE, las=1)
	text(y=xx, x=1.01, label=as.character(sprintf("%.3f", round(rev(mip$x),3))), pos=4)
	try(segments(mip$x-mip$sd, rev(xx), 
		mip$x+mip$sd, rev(xx), lwd = 1.5))
	dev.off()
}

#' @title plotMipsOverMimps: Plot MIPs over reps.
#' @description Plot model inclusion probabilities over the multiple imputation matched pair/quartets.
#' 
#' @param fdir file directory
#' @param fname file name
#' @param its number of MIMPs
#' @param ... additional arguments
#' @return make MIP plot panel and save them to file directory
#' @examples
#' ## not run
#' @export
plotMipsOverMIMPs <- function(fdir, fname, its, ...){
	mip <- read.csv(file.path(fdir, fname))
	pdf(file.path(fdir, paste(file_path_sans_ext(fname), "_over_MIMPs_fig.pdf", sep="")), 
		width=8, height=12)
	mip$X <- simplify.labels(mip$X)
	rows <- nrow(mip)
	par(mfrow=c(6,3))
	par(mar=c(2.5,2.5,2.5,1.5))
	for(i in c(1:rows)){
		plot(y=as.numeric(mip[i,c(2:(its+1))]), x=c(1:its), type="l", 
			ylab="MIP", xlab="MIMP", main=mip[i,"X"],
			ylim=c(-0.1,1.1), las=1)
	}
	dev.off()
}

#' @title averageMips: Average the MIP values
#' @description Find mean (and sd) the posterior MIP values in mip table.
#' 
#' @param filenames string specifying path and MIP files of interest
#' @param fdir specifies where the new file should be placed
#' @param ... additional arguments
#' @return saves csv, returns path and name of csv with mean MIP values
#' @examples
#' ## not run
#' @export
averageMips <- function(filenames, fdir, ...){
	## Use this to combine MIP data across MIMPs
	## Input the MIP files from the different reps
	## Output a single MIP file, output the plot

	dat <- lapply(filenames, read.csv, row.names=1)
	mip.table <- dat[[1]]

	for(i in c(2:length(dat))){
		mip.table <- cbind(mip.table, dat[[i]])
	}

	names(mip.table) <- paste("rep", c(1:length(dat)), sep="_")
	new.mip.table <- mip.table
	new.mip.table$x  <- rowMeans(mip.table, na.rm=TRUE)
	new.mip.table$sd 	<- apply(X=mip.table, MARGIN=1, FUN=sd, na.rm=TRUE)

	mip.file <- file.path(fdir, "new.mip.table.csv")
	write.csv(new.mip.table, mip.file, row.names=TRUE)
	return("new.mip.table.csv")
}

#' @title Mx1HaploClass: Convert CC founder strain letter to Mx subspecies haplotype
#' @description Convert CC founder strain letter to functional Mx1 subspecies haplotype group class.
#' 
#' @param arg A through H (shorthand for founders)
#' @param ... additional arguments
#' @return returns dom, mus, cast based on substrain origin
#' @examples
#' ## not run
#' ## Note: Ferris et al., 2013, PLoS Pathogens: 
#' "three functionally distinct classes corresponding to 
#' ## domesticus (dom: A/J, C57BL6/J, 129s1/SvImJ, NOD/HiLtJ and WSB/EiJ), 
#' ## castaneus (cast: CAST/EiJ) and musculus (mus: PWK/PhJ and NZO/ShILtJ)."
#' @export
Mx1HaploClass <- function(arg, ...) {
	ret <- NULL
	for(i in 1:length(arg)){
		ret[i] <- switch(as.character(arg[i]), 
			  	A = "dom", 
			  	B = "dom", 
			  	C = "dom", 
			  	D = "dom", 
			  	'E' = "mus", 
			  	F = "cast", 
			  	G = "mus", 
			  	H = "dom", 
			  	"NA")
  	}
	return(ret)
}

#' @title imputeMissing: Impute missing phenotypes/covariates
#' @description Impute missing phenotype(s)/covariate, with selective censoring,
#'              from specific diallel categories based on posteriors.
#' @param dat original data set
#' @param phenotype the phenotype column name
#' @param covariate the covariate column name
#' @param postPred posterior predictive data for phenotype
#' @param postPred.cov posterior predictive data for covariate
#' @param effectEst posterior effect estimates for phenotype
#' @param effectEst.cov posterior effect estimates for covariate
#' @param SigmaEst.cov posterior sigma estimates for covariate
#' @param BlockEst.cov posterior block estimates for covariate
#' @param SigmaEst posterior sigma estimates for phenotype
#' @param BlockEst posterior block estimates for phenotype
#' @param FixedEst posterior fixed estimates for phenotype
#' @param toImpute categories for imputation 
#' @param toCensor categories for censoring 
#' @param colsToCopy data columns from imputation frame to copy 
#' @param numImps number of timesteps to select from total (# imputed datasets to generate)
#' @param savedir the directory the files are saved in
#' @param ... additional arguments
#' @return generates (n=timestepLen) imputed datasets, returns vector of file names
#' @examples
#' ## not run
#' @export
imputeMissing <- function(
	dat, phenotype, covariate, postPred, postPred.cov,
	effectEst, SigmaEst, BlockEst, FixedEst,
	effectEst.cov, SigmaEst.cov, BlockEst.cov,
	toImpute, toCensor, colsToCopy, numImps, 
	savedir, ...
	){

	# Impute animals from specified categories
	ncolumns <- length(names(dat))
	bat.vec <- toImpute$batname
	cat.vec <- toImpute$catname

	# Random timesteps to include; assume length of SigmaEst and SigmaEst.cov are the same
	timesteps <- floor(seq(from=1, to=length(SigmaEst), by=length(SigmaEst)/numImps))
	timesteps.cov <- floor(seq(from=1, to=length(SigmaEst.cov), by=length(SigmaEst.cov)/numImps))
	filenames <- NULL

	for(tt in 1:length(timesteps)){
		this.timestep.cov <- timesteps.cov[tt]
		this.timestep <- timesteps[tt]
		data.this <- dat
		post.sigma.cov <- as.numeric(SigmaEst.cov[this.timestep.cov])
		post.sigma <- as.numeric(SigmaEst[this.timestep])
		post.fixed <- as.numeric(FixedEst[this.timestep])
		for(i in 1:length(bat.vec)){
			imputed <- as.data.frame(matrix(data=rep(NA, times=ncolumns), nrow=1, ncol=ncolumns))
			names(imputed) <- names(dat)
			imputed[,colsToCopy] <- toImpute[i,colsToCopy]
			post.mean.cov <- as.numeric(postPred.cov[this.timestep.cov,cat.vec[i]])
			post.batch.cov <- as.numeric(BlockEst.cov[this.timestep.cov,bat.vec[i]])
			post.mean <- as.numeric(postPred[this.timestep,cat.vec[i]])
			post.batch <- as.numeric(BlockEst[this.timestep,bat.vec[i]])
			imputed[,covariate] <- 	post.mean.cov + rnorm(n=1, mean=0, sd=sqrt(post.sigma.cov)) + 
									post.batch.cov
			imputed[,phenotype] <- 	post.mean + rnorm(n=1, mean=0, sd=sqrt(post.sigma)) + 
									post.batch + post.fixed*imputed[,covariate]
			imputed[,"ID.1"] <- paste(sprintf("%04d",i), "imp", sep="")
			imputed[,"ID"]		<- paste(imputed[,"Strain"], "_", imputed[,"Sex"], imputed[,"ID.1"], sep="")
			data.this <- rbind(data.this, imputed)
		}
		for(j in 1:nrow(toCensor)){
		  dat.toCensor <- subset(data.this, dam==toCensor[j,]$dam & sire==toCensor[j,]$sire & Block==toCensor[j,]$Block & Trt==toCensor[j,trt_colname])
		  thisID <- dat.toCensor[sample(nrow(dat.toCensor),1),]$ID
      cat("Censored animal ID: ", thisID, "\n")
		  data.this <- subset(data.this, !(data.this$ID==thisID))
		}
		dir.create(file.path(savedir, sprintf("%04d",tt)))
		imputedFile <- file.path(savedir, sprintf("%04d",tt), "ImputedData.csv")
		filenames <- c(filenames, imputedFile)
		write.csv(data.this, imputedFile, row.names=FALSE)
		cat("Impute File written:", imputedFile, "\n")
	}
	return(filenames)
}


#' @title plotBatchEffects: Plot batch effects across timecourse
#' @description Plot batch effects across the timecourse.

#' @param savedirs one or more directories containing analysis files
#' @param batchIDs the IDs for the batches to be plotted
#' @param batchNames the names of the batches to be plotted
#' @param muID the ID for the posterior mean
#' @param xlim gives limits of x-axis
#' @param ... additional arguments
#' @return generates timecourse plot and effect estimate plots (ordered, unordered)
#' @examples
#' ## not run
#' @export
plotBatchEffects <- function(savedirs, batchIDs, batchNames, muID, xlim=c(85,101), ...){
	rowdim <- length(batchIDs)
	coldim <- length(savedirs)
	dat <- matrix(data=rep(NA, rowdim*(coldim+1)), nrow=rowdim, ncol=coldim+1)
	dat[,1] <- 100
	colnames(dat) <- c("D0", paste("D", c(1:coldim), sep=""))
	rownames(dat) <- batchNames
	for(i in c(1:coldim)){
		stacked.matches <- as.mcmc(read.csv(file.path(savedirs[i], "output/stackedTreatDiMIMPs.csv"),
			check.names=FALSE))
		mu <- stacked.matches[,muID]
		mx <- stacked.matches[,batchIDs]
#		rowmeans <- rowMeans(as.matrix(mx), na.rm=TRUE) # to center
#		mx <- mx - rowmeans
		mx <- mx + mu + 100
		colmeans <- colMeans(mx, na.rm=TRUE)
		dat[,i+1] <- colmeans
		colnames(mx) <- batchNames
		mx.unordered <- mx
		mx <- mx[,order(colmeans, decreasing=TRUE)]
		
		pdf(file.path(savedirs[i], "output/batch_timecourse.pdf"), width=6, height=6)
		plot(dat[i,], ylim=c(85,101), type="b", 
			ylab="pct starting wt", xlab="DPI", xlim=c(1,6.5), xaxt="n");
		for(j in c(2:rowdim)){
			lines(dat[j,], type="b");
		}
		text(x=5.5, y=86+(rank(dat[,i+1])), labels=names(dat[,i+1]), pos=4, cex=0.9)
		axis(1, at=c(1:5), labels=c(0:coldim))
		abline(h=100, lty=2)
		dev.off()

		pdf(file.path(savedirs[i], "output", 
			paste("batch_effects_on_", colnames(dat)[i+1],".pdf", sep="")), 
			width=4, height=6)
		plot.hpd(mx, ylab="", 
			xlab=paste("pct wt loss on ", as.character(colnames(dat)[i+1]), sep=""),
			xlim=xlim)
		abline(v=100, col="grey")
		dev.off()

		pdf(file.path(savedirs[i], "output", 
			paste("batch_unordered_effects_on_", colnames(dat)[i+1],".pdf", sep="")), 
			width=4, height=6)
		plot.hpd(mx.unordered, ylab="", 
			xlab=paste("pct wt loss on ", as.character(colnames(dat)[i+1]), sep=""), 
			xlim=xlim)
		abline(v=100, col="grey")
		dev.off()

	}

}

#' @title makeRotationMatrix: Make a rotation matrix
#' @description Turn an n-column design matrix into an n-1, sum to 0, design matrix
#'              while maintaining independence. 
#'              See Appendix (Lenarcic, 2012, Genetics; Crowley, 2014, Genetics)
#' @param X original design matrix
#' @param n number of original columns of design matrix
#' @param ... additional arguments
#' @return generates rotation matrix to reduce dimensions for better sampling
#' @examples
#' ## not run
#' @export
makeRotationMatrix <- function(X, n, ...)
{
 j <- ncol(X)
 k <- (-1 + sqrt(j))*(j - 1)^(-3/2)
 m <- 1/sqrt(j - 1)
 c <- (n - 2)*k + m
 M <- diag(j - 1)
 M[0 == M] <- -k
 M <- rbind(M, rep(-m, j - 1))
 return(M)
}

#' @title diallelMatrixMakerShortname
#' @description Make design matrices for diallel, based on dam and sire columns with short names.
#' @param data data frame
#' @param dam.col.name dam column name
#' @param sire.col.name sire column name
#' @param batch.col.name name of batch/random effect column
#' @param ... additional arguments
#' @return returns diallel incidence matrices
#' @examples
#' ## not run
#' @export
diallelMatrixMakerShortname <- function(data, dam.col.name, sire.col.name, 
                               batch.col.name=NULL, ...){
  dam.mat <- incidence.matrix(data[, as.character(dam.col.name)])
  sire.mat <- incidence.matrix(data[, as.character(sire.col.name)])
  add.mat <- dam.mat + sire.mat
  mat.mat <- dam.mat - sire.mat

  strains <-  c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
  colnames(dam.mat) <- paste0("dam:", strains)
  colnames(sire.mat) <- paste0("sire:", strains)
  colnames(add.mat) <- paste0("additive:", strains)
  colnames(mat.mat) <- paste0("maternal:", strains)
  # dam.mat <- diag(data$het) %*% dam.mat 
  # sire.mat <- diag(data$het) %*% sire.mat
  
  jk <- apply(cbind.data.frame(data[, as.character(dam.col.name)], data[, as.character(sire.col.name)]), 1, 
                      function(x) paste(sort(as.character(substr(x, start=1, stop=1))), collapse=""))
  jk.asymm <- apply(cbind.data.frame(data[, as.character(dam.col.name)], data[, as.character(sire.col.name)]), 1, 
                      function(x) paste(as.character(substr(x, start=1, stop=1)), collapse=""))
  asymm <- apply(cbind.data.frame(jk, jk.asymm), 1, 
                      function(x){ifelse(x[1]==x[2], -1, 1)})
  
  jk.mat <- incidence.matrix(as.factor(jk))
  drops <- c("AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH")
  data$inbred <- ifelse(data[, as.character(dam.col.name)]==data[, as.character(sire.col.name)], 1, 0)
  inbred.mat <- diag(data$inbred) %*% jk.mat
  jk.mat <- jk.mat[, !(colnames(jk.mat) %in% drops)]
  
  inbred.mat <- inbred.mat[, colnames(inbred.mat) %in% drops]
  asymm.mat <- diag(asymm) %*% jk.mat
  
  colnames(jk.mat) <- paste0("v:", colnames(jk.mat))
  colnames(asymm.mat) <- paste0("w:", colnames(asymm.mat))
  colnames(inbred.mat) <- paste0("inbreeding:", strains)

  if(is.null(batch.col.name)){
    return(list(	dam.mat=dam.mat, sire.mat=sire.mat, 
                 add.mat=add.mat, mat.mat=mat.mat, inbred.mat=inbred.mat, 
                 jk.mat=jk.mat, asymm.mat=asymm.mat))
  }else{
    if(1==length(batch.col.name)){
      batch.mat <- incidence.matrix(as.factor(data[, as.character(batch.col.name)]))
      return(list(	dam.mat=dam.mat, sire.mat=sire.mat, 
                   add.mat=add.mat, mat.mat=mat.mat, inbred.mat=inbred.mat, 
                   jk.mat=jk.mat, asymm.mat=asymm.mat, batch.mat=batch.mat))
    }else{
      stop("Not implemented for length(batch.col.name)>1")
    }
  }
  
}

#' @title diallelMatrixMaker
#' @description Make design matrices for diallel, based on dam and sire columns with long strain names.
#' @param data data frame
#' @param dam.col.name dam column name
#' @param sire.col.name sire column name
#' @param batch.col.name name of batch/random effect column
#' @param batch.1.col.name name of additional batch/random effect column
#' @param ... additional arguments
#' @return returns diallel incidence matrices
#' @examples
#' ## not run
#' @export
diallelMatrixMaker <- function(data, dam.col.name, sire.col.name, batch.col.name = NULL, batch.1.col.name = NULL,
    ...){
  dam.mat <- incidence.matrix(data[, as.character(dam.col.name)])[,c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")]
  sire.mat <- incidence.matrix(data[, as.character(sire.col.name)])[,c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")]
  add.mat <- dam.mat + sire.mat
  mat.mat <- dam.mat - sire.mat
  strains <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", 
               "WSB")
  colnames(dam.mat) <- paste0("dam:", strains)
  colnames(sire.mat) <- paste0("sire:", strains)
  colnames(add.mat) <- paste0("additive:", strains)
  colnames(mat.mat) <- paste0("maternal:", strains)
  jk <- apply(cbind.data.frame(data[, as.character(dam.col.name)], 
                               data[, as.character(sire.col.name)]), 1, function(x) paste(sort(as.character(x)), 
                                                                                          collapse = ";"))
  jk.asymm <- apply(cbind.data.frame(data[, as.character(dam.col.name)], 
                                     data[, as.character(sire.col.name)]), 1, function(x) paste(as.character(x), 
                                                                                                collapse = ";"))
  asymm <- apply(cbind.data.frame(jk, jk.asymm), 1, function(x) {
    ifelse(x[1] == x[2], -1, 1)
  })
  jk.mat <- incidence.matrix(as.factor(jk))
  drops <- c("AJ;AJ", "B6;B6", "129;129", "NOD;NOD", "NZO;NZO", 
             "CAST;CAST", "PWK;PWK", "WSB;WSB")
  data$inbred <- ifelse(data[, as.character(dam.col.name)] == 
                          data[, as.character(sire.col.name)], 1, 0)
  inbred.mat <- diag(data$inbred) %*% jk.mat
  jk.mat <- jk.mat[, !(colnames(jk.mat) %in% drops)]
  inbred.mat <- inbred.mat[, colnames(inbred.mat) %in% drops]
  asymm.mat <- diag(asymm) %*% jk.mat
  colnames(jk.mat) <- paste0("v:", colnames(jk.mat))
  colnames(asymm.mat) <- paste0("w:", colnames(asymm.mat))
  colnames(inbred.mat) <- paste0("inbreeding:", strains)
  if (is.null(batch.col.name) & is.null(batch.1.col.name)) {
    return(list(dam.mat = dam.mat, sire.mat = sire.mat, add.mat = add.mat, 
                mat.mat = mat.mat, inbred.mat = inbred.mat, jk.mat = jk.mat, 
                asymm.mat = asymm.mat))
  }
  else {
    if(1==length(batch.col.name)){
      batch.mat <- incidence.matrix(as.factor(data[, as.character(batch.col.name)]))
      batch.1.mat <- NULL
      try(batch.1.mat <- incidence.matrix(as.factor(data[, 
                                                         as.character(batch.1.col.name)])))
      return(list(dam.mat = dam.mat, sire.mat = sire.mat, 
                  add.mat = add.mat, mat.mat = mat.mat, inbred.mat = inbred.mat, 
                  jk.mat = jk.mat, asymm.mat = asymm.mat, batch.mat = batch.mat, 
                  batch.1.mat = batch.1.mat))
    }
    else{
      stop("Not implemented for length(batch.col.name)>1")
    }
  }
}


#' @title diallelMatrixMakeAndRotate
#' @description Make design matrices for diallel, rotate to n-1 space.
#' @param data data frame
#' @param dam.col.name dam column name
#' @param sire.col.name sire column name
#' @param batch.col.name name of batch/random effect column
#' @param batch.1.col.name name of additional batch/random effect column
#' @param ... additional arguments
#' @return returns diallel incidence matrices, rotated
#' @examples
#' ## not run
#' @export
diallelMatrixMakeAndRotate <- function(data, dam.col.name, sire.col.name, 
                                       batch.col.name=NULL, batch.1.col.name=NULL, ...){
  matrices <- diallelMatrixMaker(data, dam.col.name, sire.col.name, batch.col.name, batch.1.col.name)
  dam.mat <- matrices$dam.mat
  sire.mat <- matrices$sire.mat
  add.mat <- matrices$add.mat
  mat.mat <- matrices$mat.mat
  inbred.mat <- matrices$inbred.mat
  jk.mat <- matrices$jk.mat
  asymm.mat <- matrices$asymm.mat
  batch.mat <- matrices$batch.mat
  batch.1.mat <- matrices$batch.1.mat
  M.dam <- makeRotationMatrix(X = dam.mat, n = 8)
  M.sire <- makeRotationMatrix(X = sire.mat, n = 8)
  M.add <- makeRotationMatrix(X = add.mat, n = 8)
  M.mat <- makeRotationMatrix(X = mat.mat, n = 8)
  M.inbred <- makeRotationMatrix(X = inbred.mat, n = 8)
  M.jk <- makeRotationMatrix(X = jk.mat, n = 28)
  M.asymm <- makeRotationMatrix(X = asymm.mat, n = 28)
  M.batch <- makeRotationMatrix(X = batch.mat, n = ncol(batch.mat))
  M.batch.1 <- makeRotationMatrix(X = batch.1.mat, n = ncol(batch.1.mat))
  t.dam.mat <- dam.mat %*% M.dam
  t.sire.mat <- sire.mat %*% M.sire
  t.add.mat <- add.mat %*% M.add
  t.mat.mat <- mat.mat %*% M.mat
  t.inbred.mat <- inbred.mat %*% M.inbred
  t.jk.mat <- jk.mat %*% M.jk
  t.asymm.mat <- asymm.mat %*% M.asymm
  t.batch.mat <- batch.mat %*% M.batch
  t.batch.1.mat <- batch.1.mat %*% M.batch.1
  return(list(RotMat=list(	M.dam=M.dam, M.sire=M.sire,
  							M.add=M.add, M.mat=M.mat, M.inbred=M.inbred, 
  							M.jk=M.jk, M.asymm=M.asymm,
                         	M.batch=M.batch, M.batch.1=M.batch.1),
              DesignMat=list(t.dam.mat=t.dam.mat, t.sire.mat=t.sire.mat,
              			t.add.mat = t.add.mat, t.mat.mat = t.mat.mat, t.inbred.mat = t.inbred.mat,
              			t.jk.mat = t.jk.mat, t.asymm.mat = t.asymm.mat, 
              			t.batch.mat = t.batch.mat, t.batch.1.mat = t.batch.1.mat)))
}

#' @title plotHPDccColors
#' @description Plot HPD function, allowing color specification.
#' @param coda.object The CODA MCMC object
#' @param wanted subset of variable names to plot
#' @param prob.wide fractional probability for wide interval
#' @param prob.narrow fractional probability for narrow interval
#' @param xlab x-axis label
#' @param names names of variables to use on plot
#' @param type type of plot, i.e. "p"
#' @param name.margin Margin width
#' @param plt.left left plot
#' @param plt.right right plot
#' @param plt.bottom bottom plot
#' @param plt.title title of plot
#' @param ylab y-axis label
#' @param name.line name line spacing
#' @param main main title
#' @param main.line main line spacing
#' @param col color(s) for HPD
#' @param stack logical indicating whether to stack MCMC chains first
#' @param ... additional arguments
#' @return plots data, using colored HPDs
#' @examples
#' ## not run
#' @export
plotHPDccColors <- function(coda.object, wanted = varnames(coda.object), prob.wide = 0.95,
    prob.narrow = 0.5, xlab = "HPD interval", names = NULL, type = "p",
    name.margin = 6.1, plt.left = NULL, plt.right = NULL, plt.bottom = NULL,
    plt.title = NULL, ylab = "", name.line = 3.9, main = "", main.line = 2, 
    col="black", stack=TRUE, ...){
    which.wanted = ifow(is.integer(wanted), wanted, match(wanted,
        varnames(coda.object)))
    num.wanted = length(which.wanted)
    if (!exists("name.margin") || is.null(name.margin)) {
        name.margin = 6.1
    }
    if(TRUE==stack){
    	chain <- mcmc.stack(coda.object)
    }else{
    	chain <- coda.object
    }
    mu <- colMeans(chain[, which.wanted], na.rm=TRUE)
    med <- apply(coda::HPDinterval(chain, prob = 0.01)[which.wanted,
        ], 1, mean, na.rm=TRUE)
    hpd.wide <- coda::HPDinterval(chain, prob = prob.wide)[which.wanted,
        ]
    hpd.narrow <- coda::HPDinterval(chain, prob = prob.narrow)[which.wanted,
        ]
    mid.vals <- med
    if (is.null(names))
        names <- varnames(chain)[which.wanted]
#    else names <- rep(names, length.out = length(wanted))
    strsplit.names <- strsplit(sub.labels(names), split=";")
    col <- col.switch(lapply(strsplit.names, function(x){x[1]}))
    border.col <- col.switch(lapply(strsplit.names, 
    	function(x){ifelse(1==length(x), "black", x[2])}))
    # border.col <- col.switch(lapply(strsplit.names, 
    # 	function(x){if(length(x)==1){x[1]}else{x[2]}}))
    outer.col <- col.switch(lapply(strsplit.names, 
    	function(x){ifelse(1==length(x),"white","black")}))
    ypos <- plotCIccColors(med, hpd.narrow, hpd.wide, names = names,
        xlab = xlab, col.midvals = "white", pch.midvals = "|",
        type = type, name.margin = name.margin, plt.left = plt.left,
        plt.right = plt.right, plt.bottom = plt.bottom, plt.title = plt.title,
        ylab = ylab, name.line = name.line, main = main, main.line = main.line,
        col=col, border.col=border.col, outer.col=outer.col, ...)
    if ("p" == type) {
        points(mu, ypos, pch = "|")
    }
    invisible(ypos)
}

plotCIccColors <- function (midvals, narrow.intervals, wide.intervals, names = 1:length(midvals),
    add = FALSE, main = "", main.line = 2, xlab = "Estimate",
    xlab.line = 2.5, xlim = NULL, ylab = "", yaxis = TRUE, ylim = c(0,
        length(midvals)), name.line = 4, pch.midvals = 19, col = NULL, border.col = NULL,
    outer.col=NULL, col.midvals = "white", cex.labels = 1, type = "p", name.margin = 6.1,
    title.margin = 4.1, title.line = 3.5, bottom.margin = 5.1,
    bottom.line = 4.5, right.margin = 2.1, right.line = 1.5,
    mar = sides(left = name.margin, bottom = bottom.margin, top = title.margin,
        right = right.margin), mar.update = sides(), before.data = function() {
    }, plt.left = NULL, plt.right = NULL, plt.bottom = NULL,
    plt.title = NULL, ...)
{
    nvals <- length(midvals)
    col.midvals <- rep(col.midvals, length.out = nvals)
    col <- rep(col, length.out = nvals)
    border.col <- rep(border.col, length.out = nvals)
    y.pos <- (1:nvals) - 0.5
    if (!add) {
        if (is.null(xlim)) {
            xlim <- range(c(wide.intervals, narrow.intervals,
                midvals), na.rm = TRUE)
            xlim <- c(-1, 1) * diff(xlim) * 0.1 + xlim
        }
        if (name.margin == 6.1 && !is.null(plt.left) && is.numeric(plt.left) &&
            plt.left >= 0 && plt.left <= 1) {
            name.margin = 6.1 * plt.left/0.2
        }
        if (right.margin == 2.1 && !is.null(plt.right) && is.numeric(plt.right) &&
            plt.right >= 0 && plt.right <= 1) {
            right.margin = 2.1 * (1 - plt.right)/0.05
        }
        if (title.margin == 4.1 && !is.null(plt.title) && is.numeric(plt.title) &&
            plt.title >= 0 && plt.title <= 1) {
            title.margin = 4.1 * (1 - plt.title)/0.12
        }
        if (bottom.margin == 5.1 && !is.null(plt.bottom) && is.numeric(plt.bottom) &&
            plt.bottom >= 0 && plt.bottom <= 1) {
            bottom.margin = 5.1 * plt.bottom/0.16
        }
        mar <- c(bottom.margin, name.margin, title.margin, right.margin) +
            0.1
        mar = update.sides(mar, mar.update)
        oldmar <- par(mar = mar)
        on.exit(par(mar = oldmar))
        MyD = FALSE
        AT = "plot(x = xlim, y=ylim, type=\"n\", axes=FALSE, ylab=ylab, xlim=xlim, ylim=ylim, xlab=\"\", main=\"\", \n       ...); MyD = TRUE"
        try(eval(parse(text = AT)), silent = TRUE)
        if (FALSE == MyD) {
            try(plot(x = xlim, y = ylim, type = "n", axes = FALSE,
                ylab = ylab, ylim = ylim, xlim = xlim, xlab = "",
                main = "", ))
        }
        if (!is.null(main) && is.character(main) && main[1] !=
            "") {
            try(title(main = main, line = main.line, cex.main = 1.5))
        }
        title(xlab = xlab, line = xlab.line)
        axis(1)
        axis(3, line = -0.8)
        if (yaxis) {
            axis(2, at = y.pos, labels = rev(names), las = 1,
                lty = 0, hadj = 0, line = name.line, cex.axis = cex.labels)
        }
    }
    before.data()
    if ("p" == type) {
        for (i in 1:nvals) {
            pos <- nvals - i + 0.5
            lines(wide.intervals[i, ], rep(pos, 2), lwd=4.5, col=outer.col[i])
            lines(narrow.intervals[i, ], rep(pos, 2), lwd = 6.0, col=outer.col[i])
            lines(wide.intervals[i, ], rep(pos, 2), lwd = 3.0, col=border.col[i])
            lines(narrow.intervals[i, ], rep(pos, 2), lwd = 4.5, col=border.col[i])
            lines(wide.intervals[i, ], rep(pos, 2), lwd=1.0, col=col[i])
            lines(narrow.intervals[i, ], rep(pos, 2), lwd = 2.5, col=col[i])
            points(midvals[i], pos, pch = pch.midvals, col = col.midvals[i])
        }
    }
    invisible(rev(y.pos))
}
