F_translate_headers <- function(header, standard = c("MARKER", "EFFECTALLELE", "OTHERALLELE", "BETA", "SE", "Z", "P", "LP", "EAF", "N", "NSTUDIES", "QHET", "QHETP", "I2HET", "DIRECTIONS"), alternative) {
  capitalized <- toupper(header)
  for(forI in 1:length(standard)) {
    column_no <- F_find_header(standard[forI], alternative, capitalized)
    if (column_no > 0) { 
      header[column_no] <- standard[forI]
    }
  }	
  return(list(header_N = length(header), header_h = header))
}

F_find_header <- function(std_name, alt_names, results_header) {
  relevant_alts <- alt_names[which(alt_names[, 1] ==  std_name), 2]
  if (sum(relevant_alts %in% results_header) ==  1) {
    for (ih in 1 : length(results_header)) {
      if (results_header[ih] %in% relevant_alts) { return(ih) }
    }
  } else { return(0) }
}	# The first line collects the alternative names for the default header
#	name. These are then matched to the headers of the current
#	dataset: if one matches, the number of the dataset column
#	bearing that name is passed on.

complement <- function (allele) {
  compl_allele <- allele
  compl_allele[allele == "A"] <- "T"
  compl_allele[allele == "C"] <- "G"
  compl_allele[allele == "G"] <- "C"
  compl_allele[allele == "T"] <- "A"
  return(compl_allele)
}  

GClambda <- function(P) {
  P<-P[!is.na(P)]
  lambda <- median(qchisq(P, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
  return(lambda)
}

subtract.directions <- function (meta_cohort) {
  meta_directions <- unlist(strsplit(meta_cohort$DIRECTIONS,split=""))
  meta_directions[meta_directions=="-"] <- -1
  meta_directions[meta_directions=="+"] <-  1
  meta_directions[meta_directions=="?"] <-  0
  meta_directions <- matrix(as.integer(meta_directions),nrow=nrow(meta_cohort),byrow=T)
  cohort_directions <- ifelse(meta_cohort$BETA.y < 0, -1, ifelse(meta_cohort$BETA.y > 0, 1, 0))
  cohort_directions[is.na(cohort_directions)] <- 0
  
  ic<-1
  maxr<-0
  for (i in 1:ncol(meta_directions)) {
    if (cor(cohort_directions,meta_directions[,i]) > maxr) {
	  maxr<-cor(cohort_directions,meta_directions[,i])
	  ic <- i
	}
  }
  meta_directions.adj <- meta_directions[, -ic]
  meta_directions.adj[meta_directions.adj==-1] <- "-"
  meta_directions.adj[meta_directions.adj==1] <-  "+"
  meta_directions.adj[meta_directions.adj==0] <-  "?"
  if (ncol(meta_directions.adj)>1) {
    for (i in 2:ncol(meta_directions.adj)) {
      meta_directions.adj[,1] <- paste0(meta_directions.adj[,1], meta_directions.adj[,i])
	}
  }
  return(meta_directions.adj[,1])
}

meta.min.1 <- function(meta_cohort, metamethod = "FIV", lambda.meta = 1, lambda.cohort = 1, gc_meta = TRUE, calculate_lambda.meta = TRUE, calculate_lambdas.cohort = TRUE, logfile = "MetaSubtract.log") {
  cat(paste(" - Pre-set genomic control lambda of the imported meta-results = ",format(lambda.meta, digits=4)," \n", sep=""))
  write.table(paste(" - Pre-set genomic control lambda of the imported meta-results = ",format(lambda.meta, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
  if (calculate_lambdas.cohort) {
    lambda.cohort <- GClambda(meta_cohort$P.y)
    cat(paste(" - Genomic control lambda calculated from the cohort results = ",format(lambda.cohort, digits=4)," \n", sep=""))
    write.table(paste(" - Genomic control lambda calculated from the cohort results = ",format(lambda.cohort, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
  } else {
    cat(paste(" - Pre-set genomic control lambda of the cohort results = ",format(lambda.cohort, digits=4)," \n", sep=""))
    write.table(paste(" - Pre-set genomic control lambda of the cohort results = ",format(lambda.cohort, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
  }

  if (toupper(metamethod) == "FIV") { #fixed effect inverse variance
    cat(" - Assuming a fixed-effects inverse-variance weighted meta-analysis method \n")
    write.table(" - Assuming a fixed-effects inverse-variance weighted meta-analysis method",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
    adj <- NULL
    if ("BETA.x" %in% names(meta_cohort) & "SE.x" %in% names(meta_cohort)) {
      if ("BETA.y" %in% names(meta_cohort) & "SE.y" %in% names(meta_cohort)) {
        meta_cohort$BETA.adj <- meta_cohort$BETA.x
        meta_cohort$SE.adj <- meta_cohort$SE.x/sqrt(lambda.meta)
        meta_cohort$P.adj <- 2*pnorm(-abs(meta_cohort$BETA.adj)/meta_cohort$SE.adj)
        meta_cohort$LP.adj <- -log10(meta_cohort$P.adj)
        meta_cohort$EAF.adj <- meta_cohort$EAF.x
        meta_cohort$N.adj <- meta_cohort$N.x
        meta_cohort$NSTUDIES.adj <- meta_cohort$NSTUDIES
		meta_cohort$DIRECTIONS.adj <- subtract.directions(meta_cohort)
        meta_cohort$SE.y <- meta_cohort$SE.y/sqrt(lambda.cohort)
        
        subset <- which(!is.na(meta_cohort$BETA.y) & !is.na(meta_cohort$SE.y))

        meta.weight <- 1/(meta_cohort$SE.adj[subset])^2
        cohort.weight <- 1/(meta_cohort$SE.y[subset])^2
        meta_cohort$BETA.adj[subset] <- (meta.weight*meta_cohort$BETA.adj[subset]-cohort.weight*meta_cohort$BETA.y[subset])/(meta.weight-cohort.weight)
        meta_cohort$SE.adj[subset] <- 1/sqrt(meta.weight-cohort.weight)
        meta_cohort$P.adj[subset] <- 2*pnorm(-abs(meta_cohort$BETA.adj[subset])/meta_cohort$SE.adj[subset])
        if (gc_meta) {
          if (calculate_lambda.meta) {
            lambda.meta.adj <- GClambda(meta_cohort$P.adj)
            cat(paste(" - Genomic control lambda calculated from the corrected meta-results = ",format(lambda.meta.adj, digits=4)," \n", sep=""))
            write.table(paste(" - Genomic control lambda calculated from the corrected meta-results = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
    	} else {
            lambda.meta.adj <- lambda.meta
            cat(paste(" - Pre-set genomic control lambda applied to the corrected meta-results = ",format(lambda.meta.adj, digits=4)," \n", sep=""))
            write.table(paste(" - Pre-set genomic control lambda applied to the corrected meta-results = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
          }
          meta_cohort$P.adj[subset] <- 2*pnorm(-abs(meta_cohort$BETA.adj[subset])/(meta_cohort$SE.adj[subset]*sqrt(lambda.meta.adj)))
        }
        meta_cohort$LP.adj[subset] <- -log10(meta_cohort$P.adj[subset])
        meta_cohort$EAF.adj[subset] <- (2*meta_cohort$N.adj[subset]*meta_cohort$EAF.adj[subset]-2*meta_cohort$N.y[subset]*meta_cohort$EAF.y[subset])/(2*meta_cohort$N.adj[subset]-2*meta_cohort$N.y[subset])
        meta_cohort$N.adj[subset] <- meta_cohort$N.adj[subset]-meta_cohort$N.y[subset]
        meta_cohort$NSTUDIES.adj[subset] <- meta_cohort$NSTUDIES.adj[subset]-1

        meta_cohort$QHET.adj[subset] <- meta_cohort$QHET[subset]*meta_cohort$SE.adj[subset]^2/meta_cohort$SE.x[subset]^2-(meta_cohort$BETA.adj[subset]-meta_cohort$BETA.x[subset])^2-1/meta_cohort$SE.y[subset]^2*(meta_cohort$BETA.y[subset]-meta_cohort$BETA.x[subset])^2*meta_cohort$SE.adj[subset]^2
        meta_cohort$QHETP.adj[subset] <- 1-pchisq(meta_cohort$QHET.adj[subset],meta_cohort$NSTUDIES.adj[subset]-1)
        meta_cohort$I2HET.adj[subset] <- 100*(meta_cohort$QHET.adj[subset]-meta_cohort$NSTUDIES.adj[subset]+1)/meta_cohort$QHET.adj[subset]

        adj <- meta_cohort[, c("BETA.adj", "SE.adj", "P.adj", "LP.adj", "EAF.adj", "N.adj", "NSTUDIES.adj", "DIRECTIONS.adj","QHET.adj", "QHETP.adj", "I2HET.adj")]
      } else { 
        write.table("ERROR: 'BETA' or 'SE' not found in cohort results",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
        stop("ERROR: 'BETA' or 'SE' not found in cohort results")
      }
    } else {
      write.table("ERROR: 'BETA' or 'SE' not found in meta results",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
      stop("ERROR: 'BETA' or 'SE' not found in meta results")
    }
  } else if (toupper(metamethod) == "FSSW") { #fixed effect sample size weighted
    cat(" - Assuming a fixed-effects sample size weighted meta-analysis method \n")
    write.table(" - Assuming a fixed-effects sample size weighted meta-analysis method",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
    adj <- NULL
    if ("BETA.x" %in% names(meta_cohort) & "N.x" %in% names(meta_cohort)) {
      if ("BETA.y" %in% names(meta_cohort) & "N.y" %in% names(meta_cohort)) {
        meta_cohort$BETA.adj <- meta_cohort$BETA.x
        meta_cohort$SE.adj <- meta_cohort$SE.x/sqrt(lambda.meta)
        meta_cohort$P.adj <- 2*pnorm(-abs(meta_cohort$BETA.adj)/meta_cohort$SE.adj)
        meta_cohort$LP.adj <- -log10(meta_cohort$P.adj)
        meta_cohort$EAF.adj <- meta_cohort$EAF.x
        meta_cohort$N.adj <- meta_cohort$N.x
        meta_cohort$NSTUDIES.adj <- meta_cohort$NSTUDIES
		meta_cohort$DIRECTIONS.adj <- subtract.directions(meta_cohort)

        meta_cohort$SE.y <- meta_cohort$SE.y/sqrt(lambda.cohort)

        subset <- which(!is.na(meta_cohort$BETA.y) & !is.na(meta_cohort$SE.y))
        
        meta.weight <- meta_cohort$N.adj[subset]
        cohort.weight <- meta_cohort$N.y[subset]
        meta_cohort$BETA.adj[subset] <- (meta.weight*meta_cohort$BETA.adj[subset]-cohort.weight*meta_cohort$BETA.y[subset])/(meta.weight-cohort.weight)
        meta_cohort$SE.adj[subset] <- sqrt(1/(1/(meta_cohort$SE.adj[subset])^2-1/(meta_cohort$SE.y[subset])^2))
        meta_cohort$P.adj[subset] <- 2*pnorm(-abs(meta_cohort$BETA.adj[subset])/meta_cohort$SE.adj[subset])
        if (gc_meta) {
          if (calculate_lambda.meta) {
            lambda.meta.adj <- GClambda(meta_cohort$P.adj)
            cat(paste(" - Genomic control lambda calculated from the corrected meta-results = ",format(lambda.meta.adj, digits=4)," \n", sep=""))
            write.table(paste(" - Genomic control lambda calculated from the corrected meta-results = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
          } else {
            lambda.meta.adj <- lambda.meta
            cat(paste(" - Pre-set genomic control lambda applied to the corrected meta-results = ",format(lambda.meta.adj, digits=4)," \n", sep=""))
            write.table(paste(" - Pre-set genomic control lambda applied to the corrected meta-results = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
          }
          meta_cohort$P.adj[subset] <- 2*pnorm(-abs(meta_cohort$BETA.adj[subset])/(meta_cohort$SE.adj[subset]*sqrt(lambda.meta.adj)))
        }
        meta_cohort$LP.adj[subset] <- -log10(meta_cohort$P.adj[subset])
        meta_cohort$EAF.adj[subset] <- (2*meta_cohort$N.adj[subset]*meta_cohort$EAF.adj[subset]-2*meta_cohort$N.y[subset]*meta_cohort$EAF.y[subset])/(2*meta_cohort$N.adj[subset]-2*meta_cohort$N.y[subset])
        meta_cohort$N.adj[subset] <- meta_cohort$N.adj[subset]-meta_cohort$N.y[subset]
        meta_cohort$NSTUDIES.adj[subset] <- meta_cohort$NSTUDIES.adj[subset]-1

        meta_cohort$QHET.adj[subset] <- meta_cohort$QHET[subset]*meta_cohort$SE.adj[subset]^2/meta_cohort$SE.x[subset]^2-(meta_cohort$BETA.adj[subset]-meta_cohort$BETA.x[subset])^2-1/meta_cohort$SE.y[subset]^2*(meta_cohort$BETA.y[subset]-meta_cohort$BETA.x[subset])^2*meta_cohort$SE.adj[subset]^2
        meta_cohort$QHETP.adj[subset] <- 1-pchisq(meta_cohort$QHET.adj[subset],meta_cohort$NSTUDIES.adj[subset]-1)
        meta_cohort$I2HET.adj[subset] <- 100*(meta_cohort$QHET.adj[subset]-meta_cohort$NSTUDIES.adj[subset]+1)/meta_cohort$QHET.adj[subset]
  
        adj <- meta_cohort[, c("BETA.adj", "SE.adj", "P.adj", "LP.adj", "EAF.adj", "N.adj", "NSTUDIES.adj", "DIRECTIONS.adj","QHET.adj", "QHETP.adj", "I2HET.adj")]
      } else { 
        write.table("ERROR: 'BETA' or 'N' not found in cohort results",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
        stop("ERROR: 'BETA' or 'N' not found in cohort results")
      }
    } else {
      write.table("ERROR: 'BETA' or 'N' not found in meta results",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
      stop("ERROR: 'BETA' or 'N' not found in meta results")
    }
  } else if (toupper(metamethod) == "FSZ") { #fixed effects sample size z-score method
    cat(" - Assuming a fixed-effects sqrt sample size z-score weighted meta-analysis method \n")
    write.table(" - Assuming a fixed-effects sqrt sample size z-score weighted meta-analysis method",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
    adj <- NULL
    if (!("Z.x" %in% names(meta_cohort))) {
      if ("BETA.x" %in% names(meta_cohort) & "SE.x" %in% names(meta_cohort)) {
        meta_cohort$Z.x <- meta_cohort$BETA.x/meta_cohort$SE.x
      }
    }
    if (!("Z.y" %in% names(meta_cohort))) {
      if ("BETA.y" %in% names(meta_cohort) & "SE.y" %in% names(meta_cohort)) {
        meta_cohort$Z.y <- meta_cohort$BETA.y/meta_cohort$SE.y
      }
    }
    
    if ("Z.x" %in% names(meta_cohort) & "N.x" %in% names(meta_cohort)) {
      if ("Z.y" %in% names(meta_cohort) & "N.y" %in% names(meta_cohort)) {
        meta_cohort$Z.adj <- meta_cohort$Z.x*sqrt(lambda.meta)
        meta_cohort$SE.adj <- meta_cohort$SE.x/sqrt(lambda.meta)
        meta_cohort$P.adj <- 2*pnorm(-abs(meta_cohort$Z.adj))
        meta_cohort$LP.adj <- -log10(meta_cohort$P.adj)
        meta_cohort$EAF.adj <- meta_cohort$EAF.x
        meta_cohort$N.adj <- meta_cohort$N.x
        meta_cohort$NSTUDIES.adj <- meta_cohort$NSTUDIES
		meta_cohort$DIRECTIONS.adj <- subtract.directions(meta_cohort)

        meta_cohort$Z.y <- meta_cohort$Z.y*sqrt(lambda.cohort)
        meta_cohort$SE.y <- meta_cohort$SE.y/sqrt(lambda.cohort)
        
        subset <- which(!is.na(meta_cohort$Z.y) & !is.na(meta_cohort$N.y))
        
        meta_cohort$Z.adj[subset] <- (sqrt(meta_cohort$N.adj[subset])*meta_cohort$Z.adj[subset]-sqrt(meta_cohort$N.y[subset])*meta_cohort$Z.y[subset]) / sqrt((meta_cohort$N.adj[subset]-meta_cohort$N.y[subset]))
        meta_cohort$P.adj[subset] <- 2*pnorm(-abs(meta_cohort$Z.adj[subset]))
        if (gc_meta) {
          if (calculate_lambda.meta) {
            lambda.meta.adj <- GClambda(meta_cohort$P.adj)
            cat(paste(" - Genomic control lambda calculated from the corrected meta-results = ",format(lambda.meta.adj, digits=4)," \n", sep=""))
            write.table(paste(" - Genomic control lambda calculated from the corrected meta-results = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
          } else {
            lambda.meta.adj <- lambda.meta
            cat(paste(" - Pre-set genomic control lambda applied to the corrected meta-results = ",format(lambda.meta.adj, digits=4)," \n", sep=""))
            write.table(paste(" - Pre-set genomic control lambda applied to the corrected meta-results = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
          }
          meta_cohort$P.adj[subset] <- 2*pnorm(-abs(meta_cohort$Z.adj[subset]/sqrt(lambda.meta.adj)))
        }
        meta_cohort$LP.adj[subset] <- -log10(meta_cohort$P.adj[subset])
        meta_cohort$EAF.adj[subset] <- (2*meta_cohort$N.adj[subset]*meta_cohort$EAF.adj[subset]-2*meta_cohort$N.y[subset]*meta_cohort$EAF.y[subset])/(2*meta_cohort$N.adj[subset]-2*meta_cohort$N.y[subset])
        meta_cohort$N.adj[subset] <- meta_cohort$N.adj[subset]-meta_cohort$N.y[subset]
        meta_cohort$NSTUDIES.adj[subset] <- meta_cohort$NSTUDIES.adj[subset]-1

        adj <- meta_cohort[, c("Z.adj", "P.adj", "LP.adj", "EAF.adj", "N.adj", "NSTUDIES.adj", "DIRECTIONS.adj")]
      } else { 
        write.table("ERROR: 'Z' or 'N' not found in cohort results",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
        stop("ERROR: 'Z' or 'N' not found in cohort results")
      }
    } else {
      write.table("ERROR: 'Z' or 'N' not found in meta results",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
      stop("ERROR: 'Z' or 'N' not found in meta results")
    }
  }
  lambda.meta.adj <- GClambda(adj$P.adj)
  cat(paste(" - Genomic control lambda calculated of the final corrected meta-results = ",format(lambda.meta.adj, digits=4)," \n", sep=""))  
  write.table(paste(" - Genomic control lambda calculated of the final corrected meta-results = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
  
  return(adj)
}

meta.subtract <- function(metafile, cohortfiles, metamethod = "FIV", lambda.meta = 1, lambdas.cohort = 1, gc_meta = TRUE, calculate_lambda.meta = TRUE, calculate_lambdas.cohort = TRUE, alternative = "alternative_headers.txt", logfile = "MetaSubtract.log", dir = getwd()) {
  print(date())
  write.table(date(),logfile,col.names=F,row.names=F,quote=F,append=F)
  if (!is.character(dir)) {
    warning(paste("Invalid format of dir '",dir,"'"))	
    write.table(paste("Invalid format of dir '",dir,"'"),logfile,col.names=F,row.names=F,quote=F,append=T)
  } else if (dir.exists(dir)) {
    setwd(dir)
  } else {
    warning(paste("Directory '",dir,"' does not exist. Current directory is ",getwd(), sep=""))
	write.table(paste("Directory '",dir,"' does not exist. Current directory is ",getwd(), sep=""),logfile,col.names=F,row.names=F,quote=F,append=T)
  }
  
  if (length(lambdas.cohort)<length(cohortfiles) & !calculate_lambdas.cohort) {
    if (!(length(lambdas.cohort)==1 & lambdas.cohort==1)) {
      warning(paste("Less GC lambdas have been given than cohort files.\n The last",length(cohortfiles)-length(lambdas.cohort),"cohort files will be assumed to have GC lambda=1."))
	  write.table(paste("Less GC lambdas have been given than cohort files.\n The last",length(cohortfiles)-length(lambdas.cohort),"cohort files will be assumed to have GC lambda=1."),logfile,col.names=F,row.names=F,quote=F,append=T)
    }
    lambdas.cohort <- c(lambdas.cohort,rep(1,length(cohortfiles)-length(lambdas.cohort)))
  }  

  if (any(duplicated(cohortfiles))) {
    warning("Some cohort names appear more than once in the list.")
    write.table("Some cohort names appear more than once in the list.",logfile,col.names=F,row.names=F,quote=F,append=T)
    lambdas.cohort <- lambdas.cohort[!duplicated(cohortfiles)]
    cohortfiles <- cohortfiles[!duplicated(cohortfiles)]
  }

  if (!file.exists(alternative)) { 
    if (!file.exists(system.file("extdata",alternative,package="MetaSubtract"))) { 
      write.table(paste("Alternative header file '",alternative,"' does not exist.", sep=""),logfile,col.names=F,row.names=F,quote=F,append=T)
      stop(paste("Alternative header file '",alternative,"' does not exist.", sep=""))
    } else {
      alternative<-system.file("extdata",alternative,package="MetaSubtract")
    }
  }
  header_translations <- read.table(alternative)
  
  if (!(metamethod %in% c("FIV", "FSSW", "FSZ"))) {
    write.table("ERROR: Please specify one of the following meta-analysis methods: FIV, FSSW, or FSZ",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
    stop("ERROR: Please specify one of the following meta-analysis methods: FIV, FSSW, or FSZ")
  }

  if (!file.exists(metafile)) { 
    write.table(paste("'",metafile,"' does not exist", sep=""),logfile,col.names=F,row.names=F,quote=F,append=T)
    stop(paste("'",metafile,"' does not exist", sep=""))
  }
  cat(paste("Reading",metafile,"\n"))
  write.table(paste("Reading",metafile),logfile,col.names=F,row.names=F,quote=F,append=T)
  meta <- read.table(metafile, header = TRUE, as.is = TRUE)
  header_meta <- names(meta)
  headersinfometa <- F_translate_headers(header = header_meta, alternative = header_translations)
  names(meta) <- headersinfometa$header_h
  meta$EFFECTALLELE <- toupper(meta$EFFECTALLELE)
  meta$OTHERALLELE <- toupper(meta$OTHERALLELE)
  markercol<-which("MARKER" %in% names(meta))
  meta$order<-c(1:nrow(meta))
  
  for (cohortfile in cohortfiles) {
    if (!file.exists(cohortfile)) { 
      warning(paste("'",cohortfile,"' does not exist. Meta results will not be corrected.", sep=""))
      write.table(paste("'",cohortfile,"' does not exist. Meta results will not be corrected.", sep=""),logfile,col.names=F,row.names=F,quote=F,append=T)
    } else {
      cohort <- read.table(cohortfile, header = TRUE, as.is = TRUE)
      cat(paste("Subtracting effects of",cohortfile," \n"))
      write.table(paste("Subtracting effects of",cohortfile),logfile,col.names=F,row.names=F,quote=F,append=T)
      headersinfo <- F_translate_headers(header = names(cohort), alternative = header_translations)
      names(cohort) <- headersinfo$header_h
      cohort$EFFECTALLELE <- toupper(cohort$EFFECTALLELE)
      cohort$OTHERALLELE <- toupper(cohort$OTHERALLELE)
    
      meta_cohort <- merge(meta, cohort, by = "MARKER", all.x = TRUE, all.y = FALSE)
    
      # flip alleles in cohort so that they match with those in the meta and adjust BETA and EAF
      subset = !is.na(meta_cohort$EFFECTALLELE.y) & meta_cohort$EFFECTALLELE.x != meta_cohort$EFFECTALLELE.y
      meta_cohort[subset, c("BETA.y")] <- -meta_cohort[subset, c("BETA.y")]
      meta_cohort[subset, c("EAF.y")] <- 1-meta_cohort[subset, c("EAF.y")]
      meta_cohort[subset, c("EFFECTALLELE.y", "OTHERALLELE.y")] <- meta_cohort[subset, c("OTHERALLELE.y", "EFFECTALLELE.y")]  
    
      # strand switch of alleles
      subset = !is.na(meta_cohort$EFFECTALLELE.y) & meta_cohort$EFFECTALLELE.x == complement(meta_cohort$EFFECTALLELE.y) & meta_cohort$EFFECTALLELE.x != complement(meta_cohort$OTHERALLELE.y)
      if (sum(subset>0)) { meta_cohort[subset, c("EFFECTALLELE.y", "OTHERALLELE.y")] <- cbind(complement(meta_cohort[subset, "EFFECTALLELE.y"]), complement(meta_cohort[subset, "OTHERALLELE.y"])) }
    
      # strand switch of alleles and flip alleles in cohort so that they match with thoSE in the meta and adjust BETA and EAF
      subset = !is.na(meta_cohort$EFFECTALLELE.y) & meta_cohort$EFFECTALLELE.x == complement(meta_cohort$OTHERALLELE.y) & meta_cohort$EFFECTALLELE.x != complement(meta_cohort$OTHERALLELE.y)
      meta_cohort[subset, c("BETA.y")] <-  -meta_cohort[subset, c("BETA.y")]
      meta_cohort[subset, c("EAF.y")] <- 1-meta_cohort[subset, c("EAF.y")]
      if (sum(subset>0)) { meta_cohort[subset, c("EFFECTALLELE.y", "OTHERALLELE.y")] <- cbind(complement(meta_cohort[subset, "OTHERALLELE.y"]), complement(meta_cohort[subset, "EFFECTALLELE.y"])) }
      
      # added by L.Corbin 7th August 2017
      # corrects error which was leading to assignment of adjusted betas to wrong SNPs
      meta_cohort <- meta_cohort[order(meta_cohort$order),]
      # end of edit

      adj <- meta.min.1(meta_cohort, metamethod, lambda.meta, lambdas.cohort[which(cohortfile %in% cohortfiles)], gc_meta, calculate_lambda.meta, calculate_lambdas.cohort, logfile=logfile)
      s <- " Columns that have been corrected are: "
      for (i in names(adj)) {
        orig <- sub(".adj","",i)
        if (orig %in% names(meta)) {
		  index <- which(orig==names(meta))
          meta[,orig] <-  adj[,i]
          s<-paste(s,header_meta[index], sep=", ")
        }  
      }
      s <- paste(substr(s,1,38), substr(s,42,nchar(s)))
      cat(paste(s," \n"))
      write.table(s,logfile,col.names=F,row.names=F,quote=F,append=T)

	  names(adj)<-gsub(".adj","",names(adj))
      s <- " No corrections have been made for: "
      for (i in (1:length(header_meta))) {
        orig <- header_meta[i]
		conv <- names(meta)[i]
        if (!(orig %in% names(adj)) & !(conv %in% names(adj))) {
          s<-paste(s,header_meta[i], sep=", ")
        }  
      }
      s <- paste(substr(s,1,35), substr(s,39,nchar(s)))
      cat(paste(s," \n \n"))
      write.table(paste(s," \n"),logfile,col.names=F,row.names=F,quote=F,append=T)
    }
  }
  meta <- meta[order(meta$order),]
  if (markercol==1) { 
    meta.adj <- meta[,c(1:headersinfometa$header_N)]
  } else {
    meta.adj <- meta[,c(2:markercol,1,(markercol+1):headersinfometa$header_N)]
  }
  names(meta.adj) <- header_meta
  return(meta.adj)
}

