F_translate_headers <- function(header, standard = c("MARKER", "EFFECTALLELE", "OTHERALLELE", "BETA", "SE", "Z", "P", "LP", "EAF", "N", "NSTUDIES", "QHET", "QHETP", "I2HET", "DIRECTIONS"), alternative) {
  capitalized <- toupper(header)
  for (forI in 1:length(standard)) {
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

subtract.directions <- function (meta_cohort, logfile = "MetaSubtract.log") {
  if (is.factor(meta_cohort$DIRECTIONS)) meta_cohort$DIRECTIONS <- as.character(levels(meta_cohort$DIRECTIONS))[meta_cohort$DIRECTIONS]
  meta_directions <- unlist(strsplit(meta_cohort$DIRECTIONS,split=""))
  meta_directions[meta_directions=="-"] <- -1
  meta_directions[meta_directions=="+"] <-  1
  meta_directions[meta_directions=="?"] <-  0
  meta_directions <- matrix(as.integer(meta_directions),nrow=nrow(meta_cohort),byrow=T)
  cohort_directions <- ifelse(meta_cohort$BETA.cohort < 0, -1, ifelse(meta_cohort$BETA.cohort > 0, 1, 0))
  cohort_directions[is.na(cohort_directions)] <- 0
  
  ic<-1
  maxr<-0
  for (i in 1:ncol(meta_directions)) {
    concordance_table <- prop.table(table(cohort_directions[cohort_directions!=0]==meta_directions[cohort_directions!=0,i]))
    if (length(concordance_table)==1) {
	  maxr <- 1
	  ic <- i
	} else if (prop.table(table(cohort_directions[cohort_directions!=0]==meta_directions[cohort_directions!=0,i]))[2] > maxr) {
	  maxr<-prop.table(table(cohort_directions[cohort_directions!=0]==meta_directions[cohort_directions!=0,i]))[2]
	  ic <- i
	}
  }
  if (maxr>0.99) { 
	s<-paste(" - Concordance of cohort's directions with those from column ",ic," of directions in meta summary statistics is ",signif(maxr,2),".",sep="")
    write.table(s,logfile,col.names=F,row.names=F,quote=F,append=T)
    meta_directions.adj <- meta_directions[, -ic]
    meta_directions.adj[meta_directions.adj==-1] <- "-"
    meta_directions.adj[meta_directions.adj==1] <-  "+"
    meta_directions.adj[meta_directions.adj==0] <-  "?"
    if (ncol(meta_directions.adj)>1) {
      for (i in 2:ncol(meta_directions.adj)) {
        meta_directions.adj[,1] <- paste0(meta_directions.adj[,1], meta_directions.adj[,i])
	  }
    }
	adj <- TRUE
  } else {
	meta_directions.adj <- data.frame(v1=meta_cohort$DIRECTIONS,v2=NA)
	adj <- FALSE
	s<-paste(" - WARNING: No column of directions in meta summary statistics matches well with those of the cohort (max concordance=",signif(maxr,2),"). Directions are not corrected and input data should be checked (quality filters might have been applied before meta-analysis).",sep="")
    write.table(s,logfile,col.names=F,row.names=F,quote=F,append=T)
    cat(paste(s,"\n"))
  }
  return(list(meta_directions.adj[,1],adj=adj))
}

meta.min.1 <- function(meta_cohort, metamethod = "FIV", lambda.meta = 1, lambda.cohort = 1, gc_meta = TRUE, calculate_lambda.meta = TRUE, calculate_lambdas.cohort = TRUE, logfile = "MetaSubtract.log", namesmeta, namescohort, index = 1, dir = tempdir()) {
  if (index==1) {
    cat(paste(" - Pre-set genomic control lambda of the imported meta summary statistics = ",format(lambda.meta, digits=4)," \n", sep=""))
    write.table(paste(" - Pre-set genomic control lambda of the imported meta summary statistics = ",format(lambda.meta, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
  } else {
    if ("P.meta" %in% names(meta_cohort)) {
	  lambda.meta <- GClambda(meta_cohort$P.meta)
	} else if ("P" %in% namesmeta) {
	  lambda.meta <- GClambda(meta_cohort$P)
    }	
    cat(paste(" - Computed genomic control lambda of the imported meta summary statistics = ",format(lambda.meta, digits=4)," \n", sep=""))
    write.table(paste(" - Computed genomic control lambda of the imported meta summary statistics = ",format(lambda.meta, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
  }
  
  if (calculate_lambdas.cohort) {
    if (!("P" %in% namescohort) & "BETA" %in% namescohort & "SE" %in% namescohort) meta_cohort$P.cohort <- 2*pnorm(-abs(meta_cohort$BETA.cohort)/meta_cohort$SE.cohort)
    if ("P.cohort" %in% names(meta_cohort)) {
	  lambda.cohort <- GClambda(meta_cohort$P.cohort)
      cat(paste(" - Genomic control lambda calculated from the cohort results = ",format(lambda.cohort, digits=4)," \n", sep=""))
      write.table(paste(" - Genomic control lambda calculated from the cohort results = ",format(lambda.cohort, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
	} else {
	  s <- " - No P-values are present or can be determined in the cohort results file. No genomic control correction is applied to the cohort results."
	  cat(paste(s,"\n"))
      write.table(s,logfile,col.names=F,row.names=F,quote=F,append=TRUE)
	}
  } else {
    cat(paste(" - Pre-set genomic control lambda of the cohort results = ",format(lambda.cohort, digits=4)," \n", sep=""))
    write.table(paste(" - Pre-set genomic control lambda of the cohort results = ",format(lambda.cohort, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
  }

  if (toupper(metamethod) == "FIV") { #fixed effect inverse variance
    cat(" - Assuming a fixed-effects inverse-variance weighted meta-analysis method \n")
    write.table(" - Assuming a fixed-effects inverse-variance weighted meta-analysis method",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
    adj <- NULL
    if ("BETA.meta" %in% names(meta_cohort) & "SE.meta" %in% names(meta_cohort)) {
      if ("BETA.cohort" %in% names(meta_cohort) & "SE.cohort" %in% names(meta_cohort)) {
        meta_cohort$BETA.adj <- meta_cohort$BETA.meta
        meta_cohort$SE.adj <- meta_cohort$SE.meta/sqrt(lambda.meta)
        meta_cohort$P.adj <- 2*pnorm(-abs(meta_cohort$BETA.adj)/meta_cohort$SE.adj)
		if ("LP" %in% names(meta_cohort)) meta_cohort$LP.adj <- -log10(meta_cohort$P.adj)
        if ("EAF.meta" %in% names(meta_cohort)) { meta_cohort$EAF.adj <- meta_cohort$EAF.meta } else  { meta_cohort$EAF.adj <- NA }
        if ("N.meta" %in% names(meta_cohort)) { meta_cohort$N.adj <- meta_cohort$N.meta } else { meta_cohort$N.adj <- NA }
        if ("NSTUDIES" %in% names(meta_cohort)) { meta_cohort$NSTUDIES.adj <- meta_cohort$NSTUDIES } else { meta_cohort$NSTUDIES <- NA }
        meta_cohort$DIRECTIONS.adj <- NA
		if ("DIRECTIONS" %in% names(meta_cohort)) { 
		  adj_directions <- subtract.directions(meta_cohort, logfile)
          if (adj_directions[[2]]) meta_cohort$DIRECTIONS.adj <- adj_directions[[1]]
        }
        meta_cohort$SE.cohort <- meta_cohort$SE.cohort/sqrt(lambda.cohort)
			
        subset <- which(!is.na(meta_cohort$BETA.cohort) & !is.na(meta_cohort$SE.cohort))

        meta.weight <- 1/(meta_cohort$SE.adj[subset])^2
        cohort.weight <- 1/(meta_cohort$SE.cohort[subset])^2
        meta_cohort$BETA.adj[subset] <- (meta.weight*meta_cohort$BETA.adj[subset]-cohort.weight*meta_cohort$BETA.cohort[subset])/(meta.weight-cohort.weight)
        meta_cohort$SE.adj[subset] <- 1/sqrt(meta.weight-cohort.weight)
        meta_cohort$P.adj[subset] <- 2*pnorm(-abs(meta_cohort$BETA.adj[subset])/meta_cohort$SE.adj[subset])
        if (gc_meta) {
          if (calculate_lambda.meta) {
            lambda.meta.adj <- GClambda(meta_cohort$P.adj)
            cat(paste(" - Genomic control lambda calculated from the corrected meta summary statistics = ",format(lambda.meta.adj, digits=4)," \n", sep=""))
            write.table(paste(" - Genomic control lambda calculated from the corrected meta summary statistics = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
    	} else {
            lambda.meta.adj <- lambda.meta
            cat(paste(" - Pre-set genomic control lambda applied to the corrected meta summary statistics = ",format(lambda.meta.adj, digits=4)," \n", sep=""))
            write.table(paste(" - Pre-set genomic control lambda applied to the corrected meta summary statistics = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
          }
          meta_cohort$P.adj[subset] <- 2*pnorm(-abs(meta_cohort$BETA.adj[subset])/(meta_cohort$SE.adj[subset]*sqrt(lambda.meta.adj)))
        }
        if ("LP.adj" %in% names(meta_cohort)) meta_cohort$LP.adj[subset] <- -log10(meta_cohort$P.adj[subset])
        if ("EAF.cohort" %in% names(meta_cohort)) meta_cohort$EAF.adj[subset] <- (2*meta_cohort$N.adj[subset]*meta_cohort$EAF.adj[subset]-2*meta_cohort$N.cohort[subset]*meta_cohort$EAF.cohort[subset])/(2*meta_cohort$N.adj[subset]-2*meta_cohort$N.cohort[subset])
        if ("N.cohort" %in% names(meta_cohort)) meta_cohort$N.adj[subset] <- meta_cohort$N.adj[subset]-meta_cohort$N.cohort[subset]
        if ("NSTUDIES" %in% names(meta_cohort)) meta_cohort$NSTUDIES.adj[subset] <- meta_cohort$NSTUDIES.adj[subset]-1
print("test4")

        if ("QHET" %in% names(meta_cohort)) {
		  meta_cohort$QHET.adj <- meta_cohort$QHET
		  meta_cohort$QHET.adj[subset] <- meta_cohort$QHET[subset]*meta_cohort$SE.adj[subset]^2/meta_cohort$SE.meta[subset]^2-(meta_cohort$BETA.adj[subset]-meta_cohort$BETA.meta[subset])^2-1/meta_cohort$SE.cohort[subset]^2*(meta_cohort$BETA.cohort[subset]-meta_cohort$BETA.meta[subset])^2*meta_cohort$SE.adj[subset]^2
          meta_cohort$QHETP.adj <- 1-pchisq(meta_cohort$QHET.adj,meta_cohort$NSTUDIES.adj-1)
          meta_cohort$I2HET.adj <- 100*(meta_cohort$QHET.adj-meta_cohort$NSTUDIES.adj+1)/meta_cohort$QHET.adj
        }
print("test5")

        if ("Z" %in% names(meta_cohort) | "Z.meta" %in% names(meta_cohort)) meta_cohort$Z.adj[subset] <- meta_cohort$BETA.adj[subset] / meta_cohort$SE.adj[subset] 
		
        vars <- c("LP.adj", "EAF.adj", "N.adj", "NSTUDIES.adj", "DIRECTIONS.adj","QHET.adj", "QHETP.adj", "I2HET.adj", "Z.adj")
		adj <- meta_cohort[,c("BETA.adj", "SE.adj")]
		if ("P" %in% namesmeta) adj <- data.frame(adj,P.adj=meta_cohort$P.adj)
		for (var in vars) {
		  if (var %in% names(meta_cohort)) {
		    if (any(!is.na(meta_cohort[,var]))) {
              adj <- data.frame(adj, var=meta_cohort[,var])
			  names(adj)[names(adj)=="var"] <- var
			}
          }
		}
      } else { 
        write.table("ERROR: 'BETA' or 'SE' not found in cohort results",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
        stop("ERROR: 'BETA' or 'SE' not found in cohort results")
      }
    } else {
      write.table("ERROR: 'BETA' or 'SE' not found in meta summary statistics",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
      stop("ERROR: 'BETA' or 'SE' not found in meta summary statistics")
    }
  } else if (toupper(metamethod) == "FSSW") { #fixed effect sample size weighted
    cat(" - Assuming a fixed-effects sample size weighted meta-analysis method \n")
    write.table(" - Assuming a fixed-effects sample size weighted meta-analysis method",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
    adj <- NULL
    if ("BETA.meta" %in% names(meta_cohort) & "N.meta" %in% names(meta_cohort)) {
      if ("BETA.cohort" %in% names(meta_cohort) & "N.cohort" %in% names(meta_cohort)) {
        meta_cohort$BETA.adj <- meta_cohort$BETA.meta
        meta_cohort$N.adj <- meta_cohort$N.meta
        if ("SE.meta" %in% names(meta_cohort)) { 
		  meta_cohort$SE.adj <- meta_cohort$SE.meta/sqrt(lambda.meta)
          meta_cohort$P.adj <- 2*pnorm(-abs(meta_cohort$BETA.adj)/meta_cohort$SE.adj) 
          if ("SE.cohort" %in% names(meta_cohort)) meta_cohort$SE.cohort <- meta_cohort$SE.cohort/sqrt(lambda.cohort)
		} else if ("SE" %in% namesmeta) { 
		  meta_cohort$SE.adj <- meta_cohort$SE/sqrt(lambda.meta) 
		}
		if ("LP" %in% names(meta_cohort)) meta_cohort$LP.adj <- -log10(meta_cohort$P.adj)
        if ("EAF.meta" %in% names(meta_cohort)) { meta_cohort$EAF.adj <- meta_cohort$EAF.meta } else { meta_cohort$EAF.adj <- NA }
        if ("NSTUDIES" %in% names(meta_cohort)) { meta_cohort$NSTUDIES.adj <- meta_cohort$NSTUDIES } else { meta_cohort$NSTUDIES.adj <- NA }
        meta_cohort$DIRECTIONS.adj <- NA
		if ("DIRECTIONS" %in% names(meta_cohort)) { 
		  adj_directions <- subtract.directions(meta_cohort, logfile)
          if (adj_directions[[2]]) meta_cohort$DIRECTIONS.adj <- adj_directions[[1]]
        }

        subset <- which(!is.na(meta_cohort$BETA.cohort) & !is.na(meta_cohort$SE.cohort))
        
        meta.weight <- meta_cohort$N.adj[subset]
        cohort.weight <- meta_cohort$N.cohort[subset]
        meta_cohort$BETA.adj[subset] <- (meta.weight*meta_cohort$BETA.adj[subset]-cohort.weight*meta_cohort$BETA.cohort[subset])/(meta.weight-cohort.weight)
        if ("SE.adj" %in% names(meta_cohort)) meta_cohort$SE.adj[subset] <- sqrt(1/(1/(meta_cohort$SE.adj[subset])^2-1/(meta_cohort$SE.cohort[subset])^2))
        if ("P.adj" %in% names(meta_cohort)) {
		  meta_cohort$P.adj[subset] <- 2*pnorm(-abs(meta_cohort$BETA.adj[subset])/meta_cohort$SE.adj[subset])
          if (gc_meta) {
            if (calculate_lambda.meta) {
              lambda.meta.adj <- GClambda(meta_cohort$P.adj)
              cat(paste(" - Genomic control lambda calculated from the corrected meta summary statistics = ",format(lambda.meta.adj, digits=4)," \n", sep=""))
              write.table(paste(" - Genomic control lambda calculated from the corrected meta summary statistics = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
            } else {
              lambda.meta.adj <- lambda.meta
              cat(paste(" - Pre-set genomic control lambda applied to the corrected meta summary statistics = ",format(lambda.meta.adj, digits=4)," \n", sep=""))
              write.table(paste(" - Pre-set genomic control lambda applied to the corrected meta summary statistics = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
            }
            meta_cohort$P.adj[subset] <- 2*pnorm(-abs(meta_cohort$BETA.adj[subset])/(meta_cohort$SE.adj[subset]*sqrt(lambda.meta.adj)))
          }
		} else {
          if (gc_meta) {
		    s <- "No genomic control correction is applied as there are no P-values in the meta summary statistics file."
			cat(paste(s,"\n"))
            write.table(s,logfile,col.names=F,row.names=F,quote=F,append=TRUE)
		  }
        }
		
        if ("LP.adj" %in% names(meta_cohort)) meta_cohort$LP.adj[subset] <- -log10(meta_cohort$P.adj[subset])
        if ("EAF.adj" %in% names(meta_cohort)) meta_cohort$EAF.adj[subset] <- (2*meta_cohort$N.adj[subset]*meta_cohort$EAF.adj[subset]-2*meta_cohort$N.cohort[subset]*meta_cohort$EAF.cohort[subset])/(2*meta_cohort$N.adj[subset]-2*meta_cohort$N.cohort[subset])
        meta_cohort$N.adj[subset] <- meta_cohort$N.adj[subset]-meta_cohort$N.cohort[subset]
        if ("NSTUDIES" %in% names(meta_cohort)) meta_cohort$NSTUDIES.adj[subset] <- meta_cohort$NSTUDIES.adj[subset]-1

        if ("QHET" %in% names(meta_cohort)) {
		  meta_cohort$QHET.adj <- meta_cohort$QHET
		  meta_cohort$QHET.adj[subset] <- meta_cohort$QHET[subset]*meta_cohort$SE.adj[subset]^2/meta_cohort$SE.meta[subset]^2-(meta_cohort$BETA.adj[subset]-meta_cohort$BETA.meta[subset])^2-1/meta_cohort$SE.cohort[subset]^2*(meta_cohort$BETA.cohort[subset]-meta_cohort$BETA.meta[subset])^2*meta_cohort$SE.adj[subset]^2
          meta_cohort$QHETP.adj <- 1-pchisq(meta_cohort$QHET.adj,meta_cohort$NSTUDIES.adj-1)
          meta_cohort$I2HET.adj <- 100*(meta_cohort$QHET.adj-meta_cohort$NSTUDIES.adj+1)/meta_cohort$QHET.adj
        }
		
        if ("Z.meta" %in% names(meta_cohort)) { 
		  meta_cohort$Z.adj[subset] <- meta_cohort$BETA.adj[subset] / meta_cohort$SE.adj[subset] 
		}

        vars <- c("LP.adj", "EAF.adj", "N.adj", "NSTUDIES.adj", "DIRECTIONS.adj","QHET.adj", "QHETP.adj", "I2HET.adj", "Z.adj")
		adj <- meta_cohort[,c("BETA.adj", "SE.adj")]
		if ("P" %in% namesmeta) adj <- data.frame(adj,P.adj=meta_cohort$P.adj)
		for (var in vars) {
		  if (var %in% names(meta_cohort)) {
		    if (any(!is.na(meta_cohort[,var]))) {
              adj <- data.frame(adj, var=meta_cohort[,var])
			  names(adj)[names(adj)=="var"] <- var
			}
          }
		}
      } else { 
        write.table("ERROR: 'BETA' or 'N' not found in cohort results",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
        stop("ERROR: 'BETA' or 'N' not found in cohort results")
      }
    } else {
      write.table("ERROR: 'BETA' or 'N' not found in meta summary statistics",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
      stop("ERROR: 'BETA' or 'N' not found in meta summary statistics")
    }
  } else if (toupper(metamethod) == "FSZ") { #fixed effects sample size z-score method
    cat(" - Assuming a fixed-effects sqrt sample size z-score weighted meta-analysis method \n")
    write.table(" - Assuming a fixed-effects sqrt sample size z-score weighted meta-analysis method",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
    adj <- NULL
    if (!("Z.meta" %in% names(meta_cohort))) {
      if ("BETA.meta" %in% names(meta_cohort) & "SE.meta" %in% names(meta_cohort)) {
        meta_cohort$Z.meta <- meta_cohort$BETA.meta/meta_cohort$SE.meta
      }
    }
    if (!("Z.cohort" %in% names(meta_cohort))) {
      if ("BETA.cohort" %in% names(meta_cohort) & "SE.cohort" %in% names(meta_cohort)) {
        meta_cohort$Z.cohort <- meta_cohort$BETA.cohort/meta_cohort$SE.cohort
      }
    }
    
    if ("Z.meta" %in% names(meta_cohort) & "N.meta" %in% names(meta_cohort)) {
      if ("Z.cohort" %in% names(meta_cohort) & "N.cohort" %in% names(meta_cohort)) {
        meta_cohort$Z.adj <- meta_cohort$Z.meta*sqrt(lambda.meta)
        if ("SE.meta" %in% names(meta_cohort)) { meta_cohort$SE.adj <- meta_cohort$SE.meta/sqrt(lambda.meta) } else if ("SE" %in% namesmeta) { meta_cohort$SE.adj <- meta_cohort$SE/sqrt(lambda.meta) } 
        meta_cohort$P.adj <- 2*pnorm(-abs(meta_cohort$Z.adj))
        if ("LP" %in% namesmeta) meta_cohort$LP.adj <- -log10(meta_cohort$P.adj)
		if ("EAF.meta" %in% names(meta_cohort)) { meta_cohort$EAF.adj <- meta_cohort$EAF.meta } else { meta_cohort$EAF.adj <- NA } 
        meta_cohort$N.adj <- meta_cohort$N.meta
        if ("NSTUDIES" %in% names(meta_cohort)) meta_cohort$NSTUDIES.adj <- meta_cohort$NSTUDIES
        meta_cohort$DIRECTIONS.adj <- NA
		if ("DIRECTIONS" %in% names(meta_cohort)) { 
		  adj_directions <- subtract.directions(meta_cohort, logfile)
          if (adj_directions[[2]]) meta_cohort$DIRECTIONS.adj <- adj_directions[[1]]
        }

        meta_cohort$Z.cohort <- meta_cohort$Z.cohort*sqrt(lambda.cohort)
        if ("SE" %in% namescohort) meta_cohort$SE.cohort <- meta_cohort$SE.cohort/sqrt(lambda.cohort)
        
        subset <- which(!is.na(meta_cohort$Z.cohort) & !is.na(meta_cohort$N.cohort))
        
        meta_cohort$Z.adj[subset] <- (sqrt(meta_cohort$N.adj[subset])*meta_cohort$Z.adj[subset]-sqrt(meta_cohort$N.cohort[subset])*meta_cohort$Z.cohort[subset]) / sqrt((meta_cohort$N.adj[subset]-meta_cohort$N.cohort[subset]))
        meta_cohort$P.adj[subset] <- 2*pnorm(-abs(meta_cohort$Z.adj[subset]))
        if (gc_meta) {
          if (calculate_lambda.meta) {
            lambda.meta.adj <- GClambda(meta_cohort$P.adj)
            cat(paste(" - Genomic control lambda calculated from the corrected meta summary statistics = ",format(lambda.meta.adj, digits=4)," \n", sep=""))
            write.table(paste(" - Genomic control lambda calculated from the corrected meta summary statistics = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
          } else {
            lambda.meta.adj <- lambda.meta
            cat(paste(" - Pre-set genomic control lambda applied to the corrected meta summary statistics = ",format(lambda.meta.adj, digits=4)," \n", sep=""))
            write.table(paste(" - Pre-set genomic control lambda applied to the corrected meta summary statistics = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
          }
          meta_cohort$P.adj[subset] <- 2*pnorm(-abs(meta_cohort$Z.adj[subset]/sqrt(lambda.meta.adj)))
        }
        if ("LP.adj" %in% names(meta_cohort)) meta_cohort$LP.adj[subset] <- -log10(meta_cohort$P.adj[subset])
        if ("EAF.adj" %in% names(meta_cohort)) meta_cohort$EAF.adj[subset] <- (2*meta_cohort$N.adj[subset]*meta_cohort$EAF.adj[subset]-2*meta_cohort$N.cohort[subset]*meta_cohort$EAF.cohort[subset])/(2*meta_cohort$N.adj[subset]-2*meta_cohort$N.cohort[subset])
        meta_cohort$N.adj[subset] <- meta_cohort$N.adj[subset]-meta_cohort$N.cohort[subset]
        if ("NSTUDIES" %in% names(meta_cohort)) meta_cohort$NSTUDIES.adj[subset] <- meta_cohort$NSTUDIES.adj[subset]-1

        if ("QHET" %in% names(meta_cohort) & "BETA.adj" %in% names(meta_cohort) & "SE.adj" %in% names(meta_cohort) & "BETA.meta" %in% names(meta_cohort) & "SE.meta" %in% names(meta_cohort) & "BETA.cohort" %in% names(meta_cohort) & "SE.cohort" %in% names(meta_cohort)) {
		  meta_cohort$QHET.adj <- meta_cohort$QHET
		  meta_cohort$QHET.adj[subset] <- meta_cohort$QHET[subset]*meta_cohort$SE.adj[subset]^2/meta_cohort$SE.meta[subset]^2-(meta_cohort$BETA.adj[subset]-meta_cohort$BETA.meta[subset])^2-1/meta_cohort$SE.cohort[subset]^2*(meta_cohort$BETA.cohort[subset]-meta_cohort$BETA.meta[subset])^2*meta_cohort$SE.adj[subset]^2
          meta_cohort$QHETP.adj <- 1-pchisq(meta_cohort$QHET.adj,meta_cohort$NSTUDIES.adj-1)
          meta_cohort$I2HET.adj <- 100*(meta_cohort$QHET.adj-meta_cohort$NSTUDIES.adj+1)/meta_cohort$QHET.adj
        }

        vars <- c("LP.adj", "BETA.adj", "SE.adj", "EAF.adj", "N.adj", "NSTUDIES.adj", "DIRECTIONS.adj","QHET.adj", "QHETP.adj", "I2HET.adj")
		adj <- meta_cohort[,c("Z.adj")]
		if ("P" %in% namesmeta) adj <- data.frame(adj,P.adj=meta_cohort$P.adj)
		for (var in vars) {
		  if (var %in% names(meta_cohort)) {
		    if (any(!is.na(meta_cohort[,var]))) {
              adj <- data.frame(adj, var=meta_cohort[,var])
			  names(adj)[names(adj)=="var"] <- var
			}
          }
		}

      } else { 
        write.table("ERROR: 'Z' or 'N' not found in cohort results",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
        stop("ERROR: 'Z' or 'N' not found in cohort results")
      }
    } else {
      write.table("ERROR: 'Z' or 'N' not found in meta summary statistics",logfile,col.names=F,row.names=F,quote=F,append=TRUE)
      stop("ERROR: 'Z' or 'N' not found in meta summary statistics")
    }
  }
  
  lambda.meta.adj <- GClambda(adj$P.adj)
  cat(paste(" - Genomic control lambda calculated of the final corrected meta summary statistics = ",format(lambda.meta.adj, digits=4)," \n", sep=""))  
  write.table(paste(" - Genomic control lambda calculated of the final corrected meta summary statistics = ",format(lambda.meta.adj, digits=4), sep=""),logfile,col.names=F,row.names=F,quote=F,append=TRUE)
  return(adj)
}

meta.subtract <- function(metafile, cohortfiles, metamethod = "FIV", lambda.meta = 1, lambdas.cohort = 1, gc_meta = TRUE, calculate_lambda.meta = TRUE, calculate_lambdas.cohort = TRUE, alternative = "alternative_headers.txt", save.as.data.frame = TRUE, savefile = "meta.results_corrected.with.MetaSubtract.txt.gz", logfile = "MetaSubtract.log", dir = tempdir()) {
  cat(paste("Analysis started",date(),"\n"))
  if (!is.character(dir)) {
    stop(paste("Invalid format of dir '",dir,"'. Current directory is ",getwd(), sep=""))	
    write.table(paste("Invalid format of dir '",dir,"'"),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
  } else if (!dir.exists(dir)) {
    stop(paste("Directory '",dir,"' does not exist. Current directory is ",getwd(), sep=""))
	write.table(paste("Directory '",dir,"' does not exist. Current directory is ",getwd(), sep=""),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
  }
  write.table(paste("Analysis started",date(),"\n"),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=F)

  
  if (length(lambdas.cohort)<length(cohortfiles) & !calculate_lambdas.cohort) {
    warning(paste("Less GC lambdas have been given than cohort files.\n The last",length(cohortfiles)-length(lambdas.cohort),"cohort files will be assumed to have GC lambda=1."))
    write.table(paste("Less GC lambdas have been given than cohort files.\n The last",length(cohortfiles)-length(lambdas.cohort),"cohort files will be assumed to have GC lambda=1."),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
    lambdas.cohort <- c(lambdas.cohort,rep(1,length(cohortfiles)-length(lambdas.cohort)))
  }  

  if (any(duplicated(cohortfiles))) {
    warning("Some cohort names appear more than once in the list.")
    write.table("Some cohort names appear more than once in the list.",file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
    lambdas.cohort <- lambdas.cohort[!duplicated(cohortfiles)]
    cohortfiles <- cohortfiles[!duplicated(cohortfiles)]
  }

  if (!file.exists(alternative)) { 
    if (!file.exists(system.file("extdata",alternative,package="MetaSubtract"))) { 
      write.table(paste("Alternative header file '",alternative,"' does not exist.", sep=""),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
      stop(paste("Alternative header file '",alternative,"' does not exist.", sep=""))
    } else {
      alternative<-system.file("extdata",alternative,package="MetaSubtract")
    }
  }
  header_translations <- read.table(alternative)
  
  if (!(metamethod %in% c("FIV", "FSSW", "FSZ"))) {
    write.table("ERROR: Please specify one of the following meta-analysis methods: FIV, FSSW, or FSZ",file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=TRUE)
    stop("ERROR: Please specify one of the following meta-analysis methods: FIV, FSSW, or FSZ")
  }

  if (!file.exists(metafile)) { 
    if (!file.exists(system.file("extdata",metafile,package="MetaSubtract"))) { 
      write.table(paste("'",metafile,"' does not exist", sep=""),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
      stop(paste("'",metafile,"' does not exist", sep=""))
	} else {
	  metafile <- system.file("extdata",metafile,package="MetaSubtract")
	}
  }
  cat(paste("Reading",metafile,"\n"))
  write.table(paste("Reading",metafile),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
  meta <- read.table(metafile, header = TRUE, as.is = TRUE)
  header_meta <- names(meta)
  headersinfometa <- F_translate_headers(header = header_meta, alternative = header_translations)
  names(meta) <- headersinfometa$header_h
  if ("EFFECTALLELE" %in% names(meta)) meta$EFFECTALLELE <- toupper(meta$EFFECTALLELE)
  if ("OTHERALLELE" %in% names(meta)) meta$OTHERALLELE <- toupper(meta$OTHERALLELE)
  markercol<-which("MARKER" %in% names(meta))
  if (markercol==0) {
    write.table(paste("No marker names are present in '",metafile,"'. Meta results will not be corrected.", sep=""),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
    stop(paste("No marker names are present in '",metafile,"'. Meta results will not be corrected.", sep=""))
  }
  if (any(duplicated(meta$MARKER))) {
    write.table(paste("Duplicate marker names are present in '",metafile,"'. Meta results will not be corrected.", sep=""),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
    stop(paste("Duplicate marker names are present in '",metafile,"'. Meta results will not be corrected.", sep=""))
  }
  
  meta$order<-c(1:nrow(meta))

  
  for (cohortfile in cohortfiles) {
    cohortfile.o <- cohortfile
    if (!file.exists(cohortfile) & !file.exists(system.file("extdata",cohortfile,package="MetaSubtract"))) { 
      warning(paste("'",cohortfile,"' does not exist. Meta results will not be corrected.", sep=""))
      write.table(paste("'",cohortfile,"' does not exist. Meta results will not be corrected.", sep=""),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
    } else {
      if (!file.exists(cohortfile)) { 
	  	cohortfile <- system.file("extdata",cohortfile,package="MetaSubtract")
      }
      cohort <- read.table(cohortfile, header = TRUE, as.is = TRUE)
      cat(paste("Subtracting effects of",cohortfile," \n"))
      write.table(paste("Subtracting effects of",cohortfile),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
      headersinfo <- F_translate_headers(header = names(cohort), alternative = header_translations)
      names(cohort) <- headersinfo$header_h
      if ("EFFECTALLELE" %in% names(cohort)) cohort$EFFECTALLELE <- toupper(cohort$EFFECTALLELE)
      if ("OTHERALLELE" %in% names(cohort)) cohort$OTHERALLELE <- toupper(cohort$OTHERALLELE)
      markercol2<-which("MARKER" %in% names(cohort))
      if (markercol2==0) {
        write.table(paste("No marker names are present in '",cohortfile,"'. Meta results will not be corrected for this cohort.", sep=""),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
        stop(paste("No marker names are present in '",cohortfile,"'. Meta results will not be corrected for this cohort.", sep=""))
      }
      if (any(duplicated(cohort$MARKER))) {
        write.table(paste("Duplicate marker names are present in '",cohortfile,"'. Meta results will not be corrected for this cohort.", sep=""),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
        stop(paste("Duplicate marker names are present in '",cohortfile,"'. Meta results will not be corrected for this cohort.", sep=""))
      }

      meta_cohort <- merge(meta, cohort, by = "MARKER", all.x = TRUE, all.y = FALSE, suffixes=c(".meta", ".cohort"))
 
      if ("EFFECTALLELE.meta" %in% names(meta_cohort) & "EFFECTALLELE.cohort" %in% names(meta_cohort)) {
        # flip alleles in cohort so that they match with those in the meta and adjust BETA and EAF
        subset = !is.na(meta_cohort$EFFECTALLELE.cohort) & meta_cohort$EFFECTALLELE.meta != meta_cohort$EFFECTALLELE.cohort
        if ("OTHERALLELE.cohort" %in% names(meta_cohort)) { 
          meta_cohort[subset, c("BETA.cohort")] <- -meta_cohort[subset, c("BETA.cohort")]
          if ("EAF.cohort" %in% names(meta_cohort)) meta_cohort[subset, c("EAF.cohort")] <- 1-meta_cohort[subset, c("EAF.cohort")]
          meta_cohort[subset, c("EFFECTALLELE.cohort", "OTHERALLELE.cohort")] <- meta_cohort[subset, c("OTHERALLELE.cohort", "EFFECTALLELE.cohort")]
          write.table(paste(" - For",sum(subset),"markers alleles have been flipped"),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
		} else {
          if ("OTHERALLELE.meta" %in% names(meta_cohort)) {
		    subset2 <- subset & !is.na(meta_cohort$OTHERALLELE.meta) & meta_cohort$OTHERALLELE.meta == meta_cohort$EFFECTALLELE.cohort
            meta_cohort[subset2, c("BETA.cohort")] <- -meta_cohort[subset2, c("BETA.cohort")]
            meta_cohort[subset2, c("EAF.cohort")] <- 1-meta_cohort[subset2, c("EAF.cohort")]
            meta_cohort[subset2, "EFFECTALLELE.cohort"] <- meta_cohort[subset2, "EFFECTALLELE.meta"]
            write.table(paste(" - For",sum(subset2),"markers alleles have been flipped"),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
		  } else {
            meta_cohort[subset, c("BETA.cohort")] <- -meta_cohort[subset, c("BETA.cohort")]
            if ("EAF.cohort" %in% names(meta_cohort)) meta_cohort[subset, c("EAF.cohort")] <- 1-meta_cohort[subset, c("EAF.cohort")]
            meta_cohort[subset, "EFFECTALLELE.cohort"] <- meta_cohort[subset, "EFFECTALLELE.meta"]
            write.table(paste(" - For",sum(subset),"markers alleles have been flipped"),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
		  }	
        } 
		
        # strand switch of alleles
        if ("OTHERALLELE.cohort" %in% names(meta_cohort)) { 
          subset = !is.na(meta_cohort$EFFECTALLELE.cohort) & meta_cohort$EFFECTALLELE.meta == complement(meta_cohort$EFFECTALLELE.cohort) & meta_cohort$EFFECTALLELE.meta != meta_cohort$OTHERALLELE.cohort
          write.table(paste(" - For",sum(subset),"markers alleles have been swapped for strand"),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
          if (sum(subset>0)) { meta_cohort[subset, c("EFFECTALLELE.cohort", "OTHERALLELE.cohort")] <- cbind(complement(meta_cohort[subset, "EFFECTALLELE.cohort"]), complement(meta_cohort[subset, "OTHERALLELE.cohort"])) }
		} else {
          if ("OTHERALLELE.meta" %in% names(meta_cohort)) { 
            subset = !is.na(meta_cohort$EFFECTALLELE.cohort) & meta_cohort$EFFECTALLELE.meta == complement(meta_cohort$EFFECTALLELE.cohort) & meta_cohort$OTHERALLELE.meta != meta_cohort$OTHERALLELE.cohort
            write.table(paste(" - For",sum(subset),"markers alleles have been swapped for strand"),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
            if (sum(subset>0)) { meta_cohort[subset, c("EFFECTALLELE.cohort", "OTHERALLELE.cohort")] <- cbind(complement(meta_cohort[subset, "EFFECTALLELE.cohort"]), complement(meta_cohort[subset, "OTHERALLELE.cohort"])) }
		  } 		    
		}
		
    
        # strand switch of alleles and flip alleles in cohort so that they match with those in the meta and adjust BETA and EAF
        if ("OTHERALLELE.cohort" %in% names(meta_cohort)) { 
          subset = !is.na(meta_cohort$EFFECTALLELE.cohort) & meta_cohort$EFFECTALLELE.meta == complement(meta_cohort$OTHERALLELE.cohort) & meta_cohort$EFFECTALLELE.meta != meta_cohort$EFFECTALLELE.cohort
          write.table(paste(" - For",sum(subset),"markers alleles have been flipped and swapped for strand"),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
          meta_cohort[subset, c("BETA.cohort")] <-  -meta_cohort[subset, c("BETA.cohort")]
          if ("EAF.cohort" %in% names(meta_cohort)) meta_cohort[subset, c("EAF.cohort")] <- 1-meta_cohort[subset, c("EAF.cohort")]
          if (sum(subset>0)) { meta_cohort[subset, c("EFFECTALLELE.cohort", "OTHERALLELE.cohort")] <- cbind(complement(meta_cohort[subset, "OTHERALLELE.cohort"]), complement(meta_cohort[subset, "EFFECTALLELE.cohort"])) }
		} else {
		  if ("OTHERALLELE.meta" %in% names(meta_cohort)) { 
            subset = !is.na(meta_cohort$EFFECTALLELE.cohort) & meta_cohort$OTHERALLELE.meta == complement(meta_cohort$EFFECTALLELE.cohort) & meta_cohort$EFFECTALLELE.meta != meta_cohort$EFFECTALLELE.cohort
            write.table(paste(" - For",sum(subset),"markers alleles have been flipped and swapped for strand"),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
            meta_cohort[subset, c("BETA.cohort")] <-  -meta_cohort[subset, c("BETA.cohort")]
            if ("EAF.cohort" %in% names(meta_cohort)) meta_cohort[subset, c("EAF.cohort")] <- 1-meta_cohort[subset, c("EAF.cohort")]
            if (sum(subset>0)) { meta_cohort[subset, c("EFFECTALLELE.cohort", "OTHERALLELE.cohort")] <- cbind(complement(meta_cohort[subset, "OTHERALLELE.cohort"]), complement(meta_cohort[subset, "EFFECTALLELE.cohort"])) }
          }
		}
      
        # remove genetic markers that still do not match alleles
        subset = !is.na(meta_cohort$EFFECTALLELE.cohort) & !( (meta_cohort$EFFECTALLELE.meta == meta_cohort$EFFECTALLELE.meta & meta_cohort$OTHERALLELE.meta == meta_cohort$OTHERALLELE.cohort) |
	                                                     (meta_cohort$EFFECTALLELE.meta == meta_cohort$OTHERALLELE.cohort & meta_cohort$OTHERALLELE.meta == meta_cohort$EFFECTALLELE.cohort) |
               	                                         (meta_cohort$EFFECTALLELE.meta == complement(meta_cohort$EFFECTALLELE.cohort) & meta_cohort$OTHERALLELE.meta == complement(meta_cohort$OTHERALLELE.cohort)) |
                                                         (meta_cohort$EFFECTALLELE.meta == complement(meta_cohort$OTHERALLELE.cohort) & meta_cohort$OTHERALLELE.meta == complement(meta_cohort$EFFECTALLELE.cohort)) )
        s<-paste(" - Number of genetic markers removed because of mismatch in alleles: ", sum(subset), " (", signif(100*sum(subset)/nrow(cohort),2), "%)",sep="")
	    if (sum(subset)>0) { 
          cat(paste(s," \n"))
	      meta_cohort[subset, c("BETA.cohort", "SE.cohort", "P.cohort")] <- NA 
		  write.table(meta_cohort[subset,c("MARKER", "EFFECTALLELE.meta", "OTHERALLELE.meta", "EFFECTALLELE.cohort", "OTHERALLELE.cohort")],file.path(dir,paste0(basename(cohortfile),".allele_mismatch.txt")), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t") 
	    }
	    write.table(s,file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
	  } else {
        write.table(paste(" - No alleles are given, so effects are assumed to be given for same allele as in the meta summary statistics."),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
      }
	  
      # added by L.Corbin 7th August 2017
      # corrects error which was leading to assignment of adjusted betas to wrong SNPs
      meta_cohort <- meta_cohort[order(meta_cohort$order),]
      # end of edit

	  # adjusting meta-GWAS summary statistic for current cohort's GWAD results
      adj <- meta.min.1(meta_cohort, metamethod, lambda.meta, lambdas.cohort[which(cohortfile.o==cohortfiles)], gc_meta, calculate_lambda.meta, calculate_lambdas.cohort, logfile = file.path(dir,logfile), names(meta), names(cohort), index=which(cohortfile.o==cohortfiles))

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
      write.table(s,file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)

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
      write.table(paste(s," \n"),file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
    }
  }
  meta <- meta[order(meta$order),]
  if (markercol==1) { 
    meta.adj <- meta[,c(1:headersinfometa$header_N)]
  } else {
    meta.adj <- meta[,c(2:markercol,1,(markercol+1):headersinfometa$header_N)]
  }
  names(meta.adj) <- header_meta
  
  if (!is.null(savefile) & !is.na(savefile)) {
    if (is.character(savefile)) {
	  s <- paste("Corrected meta summary statistics are saved to ",savefile,".",sep="")
	  cat(paste(s,"\n"))
	  write.table(s,file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
	  if (substr(savefile,nchar(savefile)-2,nchar(savefile)) == ".gz") {
        write.table(meta.adj, gzfile(file.path(dir,savefile)), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
	  } else {
        write.table(meta.adj, file.path(dir,savefile), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
	  }	
	} else {
	  s <- paste(savefile," is not a character string. Corrected meta summary statistics have been saved to 'meta.results_corrected.with.MetaSubtract.txt.gz'.",sep="")
	  cat(paste(s,"\n"))
	  savefile <- "meta.results_corrected.with.MetaSubtract.txt.gz"
	  write.table(s,file.path(dir,logfile),col.names=F,row.names=F,quote=F,append=T)
      write.table(meta.adj, file.path(dir,savefile), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    }
  }	
	
  if (save.as.data.frame) {
    return(meta.adj)
  } else {
    return(NULL)
  }	
}

