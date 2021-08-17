# Load Required Packages
for (lib in c(
  'dplyr',
  'pbapply',
  'lmerTest',
  'car',
  'parallel',
  'MuMIn',
  'glmmTMB',
  'MASS',
  'cplm',
  'pscl'
)) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}


#####################
## multitesting correction ##
#####################
# Combine results and multiplicity correction
rerun_multiplicity_correction <- function (infile, all_outfile, sig_outfile, effects_name, sig_level, correction) {
    effects_name <- unlist(strsplit(effects_name, ",", fixed = TRUE))
    paras <- read.table(infile, header=TRUE, sep="\t", quote="", comment.char = "", check.names = FALSE)
    paras <- subset(paras, metadata %in% effects_name)
    
    # multiplicity correction based on all metadata
    paras$qval <- as.numeric(p.adjust(as.numeric(paras$pval), correction))
    ordered_results <- paras[order(paras$qval), ]
    write.table(ordered_results, file = all_outfile, append = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE, quote=FALSE)
    
    sig <- subset(ordered_results, qval < as.numeric(sig_level))
    write.table(sig, file = sig_outfile, append = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE, quote=FALSE)
}

# Combine results and multiplicity correction
multiplicity_correction <- function (infile, outfile, effects_name, sig_level, correction) {
  print(infile)
  effects_name <- unlist(strsplit(effects_name, ",", fixed = TRUE))
  paras <- read.table(infile, header=TRUE, sep="\t", quote="", comment.char = "", check.names = FALSE)
  #paras <- subset(paras, metadata %in% effects_name)
  new_col_name <- c(colnames(paras), "qvalue")
  
  # multiplicity correction based on all metadata
  paras$qvalue <- as.numeric(p.adjust(as.numeric(paras$pval), correction))
  ordered_results <- paras[order(paras$qvalue), ]
  write.table(t(new_col_name), file = outfile, append = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
  write.table(ordered_results, file = outfile, append = TRUE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
  
  # multiplicity correction based on each variable
  paras <- subset(paras, metadata %in% effects_name)
  outfile_per <- gsub(".tsv", ".correct_per_variable.tsv", outfile)
  paras_per <- ddply(paras, .(metadata), transform, qvalue=as.numeric(p.adjust(as.numeric(pval), correction)))
  ordered_results_per <- paras_per[order(paras_per$metadata, paras_per$qvalue), ]
  write.table(t(new_col_name), file = outfile_per, append = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
  write.table(ordered_results_per, file = outfile_per, append = TRUE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
  
  # multiplicity correction based on each level
  outfile_per <- gsub(".tsv", ".correct_per_level.tsv", outfile)
  paras_per <- ddply(paras, .(metadata, value), transform, qvalue=as.numeric(p.adjust(as.numeric(pval), correction)))
  ordered_results_per <- paras_per[order(paras_per$metadata, paras_per$value, paras_per$qvalue), ]
  write.table(t(new_col_name), file = outfile_per, append = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
  write.table(ordered_results_per, file = outfile_per, append = TRUE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
}  


##################### 
## fit each feature ##
#####################
fit.each <- function(
    non_transform_features,
    non_transform_dnadata,
    features,
    metadata,
    dnadata,
    rna_dna_flt,
    min_abundance,
    min_samples,
    model,
    fixed_effects,
    effects_names,
    random_effects_formula = NULL) {

    #############################################################
    # Determine the function and summary for the model selected #
    #############################################################
  
    ################
    # Linear Model #
    ################
  
    if (model == "LM") {
      if (is.null(random_effects_formula)) {
        model_function <-
        function(formula, data, na.action) {
          return(glm(
            formula,
            data = data,
            family = 'gaussian',
            na.action = na.action
          ))
        }
        summary_function <- function(fit) {
        lm_summary <- summary(fit)$coefficients
        para <- as.data.frame(lm_summary)[-1, -3]
        para$name <- rownames(lm_summary)[-1]
        return(para)
        }
      } else {
        ranef_function <- lme4::ranef
        model_function <-
          function(formula, data, na.action) {
            return(lmerTest::lmer(
            formula, 
            data = data, 
            na.action = na.action))
          }
        summary_function <- function(fit) {
          lm_summary <- coef(summary(fit))
          para <- as.data.frame(lm_summary)[-1, -c(3:4)]
          para$name <- rownames(lm_summary)[-1]
          return(para)
        }
      }
    }
  
    ################
    # Logistic Model #
    ################
  
    if (model == "LOGIT") {
      if (is.null(random_effects_formula)) {
        model_function <-
          function(formula, data, na.action) {
            return(glm(
              formula,
              data = data,
              family = 'binomial',
              na.action = na.action
            ))
          }
        summary_function <- function(fit) {
          lm_summary <- summary(fit)$coefficients
          para <- as.data.frame(lm_summary)[-1, -3]
          para$name <- rownames(lm_summary)[-1]
          return(para)
        }
      } else {
        ranef_function <- lme4::ranef
        model_function <-
          function(formula, data, na.action) {
            return(glmer(
              formula, 
              data = data,
              family = 'binomial',
              na.action = na.action))
          }
        summary_function <- function(fit) {
          lm_summary <- coef(summary(fit))
          para <- as.data.frame(lm_summary)[-1, -3]
          para$name <- rownames(lm_summary)[-1]
          return(para)
        }
      }
    }
  
    ####################
    # Compound Poisson #
    ####################
  
    if (model == "CPLM") {
      if (is.null(random_effects_formula)) {
        model_function <- cplm::cpglm
        summary_function <- function(fit) {
          cplm_out <-
            capture.output(
              cplm_summary <- cplm::summary(fit)$coefficients)
          para <- as.data.frame(cplm_summary)[-1, -3]
          para$name <- rownames(cplm_summary)[-1]
          logging::logdebug(
            "Summary output\n%s", 
            paste(cplm_out, collapse = "\n"))
          return(para)
        }
      } else {
        ranef_function <- glmmTMB::ranef
        model_function <-
          function(formula, data, na.action) {
            return(glmmTMB::glmmTMB(
              formula,
              data = data,
              family=glmmTMB::tweedie(link = "log"),
              ziformula = ~0,
              na.action = na.action
            ))
          }
        summary_function <- function(fit) {
          glmmTMB_summary <- coef(summary(fit))
          para <- as.data.frame(glmmTMB_summary$cond)[-1, -3]
          para$name <- rownames(glmmTMB_summary$cond)[-1]
          return(para)
        }
      }
    }
  
    #####################
    # Negative Binomial #
    #####################
  
    if (model == "NEGBIN") {
      if (is.null(random_effects_formula)) {
        model_function <- MASS::glm.nb
        summary_function <- function(fit) {
          glm_summary <- summary(fit)$coefficients
          para <- as.data.frame(glm_summary)[-1, -3]
          para$name <- rownames(glm_summary)[-1]
          return(para)
        }
      } else {
        ranef_function <- glmmTMB::ranef
        model_function <-
          function(formula, data, na.action) {
            return(glmmTMB::glmmTMB(
              formula,
              data = data,
              family=glmmTMB::nbinom2(link = "log"),
              ziformula = ~0,
              na.action = na.action
            ))
          }
        summary_function <- function(fit) {
          glmmTMB_summary <- coef(summary(fit))
          para <- as.data.frame(glmmTMB_summary$cond)[-1, -3]
          para$name <- rownames(glmmTMB_summary$cond)[-1]
          return(para)
        }
      }
    }
  
    ###################################
    # Zero-inflated Negative Binomial #
    ###################################
  
    if (model == "ZINB") {
      if (is.null(random_effects_formula)) {
        model_function <-
          function(formula, data, na.action) {
            return(pscl::zeroinfl(
              formula,
              data = data,
              dist = "negbin",
              na.action = na.action))
          }
        summary_function <- function(fit) {
          pscl_summary <- summary(fit)$coefficients$count
          para <-as.data.frame(pscl_summary)[-c(1, (ncol(metadata) + 2)), -3]
          para$name <- rownames(pscl_summary)[-c(1, (ncol(metadata) + 2))]
          return(para)
        }
      } else {
        ranef_function <- glmmTMB::ranef
        model_function <-
          function(formula, data, na.action) {
            return(glmmTMB::glmmTMB(
              formula,
              data = data,
              family=glmmTMB::nbinom2(link = "log"),
              ziformula = ~1,
              na.action = na.action))
          }
        summary_function <- function(fit) {
          glmmTMB_summary <- coef(summary(fit))
          para <- as.data.frame(glmmTMB_summary$cond)[-1, -3]
          para$name <- rownames(glmmTMB_summary$cond)[-1]
          return(para)
        }
      }
    } 
    
  
    #############################################################
    # Per feature fitting #
    #############################################################
    output <- list()
    rownames <- c()
    for (i in seq_along(colnames(features))) {
        i <- colnames(features)[i]
        myname <- i
        mymeta <- metadata
        if (i %in% colnames(dnadata)) {
            flag <- 1
            mymeta <- data.frame(mymeta, i=dnadata[, i])
            names(mymeta)[names(mymeta) == "i"] <- i
        } else {
            flag <- 0
            tmp <- unlist(strsplit(i, "__", fixed = TRUE))
            if (length(tmp) > 1) { 
                if (tmp[1] %in% colnames(dnadata)) {
                    flag <- 1
                    mymeta <- data.frame(mymeta, i=dnadata[, tmp[1]])
                    names(mymeta)[names(mymeta) == "i"] <- tmp[1]
                    myname <- tmp[1]
                } else {
                    if (tmp[-1] %in% colnames(dnadata)) {
                        flag <- 1
                        mymeta <- data.frame(mymeta, i=dnadata[, tmp[-1]])
                        names(mymeta)[names(mymeta) == "i"] <- tmp[-1]
                        myname <- tmp[-1]
                    }   
                }
            }
            if (flag == 0) {
                tmp <- unlist(strsplit(i, "_", fixed = TRUE))
                if (length(tmp) > 1) {
                    if (tmp[1] %in% colnames(dnadata)) {
                        flag <- 1
                        mymeta <- data.frame(mymeta, i=dnadata[, tmp[1]])
                        names(mymeta)[names(mymeta) == "i"] <- tmp[1]
                        myname <- tmp[1]
                    } else {
                        if (tmp[-1] %in% colnames(dnadata)) {
                            flag <- 1
                            mymeta <- data.frame(mymeta, i=dnadata[, tmp[-1]])
                            names(mymeta)[names(mymeta) == "i"] <- tmp[-1]
                            myname <- tmp[-1]
                        }   
                    }
                }
                if (flag == 0) {
                    logging::logwarn(
                        paste("Feature name not found in dnadata",
                              "so not applied to formula as fixed effect: %s"), toString(colnames(mymeta)))
                    next        
                }   
            }
        }
        
        # check min(RNA, DNA) > 0 across sample
        if (rna_dna_flt != "none") {
            #rna_dna_flt <- as.numeric(rna_dna_flt)
            min_abundance <- as.numeric(min_abundance)
            mydna <- as.numeric(non_transform_dnadata[, myname])
            myrna <- as.numeric(non_transform_features[, i])
            
            # check abundance detectable
            mydna[mydna <= min_abundance] <- 0
            myrna[myrna <= min_abundance] <- 0
            mydna[is.na(mydna)] <- 0
            myrna[is.na(myrna)] <- 0
            #if (length(which(mydna > min_abundance)) <= min_samples) {
            #	next
            #}
            #if (length(which(myrna > min_abundance)) <= min_samples) {
            #	next
            #}
            
            # min(RNA, DNA) > 0 across sample
            if (rna_dna_flt == "lenient") {
                mydna_sum <- sum(mydna)
                myrna_sum <- sum(myrna)
                if (as.numeric(min(myrna_sum, mydna_sum)) <= 0) {
                    next
                }
            }
            
            # filter 0/0 cases per sample
            if (rna_dna_flt == "semi_strict") {
                mydna_sum <- sum(mydna)
                myrna_sum <- sum(myrna)
                if (as.numeric(min(myrna_sum, mydna_sum)) <= 0) {
                    next
                }
                myratio <- myrna/mydna
                indices <- which(myratio %in% NaN)
                if (length(indices) > 0) {
                    mymeta[indices, myname] <- NaN
                    features[indices, i] <- NaN
                }
            }
            
			      # filter 0/0, x/0, 0/x cases per sample
            if (rna_dna_flt == "strict") {
                mydna_sum <- sum(mydna)
                myrna_sum <- sum(myrna)
                if (as.numeric(min(myrna_sum, mydna_sum)) <= 0) {
                    next
                }
                myratio <- myrna/mydna
                indices <- which(myratio %in% NaN)
                if (length(indices) > 0) {
                    mymeta[indices, myname] <- NaN
                    features[indices, i] <- NaN
                }
                indices <- which(myratio %in% 0)
                if (length(indices) > 0) {
                    mymeta[indices, myname] <- NaN
                    features[indices, i] <- NaN
                }
                indices <- which(myratio %in% Inf)
                if (length(indices) > 0) {
                    mymeta[indices, myname] <- NaN
                    features[indices, i] <- NaN
                }
            }

        }
        
        mydf <- data.frame(expr = features[,i], mymeta)
        formula_text <-
            paste("expr ~ ", paste(colnames(mymeta), collapse = " + "))
        #logging::loginfo("Formula for fixed effects: %s", formula_text)
        formula <-
            tryCatch(
                as.formula(formula_text),
                error = function(e)
                    stop(
                        paste(
                            "Invalid formula.",
                            "Please provide a different formula: ",
                            formula_text
                        )   
                    )   
            ) 
        if (! is.null(random_effects_formula)) {
            formula <-
                paste(
                    '. ~', 
                    paste(all.vars(formula)[-1], collapse = ' + '), 
                    '.', 
                    sep = ' + ')
            formula <- update(random_effects_formula, formula)
        }
        logging::loginfo("Final formula: %s", formula)
        logging::loginfo("Metadata items: %s", toString(colnames(mymeta)))
        
        fit <- tryCatch({
            fit1 <-
                model_function(
                    formula, 
                    data = mydf, 
                    na.action = na.exclude)
        }, error = function(err) {
            fit1 <-
                try({
                    model_function(
                        formula, 
                        data = mydf, 
                        na.action = na.exclude)
                })
            return(fit1)
        })
        
        # Gather Output
        if (all(!inherits(fit, "try-error"))) {
            para <- summary_function(fit)
            residuals <- residuals(fit)
            myfitted <- fitted(fit)
            if (!(is.null(random_effects_formula))) {
              l <- ranef_function(fit)
              d<-as.vector(unlist(l))
              names(d)<-unlist(lapply(l, row.names))
              ranef<-d
            }
        }
        else{
            logging::logwarn(paste(
                "Fitting problem for feature", 
                i, 
                "returning NA"))
            para <-
                as.data.frame(matrix(NA, 
                                     nrow = ncol(mymeta), ncol = 3))
            para$name <- colnames(mymeta)
            residuals <- rep(NA, nrow(features))
            myfitted <- rep(NA, nrow(features))
            if (!(is.null(random_effects_formula))) ranef <- NA
        }
        colnames(para) <- c('coef', 'stderr' , 'pval', 'name')
        if (nrow(para) == 0) {
            para <- data.frame()
            residuals <- rep(NA, nrow(features))
            myfitted <- rep(NA, nrow(features))
            if (!(is.null(random_effects_formula))) ranef <- NA
        } else {
            para$feature <- i
        }
        output$para <- rbind(output$para, para)
        rownames <- c(rownames, i)
        output$residuals <- rbind(output$residuals, residuals) 
        output$fitted <- rbind(output$fitted, myfitted)
        if (!(is.null(random_effects_formula))) {
          output$ranef <- rbind(output$ranef, ranef) 
        }
    }
    row.names(output$residuals) <- rownames
    row.names(output$fitted) <- rownames
    if (!(is.null(random_effects_formula))) {
      row.names(output$ranef) <- rownames
    }
    return(output) 
}


#####################
## fit the data using the model based on each feature and applying the correction ##
#####################
fit.dnadata <-
    function(
        non_transform_features,
        non_transform_dnadata,
        features,
        metadata,
        dnadata,
        rna_dna_flt,
        min_abundance,
        min_samples,
        model,
        fixed_effects,
        effects_names,
        random_effects_formula = NULL,
        correction = "BH",
        cores = 1) {                
        
        ##############################
        # Apply per-feature modeling #
        ##############################
        outputs <- fit.each(
            non_transform_features,
            non_transform_dnadata,
            features,
            metadata,
            dnadata,
            rna_dna_flt,
            min_abundance,
            min_samples,
            model,
            fixed_effects,
            effects_names,
            random_effects_formula = NULL)
        
        
        # bind the results for each feature
        paras <- outputs$para
        residuals <- outputs$residuals
        colnames(residuals) <- rownames(features)
        fitted <- outputs$fitted
        #colnames(fitted) <- rownames(features)
        if (!(is.null(random_effects_formula))) {
          ranef <- outputs$ranef
          #colnames(ranef) <- rownames(features)
        }
        
        ################################
        # Apply correction to p-values #
        ################################
        paras$qval <- as.numeric(p.adjust(paras$pval, method = correction))
        
        #####################################################
        # Determine the metadata names from the model names #
        #####################################################
        metadata_names <- colnames(metadata[effects_names])
        # order the metadata names by decreasing length
        metadata_names_ordered <-
            metadata_names[order(
                nchar(metadata_names), decreasing = TRUE)]
        # find the metadata name based on the match 
        # to the beginning of the string
        extract_metadata_name <- function(name) {
            myvalue <- metadata_names_ordered[mapply(
                startsWith, 
                name, 
                metadata_names_ordered)][1]
            if(is.na(myvalue)) {
                myvalue <- name
            }
            return(myvalue)
        }
        paras$metadata <- unlist(lapply(paras$name, extract_metadata_name))
        # compute the value as the model contrast minus metadata
        paras$value <-
            mapply(function(x, y) {
                if (x == y)
                    x
                else
                    gsub(x, "", y)
            }, paras$metadata, paras$name)
        
        ##############################
        # Sort by decreasing q-value #
        ##############################
        
        if (length(paras$qval) > 0) {
            paras <- paras[order(paras$qval, decreasing = FALSE), ]
            paras <-
                dplyr::select(
                    paras,
                    c('feature', 'metadata', 'value'),
                    dplyr::everything())
            if (!(is.null(random_effects_formula))) {
              return(list("results" = paras, "residuals" = residuals, "fitted" = fitted, "ranef" = ranef))
            } else {
              return(list("results" = paras, "residuals" = residuals, "fitted" = fitted))
            }
        } else {
            return(NULL)
        }
    }


#####################
## Transformation ##
###################
# Presence/Absence Transformation #
PA <- function(x) { 
  y <- replace(x, x > 0, 1)
  return(y)
}

transformFeatures_ext <- function(features, transformation) {
  if (transformation == 'PA')     {
    features <- apply(features, 2, PA)
  }
  
  return(features)
}
