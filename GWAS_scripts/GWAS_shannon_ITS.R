################## This script runs the GWAS on the ITS Shannon diversity, using read depth as cofactor #####################

##################### We use the lme4qtl package. 

library(Matrix)
library(MASS)
library(ggplot2)
library(cowplot)
library(viridis)
library(lme4qtl) 
library(matlm) 
library(wlm) 
library(RNOmni)
library(parallel)
library(tictoc)
library(dplyr)
options(scipen = 999)

############## From R 4.3, some of the functions in the package stopped working and were not updated. Hence, we had to correct them ourselves and put them in the script explicitly.
############## The two functions below do exactly the same as in the original package.

matlm_pred3 <- function (x, ind = NULL, num_batches = 1, batch_size = NULL, 
                         path_pred = ".", ...) 
{
  stopifnot(!any(duplicated(ind)))
  out <- list(path = path_pred)
  # if (class(x) == "matrix") {
  # if (T) {
  if (any(class(x) == "matrix")) {
    out$data <- x
    oldClass(out) <- c("matlmPredMat", "matlmPred")
  }
  else if (class(x) == "big.matrix") {
    out$data <- describe(x)
    oldClass(out) <- c("matlmPredBigMat", "matlmPred")
  }
  else if (class(x) == "big.matrix.descriptor") {
    out$data <- x
    oldClass(out) <- c("matlmPredBigMat", "matlmPred")
  }
  else {
    stop("not supported class")
  }
  if (is.null(ind)) {
    ind <- seq(1, pred_nrow(out))
  }
  out$ind <- ind
  out$complete <- (length(ind) == pred_nrow(out))
  cnames <- pred_colnames(out)
  if (is.null(cnames)) {
    out$cnames <- as.character(paste0("pred", seq(1, pred_ncol(out))))
    out$null_cnames <- TRUE
  }
  else {
    out$cnames <- cnames
    out$null_cnames <- FALSE
  }
  if (num_batches > pred_ncol(out)) {
    warning("num_batches > pred_ncol: force num_batches <- pred_ncol")
    num_batches <- pred_ncol(out)
  }
  if (is.null(batch_size)) {
    beg <- seq(1, pred_ncol(out), length = num_batches) %>% 
      floor
    end <- c(beg[-1] - 1, pred_ncol(out))
  }
  else {
    beg <- seq(1, pred_ncol(out), by = batch_size)
    end <- c(beg[-1] - 1, pred_ncol(out))
  }
  stopifnot(all(beg <= pred_ncol(out)))
  stopifnot(all(end <= pred_ncol(out)))
  out$num_batches <- length(beg)
  out$batch_size <- batch_size
  out$beg <- beg
  out$end <- end
  return(out)
}

matlm2 <- function (formula, data, ..., varcov = NULL, transform = NULL, 
                    ids = NULL, pred, int, num_batches = 1, batch_size = NULL, 
                    path_pred = ".", num_perm = 0, seed_perm, returnRespOrth = FALSE, 
                    returnPredOrth = FALSE, returnPredOrthSc = FALSE, stats_full = FALSE, 
                    cores = 1, verbose = 0) 
{
  tic.clearlog()
  tic("matlm")
  tic("args")
  mc <- match.call()
  env <- parent.frame(1)
  missing_transform <- missing(transform)
  missing_varcov <- missing(varcov)
  missing_ids <- missing(ids)
  missing_seed_perm <- missing(seed_perm)
  stopifnot(class(formula) == "formula")
  data <- as.data.frame(data)
  model <- ifelse(missing(int), "marginal", "interaction")
  permutation <- ifelse(num_perm > 0, "perm_x", "none")
  weighted <- (!missing_transform | !missing_varcov)
  toc(log = TRUE, quiet = TRUE)
  tic("indices")
  if (verbose > 0) {
    cat(" - computing indices...\n")
  }
  ind_model <- matlm_ind(formula, data, ...)
  nobs_data <- nrow(data)
  nobs_model <- length(ind_model)
  nobs_omit <- nobs_data - nobs_model
  if (verbose > 1) {
    cat("  -- nobs_data", nobs_data, "/ nobs_model", nobs_model, 
        "/ nobs_omit", nobs_omit, "\n")
  }
  toc(log = TRUE, quiet = TRUE)
  tic("transform")
  if (verbose > 0) {
    cat(" - computing transform...\n")
  }
  if (weighted) {
    if (missing_transform) {
      stopifnot(!missing_varcov)
      nobs_varcov <- nrow(varcov)
      if (nobs_varcov == nobs_data) {
        if (nobs_omit) {
          varcov <- varcov[ind_model, ind_model]
        }
      }
      else {
        stopifnot(!missing_ids)
        ids_model <- ids[ind_model]
        stopifnot(!is.null(rownames(varcov)))
        ids_varcov <- rownames(varcov)
        stopifnot(all(ids_model %in% ids_varcov))
        ind <- which(ids_varcov %in% ids_model)
        varcov <- varcov[ind, ind]
      }
      evd <- eigen(varcov, symmetric = TRUE)
      vectors <- evd$vectors
      values <- evd$values
      transform <- vectors %*% diag(1/sqrt(values)) %*% 
        t(vectors)
    }
    else {
      nobs_transform <- nrow(transform)
      if (nobs_transform == nobs_data) {
        if (nobs_omit) {
          transform <- transform[ind_model, ind_model]
        }
      }
      else {
        stopifnot(!missing_ids)
        stopifnot(!is.null(rownames(transform)))
        ids_transform <- rownames(transform)
        ids_model <- ids[ind_model]
        stopifnot(all(ids_model %in% ids_transform))
        ind <- which(ids_transform %in% ids_model)
        transform <- transform[ind, ind]
      }
    }
  }
  toc(log = TRUE, quiet = TRUE)
  tic("pred")
  if (verbose > 0) {
    cat(" - creating `matlmPred`...\n")
  }
  pred <- matlm_pred3(pred, ind = ind_model, num_batches = num_batches, 
                      batch_size = batch_size, path_pred = path_pred)
  num_batches <- pred$num_batches
  stopifnot(pred_nrow(pred) == nobs_data)
  toc(log = TRUE, quiet = TRUE)
  tic("model matrices")
  if (verbose > 0) {
    cat(" - model matrices & response...\n")
  }
  y <- model.extract(model.frame(formula, data), "response")
  C <- model.matrix(formula, data)
  stopifnot(nrow(C) == nobs_model)
  stopifnot(length(y) == nobs_model)
  if (model == "interaction") {
    d_fct <- data[ind_model, int]
    stopifnot(class(d_fct) == "factor")
    stopifnot(nlevels(d_fct) == 2)
    d <- as.numeric(d_fct) - 1
    d_cname <- paste0(int, levels(d_fct)[2])
    stopifnot(d_cname %in% colnames(C))
  }
  N <- nobs_model
  toc(log = TRUE, quiet = TRUE)
  tic("apply transform")
  if (weighted) {
    y <- crossprod(transform, y)
    C <- crossprod(transform, C)
  }
  toc(log = TRUE, quiet = TRUE)
  tic("tests")
  y_orth <- matlm_orth(C, y)
  if (returnRespOrth) {
    return(as.numeric(y_orth))
  }
  y_sc <- matlm_scale(y_orth)
  matlm_batch_marginal <- function(batch, ind) {
    if (verbose > 1) {
      cat("  --  batch", batch, "/", num_batches, "\n")
    }
    X <- pred_batch(pred, batch, ind)
    if (weighted) {
      X <- crossprod(transform, X)
    }
    X_orth <- matlm_orth(C, X)
    if (returnPredOrth) {
      return(list(X_orth))
    }
    sd_X <- matlm_sd(X_orth)
    X_sc <- matlm_scale(X_orth, sd_X = sd_X)
    if (returnPredOrthSc) {
      return(list(X_sc))
    }
    k <- ncol(C) + 1
    r <- as.numeric(crossprod(X_sc, y_sc)/(N - 1))
    r2 <- r * r
    s <- sqrt((1 - r2)/(N - k))
    z <- r/s
    z2 <- z * z
    pvals <- pchisq(z2, df = 1, lower = FALSE)
    if (stats_full) {
      se <- s/sd_X
      b <- z * se
    }
    gc()
    if (!stats_full) {
      list(data_frame(predictor = colnames(X), zscore = z, 
                      pval = pvals))
    }
    else {
      list(data_frame(predictor = colnames(X), b = b, se = se, 
                      zscore = z, pval = pvals))
    }
  }
  matlm_batch_interaction <- function(batch, ind) {
    if (verbose > 1) {
      cat("  --  batch", batch, "/", num_batches, "\n")
    }
    X <- pred_batch(pred, batch, ind)
    M <- ncol(X)
    Xi <- diag(d) %*% X
    if (weighted) {
      X <- crossprod(transform, X)
      Xi <- crossprod(transform, Xi)
    }
    X_orth <- matlm_orth(C, X, Xi)
    if (returnPredOrth) {
      return(list(X_orth))
    }
    sd_X <- matlm_sd(X_orth)
    X_sc <- matlm_scale(X_orth, sd_X = sd_X)
    if (returnPredOrth) {
      return(list(X_sc))
    }
    Y <- tcrossprod(y, rep(1, M))
    Y_orth <- matlm_orth(C, X, Y)
    Y_sc <- matlm_scale(Y_orth)
    r <- apply(X_sc * Y_sc, 2, sum)/(N - 1)
    s <- r * sqrt((N - 2)/(1 - r * r))
    s2 <- s^2
    pvals <- pchisq(s2, df = 1, lower = FALSE)
    gc()
    list(data_frame(predictor = colnames(X), zscore = s, 
                    pval = pvals))
  }
  if (verbose > 0) {
    cat(" - computing association `tab`...\n")
  }
  if (cores > 1) {
    cl <- makeCluster(cores, type = "FORK")
    out <- switch(model, marginal = parSapply(cl, seq(1, 
                                                      num_batches), function(b) matlm_batch_marginal(b, 
                                                                                                     ind_model)), interaction = parSapply(cl, seq(1, num_batches), 
                                                                                                                                          function(b) matlm_batch_interaction(b, ind_model)), 
                  stop("switch by model (tab)"))
    stopCluster(cl)
  }
  else {
    out <- switch(model, marginal = sapply(seq(1, num_batches), 
                                           function(b) matlm_batch_marginal(b, ind_model)), 
                  interaction = sapply(seq(1, num_batches), function(b) matlm_batch_interaction(b, 
                                                                                                ind_model)), stop("switch by model (tab)"))
  }
  if (returnPredOrth | returnPredOrthSc) {
    mat <- do.call(cbind, out)
    return(mat)
  }
  tab <- bind_rows(out)
  if (verbose > 0) {
    cat(" - computing permutation `ptab`...\n")
  }
  pout <- switch(permutation, none = NULL, perm_x = {
    ind <- ind_model
    L <- num_perm
    N <- length(ind)
    P <- sample(ind, size = N * L, replace = TRUE) %>% matrix(nrow = N, 
                                                              ncol = L)
    pout <- sapply(seq(1:L), function(l) {
      if (verbose > 0) {
        cat("  --  permutation", l, "/", L, "\n")
      }
      if (cores > 1) {
        cl <- makeCluster(cores, type = "FORK")
        out <- switch(model, marginal = parSapply(cl, 
                                                  seq(1, num_batches), matlm_batch_marginal, 
                                                  ind = P[, l]), interaction = parSapply(cl, 
                                                                                         seq(1, num_batches), matlm_batch_interaction, 
                                                                                         ind = P[, l]), stop("switch by model (tab)"))
        stopCluster(cl)
      } else {
        out <- switch(model, marginal = sapply(seq(1, 
                                                   num_batches), matlm_batch_marginal, ind = P[, 
                                                                                               l]), interaction = sapply(seq(1, num_batches), 
                                                                                                                         matlm_batch_interaction, ind = P[, l]), stop("switch by model (tab)"))
      }
      tab <- bind_rows(out)
      tab <- mutate(tab, perm = l) %>% select(perm, everything())
      return(tab)
    }, simplify = FALSE)
  }, stop("switch by permutation"))
  ptab <- bind_rows(pout)
  toc(log = TRUE, quiet = TRUE)
  tic("return")
  out <- list(formula = formula, model = model, weighted = weighted, 
              nobs_data = nobs_data, nobs_model = nobs_model, nobs_omit = nobs_omit, 
              npred = pred_ncol(pred), permutation = permutation, tab = tab, 
              ptab = ptab, num_batches = num_batches, cores = cores)
  oldClass(out) <- c("matlmResults", "matlmList", oldClass(out))
  toc(log = TRUE, quiet = TRUE)
  toc(log = TRUE, quiet = TRUE)
  out$tictoc <- tic.log(format = FALSE)
  out$tictoc_formated <- tic.log(format = TRUE)
  tic.clearlog()
  tictoc_matlm <- tail(out$tictoc, 1)
  tictoc_elapsed <- tictoc_matlm[[1]]$toc - tictoc_matlm[[1]]$tic
  out$tictoc_elapsed <- unname(tictoc_elapsed)
  return(out)
}

#################### This is our own handcrafted function to run the GWAS using lme4qtl package. 

GWAS_cofac <- function(genotypes, trait, phenotype.name,snp.info, kinship, out.dir,categ,cofac) {
  
  dir.create(paste0(out.dir,"/frames"),recursive = T)
  dir.create(paste0(out.dir,"/plots"),recursive = T)
  dir.create(paste0(out.dir,"/herit"),recursive = T)
  phenotype <- toString(phenotype.name)
  letkin <- kinship
  usemat <- genotypes
  
  chr <- snp.info[,1]   
  pos <- snp.info[,2]
  
  use.trait <- trait
  names(use.trait) <- tolower(rownames(letkin))
  
  usemat.pruned <- usemat
  ID <- rownames(letkin) 
  mod <- lme4qtl::relmatLmer(use.trait ~ cofac + (1|ID), relmat = list(ID = letkin))
  
  herit.mod <- lme4qtl::VarProp(mod) 
  V <- lme4qtl::varcov(mod, idvar = "ID")
  V_thr <- V
  V_thr[abs(V) < 1e-10] <- 0
  decomp <- wlm::decompose_varcov(V, method = "evd", output = "all")
  W <- decomp$transform
  
  nnn <- rep(1,ncol(letkin))
  
  gassoc_gls <- matlm2(use.trait ~ cofac + nnn, nnn, pred =  usemat.pruned, ids = rownames(W), transform = W, batch_size = 2000, verbose = 2,cores = 1,stats_full = F)
  gassoc_gls$trait.name <- phenotype
  gassoc_gls$trait.noobs <- length(use.trait)
  
  gassoc_gls$tab$pval <- gassoc_gls$tab$pval*100000
  gassoc_gls$tab$zscore <- gassoc_gls$tab$zscore*100000
  gassoc_gls$chr <- chr
  gassoc_gls$pos <- pos
  save(herit.mod,file=paste(out.dir,"/herit/",phenotype,".out",sep=""))
  
  plot_frame <- data.frame(pval = gassoc_gls$tab$pval/100000,chr = gassoc_gls$chr,pos = gassoc_gls$pos/10^7,zscore = gassoc_gls$tab$zscore)
  save(plot_frame,file=paste(out.dir,"/frames/",phenotype,".out",sep=""))
  
  manh <- ggplot(plot_frame,aes(x = pos, y = -log10(pval),col = as.factor(chr)))+
    geom_point()+
    facet_grid(~as.factor(chr),space = "free_x",scale = "free_x")+
    xlab("Position Mb / 10")+
    ylab("-log10 p-value")+
    ggtitle(paste(phenotype,categ))+
    theme_minimal()+
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text.y = element_blank(),
          legend.position = "none",
          text = element_text(size = 35),
          axis.line = element_line(size = 0.2),
          axis.ticks= element_line(size = 0.2),
          panel.spacing.x = unit(6, "mm"),
          plot.title = element_text(hjust = 0.5)) +
    panel_border(color = "black",size = 0.2)+
    scale_color_viridis(discrete = T)
  
  png(filename = paste0(out.dir,"/plots/",phenotype,".png"),width = 2500,height = 750)
  print(manh)
  dev.off()
}

### Load in SNP matrix
load("ubec.out")

# Read in Shannon diversity
shannon <- read.table("shannon_fungi_final.txt")
rownames(shannon)[which(rownames(shannon)=="s215a")] <- "s215"
geno <- rownames(shannon)

div.root <- shannon[shannon$Habitat=="endosphere",1]
div.soil <- shannon[shannon$Habitat=="rhizosphere",1]

in.div.root <- as.numeric(gsub("r","",geno[grepl("r",geno)]))
in.div.soil <- as.numeric(gsub("s","",geno[grepl("s",geno)]))

### Calculate read depth per sample

OTU.gwas<-read.table("~/Documents/Lolium/LoliumSNPs/otutable_fungigwasfinal.txt")
TAX.gwas<-read.table("~/Documents/Lolium/LoliumSNPs/taxa.guilds.txt",sep = "\t",header = T)

### Remove OTU's which cannot be identified at the phylum level.
OTU.gwas <- OTU.gwas[!is.na(TAX.gwas$Phylum)&TAX.gwas$Phylum!=""&TAX.gwas$Phylum!="unidentified",]
TAX.gwas <- TAX.gwas[!is.na(TAX.gwas$Phylum)&TAX.gwas$Phylum!=""&TAX.gwas$Phylum!="unidentified",]

### Change "" to "unidentified"
TAX.gwas[TAX.gwas==""] <- "unidentified"

### Remove Saccharomyces contamination
OTU.gwas <- OTU.gwas[rownames(OTU.gwas)!="OTU422",]
TAX.gwas <- TAX.gwas[TAX.gwas$OTU!="OTU422",]
read.counts <- apply(OTU.gwas,2,sum)
names(read.counts)[which(names(read.counts)=="s215a")] <- "s215"

read.counts <- read.counts[which(names(read.counts)%in%geno)]
read.count.root <- read.counts[grepl("r",names(read.counts))]
read.count.soil <- read.counts[grepl("s",names(read.counts))]

### Process SNP matrix

snp.info <- ubec[,1:10]

ubec <- data.matrix(ubec[,11:249])
ubec[ubec == 9] <- 3
kinmat <- cov(ubec)

ubec.root <- ubec[,in.div.root]
ubec.soil <- ubec[,in.div.soil]

rm(ubec)
gc()

ubec.root <- t(ubec.root)
ubec.soil <- t(ubec.soil)

### Filter SNPs on MAF
threshold.root <- round(nrow(ubec.root)*0.05, digits = 0)
threshold.soil <- round(nrow(ubec.soil)*0.05, digits=0)

selc.root <- apply(ubec.root==1,2,sum) >= threshold.root & apply(ubec.root==3,2,sum) >= threshold.root
selc.soil <- apply(ubec.soil==1,2,sum) >= threshold.soil & apply(ubec.soil==3,2,sum) >= threshold.soil

ubec.root <- ubec.root[,selc.root]
ubec.soil <- ubec.soil[,selc.soil]

snp.info.root <- snp.info[selc.root,]
snp.info.soil <-snp.info[selc.soil,]

kinmat.root <- kinmat[in.div.root,in.div.root]
kinmat.soil <- kinmat[in.div.soil,in.div.soil]

### Call GWAS function
sapply(c(1), function(x) GWAS_cofac(genotypes = ubec.root, trait = div.root, phenotype.name = "Shannon depth", kinship=kinmat.root, snp.info = snp.info.root,out.dir = "GWAS_lolium_shannon_fungi_corrected_v2/root_depth",categ = "soil",cofac = read.count.root))
sapply(c(1), function(x) GWAS_cofac(genotypes = ubec.soil, trait = div.soil, phenotype.name = "Shannon depth", kinship=kinmat.soil, snp.info = snp.info.soil,out.dir = "GWAS_lolium_shannon_fungi_corrected_v2/soil_depth",categ = "soil",cofac = read.count.soil))




