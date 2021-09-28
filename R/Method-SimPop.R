makeSimPop <- function(n_mar, xo_freq, pos, ploidy, ...){
  object <- new("SimPop")
  object@lineage <- data.frame(gen = 1, n_family = 0, n_sibling = 0, n_ind = 0, crosstype = "founder")
  object@param <- data.frame(n_mar = n_mar, xo_freq = xo_freq, len = max(pos), ploidy = ploidy)
  if(length(pos) == 1){
    max_pos <- pos * 10^6
    physical <- sort(sample(1:max_pos, n_mar))
  } else{
    max_pos <- max(pos)
    physical <- sort(pos)
  }
  genetic <- xo_freq * diff(physical) / max_pos
  object@marker <- list(physical = physical, genetic = genetic)
  return(object)
}

setMethod("show",
          "SimPop",
          function(object){
            message('No of markers: ', getNumMar(object))
            message('Crossover frequency per chromosome: ', object@param$xo_freq)
            message('Chromosome length (Mb): ', object@param$len)
            message('Ploidy: ', getPloidy(object))

            if(length(object@lineage) == 0){
              message('\nNo population data.')
            } else {
              print(object@lineage)
            }
          })

setMethod("getHaplotype",
          "SimPop",
          function(object, gen, fam, sib, fmt, ...){
            if(gen == "last"){
              gen <- getLineage(object, gen)$gen
            }
            if(!is.numeric(gen)){
              stop('gen should be an integer or "last" to take the last generation.')
            }
            gen <- as.integer(gen)

            if(is.na(fam) & is.na(sib)){
              output <- object@haplotype[[gen]]

            } else if(is.na(sib)){
              fam <- as.integer(fam)
              output <- NULL
              for(fam_i in fam){
                output <- c(output, object@haplotype[[gen]][fam_i])
              }

            } else if(is.na(fam)){
              fam <- 1:getNumFam(object, gen)
              sib <- as.integer(sib)
              output <- NULL
              for(fam_i in fam){
                for(sib_i in sib){
                  output <- c(output, object@haplotype[[gen]][[fam_i]][sib])
                }
              }

            } else {
              if(fam == "all"){
                fam <- 1:getNumFam(object, gen)
              }
              if(sib == "all"){
                sib <- 1:getNumSib(object, gen)
              }
              fam <- as.integer(fam)
              sib <- as.integer(sib)
              output <- NULL
              for(fam_i in fam){
                for(sib_i in sib){
                  if(fmt == "array"){
                    output <- cbind(output, object@haplotype[[gen]][[fam_i]][[sib_i]])
                  } else if(fmt == "matrix"){
                    output <- rbind(output, object@haplotype[[gen]][[fam_i]][[sib_i]])
                  } else if(fmt == "list"){
                    output <- c(output, object@haplotype[[gen]][[fam_i]][sib_i])
                  }
                }
              }
              if(fmt == "array"){
                output <- array(output, c(2, getNumMar(object), length(fam) * length(sib)))
              }
            }
            return(output)
          })

setMethod("getGenotype",
          "SimPop",
          function(object, gen, fam, sib, fmt, ...){
            if(gen == "last"){
              gen <- getLineage(object, gen)$gen
            }
            if(!is.numeric(gen)){
              stop('gen should be an integer or "last" to take the last generation.')
            }
            gen <- as.integer(gen)

            if(is.na(fam) & is.na(sib)){
              output <- object@genotype[[gen]]

            } else if(is.na(sib)){
              fam <- as.integer(fam)
              output <- NULL
              for(fam_i in fam){
                output <- c(output, object@genotype[[gen]][fam_i])
              }

            } else if(is.na(fam)){
              fam <- 1:getNumFam(object, gen)
              sib <- as.integer(sib)
              output <- NULL
              for(fam_i in fam){
                for(sib_i in sib){
                  output <- c(output, object@genotype[[gen]][[fam_i]][sib])
                }
              }

            } else {
              if(fam == "all"){
                fam <- 1:getNumFam(object, gen)
              }
              if(sib == "all"){
                sib <- 1:getNumSib(object, gen)
              }
              fam <- as.integer(fam)
              sib <- as.integer(sib)
              output <- NULL
              for(fam_i in fam){
                for(sib_i in sib){
                  if(fmt == "array"){
                    output <- cbind(output, object@genotype[[gen]][[fam_i]][[sib_i]])
                  } else if(fmt == "matrix"){
                    output <- rbind(output, object@genotype[[gen]][[fam_i]][[sib_i]])
                  } else if(fmt == "list"){
                    output <- c(output, object@genotype[[gen]][[fam_i]][sib_i])
                  }
                }
              }
              if(fmt == "array"){
                output <- array(output, c(2, getNumMar(object), length(fam) * length(sib)))
              }
            }
            return(output)
          })

setMethod("getNumInd",
          "SimPop",
          function(object, gen, ...){
            if(gen == "last"){
              gen <- getLineage(object, gen)$gen
            }
            if(!is.numeric(gen)){
              stop('gen should be an integer or "last" to take the last generation.')
            }
            gen <- as.integer(gen)
            return(object@lineage$n_ind[gen])
          })

setMethod("getNumFam",
          "SimPop",
          function(object, gen, ...){
            if(gen == "last"){
              gen <- getLineage(object, gen)$gen
            }
            if(!is.numeric(gen)){
              stop('gen should be an integer or "last" to take the last generation.')
            }
            gen <- as.integer(gen)
            return(object@lineage$n_family[gen])
          })

setMethod("getNumSib",
          "SimPop",
          function(object, gen, ...){
            if(gen == "last"){
              gen <- getLineage(object, gen)$gen
            }
            if(!is.numeric(gen)){
              stop('gen should be an integer or "last" to take the last generation.')
            }
            gen <- as.integer(gen)
            return(object@lineage$n_sibling[gen])
          })

setMethod("getCrosstype",
          "SimPop",
          function(object, gen, ...){
            if(gen == "last"){
              gen <- getLineage(object, gen)$gen
            }
            if(!is.numeric(gen)){
              stop('gen should be an integer or "last" to take the last generation.')
            }
            gen <- as.integer(gen)
            return(object@lineage$crosstype[gen])
          })

setMethod("getLineage",
          "SimPop",
          function(object, gen, ...){
            if(gen == "last"){
              gen <- tail(object@lineage$gen, 1)
            }
            return(object@lineage[gen, ])
          })

setMethod("getPloidy",
          "SimPop",
          function(object, ...){
            return(object@param$ploidy)
          })

setMethod("getNumMar",
          "SimPop",
          function(object, ...){
            return(object@param$n_mar)
          })

setMethod("getXoFreq",
          "SimPop",
          function(object, ...){
            return(object@param$xo_freq)
          })

setMethod("getChrLen",
          "SimPop",
          function(object, ...){
            return(object@param$len)
          })

setMethod("getPhysPos",
          "SimPop",
          function(object, ...){
            return(object@marker$physical)
          })

setMethod("getGenPos",
          "SimPop",
          function(object, ...){
            return(object@marker$genetic)
          })

setMethod("getReads",
          "SimPop",
          function(object, gen, fmt, ...){
            gen <- getLineage(object, gen)$gen
            reads <- object@read[[gen]]
            if(fmt == "array"){
              reads <- do.call("cbind", reads)
              reads <- array(reads, c(2, getNumMar(object), getNumInd(object, gen)))
            }
            if(fmt == "matrix"){
              reads <- do.call("rbind", reads)
            }
            return(reads)
          })

setMethod("getID",
          "SimPop",
          function(object, gen, ...){
            return(object@id[[getLineage(object, gen)$gen]])
          })

setMethod("addFounder",
          "SimPop",
          function(object, allow_het, ...){
            n_founder <- getNumInd(object, gen = "last") + 1

            if(allow_het){
              ploidy <- getPloidy(object)
              chr_id <- seq(n_founder * 2 - ploidy + 1, length.out = getPloidy(object))
              haplotype <- rep(chr_id, each = getNumMar(object))
              haplotype <- matrix(haplotype, getPloidy(object), byrow = TRUE)

              if(n_founder == 1){
                genotype <- sample(2, getNumMar(object), replace = TRUE)
                genotype <- c(0, 1)[c(genotype, abs(genotype - 3))]

              } else {
                genotype <- sample(c(0, 1), getPloidy(object) * getNumMar(object), replace = TRUE)
              }
              genotype <- matrix(genotype, getPloidy(object), byrow = TRUE)

            } else {
              haplotype <- rep(n_founder, getPloidy(object) * getNumMar(object))
              haplotype <- matrix(haplotype, getPloidy(object), byrow = TRUE)
              if(n_founder == 1){
                genotype <- rep(0, getNumMar(object))
              } else if(n_founder == 2){
                genotype <- rep(1, getNumMar(object))
              } else {
                genotype <- sample(c(0, 1), getNumMar(object), replace = TRUE)
              }
              genotype <- matrix(rep(genotype, getPloidy(object)), getPloidy(object), byrow = TRUE)
            }

            if(n_founder == 1){
              object@haplotype <- list(list())
              object@genotype <- list(list())
              object@id <- list(character(0))
            }

            object@haplotype[[1]] <- c(object@haplotype[[1]], list(list(haplotype)))
            object@genotype[[1]] <- c(object@genotype[[1]], list(list(genotype)))
            object@id[[1]] <- c(object@id[[1]], paste0("Founder", n_founder))
            object@lineage$n_family <- object@lineage$n_family + 1
            object@lineage$n_sibling <- 1
            object@lineage$n_ind <- object@lineage$n_ind + 1
            object@read <- list(NULL)

            return(object)
          })

setMethod("getProgenies",
          "SimPop",
          function(object, crosstype, n_progeny, n_comb, randomComb, ...){
            crosstype <- match.arg(crosstype, c("selfing", "sibling", "pairing"))
            lastgen <- getLineage(object, gen = "last")

            if(crosstype == "selfing"){
              progenies <- .selfing(object, n_progeny)

            } else if(crosstype == "sibling"){
              if(length(n_comb) > 1){
                message('The first element of n_comb is used in "sibling" cross.')
              }

              message('n_comb is set at ', n_comb, '.\n',
                      n_progeny, ' progenies from ',
                      n_comb, ' combinations of sibling crosses for',
                      ' each sibling groups (or founders) of the last generation will be produced.')

              progenies <- .sibling(object, n_progeny, n_comb[1], randomComb)

            } else if(crosstype == "pairing"){

              if(length(n_comb) > 2){
                message('The two element of n_comb is used in "sibling" cross.')
              } else if(length(n_comb) < 2){
                stop('Please specify a vector with two integers to n_comb for "pairing" cross')
              }

              progenies <- .pairing(object, n_progeny, n_comb, randomComb)
            }
            object@haplotype <- c(object@haplotype, list(progenies$hap))
            object@genotype <- c(object@genotype, list(progenies$geno))
            object@id <- c(object@id, list(progenies$id))
            gen <- data.frame(gen = lastgen$gen + 1,
                              n_family = length(progenies$hap),
                              n_sibling = n_progeny,
                              n_ind = length(progenies$hap) * n_progeny,
                              crosstype = crosstype)
            object@lineage <- rbind(object@lineage, gen)
            object@read <- c(object@read, list(NULL))
            return(object)
          })

.selfing <- function(object, n_progeny){
  hap <- getHaplotype(object, "last", fam = "all", sib = "all", fmt = "list")
  geno <- getGenotype(object, "last", fam = "all", sib = "all", fmt = "list")
  parents_id <- expand.grid(list(1:getNumFam(object, "last"), 1:getNumSib(object, "last")))
  parents_id <- parents_id[order(parents_id[, 1]), ]
  parents_id <- apply(parents_id, 1, paste, collapse = "_")
  message(length(hap), " individuals are going to be subjected to selfing.")

  progenies_gen <- getLineage(object, "last")$gen + 1
  progenies_hap <- NULL
  progenies_geno <- NULL
  progenies_id <- NULL
  for(ind in 1:length(hap)){
    siblings_hap <- list()
    siblings_geno <- list()
    for(j in 1:n_progeny){
      paternal <- .getGamete(hap[[ind]],
                             geno[[ind]],
                             getNumMar(object),
                             getPloidy(object),
                             getGenPos(object))
      maternal <- .getGamete(hap[[ind]],
                             geno[[ind]],
                             getNumMar(object),
                             getPloidy(object),
                             getGenPos(object))
      genome_hap <- matrix(c(paternal$hap, maternal$hap), nrow = 2, byrow = TRUE)
      genome_geno <- matrix(c(paternal$geno, maternal$geno), nrow = 2, byrow = TRUE)
      siblings_hap <- c(siblings_hap, list(genome_hap))
      siblings_geno <- c(siblings_geno, list(genome_geno))
      progenies_id <- c(progenies_id,
                        paste0("G", progenies_gen, "_",
                               parents_id[ind], "x",
                               parents_id[ind], "_",
                               j))
    }
    progenies_hap <- c(progenies_hap, list(siblings_hap))
    progenies_geno <- c(progenies_geno, list(siblings_geno))
  }
  return(list(hap = progenies_hap, geno = progenies_geno, id = progenies_id))
}

.sibling <- function(object, n_progeny, n_comb, randomComb){
  n_fam <- getNumFam(object, gen = "last")
  message(n_fam, " families (groups) are going to be subjected to sibling crosses.")

  progenies_gen <- getLineage(object, "last")$gen + 1
  progenies_hap <- NULL
  progenies_geno <- NULL
  progenies_id <- NULL
  for(i_fam in 1:n_fam){
    n_sib <- getNumSib(object, gen = "last")
    if(n_sib < n_comb * 2){
      stop('Each sibling group have ',
           n_sib, " individuals.\nBut you required ",
           n_comb * 2, " combinations.")
    }
    if(randomComb){
      combi <- matrix(sample(1:n_sib, n_comb * 2), ncol = 2)
    } else {
      combi <- matrix(1:(n_comb * 2), ncol = 2)
    }

    for(i_combi in 1:nrow(combi)){
      siblings_hap <- list()
      siblings_geno <- list()
      for(j in 1:n_progeny){
        paternal <- .getGamete(getHaplotype(object,
                                            "last",
                                            fam = i_fam,
                                            sib = combi[i_combi, 1],
                                            fmt = "matrix"),
                               getGenotype(object,
                                           "last",
                                           fam = i_fam,
                                           sib = combi[i_combi, 1],
                                           fmt = "matrix"),
                               getNumMar(object),
                               getPloidy(object),
                               getGenPos(object))
        maternal <- .getGamete(getHaplotype(object,
                                            "last",
                                            fam = i_fam,
                                            sib = combi[i_combi, 2],
                                            fmt = "matrix"),
                               getGenotype(object,
                                           "last",
                                           fam = i_fam,
                                           sib = combi[i_combi, 2],
                                           fmt = "matrix"),
                               getNumMar(object),
                               getPloidy(object),
                               getGenPos(object))
        genome_hap <- matrix(c(paternal$hap, maternal$hap), nrow = 2, byrow = TRUE)
        genome_geno <- matrix(c(paternal$geno, maternal$geno), nrow = 2, byrow = TRUE)
        siblings_hap <- c(siblings_hap, list(genome_hap))
        siblings_geno <- c(siblings_geno, list(genome_geno))
        progenies_id <- c(progenies_id,
                          paste0("G", progenies_gen, "_",
                                 i_fam, "_", combi[i_combi, 1], "x",
                                 i_fam, "_", combi[i_combi, 2], "_",
                                 j))
      }
      progenies_hap <- c(progenies_hap, list(siblings_hap))
      progenies_geno <- c(progenies_geno, list(siblings_geno))
    }

  }
  return(list(hap = progenies_hap, geno = progenies_geno, id = progenies_id))
}


.pairing <- function(object, n_progeny, n_comb, randomComb){
  n_fam <- getNumFam(object, "last")
  if(n_fam < n_comb[1] * 2){
    stop('The last generation has ',
         n_fam, " sibling groups (or founders).\nBut you required ",
         n_comb[1] * 2, " combinations.")
  }
  message(n_fam, " sibling groups (or founders) are going to be subjected to pairing crosses.")

  if(randomComb){
    combi1 <- matrix(sample(1:n_fam, n_comb[1] * 2), ncol = 2)
  } else {
    combi1 <- matrix(1:(n_comb[1] * 2), ncol = 2, byrow = TRUE)
  }

  progenies_gen <- getLineage(object, "last")$gen + 1
  progenies_hap <- NULL
  progenies_geno <- NULL
  progenies_id <- NULL
  for(i_combi1 in 1:nrow(combi1)){
    n_sib <- getNumSib(object, "last")
    if(n_sib < n_comb[2]){
      stop('Each sibling group have ',
           n_sib, " individuals.\nBut you required ",
           n_comb[2], " combinations.")
    }
    if(randomComb){
      combi2 <- cbind(matrix(sample(1:n_sib, n_comb[2]), ncol = 1),
                      matrix(sample(1:n_sib, n_comb[2]), ncol = 1))
    } else {
      combi2 <- cbind(matrix(1:n_comb[2], ncol = 1),
                      matrix(1:n_comb[2], ncol = 1))
    }

    for(i_combi2 in 1:nrow(combi2)){
      siblings_hap <- list()
      siblings_geno <- list()
      for(j in 1:n_progeny){
        paternal <- .getGamete(getHaplotype(object,
                                            "last",
                                            fam = combi1[i_combi1, 1],
                                            sib = combi2[i_combi2, 1],
                                            fmt = "matrix"),
                               getGenotype(object,
                                           "last",
                                           fam = combi1[i_combi1, 1],
                                           sib = combi2[i_combi2, 1],
                                           fmt = "matrix"),
                               getNumMar(object),
                               getPloidy(object),
                               getGenPos(object))
        maternal <- .getGamete(getHaplotype(object,
                                            "last",
                                            fam = combi1[i_combi1, 2],
                                            sib = combi2[i_combi2, 2],
                                            fmt = "matrix"),
                               getGenotype(object,
                                           "last",
                                           fam = combi1[i_combi1, 2],
                                           sib = combi2[i_combi2, 2],
                                           fmt = "matrix"),
                               getNumMar(object),
                               getPloidy(object),
                               getGenPos(object))
        genome_hap <- matrix(c(paternal$hap, maternal$hap), nrow = 2, byrow = TRUE)
        genome_geno <- matrix(c(paternal$geno, maternal$geno), nrow = 2, byrow = TRUE)
        siblings_hap <- c(siblings_hap, list(genome_hap))
        siblings_geno <- c(siblings_geno, list(genome_geno))
        progenies_id <- c(progenies_id,
                          paste0("G", progenies_gen, "_",
                                 combi1[i_combi1, 1], "_", combi2[i_combi2, 1], "x",
                                 combi1[i_combi1, 2], "_", combi2[i_combi2, 2], "_",
                                 j))
      }
      progenies_hap <- c(progenies_hap, list(siblings_hap))
      progenies_geno <- c(progenies_geno, list(siblings_geno))
    }

  }
  return(list(hap = progenies_hap, geno = progenies_geno, id = progenies_id))
}

.getGamete <- function(hap, geno, n_mar, ploidy, genetic){
  xo_pos <- sapply(genetic, function(x) rpois(1, x)) != 0
  if(sum(xo_pos) != 0){
    xo_pos <- which(xo_pos)
    xo_pos <- c(0, xo_pos, n_mar)
    xo_interval <- diff(xo_pos)
    start_chr <- sample(1:ploidy, 1)
    chr_index <- rep(start_chr, xo_interval[1])
    for(k in 2:length(xo_interval)){
      prev_chr <- tail(chr_index, 1)
      possible_chr <- which(1:ploidy != prev_chr)
      if(length(possible_chr) == 1){
        next_chr <- possible_chr
      } else {
        next_chr <- sample(possible_chr, 1, replace = TRUE)
      }
      chr_index <- c(chr_index, rep(next_chr, xo_interval[k]))
    }

    i_offset <- c(0, cumsum(rep(ploidy, n_mar - 1)))
    chr_index <- chr_index + i_offset
    gamete_hap <- hap[chr_index]
    gamete_geno <- geno[chr_index]

  } else {
    start_chr <- sample(1:ploidy, 1, replace = TRUE)
    gamete_hap <- hap[start_chr, ]
    gamete_geno <- geno[start_chr, ]
  }

  return(list(hap = gamete_hap, geno = gamete_geno))
}


setMethod("selectSamples",
          "SimPop",
          function(object, n_samples, perFam, ...){
            n_fam <- getNumFam(object, "last")
            n_sib <- getNumSib(object, "last")
            if(perFam){
              sampled <- sample(1:n_sib, n_fam * n_samples, replace = T)
              sampled <- cbind(rep(1:n_fam, each = n_samples), sampled)
            } else {
              combi <- expand.grid(1:n_fam, 1:n_sib)
              sampled <- combi[sample(1:nrow(combi), n_samples), ]
            }
            hap <- NULL
            geno <- NULL
            id <- NULL
            for(i in 1:nrow(sampled)){
              hap <- c(hap, getHaplotype(object, "last", fam = sampled[i, 1], sim = sampled[i, 2]))
              geno <- c(geno, getGenotype(object, "last", fam = sampled[i, 1], sim = sampled[i, 2]))
              id <- c(id, getID(object, "last")[(sampled[i, 1] - 1)* n_sib + sampled[i, 2]])
            }
            object@haplotype <- c(object@haplotype, list(hap))
            object@genotype <- c(object@genotype, list(geno))
            object@id <- c(object@id, list(id))
            object@lineage <- rbind(object@lineage,
                                    data.frame(gen = getLineage(object, "last")$gen + 1,
                                               n_family = length(id),
                                               n_sibling = 1,
                                               n_ind = length(id),
                                               crosstype = "selected"))
            return(object)
          })

.getDistDP <- function(dp_dist, n_mar){
  if(length(dp_dist) == n_mar){
    return(dp_dist)
  } else {
    if(is.null(dp_dist)){
      lambda <- 1/4
    } else if(length(dp_dist) == 1 & is.numeric(dp_dist)){
      lambda <- dp_dist
    }
    dp_dist <- rexp(n_mar, lambda)
    dp_dist <- dp_dist / sum(dp_dist)
    return(dp_dist)
  }
}

.getDistAD <- function(ad_dist, n_mar){
  if(length(ad_dist) == n_mar){
    return(ad_dist)
  } else {
    if(is.null(ad_dist)){
      ad_mean <- 0.5
      ad_sd <- 0
    } else if(length(ad_dist) == 2 & is.numeric(ad_dist)){
      ad_mean <- ad_dist[1]
      ad_sd <- ad_dist[2]
    }
    ad_dist <- rnorm(n_mar, ad_mean, ad_sd)
    while(TRUE){
      ad_dist <- c(ad_dist, rnorm(n_mar, ad_mean, ad_sd))
      ad_dist <- ad_dist[ad_dist > 0 & ad_dist < 1]
      if(length(ad_dist) >= n_mar){
        ad_dist <- ad_dist[1:n_mar]
        break
      }
    }
    return(ad_dist)
  }
}

.getMissmap <- function(mismap, n_mar){
  if(is.data.frame(mismap) | is.matrix(mismap)){
    if(all(apply(mismap, 2, class) == "numeric")){
      check1 <- ncol(mismap) == 2
      check2 <- all(mismap[, 1] >= 0 & mismap[, 1] <= 1)
      check3 <- all(mismap[, 2] - as.integer(mismap[, 2]) == 0)
      check4 <- nrow(mismap) == n_mar
      check <- check1 & check2 & check3 & check4
      if(check){
        return(mismap)
      }
    }
  }
  if(is.vector(mismap)){
    check1 <- mismap >= 0 & mismap <= 1
    if(length(mismap) == 1 & check1){
      mismap <- 1 / mismap
      mismap <- rexp(n_mar, mismap)
      mismap[mismap > 1] <- 1
    }
    if(length(mismap) == n_mar & check1){
      allele <- sample(c(1, 2, 3), n_mar, TRUE)
      mismap <- data.frame(mismap, allele)
      return(mismap)
    }
  }
  stop('Please specify valid object to "mismap".')
}

.insertMissmap <- function(genotype, mismap, n_mar){
  mismapped <- sapply(1:n_mar, function(x) rbinom(1, 1, mismap$mismap))
  allele <- mismap$allele[mismapped == 1]
  mismapped_allele <- t(matrix(c(1, 0, 1, 1, 0, 1), 3, 2, T)[allele, ])
  mismapped_genotype <- genotype[, mismapped == 1] + mismapped_allele
  mismapped_genotype[mismapped_genotype > 1] <- 1
  genotype[, mismapped == 1] <- mismapped_genotype
  return(genotype)
}

.setProb <- function(genotype, dp_dist, ad_dist, seq_error, joint){
  genocall <- colSums(genotype) + 1
  if(joint){
    ad_dist <- rbind(ad_dist, 1 - ad_dist)
    ad_dist[, genocall == 1] <- rep(ad_dist[1, genocall == 1], each = 2)
    ad_dist[, genocall == 3] <- rep(ad_dist[2, genocall == 3], each = 2)
    ref_prob <- c((1 - seq_error), 0.5, (seq_error))[genocall] * dp_dist
    alt_prob <- c(seq_error, 0.5, (1 - seq_error))[genocall] * dp_dist
    read_prob <- rbind(ref_prob, alt_prob) * ad_dist

  } else {
    ad_dist <- rbind(ad_dist, 1 - ad_dist)
    ref_prob <- c((1 - seq_error), 0.5, (seq_error))[genocall]
    alt_prob <- c(seq_error, 0.5, (1 - seq_error))[genocall]
    read_prob <- rbind(ref_prob, alt_prob)
    read_prob[, genocall == 2] <- ad_dist[, genocall == 2]
    read_prob <- read_prob * rbind(dp_dist, dp_dist)
  }
  return(read_prob)
}

.allocateReads <- function(prob, n_mar, total_read, joint){
  if(joint){
    n_row <- nrow(prob)
    n_sample <- n_row / 2
    n_index <- 2*n_mar*n_sample
    reads <- sample(1:n_index, total_read * n_sample, TRUE, prob)
    reads <- table(factor(reads, levels = 1:n_index))
    reads <- matrix(reads, nrow = n_row)
    reads <- tapply(1:n_row, rep(1:n_sample, each = 2), function(x) return(reads[x, ]))

  } else {

    n_index <- 2*n_mar
    reads <- lapply(prob, function(x){
      tmp <- sample(1:n_index, total_read, TRUE, x)
      tmp <- table(factor(tmp, levels = 1:n_index))
      tmp <- matrix(tmp, nrow = 2)
      return(tmp)
    })
  }
  return(reads)
}

setMethod("simRead",
          "SimPop",
          function(object, gen, total_read, dp_dist, ad_dist, seq_error, mismap, joint, ...){
            n_mar <- getNumMar(object)
            dp_dist <- .getDistDP(dp_dist, n_mar)
            ad_dist <- .getDistAD(ad_dist, n_mar)
            mismap <- .getMissmap(mismap, n_mar)

            for(gen_i in gen){
              if(gen_i != "last"){
                gen_i <- as.integer(gen_i)
              }
              geno <- getGenotype(object, gen_i, fam = "all", sib = "all", fmt = "list")
              geno <- lapply(geno, .insertMissmap, mismap, n_mar)
              prob <- lapply(geno, .setProb, dp_dist, ad_dist, seq_error, joint)

              if(joint){
                reads <- .allocateReads(do.call("rbind", prob), n_mar, total_read, joint)

              } else {
                reads <- .allocateReads(prob, n_mar, total_read, joint)
              }

              object@read[[getLineage(object, gen_i)$gen]] <- reads
              object@dist <- data.frame(dp_dist = dp_dist,
                                        ad_dist = ad_dist,
                                        mismap_rate = mismap[, 1],
                                        mismap_type = mismap[, 2])
            }

            return(object)
          })

setMethod("writeVCF",
          "SimPop",
          function(object, gen, vcf_fn, ...){
            ind_name <- NULL
            reads <- NULL
            for(gen_i in gen){
              if(gen_i != "last"){
                gen_i <- as.integer(gen_i)
              }
              ind_name <- c(ind_name, getID(object, gen_i))
              read_i <- getReads(object, gen_i)
              if(is.null(read_i)){
                stop('No read data of the generation ', gen_i, '.\nRun simReads for it.')
              }
              reads <- c(reads, read_i)
            }

            geno <- sapply(reads, function(x){
              gt1 <- as.numeric(x[1, ] == 0)
              gt2 <- as.numeric(x[2, ] > 0)
              gtn <- colSums(x) == 0
              gt <- apply(cbind(gt1, gt2), 1, sort)
              gt <- apply(gt, 2, paste, collapse = "/")
              gt[colSums(x) == 0] <- "./."
              ad <- apply(x, 2, paste, collapse = ",")
              dp <- colSums(x)
              return(paste(gt, ad, dp, sep = ":"))
            })

            chr <- 1
            pos <- getPhysPos(object)
            marker_id <- 1:getNumMar(object)
            ref <- "G"
            alt <- "A"

            out_matrix <- cbind(chr,
                                pos,
                                marker_id,
                                ref,
                                alt,
                                ".",
                                "PASS",
                                ".",
                                "GT:AD:DP",
                                geno)

            header <- c("#CHROM",
                        "POS",
                        "ID",
                        "REF",
                        "ALT",
                        "QUAL",
                        "FILTER",
                        "INFO",
                        "FORMAT",
                        ind_name
            )

            out_matrix <- rbind(header, out_matrix)

            if(vcf_fn == ""){
              vcf_fn <- paste("simpop",
                              format(Sys.time(), "%Y%m%d_%H%M%S"),
                              sep = "_")
            }
            meta <- c('##fileformat=VCFv4.0',
                      '##FILTER=<ID=PASS,Description="All filters passed">',
                      '##FILTER=<ID=PASS,Description="All filters passed">',
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                      '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the reference and alternate alleles in the order listed">',
                      '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (only filtered reads used for calling)">')
            write.table(x = meta,
                        file = vcf_fn,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE,
                        quote = FALSE,
                        append = FALSE)

            write.table(x = out_matrix,
                        file = vcf_fn,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE,
                        quote = FALSE,
                        append = TRUE)
            invisible(TRUE)
          })
