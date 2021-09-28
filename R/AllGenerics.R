#' @export


setGeneric("getHaplotype", function(object, gen = "last", fam = NA, sib = NA, fmt = "list", ...) standardGeneric("getHaplotype"))
setGeneric("getGenotype", function(object, gen = "last", fam = NA, sib = NA, fmt = "list", ...) standardGeneric("getGenotype"))
setGeneric("getNumInd", function(object, gen = "last", ...) standardGeneric("getNumInd"))
setGeneric("getNumFam", function(object, gen = "last", ...) standardGeneric("getNumFam"))
setGeneric("getNumSib", function(object, gen = "last", ...) standardGeneric("getNumSib"))
setGeneric("getCrosstype", function(object, gen = "last", ...) standardGeneric("getCrosstype"))
setGeneric("getLineage", function(object, gen = "last", ...) standardGeneric("getLineage"))
setGeneric("getPloidy", function(object, ...) standardGeneric("getPloidy"))
setGeneric("getNumMar", function(object, ...) standardGeneric("getNumMar"))
setGeneric("getXoFreq", function(object, ...) standardGeneric("getXoFreq"))
setGeneric("getChrLen", function(object, ...) standardGeneric("getChrLen"))
setGeneric("getPhysPos", function(object, ...) standardGeneric("getPhysPos"))
setGeneric("getGenPos", function(object, ...) standardGeneric("getGenPos"))
setGeneric("getReads", function(object, gen, fmt = "list", ...) standardGeneric("getReads"))
setGeneric("getID", function(object, gen = "last", ...) standardGeneric("getID"))


setGeneric("addFounder", function(object, allow_het = FALSE, ...) standardGeneric("addFounder"))
setGeneric("getProgenies", function(object, crosstype = "selfing", n_progeny = 1, n_comb = 1, randomComb = FALSE, ...) standardGeneric("getProgenies"))
setGeneric("selectSamples", function(object, n_samples, perFam = FALSE, ...) standardGeneric("selectSamples"))
setGeneric("simRead", function(object, gen = "last", total_read = 1e6, dp_dist = NULL, ad_dist = NULL, seq_error = 0.0025, mismap = 0, joint = TRUE, ...) standardGeneric("simRead"))
setGeneric("writeVCF", function(object, gen, vcf_fn, ...) standardGeneric("writeVCF"))
