#' #' Parameter specification of PCM model
#' #' @details replaced by PCMSpecifyParameters
#' #' @param model a PCM model object
#' #' @return a list
#' #' @export
#' PCMSpecifyParams <- function(model, ...) {
#'   UseMethod("PCMSpecifyParams", model)
#' }
#'
#' #' @export
#' PCMSpecifyParams.PCM <- function(model, ...) {
#'   list()
#' }
#'
#' #' Number of free parameters describing a PCM
#' #' @param model a PCM object
#' #' @param countRegimeChanges logical indicating if regime changes should be counted.
#' #' If TRUE, the default implementation would add \code{PCMNumRegimes(model) - 1}.
#' #' Default FALSE.
#' #' @param countModelTypes logical indicating whether the model type should be
#' #'  counted. If TRUE the default implementation will add +1 only if there are more than
#' #'  one modelTypes (\code{length(attr(model, "modelTypes", exact = TRUE)) > 1}),
#' #'  assuming that all regimes are regimes of the same model type (e.g. OU). The
#' #'  implementation for MRG models will add +1 for every regime if there are more than
#' #'  one modelTypes. Default FALSE.
#' #' @param ... other arguments (possible future use)
#' #' @return an integer
#' #' @export
#' PCMNumParams <- function(model, countRegimeChanges = FALSE, countModelTypes = FALSE, ...) {
#'   UseMethod("PCMNumParams", model)
#' }
#'
#' #' @export
#' PCMNumParams.PCM <- function(model, countRegimeChanges = FALSE, countModelTypes = FALSE,  ...) {
#'   k <- attr(model, "k")
#'   R <- length(attr(model, "regimes"))
#'   specParams <- attr(model, "specParams")
#'
#'   vecParamSizes <- c(rep = 1, full = k, fixed = 0)
#'   matParamSizes <- c(diag1 = 1, diag = k, upper.tri = .5*k*(k-1), upper.tri.diag = .5*k*(k+1),
#'                      lower.tri = .5*k*(k-1), lower.tri.diag = .5*k*(k+1), symmetric = .5*k*(k+1),
#'                      full=k*k, fixed = 0)
#'   p <- 0
#'   for(name in names(specParams)) {
#'     type <- specParams[[name]]$type
#'     indices <- specParams[[name]]$indices
#'
#'     if(type[1] == "model") {
#'       p <- p + PCMNumParams(specParams[[name]]$default, ...)
#'     } else if(type[1] == "gvector") {
#'       if(type[2] == "custom") {
#'         ind <- indices(p, k)
#'         p <- p + length(unique(ind[ind>p]))
#'       } else {
#'         p <- p + vecParamSizes[type[2]]
#'       }
#'     } else if(type[1] == "gmatrix") {
#'       if(type[2] == "custom") {
#'         ind <- indices(p, k)
#'         p <- p + length(unique(ind[ind>p]))
#'       } else {
#'         p <- p + matParamSizes[type[2]]
#'       }
#'     } else if(type[1] == "vector") {
#'       if(type[2] == "custom") {
#'         for(r in 1:R) {
#'           ind <- indices(p, k)
#'           p <- p + length(unique(ind[ind>p]))
#'         }
#'       } else {
#'         p <- p + R*vecParamSizes[type[2]]
#'       }
#'     } else if(type[1] == "matrix") {
#'       if(type[2] == "custom") {
#'         for(r in 1:R) {
#'           ind <- indices(p, k)
#'           p <- p + length(unique(ind[ind>p]))
#'         }
#'       } else {
#'         p <- p + R*matParamSizes[type[2]]
#'       }
#'     }
#'   }
#'   if(countRegimeChanges) {
#'     # we don't count the root as a parameters, tha'ts why we substract one.
#'     p <- p + PCMNumRegimes.PCM(model) - 1
#'   }
#'   if(countModelTypes) {
#'     # assume that all regimes have the same model-type. If there is only one
#'     # model type than this is not counted as a parameter.
#'     if(length(attr(model, "modelTypes", exact=TRUE)) > 1)
#'       p <- p + 1
#'   }
#'   unname(p)
#' }
#'
#' #' Load/store a vector parameter from/to a vector of all parameters in a model.
#' #'
#' #' @details This function has both, a returned value and side effects. By default the function
#' #' loads elements from vecParams into v. This behavior is reversed if the argument load
#' #' is set to FALSE.
#' #'
#' #' @param v numeric k-vector. If load==TRUE (default) this has to be a vector-object in
#' #' the calling environment (not a temporary object such as the result from an algebraic expression).
#' #' @param vecParams numeric vector containing all parameters in the model. If load==FALSE, this
#' #' has to be a vector in the calling environment (not a temporary object such as the result from
#' #' an algebraic expression).
#' #' @param offset integer denoting the offset in vecParams (0 means start from vecParams[1] onwards).
#' #' @param k integer denoting the dimensionality (length) of the resulting vector.
#' #' @param type a character string indicating the type of loading. possible values are:
#' #' \describe{
#' #' \item{"rep"}{repeat \code{vecParams[offset+1]} k times}
#' #' \item{"full"}{copy \code{vecParams[offset + (1:k)]}}
#' #' \item{"custom"}{use a custom template vector which elements specified by mask are
#' #' to be replaced with \code{vecParams[indices(offset, k)]}}
#' #' }.
#' #' @param mask a logical vecgtor pointing to elements in template to be replaced. used only with type=="custom".
#' #' @param indices a function of the form function(offset, k) returning an integer vector of
#' #' length equal to the number of TRUE elements in mask, indicating the position in vecParams to
#' #' be copyed from. used only with type=="custom".
#' #' @param load a logical indicating if loading from or storing to vecParams should be done.
#' #' @return an integer denoting the number of elemnents read or written from/to vecParams.
#' #' In the case of type=="custom",
#' #' @export
#' PCMLoadVectorParameter <- function(
#'   v, vecParams, offset, k,
#'   type = c("rep", "full", "custom"),
#'   mask = rep(TRUE, k), indices = function(offset, k) offset + (1:k),
#'   load = TRUE) {
#'
#'   if(load) {
#'     "%op%" <- `<-`
#'     maskRep <- 1:k
#'   } else {
#'     "%op%"<- function(a,b) eval(substitute(b<-a), parent.frame())
#'     maskRep <- 1
#'   }
#'
#'   if(type[1] == "rep") {
#'     num <- 1
#'     eval(substitute(v[maskRep] %op% vecParams[offset + 1]), parent.frame())
#'   } else if(type[1] == "full") {
#'     num <- k
#'     eval(substitute(v[] %op% vecParams[offset + (1:k)]), parent.frame())
#'   } else if(type[1] == "custom") {
#'     if(is.function(indices)) {
#'       ind <- indices(offset, k)
#'       num <- length(unique(ind[ind>offset]))
#'       maskCustom <- mask
#'       eval(substitute(v[maskCustom] %op% vecParams[ind]), parent.frame())
#'     } else {
#'       stop("ERR:020a1:PCMBase:PCM.R:PCMLoadVectorParameter:: indices should be a
#'            function(offset, k) returning an integer vector.")
#'     }
#'     } else if(type[1] == "fixed") {
#'       num <- 0
#'   } else {
#'     stop(paste0("ERR:020a2:PCMBase:PCM.R:PCMLoadVectorParameter:: type ", type[1], " not recognized."))
#'   }
#'   num
#' }
#'
#' #' Load/store a matrix parameter from/to a vector of all parameters in a model.
#' #'
#' #' @details This function has both, a returned value and side effects. By default the function
#' #' loads elements from vecParams into m. This behavior is reversed if the argument load
#' #' is set to FALSE.
#' #'
#' #' @param m numeric k x k matrix. If load==TRUE (default) this has to be a matrix-object in
#' #' the calling environment (not a temporary object such as the result from an algebraic expression).
#' #' @param vecParams numeric vector containing all parameters in the model. If load==FALSE, this
#' #' has to be a vector in the calling environment (not a temporary object such as the result from
#' #' an algebraic expression).
#' #' @param k integer denoting the dimensionality (number of rows and colums) of the resulting matrix.
#' #' @param offset integer denoting the offset in vecParams (0 means start from vecParams[1] onwards).
#' #' @param type a character string indicating the type of loading. possible values are:
#' #' \describe{
#' #' \item{"diag1"}{create a diagonal matrix with all diagonal elements repeating \code{vecParams[offset+1]}.}
#' #' \item{"diag"}{create a diagonal matrix with a diagonal copied from \code{vecParams[offset + (1:k)]}.}
#' #' \item{"upper.tri"}{load an upper triangular matrix with zero diagonal.}
#' #' \item{"upper.tri.diag"}{load an upper triangular matrix including diagonal.}
#' #' \item{"lower.tri"}{load a lower triangular matrix with zero diagonal.}
#' #' \item{"lower.tri.diag"}{load a lower triangular matrix including diagonal.}
#' #' \item{"symmetric"}{load a symmetric matrix. Only the elements of the upper triangle and diagonal are
#' #' loaded from vecParams.}
#' #' \item{"full"}{load a full.}
#' #' \item{"custom"}{use a custom template matrix which's elements specified by mask are
#' #' to be replaced with \code{vecParams[offset + indices]}}.
#' #' }
#' #' @param mask a logical k x k matrix pointing to elements in template to be replaced. used only with type=="custom".
#' #' @param indices a function of the form function(offset, k) returning an integer vector of
#' #' length equal to the number of TRUE elements in mask, indicating the position in vecParams to
#' #' be copyed from. used only with type=="custom".
#' #' @param load a logical indicating if loading from or storing to vecParams should be done.
#' #'
#' #' @return an integer denoting the number of elemnents read or written from/to vecParams.
#' #' In the case of type=="custom",
#' #' the number of indices bigger than offset returned by the function indices(offset, k).
#' #' @export
#' PCMLoadMatrixParameter <- function(
#'   m, vecParams, offset, k,
#'   type = c("diag1", "diag", "upper.tri", "upper.tri.diag", "lower.tri", "lower.tri.diag", "symmetric", "full", "custom", "fixed"),
#'   mask = matrix(TRUE, k, k),
#'   indices = function(offset, k) offset + (1:(k*k)),
#'   load = TRUE) {
#'
#'   if(load) {
#'     "%op%" <- `<-`
#'     maskDiag1 <- 1:k
#'   } else {
#'     "%op%" <- function(a,b) eval(substitute(b<-a), parent.frame())
#'     maskDiag1 <- 1
#'   }
#'
#'
#'   if(length(dim(m)) == 0) {
#'     if(k == 1) {
#'       # m is a scalar
#'       if(type[1] == "diag1") {
#'         num <- 1
#'         eval(substitute(m[maskDiag1] %op% vecParams[offset + (1:num)]), parent.frame())
#'       } else if(type[1] == "diag") {
#'         num <- k
#'         eval(substitute(m[] %op% vecParams[offset + (1:num)]), parent.frame())
#'       } else if(type[1] == "upper.tri") {
#'         num <- k*(k-1)/2
#'         # nothing else to do
#'       } else if(type[1] == "upper.tri.diag") {
#'         num <- k*(k+1)/2
#'         eval(substitute(m[upper.tri(m, diag=TRUE)] %op% vecParams[offset + (1:num)]), parent.frame())
#'       } else if(type[1] == "lower.tri") {
#'         num <- k*(k-1)/2
#'         # nothing else to do
#'       } else if(type[1] == "lower.tri.diag") {
#'         num <- k*(k+1)/2
#'         eval(substitute(m[lower.tri(m, diag=TRUE)] %op% vecParams[offset + (1:num)]), parent.frame())
#'       } else if(type[1] == "symmetric") {
#'         num <- k*(k+1)/2
#'         eval(substitute({
#'           m[upper.tri(m, diag=TRUE)] %op% vecParams[offset + (1:num)]
#'         }), parent.frame())
#'         if(load) {
#'           eval(substitute({
#'             m[lower.tri(m)] <- 0
#'             m <- m+t(m)
#'             diag(m) <- 0.5 * diag(m)
#'           }), parent.frame())
#'         }
#'       } else if(type[1] == "full") {
#'         num <- k*k
#'         eval(substitute(m[] %op% vecParams[offset+(1:num)]), parent.frame())
#'       } else if(type[1] == "custom") {
#'         if(is.function(indices)) {
#'           ind <- indices(offset, k)
#'           num <- length(unique(ind[ind>offset]))
#'           maskCustom <- mask
#'           eval(substitute(m[maskCustom] %op% vecParams[ind]), parent.frame())
#'         } else {
#'           stop("ERR:020b5:PCMBase:PCM.R:PCMLoadMatrixParameter:: indices should be a
#'                function(offset, k) returning an integer vector.")
#'         }
#'         } else if(type[1] == "fixed") {
#'           num <- 0
#'       } else {
#'         stop(paste0("ERR:020b6:PCMBase:PCM.R:PCMLoadMatrixParameter:: type ", type[1], " not recognized."))
#'       }
#'     } else {
#'       stop("ERR:020b1:PCMBase:PCM.R:PCMLoadMatrixParameter:: m should be a matrix, because k is different from 1", k, " but has dim(m) = (", toString(dim(m)), ").")
#'     }
#'   } else if(length(dim(m)) == 2) {
#'     # m is a matrix
#'     if(type[1] == "diag1") {
#'       num <- 1
#'       eval(substitute(diag(m)[maskDiag1] %op% vecParams[offset + (1:num)]), parent.frame())
#'     } else if(type[1] == "diag") {
#'       num <- k
#'       eval(substitute(diag(m) %op% vecParams[offset + (1:num)]), parent.frame())
#'     } else if(type[1] == "upper.tri") {
#'       num <- k*(k-1)/2
#'       eval(substitute(m[upper.tri(m)] %op% vecParams[offset + (1:num)]), parent.frame())
#'     } else if(type[1] == "upper.tri.diag") {
#'       num <- k*(k+1)/2
#'       eval(substitute(m[upper.tri(m, diag=TRUE)] %op% vecParams[offset + (1:num)]), parent.frame())
#'     } else if(type[1] == "lower.tri") {
#'       num <- k*(k-1)/2
#'       eval(substitute(m[lower.tri(m)] %op% vecParams[offset + (1:num)]), parent.frame())
#'     } else if(type[1] == "lower.tri.diag") {
#'       num <- k*(k+1)/2
#'       eval(substitute(m[lower.tri(m, diag=TRUE)] %op% vecParams[offset + (1:num)]), parent.frame())
#'     } else if(type[1] == "symmetric") {
#'       num <- k*(k+1)/2
#'       eval(substitute({
#'         m[upper.tri(m, diag=TRUE)] %op% vecParams[offset + (1:num)]
#'       }), parent.frame())
#'       if(load) {
#'         eval(substitute({
#'           m[lower.tri(m)] <- 0
#'           m <- m+t(m)
#'           diag(m) <- 0.5 * diag(m)
#'         }), parent.frame())
#'       }
#'     } else if(type[1] == "full") {
#'       num <- k*k
#'       eval(substitute(m[,] %op% vecParams[offset+(1:num)]), parent.frame())
#'     } else if(type[1] == "custom") {
#'       if(is.function(indices)) {
#'         ind <- indices(offset, k)
#'         num <- length(unique(ind[ind>offset]))
#'         maskCustom <- mask
#'         eval(substitute(m[maskCustom] %op% vecParams[ind]), parent.frame())
#'       } else {
#'         stop("ERR:020b5:PCMBase:PCM.R:PCMLoadMatrixParameter:: indices should be a
#'              function(offset, k) returning an integer vector.")
#'       }
#'       } else if(type[1] == "fixed") {
#'         num <- 0
#'     } else {
#'       stop(paste0("ERR:020b6:PCMBase:PCM.R:PCMLoadMatrixParameter:: type ", type[1], " not recognized."))
#'     }
#'   } else {
#'     stop(paste0("ERR:020b2:PCMBase:PCM.R:PCMLoadMatrixParameter:: m should be a matrix, but has length(dim(m)) != 2, dim(m):(", toString(dim(m)), ")."))
#'   }
#'   num
#' }
#'
#'
#' #' Get a vector with all model parameters unrolled
#' #' @param model a PCM model object
#' #' @return a numerical vector
#' #' @export
#' PCMGetVecParamsFull <- function(model, ...) {
#'   UseMethod("PCMGetVecParamsFull", model)
#' }
#'
#' #' @export
#' PCMGetVecParamsFull.PCM <- function(model, ...) {
#'   specParams <- attr(model, "specParams", exact = TRUE)
#'   R <- PCMNumRegimes(model)
#'   res <- do.call(c, lapply(names(model), function(name) {
#'     if(specParams[[name]]$type[1]=="model") {
#'       PCMGetVecParamsFull(model[[name]], ...)
#'     } else if(specParams[[name]]$type[1] %in% c("gscalar", "gvector", "gmatrix") ) {
#'       rep(as.vector(model[[name]]), R)
#'     } else {
#'       as.vector(model[[name]])
#'     }
#'   }))
#'   unname(res)
#' }
#'
#'
#' #' Set model parameters from a named list
#' #' @param tree a phylo object (possible future use)
#' #' @param model a PCM model object
#' #' @param params a named list with elements among the names found in attr(model, "specParams")
#' #' @param inplace logical indicating if the parameters should be set "inplace" for the
#' #' model object in the calling environment or a new model object with the parameters set
#' #' as specified should be returned. Defaults to TRUE
#' #' @return If inplace is TRUE, the function only has a side effect of setting the
#' #' parameters of the model object in the calling environment; otherwise the function
#' #' returns a modified copy of the model object.
#' #' @export
#' PCMSetParams <- function(model, params, inplace = TRUE, ...) {
#'   UseMethod("PCMSetParams", model)
#' }
#'
#' #' @export
#' PCMSetParams.PCM <- function(model, params, inplace = TRUE, ...) {
#'   specParams <- attr(model, "specParams")
#'   for(name in names(params)) {
#'     if(! (name%in%names(specParams)) ) {
#'       stop(paste0("ERR:020c1:PCMBase:PCM.R:PCMSetParams:: ", name,
#'                   " is not a settable parameter of the model, check attr(model, 'specParams')."))
#'     }
#'     type <- specParams[[name]]$type
#'
#'     if(type[1] != "model") {
#'       if(! identical(length(model[[name]]), length(params[[name]])) ) {
#'         stop(paste0("ERR:020c2:PCMBase:PCM.R:PCMSetParams:: params[[", name,
#'                     "]] is not the same length as model[[", name, "]]; ",
#'                     "length(params[[", name, "]])=", length(params[[name]]),
#'                     ", length(model[[", name, "]])=", length(model[[name]]), "."))
#'       }
#'       if(! identical(dim(model[[name]]), dim(params[[name]])) ) {
#'         stop(paste0("ERR:020c3:PCMBase:PCM.R:PCMSetParams:: params[[", name,
#'                     "]] is not the same dimension as model[[", name, "]]; ",
#'                     "dim(params[[", name, "]])=", str(dim(params[[name]])),
#'                     ", length(model[[", name, "]])=", str(dim(model[[name]])), "."))
#'       }
#'     }
#'   }
#'   for(name in names(params)) {
#'     if(inplace) {
#'       eval(substitute(model[[name]] <- params[[name]]), parent.frame())
#'     } else {
#'       model[[name]] <- params[[name]]
#'     }
#'   }
#'
#'   if(!inplace) {
#'     model
#'   }
#' }
#'
#' #' Get a vector of the varying (non-fixed) parameters in a PCM model
#' #' @param model a PCM model object
#' #' @return a numeric vector of length PCMNumParams(model).
#' #' @seealso \code{\link{PCMSetOrGetVecParams}} \code{\link{PCMGetVecParamsFull}}
#' #' @export
#' PCMGetVecParams <- function(model, ...) {
#'   UseMethod("PCMGetVecParams", model)
#' }
#'
#' #' @export
#' PCMGetVecParams.PCM <- function(model, ...) {
#'   vec <- double(PCMNumParams(model))
#'   PCMSetOrGetVecParams(model, vec, set = FALSE, ...)
#'   vec
#' }
#'
#' #' Inplace set or get the parameters of a PCM from or into a numeric vector
#' #'
#' #' @param model a PCM model object
#' #' @param vecParams a numeric vector of length PCMNumParams(model) or longer
#' #' @param offset an integer indicating offset of the starting position for
#' #' writing/reading in vecParams. Default 0.
#' #' @param set a logical indicating whether to set the model parameters from
#' #' vecParams or to get them into vecParams. Default is TRUE which results in
#' #' setting the model parameters from vecParams. See details.
#' #' @details This is an S3 generic. Note that depending on the value of \code{set},
#' #' this function changes either the model object or the vecParams vector from the
#' #' calling environment.
#' #' @return an integer indicating the number of elements read/written from/to
#' #' vecParams.
#' #' @seealso \code{\link{PCMGetVecParams}} \code{\link{PCMGetVecParamsFull}}
#' #' @export
#' PCMSetOrGetVecParams <- function(
#'   model, vecParams, offset = 0, set = TRUE, ...) {
#'   UseMethod("PCMSetOrGetVecParams", model)
#' }
#'
#' #' @export
#' PCMSetOrGetVecParams.PCM <- function(
#'   model, vecParams, offset = 0, set = TRUE, ... ) {
#'
#'   #force(vecParams)
#'   p <- 0
#'
#'   specParams <- attr(model, "specParams")
#'   k <- attr(model, "k")
#'
#'   for(name in names(specParams)) {
#'     type <- specParams[[name]]$type
#'     mask <- specParams[[name]]$mask
#'     indices <- specParams[[name]]$indices
#'
#'
#'     if(type[1] == "vector") {
#'       for(r in 1:length(attr(model, "regimes"))) {
#'         p <- p + eval(substitute(PCMLoadVectorParameter(
#'           model[[name]][, r], vecParams, offset + p, k, type = type[2], mask = mask, indices = indices, load = set)),
#'           parent.frame())
#'       }
#'     } else if(type[1] == "matrix") {
#'       for(r in 1:length(attr(model, "regimes"))) {
#'         p <- p + eval(substitute(PCMLoadMatrixParameter(
#'           model[[name]][,, r], vecParams, offset + p, k, type = type[2], mask = mask, indices = indices, load = set)),
#'           parent.frame())
#'       }
#'     } else if(type[1] == "gvector") {
#'       # global vector for all regimes
#'       p <- p + eval(substitute(PCMLoadVectorParameter(
#'         model[[name]], vecParams, offset + p, k, type=type[2], mask = mask, indices = indices, load = set)),
#'         parent.frame())
#'     } else if(type[1] == "gmatrix") {
#'       # global vector for all regimes
#'       p <- p + eval(substitute(PCMLoadMatrixParameter(
#'         model[[name]], vecParams, offset + p, k, type=type[2], mask = mask, indices = indices, load = set)),
#'         parent.frame())
#'     } else if(type[1] == "model") {
#'       # nested PCM corresponding to a regime
#'       p <- p + eval(substitute(PCMSetOrGetVecParams(
#'         model[[name]], vecParams, offset + p, set, ...)), parent.frame())
#'     }
#'   }
#'   p
#' }
#'
#'
#' #' Numerical lower bound
#' #' @param model a PCM object
#' #' @return a PCM object of the same S3 classes as model. Calling
#' #' \code{\link{PCMSetOrGetVecParams}} on this object returns a lower
#' #' bound for that can be used, e.g. in a call to \code{\link{optim}}
#' #' @examples
#' #' model <- PCM("BM__NoX0__NoSigmae_x", k = 3)
#' #' PCMLowerBound(model)
#' #' @export
#' PCMLowerBound <- function(model, lowerBoundValue = -10, lowerBoundValuePositiveDiag = 0, namedLowerBoundValues = NULL, ...) {
#'   UseMethod("PCMLowerBound", model)
#' }
#'
#' #' @export
#' PCMLowerBound.PCM <- function(model, lowerBoundValue = -10, lowerBoundValuePositiveDiag = 0, namedLowerBoundValues = NULL, ...) {
#'   k <- attr(model, "k", exact = TRUE)
#'   if(lowerBoundValuePositiveDiag < 0 ) {
#'     stop("ERR:04000:PCMFit:PCMFit.R:PCMLowerBound.PCM:: lowerBoundValuePositiveDiag should be non-negative.")
#'   }
#'   # a vector with the actual parameters excluding repeated and fixed values
#'   par <- double(PCMNumParams(model))
#'   # we set the default upper bound value, but entries corresponding to diagonal
#'   # elements in Choleski factor upper triangular matrix parameters should be set
#'   # to lowerBoundValuePositiveDiag
#'   par[] <- lowerBoundValue
#'
#'   # all parameters unrolled in a vector including repeated parameter values as
#'   # well as fixed parameter values
#'   fullParamVector <- PCMGetVecParamsFull(model)
#'   maxFullParam <- max(fullParamVector, na.rm = TRUE)
#'   if(!is.finite(maxFullParam)) {
#'     maxFullParam <- as.double(1)
#'   }
#'
#'   # a tricky way to insert values in the model parameters that are certainly not
#'   # among the fixed non-countable parameter values. We will use match of these
#'   # values to assign the either value lowerBoundValuePositiveDiag to the right entries in the
#'   # model parameter matrices.
#'   parMask <- maxFullParam + (1:PCMNumParams(model))
#'
#'   # set the values that match unique positions in par.
#'   PCMSetOrGetVecParams(model, parMask)
#'   # Find Choleski factors of positive definite matrices. Such parameters need to
#'   # have positive diagonal elements, i.e. lowerBoundValuePositiveDiag.
#'   specParams <- attr(model, "specParams", exact = TRUE)
#'   for(name in names(specParams)) {
#'     if(specParams[[name]]$type[1] %in% c("matrix", "gmatrix") &&
#'        length(specParams[[name]]$type) >= 3 &&
#'        specParams[[name]]$type[3] == "positive.diag" ) {
#'       if(specParams[[name]]$type[1] == "gmatrix") {
#'         mi <- match(diag(model[[name]]), parMask)
#'         par[unique(mi)] <- lowerBoundValuePositiveDiag
#'       } else if(specParams[[name]]$type[1] == "matrix") {
#'         # model[[name]] is a k x k x R array
#'         R <- PCMNumRegimes(model)
#'         for(r in 1:R) {
#'           mi <- match(diag(as.matrix(model[[name]][,,r], k, k)), parMask)
#'           par[unique(mi)] <- lowerBoundValuePositiveDiag
#'         }
#'       }
#'     }
#'   }
#'   PCMSetOrGetVecParams(model, par)
#'
#'   if(is.list(namedLowerBoundValues)) {
#'     PCMSetParams(model, namedLowerBoundValues)
#'   }
#'
#'   model
#' }
#'
#' #' Numerical upper bound
#' #' @param model a PC
#' #' M object
#' #' @return a PCM object of the same S3 classes as model. Calling
#' #' \code{\link{PCMSetOrGetVecParams}} on this object returns an upper
#' #' bound for that can be used, e.g. in a call to \code{\link{optim}}
#' #' #' model <- PCM("BM__NoX0__NoSigmae_x", k = 3)
#' #' PCMLowerBound(model)
#' #' @export
#' PCMUpperBound <- function(model, upperBoundValue = 10, upperBoundValuePositiveDiag = 10, namedUpperBoundValues = NULL, ...) {
#'   UseMethod("PCMUpperBound", model)
#' }
#'
#' #' @export
#' PCMUpperBound.PCM <- function(model, upperBoundValue = 10, upperBoundValuePositiveDiag = 10, namedUpperBoundValues = NULL, ...) {
#'   if(upperBoundValuePositiveDiag <= 0 ) {
#'     stop("ERR:04010:PCMFit:PCMFit.R:PCMUpperBound.PCM:: upperBoundValuePositiveDiag should be positive.")
#'   }
#'   # a vector with the actual parameters excluding repeated and fixed values
#'   par <- double(PCMNumParams(model))
#'   # we set the default upper bound value, but entries corresponding to diagonal
#'   # elements in Choleski factor upper triangular matrix parameters should be set
#'   # to upperBoundValuePositiveDiag
#'   par[] <- upperBoundValue
#'
#'   # all parameters unrolled in a vector including repeated parameter values as
#'   # well as fixed parameter values
#'   fullParamVector <- PCMGetVecParamsFull(model)
#'   maxFullParam <- max(fullParamVector, na.rm = TRUE)
#'   if(!is.finite(maxFullParam)) {
#'     maxFullParam <- as.double(1)
#'   }
#'
#'   # a tricky way to insert values in the model parameters that are certainly not
#'   # among the fixed non-countable parameter values. We will use match of these
#'   # values to assign the either value upperBoundValuePositiveDiag to the right entries in the
#'   # model parameter matrices.
#'   parMask <- maxFullParam + (1:PCMNumParams(model))
#'
#'   # set the values that match unique positions in par.
#'   PCMSetOrGetVecParams(model, parMask)
#'   # Find Choleski factors of positive definite matrices. Such parameters need to
#'   # have positive diagonal elements, i.e. upperBoundValuePositiveDiag.
#'   specParams <- attr(model, "specParams", exact = TRUE)
#'   for(name in names(specParams)) {
#'     if(specParams[[name]]$type[1] %in% c("matrix", "gmatrix") &&
#'        length(specParams[[name]]$type) >= 3 &&
#'        specParams[[name]]$type[3] == "positive.diag" ) {
#'       if(specParams[[name]]$type[1] == "gmatrix") {
#'         mi <- match(diag(model[[name]]), parMask)
#'         par[unique(mi)] <- upperBoundValuePositiveDiag
#'       } else if(specParams[[name]]$type[1] == "matrix") {
#'         # model[[name]] is a k x k x R array
#'         R <- PCMNumRegimes(model)
#'         for(r in 1:R) {
#'           mi <- match(diag(as.matrix(model[[name]][,,r], k, k)), parMask)
#'           par[unique(mi)] <- upperBoundValuePositiveDiag
#'         }
#'       }
#'     }
#'   }
#'   PCMSetOrGetVecParams(model, par)
#'
#'   if(is.list(namedUpperBoundValues)) {
#'     PCMSetParams(model, namedUpperBoundValues)
#'   }
#'
#'   model
#' }
#'
#' #' Generate a random parameter vector for a model using uniform distribution between its lower and upper bounds.
#' #' @param model a model PCM model object
#' #' @param n an integer specifying the number of random vectors to generate
#' #' @param argsPCMLowerBound,argsPCMUpperBound named lists of arguments passed to
#' #' \code{\link{PCMLowerBound}} and \code{\link{PCMUpperBound}}.
#' #' @param argsPCMGetVecParams a named list of argument passed to PCMGetVecPaarmas
#' #' @param ... additional parameters.
#' #' @return if n = 1, a numeric vector of length \code{PCMNumParams(model)}; if n > 1,
#' #' a numeric matrix of dimension n x \code{PCMNumParams(model)}.
#' #' @seealso PCMUpperBound PCMLowerBound PCMGetVecParams
#' #' @export
#' PCMRandomVecParams <- function(model,
#'                                n = 1L,
#'                                argsPCMLowerBound = NULL,
#'                                argsPCMUpperBound = NULL,
#'                                argsPCMGetVecParams = NULL,
#'                                ...) {
#'   UseMethod("PCMRandomVecParams", model)
#' }
#'
#' #' @importFrom stats runif
#' #' @export
#' PCMRandomVecParams.PCM <- function(model,
#'                                    n = 1L,
#'                                    argsPCMLowerBound = NULL,
#'                                    argsPCMUpperBound = NULL,
#'                                    argsPCMGetVecParams = NULL,
#'                                    ...) {
#'
#'   lowerModel <- do.call(PCMLowerBound, c(list(model = model), argsPCMLowerBound))
#'   lowerVecParams <- do.call(PCMGetVecParams, c(list(model = lowerModel),
#'                                                argsPCMGetVecParams))
#'
#'   upperModel <- do.call(PCMUpperBound, c(list(model = model), argsPCMUpperBound))
#'   upperVecParams <- do.call(PCMGetVecParams, c(list(model = upperModel), argsPCMGetVecParams))
#'
#'   p <- PCMNumParams(model)
#'   res <- sapply(1:p, function(i) {
#'     runif(n, lowerVecParams[i], upperVecParams[i])
#'   })
#'   res
#' }
