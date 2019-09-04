# Copyright 2016-2019 Venelin Mitov
#
# This file is part of PCMBase.
#
# PCMBase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PCMBase is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PCMBase.  If not, see <http://www.gnu.org/licenses/>.

#' Get a list of PCM models currently implemented
#' @param pattern a character string specifying an optional for the model-names to search for.
#' @param parentClass a character string specifying an optional parent class of the models to look for.
#' @param ... additional arguments used by implementing methods.
#'
#' @return a character vector of the model classes found.
#' @details The function is using the S3 api function \code{\link{methods}} looking
#' for all registered implementations of the function \code{\link{PCMSpecify}}.
#' @importFrom utils methods
#' @export
#' @examples
#' PCMModels()
#' PCMModels("^OU")
PCMModels <- function(pattern = NULL, parentClass = NULL, ...) {
  models <- sub("PCMSpecify.", "", as.character(methods("PCMSpecify")), fixed = TRUE)
  if(!is.null(pattern)) {
    models <- models[grep(pattern, models, ...)]
  }
  if(!is.null(parentClass)) {
    models <- models[sapply(models, function(m) {
      length(intersect(class(PCM(m)), parentClass)) > 0
    })]
  }
  models
}

#' Global options for the PCMBase package
#'
#'
#' @return a named list with the currently set values of the following global
#' options:
#' \itemize{
#' \item{\code{PCMBase.Value.NA }}{NA value for the likelihood; used in GaussianPCM to
#' return this value in case of an error occuring
#' during likelihood calculation. By default, this is set to \code{as.double(NA)}.}
#' \item{\code{PCMBase.Errors.As.Warnings }}{a logical flag indicating if errors
#' (occuring, e.g. during likelihood calculation) should be treated as warnings
#' and added as an attribute "error" to returne likelihood values. Default TRUE.}
#' \item{\item{PCMBase.Raise.Lik.Errors} }{Should numerical and other sort of
#' errors occurring during likelihood calculation be raised either as errors or
#' as warnings, depending on the option \code{PCMBase.Errors.As.Warnings}.
#' Default TRUE. This option can be useful if too frequent warnings get raised
#' during a model fit procedure.}
#' \item{\code{PCMBase.Threshold.Lambda_ij }}{a 0-threshold for abs(Lambda_i + Lambda_j),
#' where Lambda_i and Lambda_j are eigenvalues of the parameter matrix H of an OU or
#' other model. Default 1e-8. See \code{\link{PCMPExpxMeanExp}}.}
#' \item{\code{PCMBase.Threshold.EV }}{A 0-threshold for the eigenvalues of the
#' matrix V for a given branch. The V matrix is considered singular if it has
#' eigenvalues smaller than \code{PCMBase.Threshold.EV } or when the ratio
#' min(svdV)/max(svdV) is below \code{PCMBase.Threshold.SV }. Default is 1e-5.
#' Treatment of branches with singular V matrix is defined by the option
#' \code{PCMBase.Skip.Singular}.}
#' \item{\code{PCMBase.Threshold.SV }}{A 0-threshold for min(svdV)/max(svdV), where
#' svdV is the vector of singular values of the matrix V for a given branch.
#' The V matrix is considered singular if it has eigenvalues smaller than
#' \code{PCMBase.Threshold.EV } or when the ratio min(svdV)/max(svdV) is below
#' PCMBase.Threshold.SV. Default is 1e-6. Treatment
#' of branches with singular V matrix is defined by the option \code{PCMBase.Skip.Singular}.}
#' \item{\code{PCMBase.Threshold.Skip.Singular }}{A double indicating if a branch of shorter
#' length with singular matrix V should be skipped during likelihood calculation. Setting this
#' option to a higher value, together with a TRUE value for the option PCMBase.Skip.Singular
#' will result in tolerating some parameter values resulting in singular variance covariance
#' matrix of the transition distribution. Default 1e-4.}
#' \item{\code{PCMBase.Skip.Singular }}{A logical value indicating whether branches with
#' singular matrix V and shorter than \code{getOption("PCMBase.Threshold.Singular.Skip")}
#'  should be skipped during likelihood calculation, adding their children
#' L,m,r values to their parent node. Default TRUE. Note, that setting this option to FALSE
#' may cause some models to stop working, e.g. the White model. Setting this option to FALSE
#' will also cause errors or NA likelihood values in the case of trees with very short or
#' 0-length branches.}
#' \item{\code{PCMBase.Tolerance.Symmetric }}{A double specifying the tolerance in tests
#' for symmetric matrices. Default 1e-8; see also \code{\link{isSymmetric}}.}
#' \item{\code{PCMBase.Lmr.mode }}{An integer code specifying the parallel likelihood calculation mode.}
#' \item{\code{PCMBase.ParamValue.LowerLimit}}{Default lower limit value for parameters, default setting is -10.0.}
#' \item{\code{PCMBase.ParamValue.LowerLimit.NonNegativeDiagonal}}{Default lower limit value for parameters corresponding to non-negative diagonal elements of matrices, default setting is 0.0.}
#' \item{\code{PCMBase.ParamValue.UpperLimit} }{Default upper limit value for parameters, default setting is 10.0.}
#' \item{\code{PCMBase.Transpose.Sigma_x} }{Should upper diagonal factors for variance-covariance rate matrices be transposed, e.g. should Sigma = t(Sigma_x) Sigma_x or, rather Sigma = Sigma_x t(Sigma_x)? Note that the two variants are not equal. The default is FALSE, meaning Sigma = Sigma_x t(Sigma_x). In this case, though Sigma_x is not the actual upper Cholesky factor of Sigma, i.e. chol(Sigma) != Sigma_x. See also \code{\link{chol}}. This option applies to parameters Sigma_x, Sigmae_x and Sigmaj_x.}
#' \item{\code{PCMBase.MaxLengthListCladePartitions} }{Maximum number of tree partitions returned by \code{\link{PCMTreeListCladePartitions}}. This option has the goal to interrupt the recursive search for new partitions in the case of calling PCMTreeListCladePartitions on a big tree with a small value of the maxCladeSize argument. By default this is set to Inf.}
#' \item{\code{PCMBase.PCMPresentCoordinatesFun} }{A function with the same synopsis as \code{\link{PCMPresentCoordinates}} that can be specified in case of custom setting for the present coordinates for specific nodes of the tree. See \code{\link{PCMPresentCoordinates}}, and \code{\link{PCMInfo}}.}
#' \item{\code{PCMBase.Use1DClasses} }{Logical indicating if 1D arithmetic operations
#' should be used instead of multi-dimensional ones. This can speed-up computations
#' in the case of a single trait. Currently, this feature is implemented only in
#' the PCMBaseCpp R-package and only for some model types, such as OU and BM.
#' Default: FALSE}
#' \item{\code{PCMBase.PrintSubscript_u} }{Logical indicating if a subscript 'u'
#' should be printed instead of a subscript 'x'. Used in \code{PCMTable}. Default: FALSE.}
#' \item{\code{PCMBase.MaxNForGuessSigma_x} }{Integer indicating the maximum
#' number of tips to use for analytical calculation of the evolutionary rate
#' matrix under a BM assumption. This option is used in the suggested PCMFit
#' R-package. Default: 1000.}
#' \item{\code{PCMBase.UsePCMVarForVCV} }{Logical (default: FALSE) indicating
#' if the function \code{\link{PCMTreeVCV}} should use \code{\link{PCMVar}}
#' instead of ape's function \code{\link{vcv}} to calculate the phylogenetic
#' variance covariance matrix under BM assumption. Note that setting this option
#' to TRUE would slow down the function PCMTreeVCV considerably but may be more
#' stable, particularly in the case of very big and deep trees, where previous
#' ape's versions of the \code{\link{vcv}} function have thrown stack-overflow
#' errors.}
#' }
#' @export
#' @examples
#' PCMOptions()
PCMOptions <- function() {
  list(PCMBase.Value.NA = getOption("PCMBase.Value.NA", as.double(NA)),
       PCMBase.Errors.As.Warnings = getOption("PCMBase.Errors.As.Warnings", TRUE),
       PCMBase.Threshold.Lambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8),
       PCMBase.Threshold.EV = getOption("PCMBase.Threshold.EV", 1e-5),
       PCMBase.Threshold.SV = getOption("PCMBase.Threshold.SV", 1e-6),
       PCMBase.Threshold.Skip.Singular = getOption("PCMBase.Threshold.Skip.Singular", 1e-4),
       PCMBase.Skip.Singular = getOption("PCMBase.Skip.Singular", TRUE),
       PCMBase.Tolerance.Symmetric = getOption("PCMBase.Tolerance.Symmetric", 1e-8),
       PCMBase.Lmr.mode = getOption("PCMBase.Lmr.mode", as.integer(11)),
       PCMBase.ParamValue.LowerLimit.NonNegativeDiagonal = getOption("PCMBase.ParamValue.LowerLimit.NonNegativeDiagonal", 0.0),
       PCMBase.ParamValue.LowerLimit = getOption("PCMBase.ParamValue.LowerLimit", -10.0),
       PCMBase.ParamValue.UpperLimit = getOption("PCMBase.ParamValue.UpperLimit", 10.0),
       PCMBase.Transpose.Sigma_x = getOption("PCMBase.Transpose.Sigma_x", FALSE),
       PCMBase.MaxLengthListCladePartitions = getOption("PCMBase.MaxLengthListCladePartitions", Inf),
       PCMBase.PCMPresentCoordinatesFun = getOption("PCMBase.PCMPresentCoordinatesFun", PCMPresentCoordinates),
       PCMBase.Use1DClasses = getOption("PCMBase.Use1DClasses", FALSE),
       PCMBase.Raise.Lik.Errors = getOption("PCMBase.Raise.Lik.Errors", TRUE),
       PCMBase.PrintSuffix_u = getOption("PCMBase.PrintSuffix_u", FALSE),
       PCMBase.MaxNForGuessSigma_x = getOption("PCMBase.MaxNForGuessSigma_x", 1000L),
       PCMBase.UsePCMVarForVCV = getOption("PCMBase.UsePCMVarForVCV", FALSE)
       )
}

#' @name PCM
#' @title Create a phylogenetic comparative model object
#'
#' @description This is the entry-point function for creating model objects
#' within the PCMBase framework representing a single model-type with one or
#' several model-regimes of this type associated with the branches of a tree.
#' For mixed Gaussian phylogenetic models, which enable multiple model-types,
#' use the \code{\link{MixedGaussian}} function.
#' @param model This argument can take one of the following forms:
#' \itemize{
#' \item a character vector of the S3-classes of the model object to be
#' created (one model object can have one or more S3-classes, with the class
#' PCM at the origin of the hierarchy);
#' \item an S3 object which's class inherits from the PCM S3 class.
#' }
#' The Details section explains how these two types of input are processed.
#' @param modelTypes a character string vector specifying a set (family) of
#' model-classes, to which the constructed model object belongs. These are used
#' for model-selection.
#' @param k integer denoting the number of traits (defaults to 1).
#' @param regimes a character or integer vector denoting the regimes.
#' @param params NULL (default) or a list of parameter values (scalars, vectors,
#' matrices, or arrays) or sub-models (S3 objects inheriting from the PCM class).
#' See details.
#' @param vecParams NULL (default) or a numeric vector the vector
#' representation of the variable parameters in the model. See details.
#' @param offset integer offset in vecParams; see Details.
#' @param spec NULL or a list specifying the model parameters (see
#' \code{\link{PCMSpecify}}). If NULL (default), the generic PCMSpecify
#' is called on the created object of class \code{model}.
#' @param ... additional parameters intended for use by sub-classes of the PCM
#' class.
#' @return an object of S3 class as defined by the argument model.
#'
#' @details This is an S3 generic. The PCMBase package defines three methods for
#' it:
#' \itemize{
#' \item{PCM.PCM:}{A default constructor for any object with a class inheriting
#' from "PCM".}
#' \item{PCM.character:}{A default PCM constructor from a character string
#' specifying the type of model.}
#' \item{PCM.default:}{A default constructor called when no other constructor is
#' found. When called this constructor raises an error message.}
#' }
#' @seealso \code{\link{MixedGaussian}}
#'
#' @importFrom ape rtree
#'
#' @export
#'
#' @examples
#' # a Brownian motion model with one regime
#' modelBM <- PCM(model = "BM", k = 2)
#' # print the model
#' modelBM
#'
#' # a BM model with two regimes
#' modelBM.ab <- PCM("BM", k = 2, regimes = c("a", "b"))
#' modelBM.ab
#'
#' # print a single parameter of the model (in this case, the root value)
#' modelBM.ab$X0
#'
#' # assign a value to this parameter (note that the brackets [] are necessary
#' # to preserve  the parameter attributes):
#' modelBM.ab$X0[] <- c(5, 2)
#'
#' PCMNumTraits(modelBM)
#' PCMNumRegimes(modelBM)
#' PCMNumRegimes(modelBM.ab)
#'
#' # number of numerical parameters in the model
#' PCMParamCount(modelBM)
#'
#' # Get a vector representation of all parameters in the model
#' PCMParamGetShortVector(modelBM)
#'
#' # Limits for the model parameters:
#' lowerLimit <- PCMParamLowerLimit(modelBM)
#' upperLimit <- PCMParamUpperLimit(modelBM)
#'
#' # assign the model parameters at random: this will use uniform distribution
#' # with boundaries specified by PCMParamLowerLimit and PCMParamUpperLimit
#' # We do this in two steps:
#' # 1. First we generate a random vector. Note the length of the vector equals PCMParamCount(modelBM)
#' randomParams <- PCMParamRandomVecParams(modelBM, PCMNumTraits(modelBM), PCMNumRegimes(modelBM))
#' randomParams
#' # 2. Then we load this random vector into the model.
#' PCMParamLoadOrStore(modelBM, randomParams, 0, PCMNumTraits(modelBM), PCMNumRegimes(modelBM), TRUE)
#'
#' print(modelBM)
#'
#' PCMParamGetShortVector(modelBM)
#'
#' # generate a random phylogenetic tree of 10 tips
#' tree <- ape::rtree(10)
#'
#' #simulate the model on the tree
#' traitValues <- PCMSim(tree, modelBM, X0 = modelBM$X0)
#'
#' # calculate the likelihood for the model parameters, given the tree and the trait values
#' PCMLik(traitValues, tree, modelBM)
#'
#' # create a likelihood function for faster processing for this specific model.
#' # This function is convenient for calling in optim because it recieves and parameter
#' # vector instead of a model object.
#' likFun <- PCMCreateLikelihood(traitValues, tree, modelBM)
#' likFun(randomParams)
PCM <- function(
  model, modelTypes = class(model)[1], k = 1L, regimes = 1L,
  params = NULL, vecParams = NULL, offset = 0L,
  spec = NULL, ...) {
  UseMethod("PCM", model)
}

#' @export
PCM.default <- function(
  model, modelTypes = class(model)[1], k = 1L, regimes = 1L,
  params = NULL, vecParams = NULL, offset = 0L,
  spec = NULL, ...) {
  stop(paste0("ERR:02091:PCMBase:PCM.R:PCM.default:: You should provide a PCM object, but provided a ", class(model)[1]))
}

#' @export
PCM.PCM <- function(
  model, modelTypes = class(model)[1], k = 1L, regimes = 1L,
  params = NULL, vecParams = NULL, offset = 0L,
  spec = NULL, ...) {

  if(is.null(spec)) {
    spec <- PCMSpecify(model, ...)
  }

  obj <- PCMDefaultObject(spec, model, ...)

  if(is.null(attr(obj, "k", exact = TRUE))) {
    attr(obj, "k") <- as.integer(k)
  }

  if(is.null(attr(obj, "regimes", exact = TRUE))) {
    attr(obj, "regimes") <- regimes
  }

  if(is.null(attr(obj, "p", exact = TRUE))) {
    attr(obj, "p") <- PCMParamCount(obj)
  }

  if(is.null(attr(obj, "modelTypes", exact = TRUE))) {
    attr(obj, "modelTypes") <- modelTypes
  }

  if(is.null(attr(obj, "mapping", exact = TRUE))) {
    attr(obj, "mapping") <- match(class(obj)[1],
                                    attr(obj, "modelTypes", exact = TRUE))
  }

  if(is.null(attr(obj, "spec", exact = TRUE))) {
    attr(obj, "spec") <- spec
  }

  if(!is.null(params)) {
    PCMParamSetByName(obj, params, ...)
  }

  if(!is.null(vecParams)) {
    PCMParamLoadOrStore(obj, vecParams, offset, load = TRUE)
  }

  obj
}

#' @export
PCM.character <- function(
  model, modelTypes = model[1], k = 1L, regimes = 1L,
  params = NULL, vecParams = NULL, offset = 0L,
  spec = NULL, ...) {

  modelObj <- list()
  class(modelObj) <- model
  if(length(model) == 1) {
    class(modelObj) <- c(class(modelObj), PCMParentClasses(modelObj))
  }
  attr(modelObj, "k") <- as.integer(k)
  attr(modelObj, "regimes") <- regimes

  PCM(modelObj, modelTypes, k, regimes, params, vecParams, offset, spec, ...)
}

#' Check if an object is a PCM.
#' @param x an object.
#' @return TRUE if `x` inherits from the S3 class "PCM".
#' @export
is.PCM <- function(x) inherits(x, "PCM")

#' Generate a default object of a given PCM model type or parameter type
#' @param spec any object having a class attribute. The value of this object is not
#' used, but its class is used for method-dispatch.
#' @param model a PCM object used to extract attributes needed for creating a
#' default object of class specified in \code{class(spec)}, such as the number of
#' traits (k) or the regimes and the number of regimes;
#' @param ... additional arguments that can be used by methdos.
#'
#' @description This is an S3 generic. See, e.g. `PCMDefaultObject.MatrixParameter`.
#' @return a parameter or a PCM object.
#' @export
PCMDefaultObject <- function(spec, model, ...) {
  UseMethod("PCMDefaultObject", spec)
}

#' @export
PCMDefaultObject.PCM <- function(spec, model, ...) {
  o <- spec
  PCMParamSetByName(o, lapply(spec, PCMDefaultObject, model = o), replaceWholeParameters = TRUE, ...)
  o[] <- o[!sapply(o, is.null)]
  attr(o, "spec") <- spec
  o
}


#' @method print PCM
#' @export
print.PCM <- function(x, ...) cat(format(x, ...), "\n")

#' @method print PCM
#' @export
format.PCM <- function(x, ...) {
  if( !is.PCM(x) ) {
    stop("format.PCM:: x must inherit from S3 class PCM.")
  }
  spec <- attr(x, "spec", exact = TRUE)

  res <- paste0(
    PCMDescribe(x, ...),
    "\nS3 class: ", toString(class(x)), "; ",
    "k=", attr(x, "k", exact = TRUE), "; p=", PCMParamCount(x, ...), "; ",
    "regimes: ", toString(attr(x, "regimes")), ". Parameters/sub-models:\n")

  for(name in names(spec)) {

    if(is.PCM(spec[[name]])) {
      objToPrint <- x[[name]]
      type <- class(objToPrint)
      description <- attr(objToPrint, "description", exact = TRUE)

      strList <- as.list(capture.output(print(objToPrint)))
      strList$sep = "\n"
    } else if(!is.Omitted(spec[[name]])) {
      objToPrint <- x[[name]]
      type <- class(objToPrint)
      description <- attr(objToPrint, "description", exact = TRUE)

      # avoid writing class and description attributes twice
      attr(objToPrint, "class") <- NULL
      attr(objToPrint, "description") <- NULL
      attr(objToPrint, "k") <- NULL
      attr(objToPrint, "regimes") <- NULL
      attr(objToPrint, "spec") <- NULL
      strList <- as.list(capture.output(print(objToPrint)))
      strList$sep = "\n"
    } else {

      # _Omitted
      type <- class(spec[[name]])
      description <- attr(spec[[name]], "description", exact = TRUE)

      strList <- list(sep = "\n")
    }
    res <- paste0(res, name, " (", toString(type),
                  if(!is.null(description)) paste0("; ", description) else "", "):\n",
                  do.call(paste, strList))
    res <- paste0(res, "\n")
  }
  res
}

#' Parent S3 classes for a model class
#' @param model an S3 object.
#' @details This S3 generic function is intended to be specified for user models.
#' This function is called by the `PCM.character` method to determine the parent
#' classes for a given model class.
#' @return a vector of character string denoting the names of the parent classes
#' @export
PCMParentClasses <- function(model) {
  UseMethod("PCMParentClasses", model)
}

#' @export
PCMParentClasses.PCM <- function(model) c()

#' Human friendly description of a PCM
#' @param model a PCM model object
#' @param ... additional arguments used by implementing methods.
#' @details This S3 generic function is intended to be specified for user models
#' @return a character string
#' @export
PCMDescribe <- function(model, ...) {
  UseMethod("PCMDescribe", model)
}

#' @export
PCMDescribe.PCM <- function(model, ...) {
  "PCM base class model with no parameters; serves as a basis for PCM model classes"
}


#' Parameter specification of PCM model
#' @param model a PCM model object.
#' @param ... additional arguments used by implementing methods.
#'
#' @description The parameter specification of a PCM model represents a named
#' list with an entry for each parameter of the model. Each entry in the list is
#' a structure defining the S3 class of the parameter and its verbal description.
#' This is an S3 generic. See `PCMSpecify.OU` for an example method.
#'
#' @return a list specifying the paramters of a PCM.
#' @export
PCMSpecify <- function(model, ...) {
  UseMethod("PCMSpecify", model)
}

#' Describe the parameters of a PCM
#' @param model a PCM object.
#' @param ... additional arguments that can be used by implementing methods.
#' @return a named list with character elements corresponding to each parameter.
#'
#' @description This is an S3 generic.
#'
#' @export
PCMDescribeParameters <- function(model, ...) {
  UseMethod("PCMDescribeParameters", model)
}

#' Specify the parameterizations for each parameter of a model
#' @param model a PCM.
#' @param ... additional arguments used by implementing methods.
#'
#' @description These are S3 generics. `PCMListParameterizations` should return
#' all possible parametrizations for the class of `model`.
#' `PCMListDefaultParameterizations` is a handy way to specify a subset of all
#' parametrizations. `PCMListDefaultParameterizations` should be used to avoid
#' generating too many model parametrizations which occupy space in the R-global
#' environment while they are not used (see \link{PCMGenerateParameterizations}).
#' It is mandatory to implement a specification for `PCMListParameterizations`
#' for each newly defined class of models.
#' `PCMListDefaultParameterizations` has a default implementation that calls
#' `PCMListParameterizations` and returns the first parametrization for each
#' parameter. Hence, implementing a method for `PCMListDefaultParameterizations`
#' for a newly defined model type is optional.
#'
#' @return a named list with list elements corresponding to each parameter in
#' model. Each list element is a list of character vectors, specifying the possible
#' S3 class attributes for the parameter in question. For an example, type
#' `PCMListParameterizations.BM` to see the possible parameterizations for the
#' BM model.
#'
#' @seealso PCMGenerateParameterizations
#' @export
PCMListParameterizations <- function(model, ...) {
  UseMethod("PCMListParameterizations", model)
}

#' @rdname PCMListParameterizations
#' @export
PCMListDefaultParameterizations <- function(model, ...) {
  UseMethod("PCMListDefaultParameterizations", model)
}

#' @export
PCMListDefaultParameterizations.default <- function(model, ...) {
  lapply(PCMListParameterizations(model, ...), function(item) item[1])
}

#' Cartesian product of possible parameterizations for the different parameters of a model
#' @param model a PCM object.
#' @param listParameterizations a list returned by a method for `PCMListParameterizations`.
#' Default: `PCMListParameterizations(model, ...)`.
#' @param ... additional arguments passed to `PCMListParameterizations(model, ...)`.
#'
#' @description This function generates a data.table in which each column corresponds to
#' one parameter of model and each row corresponds to one combination of parameterizations
#' for the model parameters, such that the whole table corresponds to the Cartesian product
#' of the lists found in `listParameterizations`. Usually, subsets of this table shold be
#' passed to `PCMGenerateParameterizations`
#' @return a data.table object.
#'
#' @importFrom data.table data.table as.data.table
#' @export
PCMTableParameterizations <- function(
  model, listParameterizations = PCMListParameterizations(model, ...), ...) {

  # produce the cartesian product of the parameterizations for each parameter
  dtClasses <- do.call(expand.grid, listParameterizations)
  dtClasses <- as.data.table(dtClasses)
  dtClasses
}

#' Generate possible parameterizations for a given type of model
#'
#' @description A parameterization of a PCM of given type, e.g. OU, is a PCM-class
#' inheriting from this type, which imposes some restrictions or transformations of
#' the parameters in the base-type. This function generates the S3 methods responsible
#' for creating such parameterizations, in particular it generates the definition
#' of the methods for the two S3 generics `PCMParentClasses` and `PCMSpecify` for
#' al parameterizations specified in the `tableParameterizations` argument.
#' @param model a PCM object.
#' @param listParameterizations a list or a sublist returned by `PCMListParameterizations`.
#' Default: `PCMListParameterizations(model)`.
#' @param tableParameterizations a data.table containing the parameterizations to
#' generate. By default this is generated from `listParameterizations` using a
#' call `PCMTableParameterizations(model, listParameterizations)`. If specified
#' by the user, this parameter takes precedence over `listParameterizations` and
#' `listParameterizations` is not used.
#' @param env an environment where the method definitions will be stored.
#' Default: `env = .GlobalEnv`.
#' @param useModelClassNameForFirstRow A logical specifying if the S3 class name of
#' `model` should be used as a S3 class for the model defined in the first row of
#' `tableParameterizations`. Default: FALSE.
#' @param sourceFile NULL or a character string indicating a .R filename, to
#' which the automatically generated code will be saved. If NULL (the default),
#' the generated source code is evaluated and the S3 methods are defined in the
#' global environment. Default: NULL.
#' @return This function does not return a value. It only has a side effect by
#' defining S3 methods in `env`.
#'
#' @export
PCMGenerateParameterizations <- function(
  model,
  listParameterizations = PCMListParameterizations(model),
  tableParameterizations = PCMTableParameterizations(model, listParameterizations),
  env = .GlobalEnv,
  useModelClassNameForFirstRow = FALSE,
  sourceFile = NULL) {

  if(!is.PCM(model)) {
    # We assume that the passed model's class is a single character not including the
    # corresponding hierarchy, which is specified by PCMParentClasses
    class(model) <- c(class(model)[1], PCMParentClasses(model))
  }

  classesModel <- class(model)
  paramNames <- names(tableParameterizations)
  paramDescriptions <- PCMDescribeParameters(model)

  for(i in seq_len(nrow(tableParameterizations))) {
    nameFunPCMParentClasses <- paste0("PCMParentClasses.", classesModel[1])
    nameFunPCMSpecify <- paste0("PCMSpecify.", classesModel[1])
    if(i != 1 || !useModelClassNameForFirstRow) {
      # add suffix to names based on the parameterizations for the row
      suffix <- ""
      for(name in paramNames) {
        qualifiersParam <- tableParameterizations[i][[name]][[1]]
        qualifiersParam <- qualifiersParam[startsWith(qualifiersParam, "_")]
        if(length(qualifiersParam) > 0) {
          suffix <- paste0(suffix, "_")
          suffix <- paste0(suffix, do.call(paste0, as.list(qualifiersParam)), "_", name)
        } else {
          suffix <- paste0(suffix, "__")
          suffix <- paste0(suffix, name)
        }
      }
      nameFunPCMParentClasses <- paste0(nameFunPCMParentClasses, suffix)
      nameFunPCMSpecify <- paste0(nameFunPCMSpecify, suffix)
      parentClasses <- classesModel
    } else {
      # this is the base class for a model, all parameterisations being its daughter classes
      parentClasses <- classesModel[-1]
    }

    sourcePCMParentClasses <- paste0(
      nameFunPCMParentClasses, " <- function(model) ",
      PCMCharacterVectorToRExpression(parentClasses))

    sourcePCMSpecify <- paste0(
      nameFunPCMSpecify, " <- function(model, ...) {\n",
      "spec <- list(\n")
    for(j in 1:length(paramNames)){
      sourcePCMSpecify <- paste0(
        sourcePCMSpecify,
        paramNames[j], " = structure(0.0, class = ",
        PCMCharacterVectorToRExpression(tableParameterizations[i][[paramNames[j]]][[1]]),
        ",\n",
        "description = '", paramDescriptions[j], "')")
      if(j < length(paramNames)) {
        sourcePCMSpecify <- paste0(
          sourcePCMSpecify, ",\n")
      } else {
        sourcePCMSpecify <- paste0(
          sourcePCMSpecify, ")\n")
      }
    }
    sourcePCMSpecify <- paste0(
      sourcePCMSpecify,
      "attributes(spec) <- attributes(model)\n",
      "if(is.null(names(spec))) names(spec) <- ",
      PCMCharacterVectorToRExpression(paramNames), "\n",
      "if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')\n",
      "spec\n",
      "}")

    if(!is.null(sourceFile)) {
      write(paste0(
        "#' @export\n", sourcePCMParentClasses, "\n\n",
        "#' @export\n", sourcePCMSpecify, "\n\n"),
        file = sourceFile, append = TRUE)
    } else {
      eval(parse(text = sourcePCMParentClasses), env)
      eval(parse(text = sourcePCMSpecify), env)
    }
  }
}


#' Generate default model types for given PCM base-classes
#' @description This function calls `PCMListParameterizations` or
#' `PCMListDefaultParameterizations` and generates the corresponding
#' `PCMParentClasses` and `PCMSpecify` methods in the global environment.
#' @param baseTypes a character vector specifying base S3-class names for which
#' the default parametrizations (sub-classes) will be generated. Defaults to
#' `c("BM", "OU")`.
#' @param parametrizations a character string specifying which one of
#' `PCMListParameterizations` or `PCMListDefaultParameterizations` should be used.
#' This argument should be one of:
#' \itemize{
#' \item{"all"}{for calling `PCMListParameterizations`}
#' \item{"default"}{for calling `PCMListDefaultParameterizations`}
#' }
#' @param sourceFile NULL or a character string indicating a .R filename, to
#' which the automatically generated code will be saved. If NULL (the default),
#' the generated source code is evaluated and the S3 methods are defined in the
#' global environment. Default: NULL.
#' @return This function has side effects only and does not return a value.
#' @seealso PCMListDefaultParameterizations
#' @export
PCMGenerateModelTypes <- function(
  baseTypes = c("BM", "OU"),
  parametrizations = c("default", "all"),
  sourceFile = NULL) {

  if( !is.null(sourceFile) ) {
    write(paste0(
      "# This file was auto-generated through a call to ",
      "PCMGenerateModelTypes()\n" ,
      "# Do not edit by hand.\n\n"),
      file = sourceFile)
  }

  for(bt in baseTypes) {
    o <- structure(0.0, class=bt)
    PCMGenerateParameterizations(
      o,
      listParameterizations = if(parametrizations[[1]] == "all") {
        PCMListParameterizations(o)
      } else {
        PCMListDefaultParameterizations(o)
      },
      sourceFile = sourceFile)
  }
}

#' Fix a parameter in a PCM model
#' @param model a PCM object
#' @param name a character string
#' @return a copy of the model with added class '_Fixed' to the class of the
#' parameter \code{name}
#' @export
PCMFixParameter <- function(model, name) {
  modelNew <- model
  if(!is.Fixed(modelNew[[name]])) {
    class(modelNew[[name]]) <- c(class(modelNew[[name]]), "_Fixed")
    class(attr(modelNew, "spec")[[name]]) <- c(class(attr(modelNew, "spec")[[name]]), "_Fixed")
  }
  modelNew
}

#' Unfix a parameter in a PCM model
#' @param model a PCM object
#' @param name a character string
#' @return a copy of the model with removed class '_Fixed' from the class of the
#' parameter \code{name}
#' @export
PCMUnfixParameter <- function(model, name) {
  modelNew <- model
  if(is.Fixed(modelNew[[name]])) {
    class(modelNew[[name]]) <- class(modelNew[[name]])[!class(modelNew[[name]])=="_Fixed"]
    class(attr(modelNew, "spec")[[name]]) <-
      class(attr(modelNew, "spec")[[name]])[!class(attr(modelNew, "spec")[[name]])=="_Fixed"]
  }
  modelNew
}



#' Number of traits modeled by a PCM
#' @param model a PCM object or an a parameter object (the name of this argument
#' could be misleading, because both, model and parameter objects are supported).
#' @return an integer
#' @export
PCMNumTraits <- function(model) {
  UseMethod("PCMNumTraits", model)
}

#' @export
PCMNumTraits.PCM <- function(model) {
  as.integer(attr(model, "k", exact = TRUE))
}

#' Get the regimes (aka colors) of a PCM or of a PCMTree object
#' @param obj a PCM or a PCMTree object
#' @return a character or an integer vector giving the regime names in the obj
#' @export
PCMRegimes <- function(obj) {
  UseMethod("PCMRegimes", obj)
}

#' @export
PCMRegimes.PCM <- function(obj) {
  attr(obj, "regimes", exact = TRUE)
}


#' Number of regimes in a obj
#' @param obj a PCM object
#' @return an integer
#' @export
PCMNumRegimes <- function(obj) {
  UseMethod("PCMNumRegimes", obj)
}

#' @export
PCMNumRegimes.PCM <- function(obj) {
  length(PCMRegimes(obj))
}

#' Get the model type(s) of a model
#'
#' @description For a regular PCM object, the model type is its S3 class. For a
#' MixedGaussian each regime is mapped to one of several possible model types.
#'
#' @param obj a PCM object
#' @return a character vector
#' @export
PCMModelTypes <- function(obj) {
  UseMethod("PCMModelTypes", obj)
}

#' @export
PCMModelTypes.PCM <- function(obj) {
  class(obj)[1]
}


#' Integer vector giving the model type index for each regime
#' @param model a PCM model
#' @param tree a phylo object with an edge.part member
#' @param ... additional parameters passed to methods
#' @return an integer vector with elements corresponding to the elements in
#' \code{PCMTreeGetPartNames(tree)}
#' @details This is a generic S3 method. The default implementation for the basic
#' class PCM returns a vector of 1's, because it assumes that a single model type
#' is associated with each regime. The implementation for mixed Gaussian models
#' returns the mapping attribute of the MixedGaussian object reordered to
#' correspond to \code{PCMTreeGetPartNames(tree)}.
#' @export
PCMMapModelTypesToRegimes <- function(model, tree, ...) {
  UseMethod("PCMMapModelTypesToRegimes", model)
}

#' @export
PCMMapModelTypesToRegimes.PCM <- function(model, tree, ...) {
  uniqueRegimes <- PCMTreeGetPartNames(tree)
  res <- rep(1, length(uniqueRegimes))
  names(res) <- as.character(uniqueRegimes)
  res
}


#' Get a vector of all parameters (real and discrete) describing a model on a
#' tree including the numerical parameters of each model regime, the integer ids
#' of the spliting nodes defining the regimes on the tree and the integer ids of
#' the model classes associated with each regime.
#'
#' @details This is an S3 generic.
#' In the default implementation, the last entry in the returned vector is the
#' number of numerical parameters. This is used to identify the starting positions
#' in the vector of the first splitting node.
#'
#' @param model a PCM model
#' @param tree a phylo object with an edge.part member.
#' @param ... additional parameters passed to methods.
#' @return a numeric vector concatenating the result
#' @export
PCMGetVecParamsRegimesAndModels <- function(model, tree, ...) {
  UseMethod("PCMGetVecParamsRegimesAndModels", model)
}

#' @export
PCMGetVecParamsRegimesAndModels.PCM <- function(model, tree, ...) {
  numericParams <- PCMParamGetShortVector(model)
  startingNodesRegimes <- PCMTreeGetPartition(tree)
  models <- PCMMapModelTypesToRegimes(model, tree, ...)
  c(numericParams, startingNodesRegimes, models)
}


#' Map a parametrization to its original form.
#' @description This is an S3 generic that transforms the passed argument by
#' applying the transformation rules for its S3 class.
#' @param o a PCM object or a parameter
#' @param ... additional arguments that can be used by implementing methods.
#'
#' @return a transformed version of o.
#'
#' @description This is an S3 generic. See `PCMApplyTransformation._CholeskyFactor`
#' for an example.
#'
#' @details This function returns the same object if it is not transformable.
#'
#' @seealso \code{\link{is.Transformable}}
#'
#' @export
PCMApplyTransformation <- function(o, ...) {
  if(is.Transformable(o)) {
    # this if-statement prevents a transformation when it is not needed, e.g.
    # in the case of a parameter, which is a CholeskyFactor but should not be
    # converted into a positive definite matrix, either, because this is not needed,
    # i.e. the parameter is defined as an upper triangular matrix with a non-negative
    # diagonal, or because the conversion is done at a lower level.
    UseMethod("PCMApplyTransformation", o)
  } else {
    o
  }
}

#' @export
PCMApplyTransformation.default <- function(o, ...) {
  # do nothing (act as an identity function)
  o
}

#' @export
PCMApplyTransformation.PCM <- function(o, ...) {
  if(is.Transformable(o)) {

    # Previously, the following call was used but was too slow, particularly
    # in the case of MixedGaussian models, because of their deeper nested
    # structure:
    # PCMParamSetByName(o, lapply(o, PCMApplyTransformation, ...), replaceWholeParameters = TRUE)
    # The for loop below is faster, but does no checks on the transformed object:
    for(i in seq_along(o)) {
      o[[i]] <- PCMApplyTransformation(o[[i]], ...)
    }
    classes <- class(o)
    classes <- classes[classes != "_Transformable"]
    classes <- c(classes, "_Transformed")

    class(o) <- classes
    o
  } else {
    o
  }
}

#' Conditional distribution of a daughter node given its parent node
#' @description An S3 generic function that has to be implemented for every
#'  model class.
#' @inheritParams PCMLik
#' @param r an integer specifying a model regime
#' @return an object of type specific to the type of model
#' @export
PCMCond <- function(tree, model, r=1, metaI=PCMInfo(NULL, tree, model, verbose = verbose), verbose = FALSE) {
  UseMethod("PCMCond", model)
}

#' Expected mean vector at each tip conditioned on a trait-value vector at the root
#' @inheritParams PCMLik
#' @param X0 a k-vector denoting the root trait
#' @param internal a logical indicating ig the per-node mean vectors should be returned (see Value). Default FALSE.
#'
#' @return If internal is FALSE (default), then a k x N matrix Mu, such that \code{Mu[, i]} equals the expected mean k-vector
#' at tip i, conditioned on \code{X0} and the tree. Otherwise, a k x M matrix Mu containing the mean vector for each node.
#' @export
#' @examples
#' # a Brownian motion model with one regime
#' modelBM <- PCM(model = "BM", k = 2)
#' # print the model
#' modelBM
#' # assign the model parameters at random: this will use uniform distribution
#' # with boundaries specified by PCMParamLowerLimit and PCMParamUpperLimit
#' # We do this in two steps:
#' # 1. First we generate a random vector. Note the length of the vector equals PCMParamCount(modelBM)
#' randomParams <- PCMParamRandomVecParams(modelBM, PCMNumTraits(modelBM), PCMNumRegimes(modelBM))
#' randomParams
#' # 2. Then we load this random vector into the model.
#' PCMParamLoadOrStore(modelBM, randomParams, 0, PCMNumTraits(modelBM), PCMNumRegimes(modelBM), TRUE)
#'
#' # create a random tree of 10 tips
#' tree <- ape::rtree(10)
#' PCMMean(tree, modelBM)
PCMMean <- function(tree, model, X0 = model$X0, metaI=PCMInfo(NULL, tree, model, verbose = verbose), internal = FALSE, verbose = FALSE)  {
  UseMethod("PCMMean", model)
}

#' Calculate the mean at time t, given X0, under a PCM model
#' @param t positive numeric denoting time
#' @param model a PCM model object
#' @param X0 a numeric vector of length k, where k is the number of traits in the model (Defaults to model$X0).
#' @param regime an integer or a character denoting the regime in model for
#' which to do the calculation; Defaults to PCMRegimes(model)[1L], meaning the
#' first regime in the model.
#' @param verbose a logical indicating if (debug) messages should be written on the console (Defaults to FALSE).
#' @return A numeric vector of length k
#' @export
#' @examples
#' # a Brownian motion model with one regime
#' modelBM <- PCM(model = "BM", k = 2)
#' # print the model
#' modelBM
#' # assign the model parameters at random: this will use uniform distribution
#' # with boundaries specified by PCMParamLowerLimit and PCMParamUpperLimit
#' # We do this in two steps:
#' # 1. First we generate a random vector. Note the length of the vector equals PCMParamCount(modelBM)
#' randomParams <- PCMParamRandomVecParams(modelBM, PCMNumTraits(modelBM), PCMNumRegimes(modelBM))
#' randomParams
#' # 2. Then we load this random vector into the model.
#' PCMParamLoadOrStore(modelBM, randomParams, 0, PCMNumTraits(modelBM), PCMNumRegimes(modelBM), TRUE)
#'
#' # PCMMeanAtTime(1, modelBM)
#'
#' # note that the variance at time 0 is not the 0 matrix because the model has a non-zero
#' # environmental deviation
#' PCMMeanAtTime(0, modelBM)
PCMMeanAtTime <- function(
  t, model, X0 = model$X0, regime = PCMRegimes(model)[1L], verbose = FALSE) {
  if(! (regime %in% PCMRegimes(model)) ) {
    stop("PCMMeanAtTime:: regime should be among PCMRegimes(model).")
  }

  cherry <- list(tip.label = c("t1", "t2"),
                 edge = rbind(c(3L, 1L), c(3L, 2L)),
                 edge.length = rep(t, 2),
                 Nnode = 1L)
  class(cherry) <- "phylo"
  cherry <- PCMTree(cherry)
  PCMTreeSetPartRegimes(cherry, c(`3` = regime[1]))

  metaI <- PCMInfo(X = NULL, tree = cherry, model = model, verbose = verbose)

  MeanCherry <- PCMMean(
    tree = cherry,
    model = model,
    X0 = X0,
    metaI = metaI,
    verbose = verbose)

  MeanCherry[, 1L]
}


#' Expected variance-covariance matrix for each couple of tips (i,j)
#' @inheritParams PCMLik
#' @param W0 a numeric matrix denoting the initial k x k variance covariance matrix at the
#'  root (default is the k x k zero matrix).
#' @param internal a logical indicating if the per-node variance-covariances matrices for
#' the internal nodes should be returned (see Value). Default FALSE.
#' @return If internal is FALSE, a (k x N) x (k x N) matrix W, such that k x k block
#' \code{W[((i-1)*k)+(1:k), ((j-1)*k)+(1:k)]} equals the expected
#' covariance matrix between tips i and j. Otherwise, a list with an element 'W' as described above and
#' a k x M matrix element 'Wii' containing the per-node variance covariance matrix for each node:
#' The k x k block \code{Wii[, (i-1)*k + (1:k)]} represents the variance covariance matrix for node i.
#' @export
#' @examples
#' # a Brownian motion model with one regime
#' modelBM <- PCM(model = "BM", k = 2)
#' # print the model
#' modelBM
#' # assign the model parameters at random: this will use uniform distribution
#' # with boundaries specified by PCMParamLowerLimit and PCMParamUpperLimit
#' # We do this in two steps:
#' # 1. First we generate a random vector. Note the length of the vector equals PCMParamCount(modelBM)
#' randomParams <- PCMParamRandomVecParams(modelBM, PCMNumTraits(modelBM), PCMNumRegimes(modelBM))
#' randomParams
#' # 2. Then we load this random vector into the model.
#' PCMParamLoadOrStore(modelBM, randomParams, 0, PCMNumTraits(modelBM), PCMNumRegimes(modelBM), TRUE)
#'
#' # create a random tree of 10 tips
#' tree <- ape::rtree(10)
#' covMat <- PCMVar(tree, modelBM)
PCMVar <- function(
  tree, model,
  W0 = matrix(0.0, PCMNumTraits(model), PCMNumTraits(model)),
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI=PCMInfo(NULL, tree, model, verbose = verbose),
  internal = FALSE, verbose = FALSE)  {
  UseMethod("PCMVar", model)
}

#' Calculate the variance covariance k x k matrix at time t, under a PCM model
#' @param t positive numeric denoting time
#' @param model a PCM model object
#' @param W0 a numeric matrix denoting the initial k x k variance covariance matrix at the
#'  root (default is the k x k zero matrix).
#' @param SE a k x k matrix specifying the upper triangular Cholesky factor of
#' the measurement error variance-covariance matrix. The product
#' t(SE) %*% SE is added to the variance calculated from the model.
#' Default: SE = matrix(0.0, PCMNumTraits(model), PCMNumTraits(model)).
#' @param regime an integer or a character denoting the regime in model for
#' which to do the calculation; Defaults to PCMRegimes(model)[1L], meaning the
#' first regime in the model.
#' @param verbose a logical indicating if (debug) messages should be written on the console (Defaults to FALSE).
#' @return A numeric k x k matrix
#' @export
#' @examples
#' # a Brownian motion model with one regime
#' modelBM <- PCM(model = "BM", k = 2)
#' # print the model
#' modelBM
#' # assign the model parameters at random: this will use uniform distribution
#' # with boundaries specified by PCMParamLowerLimit and PCMParamUpperLimit
#' # We do this in two steps:
#' # 1. First we generate a random vector. Note the length of the vector equals PCMParamCount(modelBM)
#' randomParams <- PCMParamRandomVecParams(modelBM, PCMNumTraits(modelBM), PCMNumRegimes(modelBM))
#' randomParams
#' # 2. Then we load this random vector into the model.
#' PCMParamLoadOrStore(modelBM, randomParams, 0, PCMNumTraits(modelBM), PCMNumRegimes(modelBM), TRUE)
#'
#' # PCMVarAtTime(1, modelBM)
#'
#' # note that the variance at time 0 is not the 0 matrix because the model has a non-zero
#' # environmental deviation
#' PCMVarAtTime(0, modelBM)
PCMVarAtTime <- function(
  t, model,
  W0 =  matrix(0.0, PCMNumTraits(model), PCMNumTraits(model)),
  SE = matrix(0.0, PCMNumTraits(model), PCMNumTraits(model)),
  regime = PCMRegimes(model)[1L], verbose = FALSE) {

  if(! (regime %in% PCMRegimes(model)) ) {
    stop("PCMVarAtTime:: regime should be among PCMRegimes(model).")
  }

  cherry <- list(tip.label = c("t1", "t2"),
                 edge = rbind(c(3L, 1L), c(3L, 2L)),
                 edge.length = rep(t, 2),
                 Nnode = 1L)
  class(cherry) <- "phylo"
  cherry <- PCMTree(cherry)
  PCMTreeSetPartRegimes(cherry, c(`3` = regime[1]))

  metaI <- PCMInfo(
    X = NULL, tree = cherry, model = model, verbose = verbose)

  VarCherry <- PCMVar(
    tree = cherry, model = model, W0 = W0,
    metaI = metaI, verbose = verbose)

  k <- PCMNumTraits(model)
  VarCherry[1:k, 1:k] + (SE %*% t(SE))
}


#' Generate a trajectory for the mean in one regime of a PCM
#'
#' @param model a PCM object.
#' @param regime a regime in `model`. Default is PCMRegimes(model)[1].
#' @param X0 a numeric vector specifying an initial point in the trait space.
#' Default is rep(0, PCMNumTraits(model))
#' @param W0 a numeric k x k symmetric positive definite matrix or 0 matrix,
#' specifying the initial variance covariance matrix at t0. By default, this is
#' a k x k 0 matrix.
#' @param tX,tVar numeric vectors of positive points in time sorted in
#' increasing order. tX specifies the points in time at which to calculate the
#' mean (conditional on X0). tVar specifies a subset of the points in tX at
#' which to generate random samples from the k-variate Gaussian distribution
#' with mean equal to the mean value at the corresponding time conditional on X0
#' and variance equal to the variance at this time, conditional on W0. Default
#' settings are `tX = seq(0, 100, by = 1)` and
#' `tVar = tX[seq(0, length(tX), length.out = 4)]`.
#' @param dims an integer vector specifying the traits for which samples at tVar
#' should be generated (see tX,tVar above).
#' Default: seq_len(PCMNumTraits(model)).
#' @param sizeSamp an integer specifying the number points in the random
#' samples (see tX and tVar above). Default 100.
#' @param doPlot2D Should a ggplot object be produced and returned. This is
#' possible only for two of the traits specified in dims.  Default: FALSE.
#' @param plot a ggplot object. This can be specified when doPlot2D is TRUE and
#' allows to add the plot of this trajectory as a layer in an existing ggplot.
#' Default: NULL
#'
#' @return if doPlot2D is TRUE, returns a ggplot. Otherwise a named list of two
#' elements:
#' \itemize{
#' \item{dt }{a data.table with columns 'regime', 't', 'X', 'V' and 'samp'. For
#' each row corresponding to time in tVar, the column samp represents a list of
#' sizeSamp k-vectors.}
#' \item{dtPlot }{a data.table with the same data as in dt, but with converted
#' columns X and samp into 2 x k columns denoted xi, i=1,...,k and xsi (i=1...k)
#' This is suitable for plotting with ggplot.}}
#' @importFrom ggplot2 ggplot scale_color_continuous geom_path aes stat_ellipse arrow
#' @importFrom data.table data.table rbindlist
#' @importFrom mvtnorm rmvnorm
#' @export
#' @examples
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#'
#' # a Brownian motion model with one regime
#' modelOU <- PCM(model = PCMDefaultModelTypes()['F'], k = 2)
#'
#' # assign the model parameters at random: this will use uniform distribution
#' # with boundaries specified by PCMParamLowerLimit and PCMParamUpperLimit
#' # We do this in two steps:
#' # 1. First we generate a random vector. Note the length of the vector equals
#' # PCMParamCount(modelBM).
#'
#' randomParams <- PCMParamRandomVecParams(
#'   modelOU, PCMNumTraits(modelOU), PCMNumRegimes(modelOU))
#' # 2. Then we load this random vector into the model.
#' PCMParamLoadOrStore(
#'   modelOU,
#'   randomParams,
#'   0, PCMNumTraits(modelBM), PCMNumRegimes(modelBM), load = TRUE)
#'
#' # let's plot the trajectory of the model starting from X0 = c(0,0)
#' PCMTrajectory(
#'   model = modelOU,
#'   X0 = c(0, 0),
#'   doPlot2D = TRUE)
#'
#'
#' # A faceted grid of plots for the two regimes in a mixed model:
#' pla <- PCMTrajectory(
#'   model = PCMBaseTestObjects$model_MixedGaussian_ab, regime = "a",
#'   X0 = c(0, 0, 0),
#'   doPlot2D = TRUE) +
#'   ggplot2::scale_y_continuous(limits = c(0, 10)) +
#'   ggplot2::facet_grid(.~regime)
#'
#' plb <- PCMTrajectory(
#'   model = PCMBaseTestObjects$model_MixedGaussian_ab, regime = "b",
#'   X0 = c(0, 0, 0),
#'   doPlot2D = TRUE) +
#'   ggplot2::scale_y_continuous(limits = c(0, 10)) +
#'   ggplot2::facet_grid(.~regime) +
#'   ggplot2::theme(
#'     axis.title.y = ggplot2::element_blank(),
#'     axis.text.y = ggplot2::element_blank(),
#'     axis.ticks.y = ggplot2::element_blank())
#' cowplot::plot_grid(pla, plb)
PCMTrajectory <- function(
  model,
  regime = PCMRegimes(model)[1],
  X0 = rep(0, PCMNumTraits(model)),
  W0 = matrix(0.0, nrow = PCMNumTraits(model), ncol = PCMNumTraits(model)),

  tX = seq(0, 100, by = 1),
  tVar = tX[seq(0, length(tX), length.out = 4)],
  dims = seq_len(PCMNumTraits(model)),
  sizeSamp = 100,
  doPlot2D = FALSE,
  plot = NULL) {

  dt <- data.table(
    regime = regime,
    t = tX,
    X = list(NULL),
    V = list(NULL),
    samp = list(NULL))

  for(i in seq_len(nrow(dt))) {
    dt[i, c("regime", "X", "V", "samp"):={
      X <- PCMMeanAtTime(
        t = t,
        model = model,
        X0 = X0,
        regime = regime)
      V = PCMVarAtTime(
        t = t,
        W0 = W0,
        model = model,
        regime = regime)
      if(t %in% tVar) {
        samp <- rmvnorm(sizeSamp, sigma = V)
      } else {
        samp <- NULL
      }
      list(regime = regime, X = list(X), V = list(V), samp = list(samp))
    }]
  }

  dtPlot <- rbindlist(lapply(seq_len(nrow(dt)), function(i) {
    dt[i,
       c(list(t = t, regime = regime),
         lapply(dims, function(d) {
           X[[1]][d]
         }),
         lapply(dims, function(d) {
           if(is.null(samp[[1]])) {
             NA_real_
           } else {
             samp[[1]][, d] + X[[1]][d]
           }
         }))]
  }))
  names(dtPlot) <- c("t", "regime", paste0("x", dims), paste0("xs", dims))

  # avoid warning from r check
  x1 <- x2 <- xs1 <- xs2 <- NULL

  if(doPlot2D) {
    if(is.null(plot)) {
      pl <- ggplot(NULL)
    } else {
      pl <- plot
    }
    pl <- pl + geom_path(
      data = dtPlot,
      mapping = aes(x = x1, y = x2))
    pl <- pl +
      stat_ellipse(
        data = dtPlot,
        aes(x = xs1, y = xs2, group = factor(t)),
        type="norm")
    pl
  } else {
    list(dt = dt, dtPlot = dtPlot)
  }
}

#' Calculate the likelihood of a model using the standard formula for multivariate pdf
#' @inheritParams PCMLik
#' @return a numerical value with named attributes as follows:
#' @importFrom mvtnorm dmvnorm
#' @export
PCMLikDmvNorm <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(X, tree, model, SE, verbose = verbose),
  log = TRUE,
  verbose = FALSE) {

  dmvnorm(as.vector(X[, 1:PCMTreeNumTips(tree)]),
          as.vector(PCMMean(tree, model, model$X0, metaI = metaI)),
          PCMVar(tree, model, metaI = metaI), log = log)
}


#' Simulation of a phylogenetic comparative model on a tree
#'
#' @description Generate trait data on a tree according to a multivariate stochastic
#' model with one or several regimes
#'
#' @param X0 a numeric vector of length k (the number of traits) specifying the
#' trait values at the root of the tree.
#' @param tree a phylo object specifying a rooted tree.
#' @param model an S3 object specifying the model (see Details).
#' @param SE a k x N matrix specifying the standard error for each measurement in
#' X. Alternatively, a k x k x N cube specifying an upper triangular k x k
#' Cholesky factor of the variance covariance matrix for the measurement error
#' for each node i=1, ..., N.
#' Default: \code{matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree))}.
#' @param metaI a named list containg meta-information about the data and the
#' model.
#' @param verbose a logical indicating if informative messages should be written
#' during execution.
#'
#' @details Internally, this function uses the \code{\link{PCMCond}} iimplementation
#'  for the given model class.
#'
#' @return numeric M x k matrix of values at all nodes of the tree, i.e. root,
#' internal and tip, where M is the number of nodes: \code{M=dim(tree$edge)[1]+1},
#' with indices from 1 to N=length(tree$tip.label) corresponding to tips, N+1
#' corresponding to the root and bigger than N+1 corresponding to internal nodes.
#' The function will fail in case that the length of the argument vector X0 differs
#' from the number of traits specified in \code{metaI$k}. Error message:
#' "PCMSim:: X0 must be of length ...".
#'
#' @importFrom mvtnorm rmvnorm
#' @seealso \code{\link{PCMLik}} \code{\link{PCMInfo}} \code{\link{PCMCond}}
#' @examples
#' library(data.table)
#' N <- 10
#' L <- 100.0
#' tr <- ape::stree(N)
#' tr$edge.length <- rep(L, N)
#' for(epoch in seq(1, L, by = 1.0)) {
#'   tr <- PCMTreeInsertSingletonsAtEpoch(tr, epoch)
#' }
#'
#' model <- PCMBaseTestObjects$model_MixedGaussian_ab
#'
#' PCMTreeSetPartRegimes(tr, c(`11` = 'a'), setPartition = TRUE)
#'
#' set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' X <- PCMSim(tr, model, X0 = rep(0, 3))
#'
#' dt <- NULL
#' for(epoch in seq(0, L, by = 1)) {
#'   nodes <- PCMTreeLocateEpochOnBranches(tr, epoch)$nodes
#'   dtEpoch <- as.data.table(t(X[, nodes]))
#'   dtEpoch[, t:=epoch]
#'   if(epoch == 0) {
#'     dtEpoch[, lineage:="root"]
#'   } else {
#'     dtEpoch[, lineage:=gsub("i.*x", "x", PCMTreeGetLabels(tr)[nodes], perl = TRUE)]
#'   }
#'
#'   dt <- rbindlist(list(dt, dtEpoch))
#' }
#'
#' @export
PCMSim <- function(
  tree, model, X0,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(
    X = NULL, tree = tree, model = model, SE = SE, verbose = verbose),
  verbose = FALSE) {

  UseMethod("PCMSim", model)
}

#' Likelihood of a multivariate Gaussian phylogenetic comparative model with non-interacting lineages
#'
#' @description The likelihood of a PCM represets the probability density function
#'   of observed trait values (data) at the tips of a tree given the tree and
#'   the model parameters. Seen as a function of the model parameters, the
#'   likelihood is used to fit the model to the observed trait data and the
#'   phylogenetic tree (which is typically inferred from another sort of data, such
#'   as an alignment of genetic sequences for the species at the tips of the tree).
#'   The \code{\link{PCMLik}} function
#'   provides a common interface for calculating the (log-)likelihood of different
#'   PCMs.
#'   Below we denote by N the number of tips, by M the total number of nodes in the
#'   tree including tips, internal and root node, and by k - the number of traits.
#'
#' @param X a \code{k x N} numerical matrix with possible \code{NA} and \code{NaN} entries. Each
#'   column of X contains the measured trait values for one species (tip in tree).
#'   Missing values can be either not-available (\code{NA}) or not existing (\code{NaN}).
#'   Thse two values have are treated differently when calculating
#'   likelihoods: see \code{\link{PCMPresentCoordinates}}.
#' @param tree a phylo object with N tips.
#' @param model an S3 object specifying both, the model type (class, e.g. "OU") as
#'   well as the concrete model parameter values at which the likelihood is to be
#'   calculated (see also Details).
#' @param SE a k x N matrix specifying the standard error for each measurement in
#' X. Alternatively, a k x k x N cube specifying an upper triangular k x k
#' Cholesky factor of the variance covariance matrix for the measurement error
#' for each node i=1, ..., N.
#' Default: \code{matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree))}.
#' @param metaI a list returned from a call to \code{PCMInfo(X, tree, model, SE)},
#'   containing meta-data such as N, M and k. Alternatively, this can be a
#'   function object that returns such a list, e.g. the function\code{PCMInfo}
#'   or the function \code{PCMInfoCpp} from the \code{PCMBaseCpp} package.
#' @param log logical indicating whether a log-liklehood should be calculated. Default
#'  is TRUE.
#' @param verbose logical indicating if some debug-messages should printed.
#'
#' @return a numerical value with named attributes as follows:
#' \describe{
#' \item{X0}{A numerical vector of length k specifying the value at the root for which
#' the likelihood value was calculated. If the model contains a member called X0, this
#' vector is used; otherwise the value of X0 maximizing the likelihood for the given
#' model parameters is calculated by maximizing the quadratic polynomial
#' 'X0 * L_root * X0 + m_root * X0 + r_root'.}
#' \item{error}{A named list containing error information if a numerical or other
#' logical error occured during likelihood calculation (this is a list returned by
#'  \code{\link{PCMParseErrorMessage}}.}
#'  If an error occured during likelihood calculation, the default behavior is to
#'  return NA with a non-NULL error attribute. This behavior can be changed in
#'  using global options:
#'  \item{"PCMBase.Value.NA"}{Allows to specify a different NA value such as \code{-Inf} or \code{-1e20} which can be used in compbination with \code{log = TRUE} when
#'   using \code{optim} to maximize the log-likelihood;}
#'  \item{"PCMBase.Errors.As.Warnings"}{Setting this option to FALSE will cause any
#'  error to result in calling the \code{\link{stop}} R-base function. If not caught
#'  in a \code{\link{tryCatch}}, this will cause the inference procedure to abort at the occurence of a numerical error. By default, this option is set to TRUE, which
#'  means that \code{getOption("PCMBase.Value.NA", as.double(NA))} is returned with
#'  an error attribute and a warning is issued.}
#'}
#' @details For efficiency, the argument \code{metaI}
#'   can be provided explicitly, because this is not supposed to change during a
#'   model inference procedure such as likelihood maximization.
#'
#' @seealso \code{\link{PCMInfo}} \code{\link{PCMAbCdEf}} \code{\link{PCMLmr}} \code{\link{PCMSim}} \code{\link{PCMCond}} \code{\link{PCMParseErrorMessage}}
#' @export
PCMLik <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(
    X = X, tree = tree, model = model, SE = SE, verbose = verbose),
  log = TRUE,
  verbose = FALSE) {

  UseMethod("PCMLik", model)
}

#' Tracing the log-likelihood calculation of a model over each node of the tree
#'
#' @description This is an S3 generic function providing tracing information
#' for the likelihood calculation for a given tree, data and model parameters.
#' Useful for illustration or for debugging purpose.
#'
#' @inheritParams PCMLik
#'
#' @return The returned object will, in general, depend on the type of model and
#' the algorithm used for likelihood calculation. For a G_LInv model and
#' pruning-wise likelihood calculation, the returned object will be a data.table
#' with columns corresponding to the node-state variables, e.g. the quadratic
#' polynomial coefficients associated with each node in the tree.
#' @seealso \code{\link{PCMInfo}} \code{\link{PCMAbCdEf}} \code{\link{PCMLmr}} \code{\link{PCMSim}} \code{\link{PCMCond}} \code{\link{PCMParseErrorMessage}}
#' @export
PCMLikTrace <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(
    X = X, tree = tree, model = model, SE = SE, verbose = verbose),
  log = TRUE,
  verbose = FALSE) {

  UseMethod("PCMLikTrace", model)
}

#' @export
logLik.PCM <- function(object, ...) {
  if(!is.PCM(object)) {
    stop("logLik.PCM:: object must inherit from class PCM.")
  }

  X <- attr(object, "X", exact = TRUE)
  if( !is.matrix(X) ) {
    stop("logLik.PCM:: When calling logLik.PCM on a model object, it should have a k x N numeric matrix attribute called 'X'.")
  }
  SE <- attr(object, "SE", exact = TRUE)
  if( !(is.matrix(SE) || length(dim(SE)) == 3) ) {
    warning("logLik.PCM:: When calling logLik.PCM on a model object, it should have a (k x N) numeric matrix or a (k x k x N) numeric cube attribute called 'SE'. Setting SE to a (k x N) zero matrix.")
    SE <- X * 0.0
  }
  tree <- attr(object, "tree", exact = TRUE)
  if( !inherits(tree, "phylo") ) {
    stop("logLik.PCM:: When calling logLik.PCM on a model object should have an attribute called 'tree' of class phylo.")
  }

  if(is.function(attr(object, "PCMInfoFun", exact = TRUE)) ) {
    value <- PCMLik(
      X, tree, object, SE, metaI = attr(object, "PCMInfoFun", exact = TRUE)(X, tree, object, SE), log = TRUE)
  } else if(is.list(attr(object, "PCMInfoFun", exact = TRUE))) {
    # In this case, it is assumed that attr(object, "PCMInfoFun", exact = TRUE) is
    # the result from calling the PCMInfoFun on the model (object) X, tree and data.
    value <- PCMLik(
      X, tree, object, SE, metaI = attr(object, "PCMInfoFun", exact = TRUE), log = TRUE)
  } else {
    value <- PCMLik(X, tree, object, SE, log = TRUE)
  }

  attr(value, "df") <- PCMParamCount(object, countRegimeChanges = TRUE, countModelTypes = TRUE)
  attr(value, "nobs") <- PCMTreeNumTips(tree)

  value
}


#' Determine which traits are present (active) on each node of the tree
#'
#' @description For every node (root, internal or tip) in tree, build a logical
#' vector of length k with TRUE values for every present coordinate. Non-present
#' coordinates arize from NA-values in the trait data. These can occur in two cases:
#' \describe{
#' \item{Missing measurements for some traits at some tips:}{the present coordinates
#' are FALSE for the corresponding tip and trait, but are full for all traits
#' at all internal and root nodes.}
#' \item{non-existent traits for some species:}{the FALSE present coordinates
#' propagate towards the parent
#' nodes - an internal or root node will have a present coordinate set to FALSE
#' for a given trait, if all of its descendants have this coordinate set to FALSE.}
#' }
#' These two cases have different effect on the likelihood calculation: missing
#' measurements (NA) are integrated out at the parent nodes; while non-existent
#' traits (NaN) are treated as reduced dimensionality of the vector at the
#' parent node.
#'
#' @param X numeric k x N matrix of observed values, with possible NA entries. The
#' columns in X are in the order of tree$tip.label
#' @param tree a phylo object
#' @param metaI  The result of calling PCMInfo.
#' @return a k x M logical matrix. The function fails in case when all traits are NAs for some of the tips.
#' In that case an error message is issued
#' "PCMPresentCoordinates:: Some tips have 0 present coordinates. Consider
#' removing these tips.".
#' @seealso \code{\link{PCMLik}}
#' @export
PCMPresentCoordinates <- function(X, tree, metaI) {

  if(is.null(metaI)) {
    N <- PCMTreeNumTips(tree)
    M <- PCMTreeNumNodes(tree)
    k <- nrow(X)
    postorder <- PCMTreePostorder(tree)
  } else {
    N <- metaI$N
    M <- metaI$M
    k <- metaI$k
    postorder <- rev(metaI$preorder)
  }

  edge <- tree$edge

  if(is.null(X)) {
    pc <- matrix(TRUE, k, M)
  } else {
    pc <- matrix(FALSE, k, M)
    for(ei in postorder) {
      i <- edge[ei, 2]
      j <- edge[ei, 1]

      if(i <= N) {
        # all es pointing to tips
        # here we count NAs as present because they are integrated out at the
        # parent nodes.
        pc[, i] <- !is.nan(X[, i])
      } else {
        # edges pointing to internal nodes, for which all children nodes have
        # been visited here we do nothing.
      }

      #update parent pc
      pc[, j] <- pc[, j] | pc[, i]
    }
    # at the end correct the present coordinates for the traits at the tips which
    # are NA but not NaN: Both NAs and NaNs result in FALSE at a tip;
    # only coordinates for which all descending tips have NaN are FALSE at
    # internal nodes.
    pc[, seq_len(N)] <- !is.na(X[, seq_len(N)])

    if(any(rowSums(pc) == 0)) {
      stop("PCMPresentCoordinates:: Some tips have 0 present coordinates.
           Consider removing these tips.")
    }
  }
  pc
}

#' Meta-information about a tree and trait data associated with a PCM
#'
#' @description This function pre-processes the given tree and data in order to
#' create meta-information used during likelihood calculaiton.
#' @inheritParams PCMLik
#' @param preorder an integer vector of row-indices in tree$edge matrix as returned
#' by PCMTreePreorder. This can be given for performance speed-up when several
#' operations needing preorder are executed on the tree. Default : \code{NULL}.
#' @param ... additional arguments used by implementing methods.
#'
#' @return a named list with the following elements:
#' \item{X}{k x N matrix denoting the trait data;}
#' \item{VE}{k x k x N array denoting the measurement error variance covariance
#' matrix for each for each tip i = 1,...,N. See the parameter \code{SE} in
#' \code{\link{PCMLik}.}}
#' \item{M}{total number of nodes in the tree;}
#' \item{N}{number of tips;}
#' \item{k}{number of traits;}
#' \item{RTree}{number of parts on the tree (distinct elements of tree$edge.part);}
#' \item{RModel}{number of regimes in the model (elements of attr(model, regimes));}
#' \item{p}{number of free parameters describing the model;}
#' \item{r}{an integer vector corresponding to tree$edge with the regime for each
#' branch in tree;}
#' \item{xi}{an integer vector of 0's and 1's corresponding to the rows in tree$edge
#' indicating the presence of a jump at the corresponding branch;}
#' \item{pc}{a logical matrix of dimension k x M denoting the present coordinates
#' for each node; in special cases this matrix can be edited by hand after calling
#' PCMInfo and before passing the returned list to PCMLik. Otherwise, this matrix
#' can be calculated in a custom way by specifying the option PCMBase.PCMPresentCoordinatesFun.
#' See also \code{\link{PCMPresentCoordinates}} and \code{\link{PCMOptions}}. }
#' This list is passed to \code{\link{PCMLik}}.
#'
#' @export
PCMInfo <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  verbose = FALSE, preorder = NULL, ...) {
  UseMethod("PCMInfo", model)
}

#' @export
PCMInfo.PCM <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  verbose = FALSE, preorder = NULL, ...) {

  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }

  tree <- PCMTree(tree)

  if(is.null(preorder)) {
    preorder <- PCMTreePreorder(tree)
  }

  k <- PCMNumTraits(model)
  N <- PCMTreeNumTips(tree)
  M <- PCMTreeNumNodes(tree)

  VE <- array(0.0, dim = c(k, k, N))
  if(is.matrix(SE) && identical(unname(dim(SE)), c(k, N)) ) {
    # SE is k x N matrix
    for(i in seq_len(N)) {
      VE[, , i] <- diag(SE[, i]*SE[, i], nrow = k)
    }
  } else if(is.array(SE) && identical(unname(dim(SE)), c(k, k, N)) ) {
    # SE is k x k x N array
    for(i in seq_len(N)) {
      VE[, , i] <- t(SE[, , i]) %*% SE[, , i]
    }
  } else {
    stop("SE should be a k x N matrix or a k x k x N array.")
  }

  res <- list(
    X = X,
    VE = VE,
    M = M,
    N = N,
    k = k,
    RTree = PCMTreeNumParts(tree),
    RModel = PCMNumRegimes(model),
    r = {
      regimeIndices <- match(PCMTreeGetRegimesForEdges(tree), PCMRegimes(model))
      if(any(is.na(regimeIndices))) {
        stop(paste0(
          "PCMInfo:: Some of the regimes for the edges in tree could not be",
          "matched with regimes in PCMRegimes(model)."))
      }
      regimeIndices
    },
    p = PCMParamCount(model),
    xi = PCMTreeJumps(tree),
    edge = tree$edge,
    edge.length = tree$edge.length,
    preorder = preorder
  )
  res <- c(res, PCMOptions())

  PCMPresentCoordinatesFun <- getOption(
    "PCMBase.PCMPresentCoordinatesFun", PCMPresentCoordinates)
  res$pc <- PCMPresentCoordinatesFun(X, tree, res)
  res$NA_double_ <- NA_real_

  res
}

#' Create a likelhood function of a numerical vector parameter
#' @inheritParams PCMLik
#' @param positiveValueGuard positive numerical value (default Inf), which
#' serves as a guard for numerical error. Values exceeding
#' this positiveGuard are most likely due to numerical error and
#' PCMOptions()$PCMBase.Value.NA is returned instead.
#' @return a function of a numerical vector parameter called p returning the
#' likelihood of X given the tree and the model with parameter values specified
#' by p.
#' @details It is possible to specify a function for the argument metaI. This
#' function should have three parameters (X, tree, model) and should return a
#' metaInfo object. (see \code{\link{PCMInfo}}).
#'
#' @export
PCMCreateLikelihood <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(X, tree, model, SE),
  positiveValueGuard = Inf) {

  if(is.function(metaI)) {
    metaI <- metaI(X, tree, model, SE)
  }
  value.NA <- PCMOptions()$PCMBase.Value.NA

  function(p, log = TRUE) {
    PCMParamLoadOrStore(model, p, offset = 0L,
                        k = PCMNumTraits(model),
                        R = PCMNumRegimes(model),
                        load = TRUE)
    value <- PCMLik(X, tree, model, SE, metaI, log = log)
    if(is.na(value) || value > positiveValueGuard) {
      value <- value.NA
    }
    value
  }
}

#' Find the S3 method for a given PCM object or class-name and an S3 generic
#' @param x a character string denoting a PCM S3 class name (e.g. "OU"), or a
#' PCM object.
#' @param method a character string denoting the name of an S3 generic function.
#' Default: "PCMCond".
#' @return a function object corresponding to the S3 method found or an error is
#' raised if no such function is found for the specified object and method.
#'
#' @importFrom utils getS3method
#' @export
PCMFindMethod <- function(x, method = "PCMCond") {
  if(is.character(x)) {
    o <- try(PCM(x), silent = TRUE)
    if(inherits(o, "try-error")) {
      stop(
        paste0(
          "PCMFindMethod:: ", toString(o), " PCM constructor called on x",
          "='", x, "'")
      )
    }
  } else if(is.PCM(x)) {
    o <- x
  } else {
    stop(
      paste0(
        "PCMFindMethod:: x should be a character string denoting a",
        "PCM class name or a PCM object."))
  }

  cls <- c(class(o), 'default')
  results <- lapply(cls, function(y) try(getS3method(method, y), silent = TRUE))
  Find(function (x) class(x) != 'try-error', results)
}

#' Given a PCM or a parameter object, extract an analogical object for a subset
#' of the dimensions (traits) in the original object.
#'
#' @details This is an S3 generic
#' @param obj a PCM or a parameter object.
#' @param dims an integer vector; should be a subset or equal to
#' \code{seq_len(PCMNumTraits(obj))} (the default).
#' @param nRepBlocks a positive integer specifying if the specified dimensions
#' should be replicated to obtain a higher dimensional model, where the parameter
#' matrices are block-diagonal with blocks corresponding to dims. Default: 1L.
#' @return an object of the same class as obj with a subset of obj's dimensions
#' multiplied \code{nRepBlocks} times.
#' @export
PCMExtractDimensions <- function(
  obj,
  dims = seq_len(PCMNumTraits(obj)),
  nRepBlocks = 1L) {
  UseMethod("PCMExtractDimensions", obj)
}

#' @export
PCMExtractDimensions.PCM <- function(
  obj,
  dims = seq_len(PCMNumTraits(obj)),
  nRepBlocks = 1L) {
  dims <- unique(dims)
  if( !isTRUE(all(dims %in% seq_len(PCMNumTraits(obj)))) ) {
    stop("PCMExtractDimensions.PCM:: Some dims are outside 1:PCMNumTraits(obj).")
  }
  obj2 <- lapply(obj, PCMExtractDimensions, dims = dims, nRepBlocks = nRepBlocks)
  attributes(obj2) <- attributes(obj)
  attr(obj2, "k") <- as.integer(length(dims) * nRepBlocks)
  attr(attr(obj2, "spec"), "k") <- as.integer(length(dims) * nRepBlocks)
  for(name in names(attr(obj2, "spec"))) {
    if(is.PCM(attr(obj2, "spec")[[name]])) {
      attr(attr(obj2, "spec")[[name]], "k") <- as.integer(length(dims) * nRepBlocks)
    }
  }
  attr(obj2, "p") <- PCMParamCount(obj2)
  if(!is.null(attr(obj, "X"))) {
    attr(obj2, "X") <- attr(obj, "X")[rep(dims, nRepBlocks), , drop = FALSE]
  }
  if(!is.null(attr(obj, "SE"))) {
    attr(obj2, "SE") <- if(is.matrix(attr(obj, "SE"))) {
      attr(obj, "SE")[rep(dims, nRepBlocks), , drop = FALSE]
    } else {
      kronecker(diag(1, nRepBlocks), attr(obj, "SE"))
    }
  }
  if(!is.null(attr(obj, "tree"))) {
    attr(obj2, "tree") <- attr(obj, "tree")
  }
  obj2
}

#' Given a PCM or a parameter object, extract an analogical object for a subset
#' of the regimes in the original object.
#'
#' @details This is an S3 generic
#' @param obj a PCM or a parameter object.
#' @param regimes an integer vector; should be a subset or equal to
#' \code{seq_len(PCMNumRegimes(obj))} (the default).
#' @return an object of the same class as obj with a subset of obj's regimes
#' @export
PCMExtractRegimes <- function(obj, regimes = seq_len(PCMNumRegimes(obj))) {
  UseMethod("PCMExtractRegimes", obj)
}

#' @export
PCMExtractRegimes.PCM <- function(obj, regimes = seq_len(PCMNumRegimes(obj))) {
  obj2 <- lapply(obj, PCMExtractRegimes, regimes = regimes)
  attributes(obj2) <- attributes(obj)
  attr(obj2, "regimes") <- attr(obj2, "regimes")[regimes]
  attr(obj2, "p") <- PCMParamCount(obj2)
  attr(attr(obj2, "spec"), "regimes") <- attr(attr(obj2, "spec"), "regimes")[regimes]

  obj2
}
