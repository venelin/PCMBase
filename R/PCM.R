# Copyright 2018 Venelin Mitov
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

#' @name PCMBase
#'
#' @title A General Framework for Gaussian Phylogenetic Comparative Models
#'
#' @description
#'
#'
NULL

#' Get a list of PCM models currently implemented
#' @param pattern a character string specifying an optional for the model-names to search for.
#' @param parentClass a character string specifying an optional parent class of the models to look for.
#' @return a character vector of the model classes found.
#' @details The function is using the S3 api function \code{\link{methods}} looking for all registered
#' implementations of the function \code{\link{PCMSpecify}}
#' @export
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
#' \item{\code{PCMBase.Threshold.Lambda_ij }}{a 0-threshold for abs(Lambda_i + Lambda_j),
#' where Lambda_i and Lambda_j are eigenvalues of the parameter matrix H of an OU or
#' other model. Default 1e-8. See \code{\link{PCMPExpxMeanExp}}.}
#' \item{\code{PCMBase.Threshold.SV }}{A 0-threshold for min(svdV)/max(svdV), where
#' svdV is the vector of eigenvalues of the matrix V for a given branch. The V matrix
#' is considered singular if it has eigenvalues equal to 0 or when the ratio
#' min(svdV)/max(svdV) is below PCMBase.Threshold.SV. Default is 1e-6. Treatment
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
#' \item{\code{PCMBase.ParamValue.UpperLimit}}{Default upper limit value for parameters, default setting is 10.0.}
#' }
#' @export
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
       PCMBase.ParamValue.UpperLimit = getOption("PCMBase.ParamValue.UpperLimit", 10.0)
       )
}

#' @name PCM
#' @title Create a phylogenetic comparative model object
#'
#' @description This is the entry-point function for creating model objects within
#' the PCMBase framework representing a single model-type with one or several
#' model-regimes of this type associated with the branches of a tree. For mixed
#' Gaussian phylogenetic models, which enable multiple model-types, use the
#' \code{\link{MRG}} function.
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
#' @param spec NULL or a list specifying the model parameters (see \code{\link{PCMSpecify}}). If NULL (default), the generic PCMSpecify
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
#' @seealso \code{\link{MRG}}
#' @export
PCM <- function(
  model, modelTypes = class(model)[1], k = 1, regimes = 1L,
  params = NULL, vecParams = NULL, offset = 0,
  spec = NULL, ...) {
  UseMethod("PCM", model)
}

#' @export
PCM.default <- function(
  model, modelTypes = class(model)[1], k = 1, regimes = 1L,
  params = NULL, vecParams = NULL, offset = 0,
  spec = NULL, ...) {
  stop(paste0("ERR:02091:PCMBase:PCM.R:PCM.default:: You should provide a PCM object, but provided a ", class(model)[1]))
}

#' @export
PCM.PCM <- function(
  model, modelTypes = class(model)[1], k = 1, regimes = 1L,
  params = NULL, vecParams = NULL, offset = 0, spec = NULL, ...) {

  if(is.null(spec)) {
    spec <- PCMSpecify(model, ...)
  }

  obj <- PCMDefaultObject(spec, model, ...)

  if(is.null(attr(obj, "k", exact = TRUE))) {
    attr(obj, "k") <- k
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
  model, modelTypes = model[1], k = 1, regimes = 1L,
  params = NULL, vecParams = NULL, offset = 0, spec = NULL, ...) {

  modelObj <- list()
  class(modelObj) <- model
  if(length(model) == 1) {
    class(modelObj) <- c(class(modelObj), PCMParentClasses(modelObj))
  }
  attr(modelObj, "k") <- k
  attr(modelObj, "regimes") <- regimes

  PCM(modelObj, modelTypes, k, regimes, params, vecParams, offset, spec, ...)
}

#' @export
is.PCM <- function(x) inherits(x, "PCM")

#' @export
print.PCM <- function(x, ...) cat (format(x, ...), "\n")

#' @export
format.PCM <- function(x, ...) {
  if( !is.PCM(x) ) {
    stop("ERR:02061:PCMBase:PCM.R:format.PCM:: x must inherit from S3 class PCM.")
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
#' This function is called by PCM.character to determine the parent classes for
#' a given model class.
#' @return a vector of character string denoting the names of the parent classes
#' @seealso \code{\link{PCM.character}}
#' @export
PCMParentClasses <- function(model) {
  UseMethod("PCMParentClasses", model)
}

#' @export
PCMParentClasses.PCM <- function(model) c()

#' Human friendly description of a PCM
#' @param model a PCM model object
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
#' @param model a PCM model object
#' @return a list
#' @export
PCMSpecify <- function(model, ...) {
  UseMethod("PCMSpecify", model)
}

#' @export
PCMDescribeParameters <- function(model, ...) {
  UseMethod("PCMDescribeParameters", model)
}

#' Specify the parameterizations for each parameter of a model
#' @export
PCMListParameterizations <- function(model, ...) {
  UseMethod("PCMListParameterizations", model)
}

#' @importFrom data.table data.table as.data.table
#' @export
PCMTableParameterizations <- function(
  model, listParameterizations = PCMListParameterizations(model, ...), ...) {

  # produce the cartesian product of the parameterizations for each parameter
  dtClasses <- do.call(expand.grid, listParameterizations)
  dtClasses <- as.data.table(dtClasses)
  dtClasses
}


#' @export
PCMGenerateParameterizations <- function(
  model,
  listParameterizations = PCMListParameterizations(model),
  tableParameterizations = PCMTableParameterizations(model, listParameterizations),
  env = .GlobalEnv,
  useModelClassNameForFirstRow = FALSE) {

  if(!is.PCM(model)) {
    # We assume that the passed model's class is a single character not including the
    # corresponding hierarchy, which is specified by PCMParentClasses
    class(model) <- c(class(model)[1], PCMParentClasses(model))
  }

  classesModel <- class(model)
  paramNames <- names(tableParameterizations)
  paramDescriptions <- PCMDescribeParameters(model)

  for(i in 1:nrow(tableParameterizations)) {
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
      "if(is.null(names(spec))) names(spec) <- ", PCMCharacterVectorToRExpression(paramNames), "\n",
      "if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')\n",
      "spec\n",
      "}")

    eval(parse(text = sourcePCMParentClasses), env)
    eval(parse(text = sourcePCMSpecify), env)
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
#' @param model a PCM object
#' @return an integer
#' @export
PCMNumTraits <- function(model) {
  UseMethod("PCMNumTraits", model)
}

#' @export
PCMNumTraits.PCM <- function(model) {
  attr(model, "k")
}

#' Regimes in a model
#' @param model a PCM object
#' @return a character or an integer vector giving the regime names in the model
#' @export
PCMRegimes <- function(model) {
  UseMethod("PCMRegimes", model)
}

#' @export
PCMRegimes.PCM <- function(model) {
  attr(model, "regimes")
}

#' #' Regimes in a model
#' #' @param model a PCM object
#' #' @param tree a phylo object or NULL. If the regimes in the model are integers and tree is not NULL,
#' #' then these integers are used as indexes in PCMTreeUniqueRegimes(tree). Default NULL.
#' #' @return a character or an integer vector giving the regime names of the models
#' #' @export
#' PCMRegimes <- function(model, tree = NULL, preorder = if(is.null(tree)) NULL else PCMTreePreorder(tree)) {
#'   UseMethod("PCMRegimes", model)
#' }
#'
#' #' @export
#' PCMRegimes.PCM <- function(model, tree = NULL, preorder = if(is.null(tree)) NULL else PCMTreePreorder(tree)) {
#'   r <- attr(model, "regimes")
#'   if(is.integer(r) && !is.null(tree)) {
#'     PCMTreeUniqueRegimes(tree, preorder)[r]
#'   } else {
#'     r
#'   }
#' }


#' Integer vector giving the model type index for each regime
#' @param model a PCM model
#' @param tree a phylo object with an edge.regime member
#' @param ... additional parameters passed to methods
#' @return an integer vector with elements corresponding to the elements in
#' \code{PCMTreeUniqueRegimes(tree)}
#' @details This is a generic S3 method. The default implementation for the basic
#' class PCM returns a vector of 1's, because it assumes that a single model type
#' is associated with each regime. The implementation for multi-regime models (MRG)
#' returns the mapping attribute of the MRG object reordered to correspond to
#' \code{PCMTreeUniqueRegimes(tree)}.
#' @export
PCMMapModelTypesToRegimes <- function(model, tree, ...) {
  UseMethod("PCMMapModelTypesToRegimes", model)
}

#' @export
PCMMapModelTypesToRegimes.PCM <- function(model, tree, ...) {
  uniqueRegimes <- PCMTreeUniqueRegimes(tree)
  res <- rep(1, length(uniqueRegimes))
  names(res) <- as.character(uniqueRegimes)
  res
}

#' Number of regimes in a model
#' @param model a PCM object
#' @return an integer
#' @export
PCMNumRegimes <- function(model) {
  UseMethod("PCMNumRegimes", model)
}

#' @export
PCMNumRegimes.PCM <- function(model) {
  length(PCMRegimes(model))
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
#' @param tree a phylo object with an edge.regime member.
#' @param ... additional parameters passed to methods.
#' @return a numeric vector concatenating the result
#' @export
PCMGetVecParamsRegimesAndModels <- function(model, tree, ...) {
  UseMethod("PCMGetVecParamsRegimesAndModels", model)
}

#' @export
PCMGetVecParamsRegimesAndModels.PCM <- function(model, tree, ...) {
  numericParams <- PCMParamGetShortVector(model)
  startingNodesRegimes <- PCMTreeGetStartingNodesRegimes(tree)
  models <- PCMMapModelTypesToRegimes(model, tree, ...)
  c(numericParams, startingNodesRegimes, models)
}


#' Map a parametrization to its original form.
#' @description This is an S3 generic that
#' @export
PCMApplyTransformation <- function(o, ...) {
  if(is.Transformable(o)) {
    # this if-statement prevents a transformation when it is not needed, e.g.
    # in the case of a parameter, which is a CholeskiFactor but should not be
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
    PCMParamSetByName(o, lapply(o, PCMApplyTransformation, ...), replaceWholeParameters = TRUE)
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
PCMCond <- function(tree, model, r=1, metaI=PCMInfo(NULL, tree, model, verbose), verbose = FALSE) {
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
PCMMean <- function(tree, model, X0 = model$X0, metaI=PCMInfo(NULL, tree, model, verbose), internal = FALSE, verbose = FALSE)  {
  UseMethod("PCMMean", model)
}

#' Calculate the mean at time t, given X0, under a PCM model
#' @param t positive numeric denoting time
#' @param model a PCM model object
#' @param X0 a numeric vector of length k, where k is the number of traits in the model (Defaults to model$X0).
#' @param regime an integer or a character denoting the regime in model for which to do the calculation;
#' (Defaults to 1L meaning the first regime in the model)
#' @param verbose a logical indicating if (debug) messages should be written on the console (Defaults to FALSE).
#' @return A numeric vector of length k
#' @export
PCMMeanAtTime <- function(t, model, X0 = model$X0, regime = 1L, verbose = FALSE) {
  if(is.character(regime)) {
    regime <- match(regime, PCMRegimes(model))
  }

  cherry <- list(tip.labels=c("t1", "t2"),
                 edge=rbind(c(3L, 1L), c(3L, 2L)),
                 edge.length=rep(t, 2),
                 Nnode = 1L,
                 edge.regime = rep(regime, 2))
  class(cherry) <- "phylo"

  metaI <- PCMInfo(NULL, cherry, model, verbose = verbose)

  metaI$r <- rep(regime, 2)

  MeanCherry <- PCMMean(cherry, model, X0, metaI = metaI, verbose = verbose)
  MeanCherry[, 1]
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
PCMVar <- function(tree, model, W0 = matrix(0.0, PCMNumTraits(model), PCMNumTraits(model)),
                   metaI=PCMInfo(NULL, tree, model, verbose), internal = FALSE, verbose = FALSE)  {
  UseMethod("PCMVar", model)
}

#' Calculate the variance covariance k x k matrix at time t, under a PCM model
#' @param t positive numeric denoting time
#' @param model a PCM model object
#' @param W0 a numeric matrix denoting the initial k x k variance covariance matrix at the
#'  root (default is the k x k zero matrix).
#' @param regime an integer or a character denoting the regime in model for which to do the calculation;
#' (Defaults to 1L meaning the first regime in the model)
#' @param verbose a logical indicating if (debug) messages should be written on the console (Defaults to FALSE).
#' @return A numeric k x k matrix
#' @export
PCMVarAtTime <- function(t, model, W0 =  matrix(0.0, PCMNumTraits(model), PCMNumTraits(model)),
                         regime = 1L, verbose = FALSE) {
  if(is.character(regime)) {
    regime <- match(regime, PCMRegimes(model))
  }

  cherry <- list(tip.labels=c("t1", "t2"),
                 edge=rbind(c(3L, 1L), c(3L, 2L)),
                 edge.length=rep(t, 2),
                 Nnode = 1L,
                 edge.regime = rep(regime, 2))
  class(cherry) <- "phylo"


  metaI <- PCMInfo(NULL, cherry, model, verbose = verbose)

  # need to manually set the r member in metaI to the regime index:
  metaI$r <- rep(regime, 2)

  VarCherry <- PCMVar(cherry, model, W0, metaI = metaI, verbose = verbose)

  k <- PCMNumTraits(model)
  VarCherry[1:k, 1:k]
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
#' @param metaI a named list containg meta-information about the data and the
#' model.
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
#' "ERR:02002:PCMBase:PCM.R:PCMSim:: X0 must be of length ...".
#'}
#' @importFrom mvtnorm rmvnorm
#' @seealso \code{\link{PCMLik}} \code{\link{PCMInfo}} \code{\link{PCMCond}}
#' @export
PCMSim <- function(
  tree, model, X0,
  metaI = PCMInfo(X = NULL, tree = tree, model = model, verbose = verbose),
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
#' @param metaI a list returned from a call to \code{PCMInfo(X, tree, model)},
#'   containing meta-data such as N, M and k.
#' @param pruneI a named list containing cached preprocessing data for the tree used
#'   to perform post-order traversal (pruning). By default, this is created
#'   using \code{PCMTreePruningOrder(tree)}. This will use the default R-implementation of the
#'   likelihood calculation, which is based on the default R-implementation of the
#'   function \code{\link{PCMLmr}} (\code{\link{PCMLmr.default}}) and the S3 specification
#'   of the function
#'   \code{\link{PCMAbCdEf}} for the given model (a function called \code{PCMAbCdEf.MODEL},
#'   where model is the name of the mode and the class attribute of the \code{model}
#'   argument. For a different implementation of the function \code{\link{PCMLmr}},
#'   provide an S3 object for which an S3 specification has been written (see Details
#'   and example section).
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
  metaI = PCMInfo(X, tree, model, verbose = verbose),
  log = TRUE,
  verbose = FALSE) {

  UseMethod("PCMLik", model)
}

#' @export
logLik.PCM <- function(object, ...) {
  if(!is.PCM(object)) {
    stop("ERR:02031:PCMBase:PCM.R:logLik.PCM:: object must inherit from class PCM.")
  }

  X <- attr(object, "X", exact = TRUE)
  if( !is.matrix(X) ) {
    stop("ERR:02032:PCMBase:PCM.R:logLik.PCM:: When calling logLik.PCM on a model object, it should have a k x N numeric matrix attribute called 'X'.")
  }
  tree <- attr(object, "tree", exact = TRUE)
  if( !inherits(tree, "phylo") ) {
    stop("ERR:02033:PCMBase:PCM.R:logLik.PCM:: When calling logLik.PCM on a model object should have an attribute called 'tree' of class phylo.")
  }

  if(is.function(attr(object, "PCMInfoFun", exact = TRUE))) {
    value <- PCMLik(X, tree, object, log = TRUE, attr(object, "PCMInfoFun", exact = TRUE)(X, tree, object))
  } else {
    value <- PCMLik(X, tree, object, log = TRUE)
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
#' measurements (NA) are integrated out at the parent nodes; while non-existent traits (NaN)
#' are treated as reduced dimensionality of the vector at the parent node.
#'
#' @param X numeric k x N matrix of observed values, with possible NA entries. The
#' columns in X are in the order of tree$tip.label
#' @param tree a phylo object
#' @param metaI  The result of calling PCMInfo.
#' @return a k x M logical matrix which can be passed as a pc argument to the PCMLik
#' function. The function fails in case when all traits are NAs for some of the tips.
#' In that case an error message is issued
#' "ERR:02001:PCMBase:PCM.R:PCMPresentCoordinates:: Some tips have 0
#' present coordinates. Consider removing these tips.".
#' @seealso \code{\link{PCMLik}}
#' @export
PCMPresentCoordinates <- function(X, tree, metaI) {

  N <- metaI$N
  M <- metaI$M
  k <- metaI$k
  postorder <- rev(metaI$preorder)

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
        # here we count NAs as present because they are integrated out at the parent nodes
        pc[, i] <- !is.nan(X[, i])
      } else {
        # edges pointing to internal nodes, for which all children nodes have been
        # visited
        # here we do nothing
      }

      #update parent pc
      pc[, j] <- pc[, j] | pc[, i]
    }
    # at the end correct the present coordinates for the traits at the tips which
    # are NA but not NaN : Both NAs and NaNs result in FALSE at a tip;
    # only coordinates for which all descending tips have NaN are FALSE at internal nodes
    pc[, 1:N] <- !is.na(X[, 1:N])

    if(any(rowSums(pc) == 0)) {
      stop("ERR:02001:PCMBase:PCM.R:PCMPresentCoordinates:: Some tips
           have 0 present coordinates. Consider removing these tips.")
    }
  }
  pc
}

#' Meta-information about a tree associated with a PCM
#'
#' @description
#' @inheritParams PCMLik
#' @param preorder an integer vector of row-indices in tree$edge matrix as returned
#' by PCMTreePreorder. This can be given for performance speed-up when several
#' operations needing preorder are executed on the tree. Default : \code{NULL}.
#' @return a named list with the following elements:
#' \item{M}{total number of nodes in the tree;}
#' \item{N}{number of tips;}
#' \item{k}{number of traits;}
#' \item{RTree}{number of regimes on the tree (distinct elements of tree$edge.regime);}
#' \item{RModel}{number of regimes in the model (distinct elements of attr(model, regimes));}
#' \item{p}{number of free parameters describing the model;}
#' \item{r}{an integer vector corresponding to tree$edge with the regime for each
#' branch in tree;}
#' \item{xi}{an integer vector of 0's and 1's corresponding to the rows in tree$edge
#' indicating the presence of a jump at the corresponding branch;}
#' \item{pc}{a logical matrix of dimension k x M denoting the present coordinates
#' for each node;}
#'
#' @export
PCMInfo <- function(X, tree, model, verbose = FALSE, preorder = NULL, ...) {
  UseMethod("PCMInfo", model)
}

#' @export
PCMInfo.PCM <- function(X, tree, model, verbose = FALSE, preorder = NULL, ...) {

  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }

  if(is.null(tree$edge.regime)) {
    PCMTreeSetDefaultRegime(tree, model)
  }

  if(is.null(preorder)) {
    preorder <- PCMTreePreorder(tree)
  }

  res <- list(
    M = PCMTreeNumNodes(tree),
    N = PCMTreeNumTips(tree),
    k = PCMNumTraits(model),
    RTree = PCMTreeNumUniqueRegimes(tree),
    RModel = PCMNumRegimes(model),
    r = PCMTreeMatchRegimesWithModel(tree, model, preorder),
    p = PCMParamCount(model),
    xi = PCMTreeJumps(tree),
    edge = tree$edge,
    edge.length = tree$edge.length,
    preorder = preorder
  )
  res <- c(res, PCMOptions())

  res$pc <- PCMPresentCoordinates(X, tree, res)

  res
}

#' Create a likelhood function of a numerical vector parameter
#' @inheritParams PCMLik
#' @param positiveValueGuard positive numerical value (default Inf), which serves as a guard for numerical error. Values exceeding
#' this positiveGuard are most likely due to numerical error and PCMOptions()$PCMBase.Value.NA is returned instead.
#' @return a function of a numerical vector parameter called p returning the likelihood
#' of X given the tree and the model with parameter values specified by p.
#' @details It is possible to specify a function for the argument metaI. This function should
#' have three parameters (X, tree, model) and should return a metaInfo object. (see \code{\link{PCMInfo}}).
#'
#' @export
PCMCreateLikelihood <- function(X, tree, model, metaI = PCMInfo(X, tree, model), positiveValueGuard = Inf) {
  if(is.function(metaI)) {
    metaI <- metaI(X, tree, model)
  }
  value.NA <- PCMOptions()$PCMBase.Value.NA

  function(p, log = TRUE) {
    PCMParamLoadOrStore(model, p, offset = 0L, k = PCMNumTraits(model), R = PCMNumRegimes(model), load = TRUE)
    value <- PCMLik(X, tree, model, metaI, log = log)
    if(value > positiveValueGuard) {
      value <- value.NA
    }
    value
  }
}

