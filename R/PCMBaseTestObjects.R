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

#' Test objects for the PCMBase package
#'
#' A list containing simulated trees, trait-values and model objects for tests
#' and examples of the PCMBase package
#'
#'
#' @format This is a list containing the following named elements representing
#'  parameters of BM, OU and MixedGaussian models with up to three traits and
#'  up to two regimes, model objects, simulated trees with partition of their
#'  nodes in up to two parts (corresponding to the two regimes), and trait data
#'  simulated on these trees.
#' \describe{
#'   \item{a.H, b.H }{ H matrices for OU models for regimes 'a' and 'b'. }
#'   \item{a.Theta, b.Theta }{ Theta vectors for OU models for regimes 'a' and 'b'. }
#'   \item{a.Sigma_x, b.Sigma_x }{ Sigma_x matrices for BM and OU models for regimes 'a' and 'b'. }
#'   \item{a.Sigmae_x, b.Sigmae_x }{ Sigmae_x matrices regimes 'a' and 'b'. }
#'   \item{a.X0, b.X0 }{ X0 vectors for regimes 'a' and 'b'. }
#'   \item{H }{ an array resulting from abind(a.H, b.H). }
#'   \item{Theta }{ a matrix resulting from cbind(Theta.a, Theta.b). }
#'   \item{Sigma_x }{ an array resulting from abind(a.Sigma_x, b.Sigma_x). }
#'   \item{Sigmae_x }{ an array resulting from abind(a.Sigmae_x, b.Sigmae_x). }
#'   \item{model.a.1, model.a.2, model.a.3 }{ univariate models with a single regime for each of 3 traits. }
#'   \item{model.a.1.Omitted_X0 }{ same as model.a.1 but omitting X0; suitable for nesting in an MGPM model. }
#'   \item{model.a.123, model.b.123 }{ single-regime 3-variate models. }
#'   \item{model.a.123.Omitted_X0 }{ single-regime 3-variate model with omitted X0 (suitable for nesting in an MGPM. }
#'   \item{model.a.123.Omitted_X0__bSigmae_x }{ same as model.a.123.Omitted_X0 but with the value of Sigmae_x copied from model.b.123. }
#'   \item{model.a.123.Omitted_X0__Omitted_Sigmae_x }{ same as model.a.123 but omitting X0 and Sigmae_x. }
#'   \item{model.b.123.Omitted_X0, model.b.123.Omitted_X0__Omitted_Sigmae_x }{ analogical to corresponding model.a.123... }
#'   \item{model.ab.123 }{ a two-regime 3-variate model. }
#'   \item{model.ab.123.bSigmae_x }{ a two-regime 3-variate model having Sigmae_x from b.Sigmae_x. }
#'   \item{model_MixedGaussian_ab }{ a two-regime MGPM model with a local Sigmae_x for each regime. }
#'   \item{model_MixedGaussian_ab_globalSigmae_x }{  a two-regime MGPM model with a global Sigmae_x. }
#'   \item{N }{ number of tips in simulated trees }
#'   \item{tree_15_tips }{ a tree of 15 tips used for testing clade extraction. }
#'   \item{tree.a }{ a tree with one part only (one regime) }
#'   \item{tree.ab }{ a tree partitioned in two parts (two regimes) }
#'   \item{traits.a.1 }{ trait values simulated with model.a.1. }
#'   \item{traits.a.123 }{ trait values simulated with model.a.123. }
#'   \item{traits.a.2 }{ trait values simulated with model.a.2. }
#'   \item{traits.a.3 }{ trait values simulated with model.a.3. }
#'   \item{traits.ab.123 }{ trait values simulated with model.ab.123 on tree.ab. }
#'   \item{tree }{a tree of 5 tips used for examples.}
#'   \item{X }{3-trait data for 5 tips used together with tree for examples. }
#'   \item{model.OU.BM }{a mixed Gaussian phylogenetic model for 3 traits and an OU and BM regime used in examples. }
#' }
"PCMBaseTestObjects"

#' Data for Fig3 in the TPB manuscript
#'
#' A list containing simulated tree, models and data used in Fig. 3
#'
#' @format This is a list containing the following named elements representing
#' simulation parameters, a simulated tree and PCM objects, used in Fig. 3. For
#' details on all these objects, read the file data-raw/Fig3.Rmd.
"dataFig3"
