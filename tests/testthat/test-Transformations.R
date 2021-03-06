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

library(PCMBase)
library(testthat)

context("Test parameter transformaitons")

if(PCMBaseIsADevRelease()) {

  # number of traits
  k <- 3

  listParameterizationsOU <- list(
    X0 = list(c("VectorParameter", "_Global")),
    H = list(c("MatrixParameter", "_Zeros"),
             c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
             c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable")),
    Theta = list(c("VectorParameter", "_Zeros")),
    Sigma_x = list(c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
                   c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
                   c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal")),

    Sigmae_x = list(c("MatrixParameter", "_Zeros"))
  )

  PCMGenerateParameterizations(structure(0.0, class="OU"), listParameterizations = listParameterizationsOU)

  PCMModels("^OU")
  model <- PCM("OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Zeros_Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Zeros_Sigmae_x",
               modelTypes = PCMModels("^OU"),
               k = k,
               regimes = letters[1:3])

  vecModelRandom <- round(PCMParamRandomVecParams(model), 1)
  modelRandom <- model
  PCMParamLoadOrStore(modelRandom, vecModelRandom, offset = 0, load=TRUE)
  model <- modelRandom
  test_that(
    "Before applying a transformation", {
      expect_true(is.Transformable(model))
    })

  model2 <- PCMApplyTransformation(model)

  test_that(
    "After applying a transformation", {
      expect_false(is.Transformable(model2))
    })

}
