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

library(testthat)
context("PCMSim")

library(PCMBase)

if(PCMBaseIsADevRelease()) {


  list2env(PCMBaseTestObjects, globalenv())

  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

  # generate traits
  traits.a.1 <- PCMSim(tree.a, model.a.1, 0, verbose=TRUE)

  test_that(
    "Simulated single-trait data on a single regime tree...", {
      expect_identical(dim(traits.a.1), c(1L, PCMTreeNumNodes(tree.a)))
    }
  )

  traits.a.2 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)
  traits.a.3 <- PCMSim(tree.a, model.a.2, 0, verbose=TRUE)

  traits.a.123 <- PCMSim(tree.a, model.a.123, c(0,0,0), verbose=TRUE)

  test_that(
    "Simulated triple-trait data on a single regime tree...", {
      expect_identical(dim(traits.a.123), c(3L, PCMTreeNumNodes(tree.a)))
    }
  )

  traits.ab.123 <- PCMSim(tree.ab, model.ab.123, c(0,0,0), verbose=TRUE)
  test_that(
    "Simulated triple-trait data on a two-regime tree...", {
      expect_identical(dim(traits.ab.123), c(3L, PCMTreeNumNodes(tree.ab)))
    }
  )

}
