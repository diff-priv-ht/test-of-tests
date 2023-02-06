# Implementation of the Test of Tests

This repository contains three files. The file `implementation.R` contains four main functions and four helper functions, as follows.
* `test.of.tests`: Run the test of tests for a given dataset, given public hypothesis test, and fixed input parameters $m, \alpha_0$.
* `power.test.of.tests`: Compute the theoretical power of the test of tests for a given public hypothesis test and fixed input parameters $m,\alpha_0$.
* `optimal.m.alpha0`: Find the optimal input parameters for a given public hypothesis test and *known* effect size.
* `practical.m.alpha0`: Find practical values of the input parameters for *unknown* effect size by searching for the minimum effect size required to achieve $\rho$ power for a given public hypothesis test.
* Helper functions `rtulap`, `ptulap`, `dp.binom.p.val`, and `dp.binom.alt.prob`: implementation of the binomial test described in Awan and Slavkovic's *Differentially Private Uniformly Most Powerful Tests for Binomial Data*.

The file `examples.Rmd` and corresponding compiled `examples.pdf` walk through an example of the use of these functions when the public test is a one-sample t-test.
