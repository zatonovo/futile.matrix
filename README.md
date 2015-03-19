Tools for studying random matrices

This package provides tools for working with random matrices.
Ensembles of random matrices have distinct eigenvalue distributions,
which can be used as a null hypothesis of sorts. Depending on the type of
matrix you have, results from random matrix theory can help you determine
whether your data is purely random or has some signal. In finance, this
information can be used to clean statistical noise from covariance and
correlation matrices.

Usage
=====
The practitioner will likely find the cutoff function most useful. This
function identifies the upper limit of the noise spectrum based on a
fit to a theoretical distribution. For example, stock market returns
are known to be non-normal. The correlation matrix of returns falls
under the class of Wishart matrices. These matrices have an eigenvalue
distribution governed by the Marcenko-Pastur distribution. However, since
stock market returns are non-normal, some eigenvalues will appear outside
of the spectrum predicted by the MP distribution. This deviation can
indicate the extent of normality and randomness within the correlations 
of a set of stock returns.
