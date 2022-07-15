# LocalPoly.jl

`LocalPoly.jl` is a Julia implementation of the local polynomial regression methods outlined in [Fan and Gijbels (1996)](https://doi.org/10.1201/9780203748725). This package is still experimental, and the API is subject to change.

## Local Polynomial Regression

Local polynomial regression is a non-parametric estimation technique that can be used to estimate both the conditional mean function and its derivatives.

Let $(X_i, Y_i)_{i=1}^N$ be observations (for sake of exposition assumed identically and independently distributed) of the random variables $(X,Y)$. Let $m(x)$ be the conditional mean function:

$$m(x) = E[Y|X=x]$$

The conditional mean function $m(x)$ can be approximated in a neighborhood of any point $x_0$ by a Taylor expansion of degree $p$:

$$m(x) \approx \sum_{j=0}^p \frac{m^{(p)}(x_0)}{j!}(x - x_0)^j$$

This suggests an estimator using the Taylor approximation with weighted data. Let $K(\cdot)$ be a valid kernel function, and $h$ the bandwidth or smoothing parameter. Denote $K_h(\cdot) = K(\cdot/h)/h$. The locally weighted sum of squared errors is:

$$\sum_{i=1}^N\left[ Y_i - \sum_{j=0}^p \beta_j \left(X_i - x_0\right)^j\right]^2 K_h(X_i - x_0)$$

Let $\widehat\beta \: (j=0,\ldots,p)$ be the $\beta$ minimizing the above expression. Then the $\nu$-th derivative of the conditional mean function evaluated at $x_0$ is:

$$\widehat m_\nu(x_0) = \nu! \widehat\beta_\nu$$

### Matrix Notation

The local polynomial estimator can be conveniently expressed using matrix notation. Define the matrices:

$$\mathbf X = \left( \begin{matrix} 1 & (X_1 - x_0) & \cdots & (X_1 - x_0)^p \\ \vdots & \vdots & & \vdots \\ 1 & (X_N - x_0) & \cdots & (X_N - x_0)^p \end{matrix} \right)$$

$$\mathbf y = \left(
    \begin{matrix}
        Y_1 \\
        \vdots \\
        Y_N
    \end{matrix}
\right)$$

$$\mathbf W = \text{diag} \left\{ K_h(X_i - x_0) \right\}$$

Then the weighted sum of squared errors is given by:

$$\left(\mathbf y - \mathbf X \beta \right)^\prime \mathbf W \left(\mathbf y - \mathbf X \beta \right)$$

The unique minimizer $\widehat\beta$ is then:

$$\widehat \beta = \left( \mathbf X^\prime \mathbf W \mathbf X \right)^{-1} \mathbf X^\prime \mathbf W \mathbf y$$
