## DGP

$$Y = m(X) + \sigma(X)\varepsilon$$

Objective function

$$\min_\beta \sum_{i=1}^N \left[ Y_i - \sum_{j=0}^p \beta_j \left(X_i - x_0\right)^j \right]^2 K_h(X_i - x_0)$$

## Matrices

$$\mathbf X = \left(
    \begin{matrix}
        1 & (X_1 - x_0) & \cdots & (X_1 - x_0)^p \\
        \vdots & \vdots & & \vdots \\
        1 & (X_N - x_0) & \cdots & (X_N - x_0)^p
    \end{matrix}
\right)$$

$$\mathbf y = \left(
    \begin{matrix}
        Y_1 \\
        \vdots \\
        Y_N
    \end{matrix}
\right)$$

$$\mathbf W = \text{diag} \left\{ K_h(X_i - x_0) \right\}$$

$$\min_\beta \left(\mathbf y - \mathbf X \beta \right)^\prime \mathbf W \left(\mathbf y - \mathbf X \beta \right)$$

$$\widehat \beta = \left( \mathbf X^\prime \mathbf W \mathbf X \right)^{-1} \mathbf X^\prime \mathbf W \mathbf y$$

## Equivalent Kernels

$$S_{n, j} = \sum_{i=1}^n K_h (X_i - x_0)(X_i-x_0)^j$$

$$S_n \equiv \mathbf X^\prime \mathbf W \mathbf X = \left( S_{n,j+l} \right)_{0 \leq j, l \leq p}$$

$$\widehat \beta_\nu = e^\prime_{\nu+1} \widehat\beta = e^\prime_{\nu+1} S_n^{-1} \mathbf X^\prime \mathbf W \mathbf y = \sum_{i=1}^n W^n_\nu \left( \frac{X_i-x_0}{h} \right) Y_i$$

$$W^n_\nu (t) = e^\prime_{\nu+1} S_n^{-1} \left( \begin{matrix}1 \\ th \\ \vdots \\ (th)^p \end{matrix}\right) \frac{K(t)}{h}$$

$$K^*_\nu (t) = e^\prime_{\nu+1} S^{-1} \left( \begin{matrix}1 \\ t \\ \vdots \\ t^p \end{matrix}\right) K(t) = \left( \sum_{l=0}^p S^{\nu l} t^l\right) K(t)$$

$$S^{-1} = \left( S^{jl} \right)_{0 \leq j, l \leq p}$$

$$S = (\mu_{j+l})_{0 \leq j, l \leq p}$$

$$\mu_j = \int \! u^j K(u) \: du$$