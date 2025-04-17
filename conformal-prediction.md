# Conformal Prediction

## Introduction

Prediction interval $C$

$$\mathbb{P}( y\in C(x)) \geq 1-\alpha$$


**Without conformal prediction**<br>
Given some training dataset $\{x_i,y_i\}_{1\leq i\leq N}$ to learn the parameter $\theta$ of the model $f_\theta$. Then given a new point $x_{\rm test}$, we have only access to the prediction $\hat{y}_{\rm test}=f_\theta(x_{\rm test}))$.

**With conformal prediction**<br>
Given some dataset $\{x_i,y_i\}_{1\leq i\leq N}$, we split into a training set $\{x_i,y_i\}_{1\leq i\leq N_{\rm train}}$, to learn the parameter $\theta$ of the model $f_\theta$, and a *calibration set* $\{x_i,y_i\}_{1\leq i\leq N_{\rm cal}}$, with $N=N_{\rm train}+N_{\rm cal}$, used to compute a conformity score, e.g., $S(x_j,y_j)=| f_\theta(x_j)-y_j|$ from which we can extract the quantile $q_{1-\alpha}$ for some small $\alpha>0$, from its distribution.
Then, given a new point $x_{\rm test}$, we have access to

$$
C(x_{\rm test}) = \{ y \;|\; S(x_{\rm test},y)\leq q_{1-\alpha} \}
$$


Note that it does assume any distribution on the data. We only assume *exchangeability*, which is a looser assumption than iid, on the calibration and test data. 


```{admonition} Example: Non-exchangeable data
:class: dropdown

Here is an example with of random variables with a trend. Consider two random variables $X_1$ and $X_2$ where:
$
X_1 \sim \mathcal{N}(0, 1)$ and $X_2 \sim \mathcal{N}(1, 1)$
and they are independent.

The joint distribution is:
$
(X_1, X_2) \sim \mathcal{N} \left( \begin{bmatrix} 0 \\ 1 \end{bmatrix}, \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix} \right)
$

If you permute them, you get $(X_2, X_1)$, whose distribution is:

$$
(X_2, X_1) \sim \mathcal{N} \left( \begin{bmatrix} 1 \\ 0 \end{bmatrix}, \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix} \right)
$$

Since the joint distributions differ, $X_1$ and $X_2$ are **non-exchangeable**.
```

## Coverage

**Statistical power vs. coverage.** The choice of the splitting between training and calibration has an important impact on the trade-off between statistical power and coverage.

**Conditional coverage.** We have seen that the prediction interval reads 
$C(x_{\rm test}) = \{ y \;|\; S(x_{\rm test},y)\leq q_{1-\alpha} \}$.
Ideally, we would like a conditional coverage, i.e.,

$$
\mathbb{P}(y\in C(x_{\rm test}) | x_{\rm test}) \geq 1-\alpha
$$
In order to have a coverage dependant of $x$. There are multiple ways to achieve this.



