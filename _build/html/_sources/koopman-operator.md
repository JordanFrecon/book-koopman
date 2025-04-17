# Deterministic Dynamics and Koopman Operators

In the previous chapter, we developed a functional-analytic perspective on stochastic processes using **Markov operators**, which evolve probability distributions and observables under random dynamics. We saw how this linear operator framework allowed us to study convergence, invariant measures, and long-term statistical behavior in a unified way.

In this chapter, we turn to **deterministic dynamical systems**, where the evolution of the system is governed by a fixed rule without randomness. Despite this fundamental shift, the operator-theoretic viewpoint remains powerful—and in fact, becomes even more striking. The central object of study is the **Koopman operator**: a linear operator that describes the evolution of observables by composition with the system’s dynamics.

Unlike the Markov operator, which averages over future possibilities, the Koopman operator captures the precise evolution of deterministic systems. This allows us to analyze nonlinear dynamics using linear tools, such as spectral decomposition, eigenfunctions, and mode analysis. Through this lens, we will explore the long-term behavior, invariant sets, and structure of complex dynamical systems.

This chapter builds a bridge between probability, functional analysis, and dynamical systems—laying the foundation for data-driven methods such as dynamic mode decomposition and Koopman learning, discussed in the next chapter.



## From Markov Operators to Koopman Operators

We now shift from stochastic to deterministic dynamics. In a deterministic system, the state evolves according to a rule without randomness. That is, the next state is entirely determined by the current one.

### Deterministic Dynamical Systems

Let $(\mathcal{X}, \Sigma_\mathcal{X})$ be a measurable space, and let $T : \mathcal{X} \to \mathcal{X}$ be a measurable, deterministic map. A **discrete-time deterministic dynamical system** describes the evolution of a system through a sequence of states determined by a rule. Such a system may take one of two forms:

- **Non-autonomous**, where the update rule varies with time:

  $$
  x_{t+1} = T_t(x_t),
  $$

  with each $T_t : \mathcal{X} \to \mathcal{X}$ depending explicitly on the time index $t$.

- **Autonomous**, where the same rule applies at every time step:

  $$
  x_{t+1} = T(x_t),
  $$

  for a fixed, time-invariant map $T$.

In this chapter, we focus exclusively on **autonomous deterministic systems**, which form the classical foundation of Koopman operator theory. Given an initial state $x_0 \in \mathcal{X}$, the system generates a deterministic trajectory $\{x_t\}_{t \geq 0}$ by iterating the map:

$$
x_t = T^t(x_0),
$$

where $T^t$ denotes $t$-fold composition of $T$. There is no randomness involved: the system follows a unique, fully determined path from any starting point.



### Primal Perspective

In the Markovian setting, we introduced the **Markov transfer operator** $\mathcal{A}_F$, which evolves observables through expectations under a transition kernel:

$$
\mathcal{A}_F f(x) = \int_{\mathcal{X}} f(y)\, p(x, dy).
$$

For deterministic systems, this simplifies: there is no averaging over possible outcomes. The evolution of an observable is obtained by composing it with the deterministic map $T$:

```{admonition} Definition: *Koopman Operator*  
:class: note

Let $T : \mathcal{X} \to \mathcal{X}$ be an autonomous deterministic map. The **Koopman operator** $\mathcal{K}$ acts on observables $f : \mathcal{X} \to \mathbb{R}$ by

$$
\mathcal{K} f(x) = f(T(x)).
$$
```

Despite the fact that $T$ may be nonlinear, the Koopman operator $\mathcal{K}$ is **linear**:

$$
\mathcal{K}(a f + b g) = a \mathcal{K} f + b \mathcal{K} g.
$$

This linearity is powerful: it allows us to study nonlinear dynamical systems using the tools of linear operator theory, including spectral analysis, eigenfunctions, and mode decompositions.



### Dual Perspective

In the stochastic setting, we studied two complementary operators:
- The **Markov transfer operator** $\mathcal{A}_F$, which acts on observables,
- Its **dual** $\mathcal{A}^*$, which acts on probability measures.

These were linked by the duality relation:

$$
\int_{\mathcal{X}} (\mathcal{A}_F f)(x)\, d\mu(x) = \int_{\mathcal{X}} f(x)\, d(\mathcal{A}^* \mu)(x).
$$

In the deterministic setting, the same dual structure arises. The Koopman operator $\mathcal{K}$ acts on observables by composition with the map $T$, while its dual is the **pushforward** operator:

$$
T_* \mu(B) := \mu(T^{-1}(B)), \quad \forall B \in \Sigma_{\mathcal{X}}.
$$

This operator evolves probability measures forward in time according to the dynamics. It is also known as the **Perron–Frobenius operator**, often denoted $\mathcal{P}$.

The Koopman and Perron–Frobenius operators satisfy the duality:

$$
\int_{\mathcal{X}} (\mathcal{K} f)(x)\, d\mu(x) = \int_{\mathcal{X}} f(x)\, d(\mathcal{P} \mu)(x).
$$

So just as $(\mathcal{A}_F, \mathcal{A}^*)$ describe the evolution of functions and distributions in the Markovian setting, the pair $(\mathcal{K}, \mathcal{P})$ plays the same role in deterministic dynamics.

This structural similarity highlights one of the key themes of this course: **both stochastic and deterministic systems admit linear operator frameworks that reveal their long-term behavior**—through fixed points, spectral properties, and invariant subspaces.


### Invariant Measures and Operator Duality

As in the Markov setting, **invariant measures** play a fundamental role in the analysis of deterministic systems through the Koopman and Perron–Frobenius operators. They serve as a bridge between the primal and dual perspectives, enabling both the definition of the Koopman operator on $L^2$ spaces and the interpretation of long-term statistical behavior.

```{admonition} Definition: *Invariant Measure for a Deterministic Map*
:class: note

Let $T : \mathcal{X} \to \mathcal{X}$ be a measurable map. A probability measure $\mu$ on $(\mathcal{X}, \Sigma_\mathcal{X})$ is **invariant** under $T$ if

$$
\mu(T^{-1}(B)) = \mu(B), \quad \forall B \in \Sigma_{\mathcal{X}}.
$$

Equivalently, the pushforward measure satisfies $T_* \mu = \mu$.
```

This means that if $x_0 \sim \mu$, then $x_t \sim \mu$ for all $t \in \mathbb{N}$. Invariant measures ensure that distributions remain stable under the dynamics and provide a natural reference measure for studying long-term averages and ergodic properties.

In the deterministic setting, invariant measures allow us to:
- Define the Koopman operator on $L^2_\mu(\mathcal{X})$ and interpret it as an isometry,
- Give statistical meaning to time-averaged observables,
- Formulate duality between $\mathcal{K}$ and the Perron–Frobenius operator $\mathcal{P} = T_*$.

> These roles mirror those of invariant measures in the Markovian setting, where the condition $\mathcal{A}^* \pi = \pi$ played the same stabilizing role.

At this stage, we will **not yet assume** the existence of an invariant measure, nor prove when one exists. That discussion will come later, once we introduce the functional spaces where $\mathcal{K}$ acts and explore ergodic and asymptotic properties of deterministic systems.

### Summary

To conclude this section, we summarize the key parallels between stochastic and deterministic dynamics through their associated operators.

These two settings share a common structure: a pair of dual linear operators describing the evolution of observables and distributions. However, they differ in how randomness enters the system and how observables evolve.

| Concept                          | Markov Processes (Stochastic)                          | Dynamical Systems (Deterministic)                        |
|----------------------------------|--------------------------------------------------------|----------------------------------------------------------|
| Evolution Rule                  | Transition kernel $p(x, dy)$                           | Map $T : \mathcal{X} \to \mathcal{X}$                    |
| Observable Operator             | Markov operator $\mathcal{A}_F f(x) = \int f(y)\, p(x, dy)$ | Koopman operator $\mathcal{K} f(x) = f(T(x))$           |
| Dual Operator (Measures)        | $\mathcal{A}^* \mu(B) = \int p(x, B)\, d\mu(x)$         | Pushforward $T_* \mu(B) = \mu(T^{-1}(B))$               |
| Dual Operator Name              | —                                                      | Perron–Frobenius operator $\mathcal{P} = T_*$           |
| Randomness                      | Probabilistic transitions                              | None (purely deterministic)                             |
| Linearity of Dynamics           | Generally nonlinear                                    | Can be nonlinear                                        |
| Linearity of Operator           | Linear                                                 | Linear                                                   |
| Common Goal                     | Study long-term distributions & statistics             | Study spectral structure & invariant sets               |

This parallel structure motivates the use of operator-theoretic tools for analyzing a wide variety of systems—stochastic or deterministic—through a unified mathematical lens.


## Functional Setting and Invariant Measures

In the previous section, we defined the Koopman operator $\mathcal{K}$ as a linear operator that evolves observables through composition with a deterministic map $T$:

$$
\mathcal{K}f(x) = f(T(x)).
$$

However, we did not specify the space to which the observable $f : \mathcal{X} \to \mathbb{C}$ belongs. This choice is not merely technical—it fundamentally shapes the **boundedness**, **continuity**, and **spectral behavior** of $\mathcal{K}$.

The function space determines not only whether $\mathcal{K}$ is well-defined and bounded, but also which dynamical features it can reveal: resonances, coherent structures, long-term averages, or transient modes.

In fact, the same dynamical system can exhibit very different spectral signatures depending on the function space in which $\mathcal{K}$ is analyzed.

In this section, we explore these functional settings in detail, and discuss how they influence the Koopman operator's analytical properties and practical usefulness.

### Common Function Spaces for Observables

Before analyzing the Koopman operator, it's important to clarify the typical **function spaces** used to represent observables in dynamical systems. Each space imposes different regularity, integrability, or boundedness constraints on functions, and these choices directly influence how the Koopman operator behaves—especially in terms of its boundedness, spectral properties, and analytic structure.

The table below summarizes the most commonly used spaces, their notation, type (Hilbert or Banach), and typical use cases.

| Function Space              | Notation                  | Type     | Description                                     | Common Use Cases                              |
|----------------------------|---------------------------|----------|-------------------------------------------------|-----------------------------------------------|
| Square-integrable          | $L^2_\mu(\mathcal{X})$     | Hilbert  | Functions with finite variance under $\mu$      | Spectral theory, ergodicity, inner products   |
| Bounded measurable         | $L^\infty(\mathcal{X})$    | Banach   | All functions with bounded essential sup norm   | General dynamics, growth/decay detection      |
| $p$-integrable             | $L^p_\mu(\mathcal{X})$     | Banach   | For $1 \leq p < \infty$, allows unbounded observables | Ergodic theorems, averages, indicator functions |
| Continuous (compact domain)| $C^0(\mathcal{X})$         | Banach   | Continuous functions with sup norm              | Topological dynamics, symbolic systems         |
| Hardy           | $H^2$                      | Hilbert  | Analytic functions with square-summable coeffs  | Resonance, control theory, analytic maps       |
| Sobolev       | $W^{k,p}$ or $H^s$         | Banach/Hilbert | Functions with $k$ derivatives in $L^p$       | Control, PDE systems, stability analysis       |



The Koopman operator is always linear, but its properties—boundedness, compactness, spectral radius, norm—depend **entirely** on the function space it acts on. **The more structure** you require of your observables (e.g., smoothness, analyticity), the **more information and spectral regularity** you gain,  but at the cost of **less generality**.


We will focus primarily on the $L^p_\mu$ family, but briefly mention others when relevant.



### Invariant Measures and Motivation

To study the Koopman operator on spaces like $L^2_\mu(\mathcal{X})$, we need to choose a measure $\mu$. Not all choices lead to well-behaved operators.

```{admonition} Definition: *Invariant Measure*
:class: note

Let $T : \mathcal{X} \to \mathcal{X}$ be a measurable map. A probability measure $\mu$ is **invariant** under $T$ if

$$
\mu(T^{-1}(B)) = \mu(B), \quad \forall B \in \Sigma_{\mathcal{X}},
$$

or equivalently, $T_* \mu = \mu$.
```

This condition ensures that the measure is preserved under the dynamics—i.e., if $X_0 \sim \mu$, then $X_t \sim \mu$ for all $t \in \mathbb{N}$.

```{admonition} Remark: *Why invariant measures matter*
:class: tip
We will see that, if $\mu$ is an invariant measure, then
- The Koopman operator **preserves the $L^p_\mu$ norm**.
- In particular, it becomes an **isometry** on $L^2_\mu$, meaning that $\|\mathcal{K} f\|_{L^2_\mu} = \|f\|_{L^2_\mu}$.
```

At this point, a natural question arises: **do invariant measures always exist**?

In general, the answer depends on both the structure of the map $T$ and the topological properties of the space $\mathcal{X}$. Fortunately, under very mild assumptions, we can guarantee the existence of at least one invariant measure.


This mirrors the situation encountered in the **Markovian setting**, where we used the **Krylov–Bogoliubov theorem** to construct invariant measures under mild compactness and continuity conditions. The key idea—averaging the empirical occupation of the system over time—can be adapted almost directly to the deterministic case.

```{admonition} Theorem: *Krylov–Bogoliubov (deterministic version)*
:class: important

Let $\mathcal{X}$ be a compact metric space, and let $T : \mathcal{X} \to \mathcal{X}$ be continuous. Then there exists at least one Borel probability measure $\mu$ such that $T_* \mu = \mu$.
```

```{admonition} Proof Sketch
:class: dropdown

Start from any initial point $x_0 \in \mathcal{X}$ and define the sequence of empirical measures:

$$
\mu_n := \frac{1}{n} \sum_{t=0}^{n-1} \delta_{T^t(x_0)}.
$$

Since $\mathcal{X}$ is compact, the space of probability measures on $\mathcal{X}$ (with the weak-* topology) is also compact by **Prokhorov's theorem**. The sequence $(\mu_n)$ is tight and admits at least one weak-* accumulation point. It can be shown that any such limit $\mu$ satisfies $T_* \mu = \mu$.
```

This construction is conceptually identical to what we did for Markov chains, where the transition kernel was used to evolve a distribution forward. In the deterministic case, the transition kernel is a Dirac measure: $p(x, \cdot) = \delta_{T(x)}$, making the deterministic case a special instance of the Markov framework.



We now present a few important examples of invariant measures in deterministic systems. These illustrate different behaviors: fixed points, periodic motion, and volume-preserving flows.

```{admonition} Examples of Invariant Measures
:class: note

Here are several classical deterministic systems that admit invariant measures:


**Fixed Point**: If $T(x^*) = x^*$, then the Dirac measure $\mu = \delta_{x^*}$ is invariant.


**Rotation on the Circle**: Let $T(x) = x + \alpha \mod 1$ on $\mathbb{T} = \mathbb{R}/\mathbb{Z}$. Then the Lebesgue probability measure on $[0,1)$ is invariant under $T$.


**Volume-Preserving Map on $\mathbb{R}^n$**: If $T : \mathbb{R}^n \to \mathbb{R}^n$ is a smooth bijection with $\det(DT(x)) = 1$ for all $x$, then the (possibly infinite) Lebesgue measure on $\mathbb{R}^n$ is invariant. When $T$ preserves a compact set $A \subset \mathbb{R}^n$, one often uses the normalized Lebesgue measure on $A$ to define a probability invariant measure.


**Expanding Map (Doubling)**: Let $T(x) = 2x \mod 1$ on $[0,1]$. The Lebesgue probability measure on $[0,1]$ is invariant under $T$.

```


```{admonition} Proofs
:class: dropdown

**Fixed Point**  
Let $T(x^*) = x^*$. Then for any measurable set $B$:

$$
T_* \delta_{x^*}(B) = \delta_{x^*}(T^{-1}(B)) = 
\begin{cases}
1 & \text{if } x^* \in T^{-1}(B) \\
0 & \text{otherwise}.
\end{cases}
$$
Since $T(x^*) = x^*$, we have $x^* \in T^{-1}(B)$ if and only if $x^* \in B$, so $T_* \delta_{x^*} = \delta_{x^*}$.

---

**Rotation on the Circle**  
Let $T(x) = x + \alpha \mod 1$. Then for any measurable set $B \subset \mathbb{T}$:

$$
T^{-1}(B) = B - \alpha.
$$
Since the Lebesgue measure $\mu$ on $\mathbb{T}$ is translation-invariant:

$$
\mu(T^{-1}(B)) = \mu(B - \alpha) = \mu(B).
$$
Hence, $T_* \mu = \mu$.

---

**Volume-Preserving Map**  
Let $T : \mathbb{R}^n \to \mathbb{R}^n$ be a smooth bijection with $\det(DT(x)) = 1$. Let $\mu$ denote the Lebesgue measure on $\mathbb{R}^n$. Then by the change of variables formula:

$$
\mu(T^{-1}(B)) = \int_{T^{-1}(B)} dx = \int_B |\det(DT^{-1}(y))| \, dy = \int_B dy = \mu(B),
$$
so $T_* \mu = \mu$.

If $T$ preserves a compact measurable subset $A \subset \mathbb{R}^n$ (i.e., $T(A) = A$), one may define a **probability measure** by normalizing the restricted Lebesgue measure:

$$
\mu_A(B) := \frac{\text{Leb}(B \cap A)}{\text{Leb}(A)}.
$$
This $\mu_A$ is then $T$-invariant on $A$.
---

**Expanding Map: Doubling on [0,1]**  
Define $T(x) = 2x \mod 1$ on $[0,1]$. Then for any measurable set $B \subset [0,1]$:

$$
T^{-1}(B) = \left\{ \frac{x}{2} \in [0,1/2] : x \in B \right\} \cup \left\{ \frac{x+1}{2} \in [1/2,1] : x \in B \right\}.
$$

That is,

$$
T^{-1}(B) = \left\{ \frac{x}{2} : x \in B \right\} \cup \left\{ \frac{x+1}{2} : x \in B \right\}.
$$

Because Lebesgue measure scales linearly under affine transformations, we get:

$$
\mu(T^{-1}(B)) = \frac{1}{2} \mu(B) + \frac{1}{2} \mu(B) = \mu(B).
$$

Thus, $T_* \mu = \mu$ and the Lebesgue probability measure on $[0,1]$ is invariant.
```


## Boundedness of the Koopman Operator

### Why Boundedness Matters

Boundedness is a central concept in functional analysis. When the Koopman operator $\mathcal{K}$ is **bounded** on a given function space (like $L^2_\mu$, $L^p$, or $C^0$), it satisfies the estimate

$$
\|\mathcal{K} f\| \leq C \|f\|
$$

for some constant $C > 0$ and all $f$ in the space. This property has powerful consequences:

- $\mathcal{K}$ is **continuous**, so it interacts well with limits and approximations;
- It can be iterated stably (i.e., $\mathcal{K}^n f$ is meaningful and uniformly controlled);
- We can apply **spectral theory** and **functional calculus**, which rely on boundedness to define the spectrum, compute projections, and analyze long-term dynamics;
- Numerically, a bounded operator allows stable computations and reliable approximations.

Boundedness is thus a prerequisite for developing most of the analytic and numerical tools used in Koopman theory.


```{admonition} Why Unbounded Koopman Operators Are Problematic
:class: warning

In infinite-dimensional spaces, an operator that is **unbounded** cannot be defined on the whole space. This is not just a technicality—it is a fundamental consequence of the **Closed Graph Theorem**:

> If a linear operator on a Banach space is defined on all of the space and is closed, then it must be bounded.

This implies the contrapositive:

> If an operator is **unbounded**, then it **cannot** be defined everywhere.

In the Koopman setting, this means:

- If $T$ is highly nonlinear or $f$ is unbounded, the composition $f \circ T$ may not lie in the same function space.
- As a result, $\mathcal{K} f$ may not be defined for all $f$ in the space, and the operator must be treated as **partial**, with domain $\mathcal{D}(\mathcal{K}) \subsetneq \mathcal{F}$.
- Such operators require special care: their **domains must be specified**, and standard spectral tools apply only if the operator is **densely defined** and **closable**.

For this reason, Koopman theory is usually developed in settings where $\mathcal{K}$ is bounded—either by assuming an invariant measure (as in $L^2_\mu$) or by choosing function spaces (like $L^\infty$ or $C^0$) where boundedness is automatic.

```

In the following subsections, we will explore in detail **when and where** the Koopman operator is bounded—and how this depends on the interplay between the dynamics $T$, the function space, and the presence of an invariant measure.


### Boundedness on $ L^p_\mu(\mathcal{X}) $

Let $(\mathcal{X}, \Sigma_{\mathcal{X}}, \mu)$ be a measure space, and suppose $ T : \mathcal{X} \to \mathcal{X} $ is measurable. We now investigate under what conditions the Koopman operator

$$
\mathcal{K} f(x) := f(T(x))
$$

is bounded on $ L^p_\mu(\mathcal{X}) $ for $ 1 \leq p < \infty $.


#### When is the operator bounded?

The key condition is that the measure $ \mu $ is **invariant** under $ T $. When this holds, the Koopman operator is not only bounded — it is an **isometry**.

```{admonition} Proposition: *Isometry on $ L^p_\mu $ when $ \mu $ is invariant*
:class: important

Let $ \mu $ be a $ T $-invariant measure. Then the Koopman operator $ \mathcal{K} $ is an **isometry** on $ L^p_\mu(\mathcal{X}) $ for all $ 1 \leq p < \infty $. That is,

$$
\|\mathcal{K} f\|_{L^p_\mu} = \|f\|_{L^p_\mu}, \quad \forall f \in L^p_\mu.
$$
```

```{admonition} Proof
:class: dropdown

We compute:

$$
\begin{align*}
\|\mathcal{K} f\|_{L^p_\mu}^p 
&= \int_{\mathcal{X}} |f(T(x))|^p \, d\mu(x) \\
&= \int_{\mathcal{X}} |f(y)|^p \, d\mu(T^{-1}(y)) \quad \text{(change of variable)} \\
&= \int_{\mathcal{X}} |f(y)|^p \, d\mu(y) \quad \text{(invariance)} \\
&= \|f\|_{L^p_\mu}^p.
\end{align*}
$$
```

#### What if $ \mu $ is not invariant?

If $ \mu $ is **not** invariant, then in general:

- $ \mathcal{K} f $ may not lie in $ L^p_\mu $, even if $ f \in L^p_\mu $;
- Or it may be defined but **not bounded**;
- Even worse, $ \mathcal{K} $ may be **unbounded**, and hence not defined on all of $ L^p_\mu $.

For example, if $ T $ expands volumes and $ f $ grows quickly, then $ f \circ T $ may lie outside $ L^p_\mu $.

This is why invariant measures play such a fundamental role: they guarantee the Koopman operator is well-behaved in $ L^p $ spaces and enable a stable spectral theory.



### Boundedness on $ L^\infty(\mathcal{X}) $

In contrast to the case of $ L^p_\mu $, where boundedness of the Koopman operator may depend on the invariance of the measure $ \mu $, the space $ L^\infty(\mathcal{X}) $ always provides a safe environment: the Koopman operator is **automatically bounded**.

```{admonition} Proposition: *Contractivity on $ L^\infty $*
:class: important

Let $ T : \mathcal{X} \to \mathcal{X} $ be a measurable map. Then the Koopman operator $ \mathcal{K} f = f \circ T $ is bounded on $ L^\infty(\mathcal{X}) $, with

$$
\|\mathcal{K} f\|_{L^\infty} \leq \|f\|_{L^\infty}.
$$

Hence, $ \mathcal{K} $ is a **contraction**: its operator norm satisfies $ \|\mathcal{K}\| \leq 1 $.
```

```{admonition} Proof
:class: dropdown

We compute:

\begin{align*}
\|\mathcal{K} f\|_{L^\infty} &= \operatorname*{ess\,sup}_{x \in \mathcal{X}} |f(T(x))| \\
&\leq \operatorname*{ess\,sup}_{y \in \mathcal{X}} |f(y)| = \|f\|_{L^\infty}.
\end{align*}

This uses the fact that the essential supremum is invariant under composition with a measurable map.
```


#### Why this matters

- **No measure required**: Unlike in $ L^p $ spaces, boundedness here does **not** rely on an invariant measure.
- **Safe default**: This makes $ L^\infty $ a **safe and general choice**, especially when the system is non-conservative, dissipative, or the invariant measure is unknown.
- **Natural space**: Observables in $ L^\infty $ include indicators, step functions, and all bounded real-valued measurements of the state.

However, because $ L^\infty $ is not a Hilbert space, we **lose the inner product structure** that underpins much of spectral theory (e.g., orthogonal decompositions). We gain robustness, but at the expense of spectral richness.


### Boundedness on $ C^0(\mathcal{X}) $

Another natural setting for the Koopman operator is the space of continuous functions, denoted $ C^0(\mathcal{X}) $, where $ \mathcal{X} $ is a compact metric space. This context is common in **topological dynamics**, where continuity rather than integrability is the main concern.

In this setting, the Koopman operator is again **always bounded**, provided that the map $ T $ itself is continuous.

```{admonition} Proposition: *Boundedness on $ C^0(\mathcal{X}) $*
:class: important

Let $ \mathcal{X} $ be a compact metric space, and let $ T : \mathcal{X} \to \mathcal{X} $ be continuous. Then the Koopman operator $ \mathcal{K} f = f \circ T $ maps $ C^0(\mathcal{X}) $ to itself and satisfies

$$
\|\mathcal{K} f\|_\infty \leq \|f\|_\infty,
$$

for all $ f \in C^0(\mathcal{X}) $. In particular, $ \mathcal{K} $ is a contraction on $ C^0(\mathcal{X}) $.
```

```{admonition} Proof
:class: dropdown

Since $ T $ is continuous and $ f \in C^0(\mathcal{X}) $, the composition $ f \circ T $ is also continuous. Moreover,

$$
\|\mathcal{K} f\|_\infty = \sup_{x \in \mathcal{X}} |f(T(x))| \leq \sup_{y \in \mathcal{X}} |f(y)| = \|f\|_\infty.
$$
```


#### Comments

- This space does **not require a measure** at all.
- The Koopman operator acts naturally on observables like distance functions, potentials, or continuous physical measurements.
- The compactness of $ \mathcal{X} $ ensures that the supremum norm is finite and well-defined.
- Just as in $ L^\infty $, we again lack a Hilbert structure, so many spectral tools from $ L^2 $ do not carry over directly.



## Spectral Properties of Koopman Operators


One of the most striking features of the Koopman operator is that even though the underlying system is often nonlinear, its action on observables is **linear**. This means we can meaningfully speak about spectral properties such as **eigenvalues** and **eigenfunctions**, as in linear algebra or Fourier analysis. The **structure of the spectrum** reveal whether the system has:

- Fixed points, periodic or quasiperiodic motion;
- Stable or unstable modes;
- Decaying transients or persistent oscillations.

But a critical subtlety is that the **spectrum of the Koopman operator is not intrinsic to the system** — it depends heavily on the **function space** where the operator acts. Whether we're working in $L^2$, $L^\infty$, or $C^0$, the spectrum can look entirely different.

We begin by exploring how this choice of functional space impacts the spectral theory.


### Functional Spaces and Their Spectral Implications

Different function spaces grant different analytical tools. In particular, the **presence or absence of a scalar product** (inner product) determines which kinds of spectral analysis can be meaningfully applied.

#### $L^2_\mu(\mathcal{X})$: The ideal setting

When $\mu$ is an invariant measure and we work in $L^2_\mu$, the Koopman operator is:

- **Linear**, **bounded**, and **isometric**;
- A **Hilbert space** structure is available (inner product, projections, orthogonality);
- The **spectral theorem** applies for unitary operators;
- Eigenfunctions can be used to build **orthogonal decompositions**;
- Spectral radius is bounded by $1$.

This makes $L^2_\mu$ the canonical setting for theoretical and computational Koopman analysis.

#### $L^\infty(\mathcal{X})$: Bounded, but no inner product

Koopman is always bounded on $L^\infty$, even without an invariant measure. However:

- There is **no scalar product**;
- We cannot define **adjoints**, **self-adjointness**, or **orthogonal decompositions**;
- Spectral theory exists, but is more limited (point spectrum only, no projections);
- Useful in general settings where no measure is known.

#### $C^0(\mathcal{X})$: Continuous dynamics

When $T$ is continuous and $\mathcal{X}$ is compact, Koopman acts boundedly on the space of continuous functions:

- Useful in **topological dynamics**;
- Still no scalar product, but strong continuity structure;
- Spectrum reflects topological, not measure-theoretic, behavior.

#### Sobolev, Hardy, and analytic spaces

More specialized spaces give additional structure:
- **Hardy spaces** (like $H^2$) are analytic and yield **discrete spectrum**;
- **Sobolev spaces** add smoothness and enable **PDE connections**;
- These spaces often lead to **compact Koopman operators** and **eigenvalue decay**.


#### Summary: Spectral Behavior by Function Space

Below is a comparison of the Koopman operator's behavior across different functional settings. While spectral theory can be applied in all of them, the tools available and the structure of the spectrum depend critically on whether the space has a scalar product, how smooth the observables are, and whether an invariant measure exists.

| Space                    | Spectral Theory   | Spectral Tools Available           | Pros                                         | Cons                                           |
|--------------------------|-------------------|------------------------------------|----------------------------------------------|------------------------------------------------|
| $ L^2_\mu(\mathcal{X}) $ | Full | Spectral theorem, adjoint, orthogonality | Inner product, isometry, modal decompositions | Requires invariant measure                    |
| $ L^\infty(\mathcal{X}) $ | Weaker         | Point spectrum, spectral radius     | Bounded always, no measure needed            | No scalar product or projections              |
| $ C^0(\mathcal{X}) $     | Weaker         | Continuous spectrum, topological tools | Suitable for topological dynamics            | Spectrum may be more abstract                 |
| Sobolev / Hardy spaces     | Specialized     | Discrete spectrum, functional calculus | Captures smoothness or analyticity           | May not be preserved under dynamics           |





### Koopman Eigenfunctions and the Point Spectrum

Eigenfunctions of the Koopman operator encode **structured, predictable modes** of evolution. Each eigenfunction evolves **multiplicatively** over time, allowing us to interpret it as a decoupled dynamical mode.

#### Koopman Eigenvalue Equation

Let $ \mathcal{K} $ be the Koopman operator acting on a space of observables $ \mathcal{F} $, such as $ L^2_\mu(\mathcal{X}) $ or $ C^0(\mathcal{X}) $. A function $ \varphi \in \mathcal{F} $ is a **Koopman eigenfunction** with eigenvalue $ \lambda \in \mathbb{C} $ if

$$
\mathcal{K} \varphi = \lambda \varphi.
$$

Since $ \mathcal{K} \varphi(x) = \varphi(T(x)) $, the eigenvalue equation becomes:

$$
\varphi(T(x)) = \lambda \varphi(x), \quad \forall x \in \mathcal{X}.
$$

This means that as the system evolves under $ T $, the value of the eigenfunction evolves geometrically:

$$
\varphi(T^t(x)) = \lambda^t \varphi(x), \quad \text{for all } t \in \mathbb{N}.
$$

Thus, Koopman eigenfunctions evolve **exponentially** (in modulus) or **periodically** (in phase) over time.


```{admonition} Remark: *Complex Koopman Eigenvalues*
:class: tip

Koopman eigenvalues are generally complex. This allows the operator to capture **oscillatory behaviors**, just as complex exponentials appear in Fourier analysis. If $ \lambda = e^{i\omega} $, then

$$
\varphi(x_t) = \lambda^t \varphi(x_0) = e^{i\omega t} \varphi(x_0)
$$

describes **pure periodic motion** at frequency $ \omega $.
```

The **modulus** of the eigenvalue $ \lambda $ reveals the **stability or amplification** of the corresponding mode:

- $ |\lambda| = 1 $: **Neutral mode**  
  The observable persists over time without growth or decay (e.g., rotations, conserved quantities).
- $ |\lambda| < 1 $: **Decaying mode**  
  The observable decays exponentially; often associated with **attractors** or **transients**.
- $ |\lambda| > 1 $: **Growing mode**  
  The observable grows exponentially — this only occurs if $ \mathcal{K} $ is **unbounded**.

#### Why Eigenfunctions Matter

Koopman eigenfunctions provide a **coordinate system** in which the nonlinear evolution becomes diagonal:

$$
\varphi(x_t) = \lambda^t \varphi(x_0).
$$

This allows:

- **Modal decomposition** of observables;
- Extraction of **frequencies** and **growth rates** in time series;
- Reduction of complex dynamics to simple exponential behaviors in the right basis.

In data-driven contexts (e.g., DMD), the goal is often to **approximate or recover** Koopman eigenfunctions and eigenvalues from observations.

### Beyond Eigenvalues: The Full Koopman Spectrum

While Koopman eigenfunctions provide valuable insight into the system’s modal structure, they are not always sufficient to describe long-term dynamics. Many systems, particularly those exhibiting chaos or mixing, require a more general spectral analysis.

The Koopman operator is linear and bounded (on appropriate spaces), so its **spectrum** can be studied using tools from spectral theory — even when no eigenfunctions exist.


```{admonition} Definition: *Koopman Spectrum*
:class: note

Let $ \mathcal{K} $ be a bounded linear operator on a Banach space $ \mathcal{F} $. The **spectrum** $ \sigma(\mathcal{K}) $ is the set of complex numbers $ \lambda \in \mathbb{C} $ for which the operator $ (\mathcal{K} - \lambda I) $ is not invertible.
```

The spectrum splits into three categories:

- **Point spectrum**: values of $ \lambda $ for which $ \mathcal{K} \varphi = \lambda \varphi $ (eigenfunctions exist);
- **Continuous spectrum**: $ \mathcal{K} - \lambda I $ is injective and has dense range, but no bounded inverse;
- **Residual spectrum**: $ \mathcal{K} - \lambda I $ is injective, but the range is not dense.

This trichotomy reflects the only ways in which invertibility can fail in infinite-dimensional spaces. Every $ \lambda \in \sigma(\mathcal{K}) $ must fall into **one and only one** of these categories.

```{admonition} More detail: How each part fits the spectrum definition
:class: dropdown

Recall that a complex number $ \lambda \in \mathbb{C} $ belongs to the **spectrum** $ \sigma(\mathcal{K}) $ if the operator $ \mathcal{K} - \lambda I $ is **not invertible** as a bounded operator from $ \mathcal{F} \to \mathcal{F} $. This can occur in exactly one of the following ways:

| Condition                                              | Consequence                                                         |
|--------------------------------------------------------|----------------------------------------------------------------------|
| $ \mathcal{K} - \lambda I $ is **not injective**      | There exists a nonzero $ \varphi $ such that $ \mathcal{K} \varphi = \lambda \varphi $ |
| $ \mathcal{K} - \lambda I $ is injective, has **dense range**, but **no bounded inverse** | Inverse exists algebraically, but solutions are unstable or diverge |
| $ \mathcal{K} - \lambda I $ is injective, but range is **not dense** | Some $ g \in \mathcal{F} $ have no approximating preimage          |


We now examine each spectral component and prove why invertibility fails.

---

- **Point spectrum** (eigenvalues):  

  There exists a nonzero $ \varphi \in \mathcal{F} $ such that
  
  $$
  \mathcal{K} \varphi = \lambda \varphi \quad \Rightarrow \quad (\mathcal{K} - \lambda I)\varphi = 0.
  $$
  Hence, $ \mathcal{K} - \lambda I $ is **not injective**, so it is not invertible.  
  ✅ Not invertible because the kernel is nontrivial.

---

- **Continuous spectrum**:  

  The operator $ \mathcal{K} - \lambda I $ is **injective**, and its range is **dense**, but it has **no bounded inverse**.

  That is, although for every $ g \in \mathcal{F} $ and every $ \varepsilon > 0 $, we can find some $ f \in \mathcal{F} $ such that  
  
  $$
  \|(\mathcal{K} - \lambda I)f - g\| < \varepsilon,
  $$
  there exists **no continuous inverse** mapping $ g \mapsto f $.

  ✅ Not invertible because the inverse (if it exists) is **unbounded**.

---

- **Residual spectrum**:  

  The operator $ \mathcal{K} - \lambda I $ is **injective**, but its range is **not dense** in $ \mathcal{F} $.  
  This means there exist functions $ g \in \mathcal{F} $ that cannot even be approximated arbitrarily well by elements in the range.

  So $ \mathcal{K} - \lambda I $ is not surjective, nor does it have dense range.  
  ✅ Not invertible because the range is not dense — no inverse can be defined on all of $ \mathcal{F} $.


```

We now define the spectral radius

```{admonition} Definition: *Spectral radius*
:class: note
The **spectral radius** is given by

$$
\rho(\mathcal{K}) := \sup \{ |\lambda| : \lambda \in \sigma(\mathcal{K}) \}.
$$
```

If $ \mathcal{K} $ is a contraction (e.g., when defined on $ L^2_\mu $ with $ \mu $ invariant), then $ \rho(\mathcal{K}) \leq 1 $, so the spectrum lies inside or on the unit disk in the complex plane.


```{admonition} Remark: *Why the Continuous Spectrum Matters*
:class: tip

Even when Koopman eigenfunctions do not exist or are insufficient, the spectrum still encodes long-term features such as mixing, statistical decorrelation, and broadband frequency content. These are typically associated with the **continuous** part of the spectrum.
```

For example, in a strongly mixing system, the entire unit circle may belong to the continuous spectrum, with no nontrivial eigenfunctions. Yet observables still converge in distribution or decorrelate — phenomena fully captured by the operator's spectrum.


## Canonical Examples

We conclude this chapter with concrete examples illustrating the spectral behavior of the Koopman operator for well-understood dynamical systems. These help bridge the gap between abstract theory and explicit computation.


### Example 1: Rotation on the Circle

Consider the circle $ \mathbb{T} = \mathbb{R}/\mathbb{Z} $, and define the rotation map:
$$
T(x) = x + \alpha \mod 1,
$$
where $ \alpha \in \mathbb{R} $ is a fixed rotation angle.

This is a classic example of a **deterministic, measure-preserving**, and **non-mixing** system. Its dynamics are simple but exhibit rich spectral structure depending on whether $ \alpha $ is rational or irrational.

#### Invariant Measure

The **Lebesgue probability measure** on $ \mathbb{T} $ is invariant:

$$
\mu(B) = \int_B dx, \quad \text{for } B \subset \mathbb{T}.
$$
Indeed, $ T $ is a translation, and Lebesgue measure is translation-invariant.

#### Koopman Operator and Eigenfunctions

The Koopman operator $ \mathcal{K} $ acts on observables $ f \in L^2(\mathbb{T}) $ via composition:

$$
\mathcal{K} f(x) = f(T(x)) = f(x + \alpha).
$$

Let us consider the Fourier basis:

$$
\varphi_n(x) := e^{2\pi i n x}, \quad n \in \mathbb{Z}.
$$

Then,

$$
\mathcal{K} \varphi_n(x) = \varphi_n(x + \alpha) = e^{2\pi i n (x + \alpha)} = e^{2\pi i n \alpha} \varphi_n(x).
$$

Hence, each $ \varphi_n $ is a Koopman eigenfunction with eigenvalue $ \lambda_n = e^{2\pi i n \alpha} $. The Koopman spectrum is **pure point** and lies on the unit circle.


```{admonition} Interpretation
:class: tip

Each Koopman eigenfunction $ \varphi_n $ represents an oscillatory mode at frequency $ n $, and evolves as

$$
\varphi_n(T^t(x)) = e^{2\pi i n \alpha t} \varphi_n(x),
$$

which is a **quasiperiodic signal** unless $ \alpha \in \mathbb{Q} $, in which case the orbit is periodic. The Koopman operator captures this perfectly via its discrete spectrum.
```


### Example 2: Doubling Map on $[0,1]$

We now consider a simple yet chaotic system: the **doubling map** defined on the interval $[0,1]$ by

$$
T(x) = 2x \mod 1.
$$

This map stretches and folds the unit interval, and is a classic example of a **strongly mixing**, **expanding**, and **ergodic** transformation.


#### Invariant Measure

The **Lebesgue probability measure** on $[0,1]$ is invariant:

$$
\mu(B) = \int_B dx, \quad \text{for } B \subset [0,1].
$$

This can be verified directly via change of variables or by noting that the preimages under $ T $ evenly divide the interval.


#### Koopman Operator and Spectrum

Let $ \mathcal{K} f(x) = f(T(x)) = f(2x \mod 1) $ act on $ L^2([0,1]) $.

We ask whether $ \mathcal{K} $ has nontrivial eigenfunctions. Suppose $ \mathcal{K} \varphi = \lambda \varphi $. Then

$$
\varphi(2x \mod 1) = \lambda \varphi(x).
$$

It turns out that the **only** $ L^2 $ eigenfunction is the constant function $ \varphi(x) = 1 $, with eigenvalue $ \lambda = 1 $. All other candidate eigenfunctions either:

- Are not measurable or square-integrable,
- Or correspond to frequencies that are destroyed by the folding and mixing of the system.

Thus, the Koopman operator has **no nontrivial point spectrum** in $ L^2([0,1]) $, and instead exhibits **purely continuous spectrum**.


```{admonition} Interpretation
:class: tip

The absence of nontrivial eigenfunctions reflects the **chaotic, mixing nature** of the dynamics. While the system is deterministic and measure-preserving, observables rapidly lose correlation over time. Koopman analysis still applies, but not through modes — instead, the spectrum is continuous, reflecting a spread of frequencies and strong statistical mixing.
```


### Example 3: Squaring Map on $[0,1]$

Let $ \mathcal{X} = [0,1] $, and consider the map

$$
T(x) = x^2.
$$

This transformation is continuous and strictly increasing on $[0,1]$, with image $ T([0,1]) = [0,1] $. Hence, $ T $ is **invertible** on this domain. However, the map **compresses values toward 0**, and does so **very slowly near 1**. This distortion plays a crucial role in the spectral properties of the associated Koopman operator.



#### Invariant Measure

We define the Koopman operator on $ L^2([0,1]) $ as

$$
\mathcal{K} f(x) = f(x^2).
$$

Although Lebesgue measure is not invariant under $ T $, the Koopman operator is still **bounded** on $ L^2([0,1]) $, since composition with a measurable function preserves square integrability.


#### Koopman Operator and Spectrum

We now analyze the spectrum of $ \mathcal{K} $ acting on $ L^2([0,1]) $.

1. **Point spectrum**:  
   The constant function $ \varphi(x) = 1 $ satisfies $ \mathcal{K} \varphi = \varphi $, so
   
   $$
   \lambda = 1 \in \sigma_p(\mathcal{K}).
   $$
   No other $ L^2 $ eigenfunctions exist — that is, for $ \lambda \ne 1 $, there is no nontrivial $ \varphi \in L^2 $ such that $ \mathcal{K} \varphi = \lambda \varphi $.

2. **Injectivity**:  
   The operator $ \mathcal{K} - \lambda I $ is **injective** for all $ \lambda \in \mathbb{C} $, since eigenfunctions do not exist except for $ \lambda = 1 $.

3. **Range not dense**:  
   Although $ \mathcal{K} $ is injective, its **range is not dense** in $ L^2([0,1]) $.  
   This is because the image of any function $ f \mapsto f(x^2) $ is distorted: since $ x^2 \leq x $, the composition shifts the “information” of $ f $ toward smaller values of $ x $. This makes it impossible for $ \mathcal{K} f $ to approximate functions that vary significantly near $ x = 1 $, even in norm. Consequently, $ \mathcal{K} - \lambda I $ is not surjective and does not have dense range.



```{admonition} Interpretation
:class: tip

The Koopman operator $ \mathcal{K} f(x) = f(x^2) $ is injective but has **non-dense range**, so for any $ \lambda \ne 1 $, the operator $ \mathcal{K} - \lambda I $ is not invertible. This places $ \lambda \in \sigma_r(\mathcal{K}) $, the **residual spectrum**.

The residual spectrum arises because the composition $ f(x^2) $ cannot capture or approximate variations in $ L^2 $ functions near $ x = 1 $. This is due to the fact that $ x^2 \to 1 $ **very slowly** as $ x \to 1 $, causing the Koopman image to lack flexibility in that region.
```


The spectrum of the Koopman operator on $ L^2([0,1]) $ is:

$$
\sigma(\mathcal{K}) = \{1\} \cup \sigma_r(\mathcal{K}),
$$
where:
- $ \lambda = 1 \in \sigma_p(\mathcal{K}) $ is the only eigenvalue,
- All other $ \lambda \in \mathbb{C} \setminus \{1\} $ lie in the **residual spectrum** $ \sigma_r(\mathcal{K}) $,


## Conclusion

This chapter developed a functional-analytic perspective on deterministic dynamics through the **Koopman operator**. By lifting nonlinear systems into a linear space of observables, we gained access to powerful tools from operator theory and spectral analysis.

We explored:
- The duality between Koopman and Perron–Frobenius operators;
- The role of invariant measures in defining and analyzing Koopman operators;
- The impact of function space choice on boundedness and spectral properties;
- The classification of the spectrum into point, continuous, and residual parts;
- Concrete examples that illustrated each spectral type in action.

This spectral viewpoint provides deep insight into long-term behavior, even when trajectories are complex or chaotic. In the next chapter, we shift to the **data-driven approximation** of Koopman operators, and how these theoretical ideas can be applied directly from observations.

