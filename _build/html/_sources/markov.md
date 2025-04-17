# Markov Chains

Markov chains are one of the most fundamental models in probability theory and stochastic processes. They describe systems that evolve over time where the future state depends only on the current state—not the full history. In this course, we build a rigorous foundation for Markov chains using tools from measure theory and functional analysis. Our goal is to understand both the qualitative behavior (e.g., stability, long-term distribution) and the quantitative tools (e.g., operators, spectral properties) that govern their dynamics.

This course builds up the mathematical foundations needed to study Markov processes and their long-term behavior. We begin by formalizing probability through measurable and measure spaces, define stochastic processes and Markov chains, and introduce operator-theoretic tools such as the Markov transfer operator. These concepts culminate in a study of ergodicity, mixing, and invariant measures. Finally, we will generalize the operator viewpoint through the Koopman framework, connecting deterministic and stochastic dynamics in a unified way.


## Mathematical Foundations

To describe systems that evolve with uncertainty—such as populations, queues, or weather—we need a precise mathematical language for probability. This is provided by *measure theory*, which offers the tools to rigorously define events, random variables, and integration. In this section, we introduce the essential objects of measure theory: measurable spaces, measures, and probability spaces. These form the foundation for defining random variables and, ultimately, stochastic processes.



### Measurable and Measure Spaces

To rigorously define probability, we first need a mathematical structure that determines which subsets of a space can be assigned probabilities. This leads to the concept of a *measurable space*, which specifies the sets on which a measure—such as a probability measure—can be meaningfully defined. We then extend this structure to a *measure space*, where a measure quantifies these sets.

```{admonition} Definition: *Measurable space*
:class: note

A measurable space, also called *Borel space*, is a tuple $(\mathcal{X}, \Sigma_\mathcal{X})$, where:  
- $\mathcal{X}$ is a set, called the *state space*,  
- $\Sigma_\mathcal{X}$ is a $\sigma$-algebra on $\mathcal{X}$, i.e., a collection of subsets of $\mathcal{X}$ satisfying:  
  1. $\mathcal{X} \in \Sigma_\mathcal{X}$,  
  2. (closure under complementation) If $A \in \Sigma_\mathcal{X}$, then $\mathcal{X} \setminus A \in \Sigma_\mathcal{X}$,  
  3. (closure under countable unions) If $\{A_n\}_{n \in \mathbb{N}}$ is a countable collection of sets in $\Sigma_\mathcal{X}$, then   $\bigcup_{n=1}^{\infty} A_n \in \Sigma_\mathcal{X}$.  

The elements of $\Sigma_\mathcal{X}$ are called *measurable sets*, and the pair $(\mathcal{X}, \Sigma_\mathcal{X})$ forms a *measurable space*.  
```

A measure space extends a measurable space by introducing a measure that assigns sizes to measurable sets

```{admonition} Definition: *Measure space*
:class: note
A *measure space* is a tuple $(\mathcal{X}, \Sigma_\mathcal{X}, \mu)$, 
where:  
- $(\mathcal{X}, \Sigma_\mathcal{X})$ is a *measurable space*,  
- $\mu: \Sigma_\mathcal{X} \to [0, \infty]$ is a *measure*, meaning it satisfies:  
  1. $\mu(\emptyset) = 0$,  
  2. *(countable additivity)* If $\{A_n\}_{n \in \mathbb{N}}$ is a countable collection of disjoint sets in $\Sigma_\mathcal{X}$, then  
     $
     \mu\left(\bigcup_{n=1}^{\infty} A_n\right) = \sum_{n=1}^{\infty} \mu(A_n).  
     $
```

A measurable space $(\mathcal{X}, \Sigma_\mathcal{X})$ provides the foundation for defining measurable sets but does not assign numerical values to them. In contrast, a measure space $(\mathcal{X}, \Sigma_\mathcal{X}, \mu)$ introduces a measure $\mu$, which assigns a non-negative size to each measurable set, enabling integration and probability calculations.

### Probability Spaces and Random Variables

A probability space is a special case of a measure space where the measure is normalized to 1, allowing us to model uncertainty and compute probabilities of events.

```{admonition} Definition: *Probability space*
:class: note
A *probability space* is a tuple $(\Omega, \mathcal{F}, \mathbb{P})$, where:  
- $(\Omega, \mathcal{F})$ is a *measurable space* where $\Omega$ is called the *sample space*,  
- $\mathbb{P}: \mathcal{F} \to [0,1]$ is a *probability measure*, meaning it satisfies:  
  1. $\mathbb{P}(\emptyset) = 0$ and $\mathbb{P}(\Omega) = 1$,  
  2. *(countable additivity)* If $\{A_n\}_{n \in \mathbb{N}}$ is a countable collection of disjoint sets in $\mathcal{F}$, then $\mathbb{P}\left(\bigcup_{n=1}^{\infty} A_n\right) = \sum_{n=1}^{\infty} \mathbb{P}(A_n).$  
```

A probability space provides the foundation for defining randomness mathematically, allowing us to assign probabilities to events.

```{admonition} Definition: *Random variable*
:class: note
Let $(\Omega, \mathcal{F}, \mathbb{P})$ be a probability space and $(\mathcal{X}, \Sigma_\mathcal{X})$ a measurable space. Then an $(\mathcal{X}, \Sigma_\mathcal{X})$-valued *random variable* is a measurable function $X: \Omega \to \mathcal{X}$. In other words, $X$ must satisfy the following measurability condition.  

- For every measurable set $B \in \Sigma_\mathcal{X}$, the preimage   $X^{-1}(B) = \{ \omega \in \Omega : X(\omega) \in B \}\in\mathcal{F}$.  
```


```{admonition} Example: Rolling a fair die
:class: dropdown
Let's take the example of a die. We would like to properly define a random variable assigning 0 when the result of the die is odd and 1 if it is even.

**Step 1: Define the probability space**

A probability space consists of three components: a sample space, a sigma-algebra, and a probability measure.

- Sample space $(\Omega)$: The set of all possible outcomes when rolling a fair six-sided die: $\Omega = \{1, 2, 3, 4, 5, 6\}.$
- Sigma-algebra $(\mathcal{F})$: The power set $ \mathcal{F} = 2^\Omega $, i.e., all subsets of $\Omega$.  
- Probability measure $(P)$: Since the die is fair, the probability of each outcome is:  

  $$
  P(\{k\}) = \frac{1}{6}, \quad \text{for } k \in \{1,2,3,4,5,6\}.
  $$  
  The probability of any event (subset of $\Omega$) is the sum of the probabilities of its elements.

**Step 2: Define measure space**

Now, we introduce the state space where the random variable takes values.

- State space $(\mathcal{X})$: The possible values the random variable can take. Here, we set $\mathcal{X} = \{0,1\}$.  
- Sigma-algebra on $\mathcal{X}$ $(\Sigma_\mathcal{X})$: The smallest sigma-algebra making $\mathcal{X}$ measurable,  
  
  $$
  \Sigma_\mathcal{X} = \{\emptyset, \{0\}, \{1\}, \{0,1\} \}.
  $$
  
**Step 3: Define the random variable and verify measurability**

A random variable is a function that maps outcomes from $\Omega$ to $\mathcal{X}$ while preserving measurability.

- Random variable $(X: \Omega \to \mathcal{X})$: We define $X$ as an indicator function:  
  
  $$
  X(\omega) =  
  \begin{cases}  
  1, & \text{if } \omega \text{ is even (i.e., } \omega \in \{2,4,6\}) \\  
  0, & \text{otherwise.}  
  \end{cases}  
  $$  
- Measurability: For $X$ to be measurable, the preimage of every measurable set in $\Sigma_\mathcal{X}$ must be in $\mathcal{F}$:  
  - $X^{-1}(\{1\}) = \{2,4,6\} \in \mathcal{F}$ (since it's a subset of $\Omega$),  
  - $X^{-1}(\{0\}) = \{1,3,5\} \in \mathcal{F}$.  

Thus, $X$ is a measurable function and hence a valid random variable.

```

### Stochastic Processes

Now that we understand what a random variable is, we can extend this concept to define a stochastic process.  

```{admonition} Definition: *Stochastic process*
:class: note

A *stochastic process* is a collection of random variables defined on a common probability space $(\Omega, \mathcal{F}, P)$ and indexed by a set $T$ (interpreted as time). Formally, a stochastic process is a family of random variables  

$$
X = \{ X_t \}_{t \in T}
$$  

where each $X_t$ is an $(\mathcal{X}, \Sigma_\mathcal{X})$-valued random variable, meaning that $X_t: \Omega \to \mathcal{X}$ is measurable with respect to $\mathcal{F}$ and $\Sigma_\mathcal{X}$.

```

For any $t \in \mathbb{N}$ and any measurable set $B \in \Sigma_{\mathcal{X}}$, we consider the probability that the process is in $B$ at time $t+1$.  

Stochastic processes provide a flexible framework for modeling systems evolving under uncertainty. By indexing a family of random variables over time, they allow us to study how randomness propagates and how dependencies arise. In the following sections, we will focus on a particularly important class: Markov chains, where the future depends only on the present. This structure will enable us to introduce operator-theoretic tools that describe how probability and information evolve over time.


## Markov Chains

Markov chains model stochastic systems where the future state depends only on the present state. In this section, we introduce the Markov property, define time-homogeneity, and interpret the transition kernel.

### Markov Property

```{admonition} Definition: *Markov chain*
:class: note

Let $\{X_t\}_{t \in \mathbb{N}}$ be a stochastic process taking values in a measurable space $(\mathcal{X}, \Sigma_{\mathcal{X}})$ and defined on a probability space $(\Omega, \mathcal{F}, \mathbb{P})$. We say that $X$ is a *Markov chain* if, for all $t \in \mathbb{N}$ and all measurable sets $B \in \Sigma_{\mathcal{X}}$,  

$$
\mathbb{P}(X_{t+1} \in B \mid X_{[t]}) = \mathbb{P}(X_{t+1} \in B \mid X_t),
$$

where $X_{[t]} = \{ X_s : s \leq t \}$ represents the history of the process up to time $t$.  

This condition, known as the *Markov property*, states that the future state $X_{t+1}$ depends only on the present state $X_t$ and not on past states. In other words, the process has no memory beyond its current state.  

```

This Markov property expresses the memoryless nature of the process: the future depends only on the present, not the full past.


### Time-Homogeneity

A Markov chain $X := \{X_t : t \in \mathbb{N}\}$ is said to be *time-homogeneous* if there exists a function  

$$
p: \mathcal{X} \times \Sigma_{\mathcal{X}} \to [0,1],
$$  

called the *transition kernel*, such that for every state $x \in \mathcal{X}$, every measurable set $B \in \Sigma_{\mathcal{X}}$, and every time step $t \in \mathbb{N}$,  

$$
P(X_{t+1} \in B \mid X_t = x) = p(x, B).
$$  

This condition means that the probability of transitioning from state $x$ to a set of states $B$ depends only on $x$ and not explicitly on time $t$. In other words, the transition probabilities remain the same at all times.  


```{admonition} Remark: *Time-homogeneity*
:class: tip
Time-homogeneity simplifies the analysis of Markov chains, as the transition probabilities do not change over time. This allows for:  
- The use of the *Chapman-Kolmogorov equations* to describe multi-step transitions,  
- The definition of an associated transition operator or stochastic matrix in the discrete case,  
- The study of long-term behavior such as stationary distributions and ergodicity.  

If the Markov chain is not time-homogeneous, the transition probabilities depend explicitly on $t$, so we would need to work with a time-dependent transition kernel $p_t(x, B)$. This significantly complicates analysis and requires more general techniques.
```

```{admonition} Example: *Simple Weather Model*
:class: dropdown

Suppose $\mathcal{X} = \{\text{Sunny}, \text{Rainy}\}$, and define a Markov chain where:

- If today is Sunny, tomorrow is Sunny with probability 0.9 and Rainy with probability 0.1.
- If today is Rainy, tomorrow is Sunny with probability 0.5 and Rainy with probability 0.5.

Then the transition kernel $p(x, \cdot)$ is defined by:

$$
p(\text{Sunny}, \cdot) = \begin{cases}
0.9 & \text{Sunny} \\
0.1 & \text{Rainy}
\end{cases}, \quad
p(\text{Rainy}, \cdot) = \begin{cases}
0.5 & \text{Sunny} \\
0.5 & \text{Rainy}
\end{cases}
$$

This defines a time-homogeneous Markov chain on a finite state space.
```

A special case of interest is when the Markov chain satisfies **detailed balance**, meaning that the process is reversible with respect to $\pi$:  

```{admonition} Definition: *Reversible Markov Chain*  
:class: note
A Markov chain is **reversible** with respect to a probability measure $\pi$ if  

$$
\pi(dx) p(x, dy) = \pi(dy) p(y, dx), \quad \forall x, y \in \mathcal{X}.
$$

```

### Structure of the Transition Kernel  

The function $p(x, B)$ represents the probability that the Markov chain moves to a state in $B$ given that the current state is $x$. More formally, for each fixed $x \in \mathcal{X}$, the mapping  

$$
B \mapsto p(x, B)
$$

defines a probability measure on the measurable space $(\mathcal{X}, \Sigma_{\mathcal{X}})$. This ensures that:  
1. $p(x, B) \geq 0$ for all $x \in \mathcal{X}$ and $B \in \Sigma_{\mathcal{X}}$,  
2. $p(x, \mathcal{X}) = 1$ for all $x \in \mathcal{X}$, ensuring total probability is preserved,  
3. If $\{B_n\}_{n \in \mathbb{N}}$ is a countable collection of disjoint sets in $\Sigma_{\mathcal{X}}$, then $p(x, \bigcup_{n=1}^{\infty} B_n) = \sum_{n=1}^{\infty} p(x, B_n)$ (This follows from the definition of a probability measure).  


Beyond describing transitions between states, we often want to understand how entire distributions or observable functions evolve. To formalize this, we turn to Markov operators—linear operators acting on functions or measures—that encode the dynamics of a Markov chain in a functional framework.

## Markov Operators  

Markov processes are often studied through the evolution of sample paths or distributions. However, a powerful alternative approach is to consider how these processes act on functions and measures using linear operators. This operator-theoretic viewpoint not only simplifies the analysis of convergence and equilibrium but also connects to spectral theory, functional analysis, and numerical approximation. It forms the basis of modern perspectives on dynamical systems, such as the Koopman operator, which we will explore in the next section.


In many applications, we are interested in how functions of the state evolve over time rather than just tracking the states themselves. The **Markov transfer operator** provides a way to describe this evolution. However, beyond functions, measures also evolve under the transition dynamics, leading to a dual operator that governs their behavior. Understanding both perspectives is key to analyzing Markov processes.



### Evolution of Functions

```{admonition} Definition: *Markov Transfer Operator*  
:class: note
Let $\mathcal{F}$ be a set of real-valued measurable functions on $\mathcal{X}$. The *Markov transfer operator* (or *Markov operator*) $\mathcal{A}_F : \mathcal{F} \to \mathcal{F}$ is defined by  

$$
\mathcal{A}_F f (x) := \int_{\mathcal{X}} p(x, dy) f(y) = \mathbb{E} [ f(X_{t+1}) \mid X_t = x ], \quad f \in \mathcal{F}, \quad x \in \mathcal{X}.
$$

```  

This operator describes how a function $ f $ is transformed under the action of the Markov process. It effectively propagates functions forward in time according to the transition dynamics of the Markov chain.

```{admonition} Remark: *Choice of function space $\mathcal{F}$*
:class: tip
A common choice of function space is $\mathcal{F} = L^\infty(\mathcal{X})$, the space of bounded functions on $\mathcal{X}$. However, in many cases, we are interested in function spaces related to the existence of an *invariant measure* $\pi$ (see after).  

In this setting, we will that that it is natural to consider $\mathcal{F} = L^2_\pi(\mathcal{X})$, the space of square-integrable functions with respect to $\pi$. 
```

In dynamical systems, $\mathcal{A}_F$ is known as the (stochastic) **Koopman operator** on the space of observables $\mathcal{F}$.  




### Evolution of Measures 

Beyond acting on functions, the Markov transfer operator induces an evolution on measures. Given a measure $\mu$ on $\mathcal{X}$, we define the **dual operator** $\mathcal{A}^*$, which acts on measures as follows:  

```{admonition} Definition: *Markov Dual (Pushforward) Operator*  
:class: note
Given a measure $\mu$ on $\mathcal{X}$, the **Markov pushforward operator** $\mathcal{A}^*$ is defined by  

$$
(\mathcal{A}^* \mu)(B) := \int_{\mathcal{X}} p(x, B) \, d\mu(x), \quad B \in \Sigma_{\mathcal{X}}.
$$

```  

This operator describes how a measure $\mu$ evolves under the dynamics of the Markov chain. In particular, if $\mu_t$ represents the distribution of $X_t$, then  

$$
\mu_{t+1} = \mathcal{A}^* \mu_t.
$$

Thus, the Markov transition kernel induces both an evolution on functions (through $\mathcal{A}_F$) and an evolution on measures (through $\mathcal{A}^*$).  


### Link Between the Two Operators

The relationship between $\mathcal{A}_F$ and $\mathcal{A}^*$ is fundamental, as they offer dual perspectives on the evolution of a Markov process. Formally, they are dual with respect to integration against a measure $\mu$, meaning that for suitable functions $f$ and measures $\mu$,  

$$
\int_{\mathcal{X}} (\mathcal{A}_F f)(x) d\mu(x) = \int_{\mathcal{X}} f(x) d(\mathcal{A}^* \mu)(x).
$$  

This duality is crucial: while $\mathcal{A}_F$ describes the evolution of observables, $\mathcal{A}^*$ governs the evolution of probability measures. A key consequence is that if $\mathcal{A}^*$ admits a fixed point—an **invariant measure** $\pi$ satisfying  

$$
\mathcal{A}^* \pi = \pi,
$$

then applying the duality relation shows that, under $\pi$, expectations remain unchanged under $\mathcal{A}_F$:  

$$
\int_{\mathcal{X}} (\mathcal{A}_F f)(x) d\pi(x) = \int_{\mathcal{X}} f(x) d\pi(x).
$$  

This means that statistical averages computed from $\pi$ remain stable over time, motivating the study of invariant measures as they characterize the **stationary statistical behavior** of the system.  



## Invariant Measure and Stationarity

To understand the long-term behavior of a Markov chain, we must examine whether the system stabilizes over time. This leads to the notion of an *invariant measure*—a distribution that remains unchanged as the process evolves. Invariant measures capture the idea of statistical equilibrium and are key to understanding ergodicity, convergence, and time-averaged behavior in stochastic systems.


### Definition and Motivation

A central question in Markov processes is whether a **stationary distribution** exists, meaning a probability measure $\pi$ that remains unchanged under the action of $\mathcal{A}^*$. This corresponds to solving  

$$
\mathcal{A}^* \pi = \pi.
$$  

```{admonition} Definition: *Invariant Measure*  
:class: note
A measure $\pi$ on $\mathcal{X}$ is called **invariant** (or **stationary**) for the Markov transition kernel $ p(x, dy) $ if it satisfies  

$$
\pi(B) = \int_{\mathcal{X}} \pi(dx) p(x, B), \quad \forall B \in \Sigma_{\mathcal{X}}.
$$  
```  

So if $X_0$ is distributed according to $\pi$, then $X_t$ remains distributed according to $\pi$ for all $t$. 

```{admonition} Remark: *Reversibility implies invariance*
:class: tip

Suppose a probability measure $\pi$ satisfies the detailed balance condition

$$
\pi(dx) p(x, dy) = \pi(dy) p(y, dx).
$$

Then $\pi$ is automatically invariant, i.e.,

$$
\pi(B) = \int_{\mathcal{X}} \pi(dx) p(x, B), \quad \forall B \in \Sigma_{\mathcal{X}}.
$$

This is why reversibility is usually studied with respect to an invariant measure. Reversibility not only implies invariance, but also leads to useful analytical properties, such as self-adjointness of the Markov operator in $L^2_\pi$.
```


The existence and uniqueness of an invariant measure are central to understanding the **long-term behavior** of the Markov process. If such a measure exists and is unique, it plays a fundamental role in describing **statistical equilibrium**.  

Why study invariant measures?

1. **Ergodicity**: If an invariant measure $\pi$ exists and is unique, it describes the long-term statistical behavior of the process. Specifically, for an ergodic Markov process, the time average of any integrable function $f$ along a trajectory converges to its expectation under $\pi$:  

   $$
   \frac{1}{T} \sum_{t=0}^{T-1} f(X_t) \xrightarrow[T \to \infty]{} \int_{\mathcal{X}} f(x) d\pi(x) \quad \text{almost surely}.
   $$  

   This ensures that the system does not depend on initial conditions in the long run and justifies using $\pi$ as a statistical description of the system.  

2. **Spectral Properties in the Ergodic Setting**: The existence of an invariant measure is often related to the spectral properties of $\mathcal{A}_F$, particularly when analyzing convergence rates to equilibrium and mixing behavior. If the Markov process is ergodic, then $\mathcal{A}_F$ has a leading eigenvalue $\lambda_1 = 1$, corresponding to the invariant measure $\pi$ (i.e., $\mathcal{A}^* \pi = \pi$), and the rest of the spectrum lies strictly inside the unit disk:  

   $$
   1 = \lambda_1 > |\lambda_2| \geq |\lambda_3| \geq \dots
   $$  

   The spectral gap, defined as  

   $$
   \gamma = 1 - |\lambda_2|,
   $$  

   controls the rate of convergence to the invariant measure. A large spectral gap implies fast mixing, while a small gap leads to slow convergence. If the process is non-ergodic, the spectrum may contain multiple eigenvalues of modulus 1, corresponding to multiple invariant measures or ergodic components.  

3. **Physical Relevance**: Many Markovian systems in physics, such as Langevin dynamics, Monte Carlo methods, and stochastic differential equations, admit an invariant measure that represents equilibrium states. Understanding the properties of $\pi$ is crucial in applications such as molecular dynamics, statistical physics, and Bayesian inference, where algorithms like MCMC rely on ergodic Markov chains to sample from $\pi$.  


### Existence

A fundamental question in Markov process theory is whether an invariant measure $\pi$ exists. {cite}`1996_Prato_D_book_ergodicity` We first recall a general existence result:  

```{admonition} **Theorem (Krylov-Bogoliubov)**  
:class: important  
Let $\mathcal{X}$ be a Polish space (complete, separable, and metric), and let $\mathcal{A}^*$ be the Markov dual operator. If there exists a sequence of probability measures $\mu_n$ such that  

$$
\mu_{n+1} = \mathcal{A}^* \mu_n
$$  

and if the sequence $ (\mu_n)$ is **tight** (i.e., does not escape to infinity), then there exists a subsequence $ (\mu_{n_k}) $ that converges weakly to an invariant measure $\pi$.
```  

This result guarantees the existence of an invariant measure under mild conditions but does not ensure uniqueness. We now derive a useful corollary when the state space is compact.  

```{admonition} **Corollary**  
:class: important  
Let $\mathcal{X}$ be a **compact** metric space, and suppose the Markov dual operator $\mathcal{A}^*$ is **Feller** (i.e., it maps continuous functions to continuous functions). Then an invariant probability measure $\pi$ exists.  
```  

```{admonition} Proof
:class: dropdown 
Since $\mathcal{X}$ is compact, the space of probability measures on $\mathcal{X}$ (endowed with the weak topology) is also compact by **Prohorov’s theorem**. If we construct a sequence of measures $\mu_n$ using the Krylov-Bogoliubov procedure, compactness guarantees that a subsequence has a weak limit. The Feller property ensures that the limiting measure $\pi$ satisfies the invariance equation

$$
\mathcal{A}^* \pi = \pi,
$$  

completing the proof.  
```

This corollary is frequently used in applications where compactness holds, such as in stochastic models with bounded state spaces.  


### Uniqueness

When is the invariant measure unique? This depends on additional conditions:

```{admonition} Theorem: *Uniqueness of the Invariant Measure*  
:class: important  
If a Markov chain is **irreducible** and **aperiodic**, then the invariant measure $\pi$ is unique, and for any initial distribution $\mu_0$, we have:  

$$
\mathcal{A}^{*t} \mu_0 \to \pi \quad \text{as} \quad t \to \infty.
$$
```

This follows from the **Perron-Frobenius theorem**, ensuring that $\mathcal{A}^*$ has a unique leading eigenvalue $1$, with all other eigenvalues strictly less than $1$ in absolute value.  

A stronger condition guaranteeing uniqueness is **Doeblin’s condition**:

```{admonition} Proposition: *Doeblin’s Condition and Exponential Convergence*  
:class: important  
If there exists $\delta > 0$ and a probability measure $\nu$ such that for some $t_0 > 0$, the transition kernel satisfies  

$$
 p^{t_0}(x, dy) \geq \delta \nu(dy), \quad \forall x \in \mathcal{X},
$$  

then the Markov chain converges to $\pi$ exponentially fast:  

$$
\|\mathcal{A}^{*t} \mu_0 - \pi\|_{\text{TV}} = O(e^{-\gamma t}),
$$  

for some rate $\gamma > 0$, where $\|\cdot\|_{\text{TV}}$ is the total variation norm.
```


### Implications

When an invariant measure $\pi$ exists and is unique, it has profound implications on the choice of function space and on the properties of the Markov transfer operator. While a common choice is to work with the space of bounded functions $\mathcal{F} = L^\infty(\mathcal{X})$, the existence of $\pi$ allows us to consider the Hilbert space $\mathcal{F} = L^2_\pi(\mathcal{X})$, which consists of functions that are square-integrable with respect to $\pi$. The corresponding Markov transfer operator, denoted $\mathcal{A}_\pi$, which acts on $L^2_\pi(\mathcal{X})$ and is known to satisfy

$$
\|\mathcal{A}_\pi\| \leq 1.
$$
```{admonition} Proof
:class: dropdown

Let $\mathcal{A}_\pi$ be the Markov transfer operator defined on $L^2_\pi(\mathcal{X})$ by  

$$
(\mathcal{A}_\pi f)(x) = \int_\mathcal{X} p(x, dy) f(y).
$$

We aim to show that $\|\mathcal{A}_\pi\| \leq 1$ where the norm is the operator norm induced by the $L^2_\pi$ norm.
Consider the $L^2_\pi$-inner product:  

$$
\langle f, g \rangle_\pi = \int_\mathcal{X} f(x) g(x) \pi(dx).
$$  

We compute the inner product of $\mathcal{A}_\pi f$ and $\mathcal{A}_\pi g$:  

$$
\langle \mathcal{A}_\pi f, \mathcal{A}_\pi g \rangle_\pi
= \int_\mathcal{X} (\mathcal{A}_\pi f)(x) (\mathcal{A}_\pi g)(x) \pi(dx).
$$  

Substituting the definition of $\mathcal{A}_\pi$,  

$$
\int_\mathcal{X} \left( \int_\mathcal{X} p(x, dy) f(y) \right) \left( \int_\mathcal{X} p(x, dz) g(z) \right) \pi(dx).
$$  

Using Fubini's theorem and the definition of an invariant measure, we obtain  

$$
\int_\mathcal{X} f(y) g(z) \left( \int_\mathcal{X} p(x, dy) p(x, dz) \pi(dx) \right).
$$  

By Jensen’s inequality and the fact that $p(x, dy)$ is a probability kernel, it follows that  

$$
\|\mathcal{A}_\pi f\|_{L^2_\pi}^2 \leq \|f\|_{L^2_\pi}^2.
$$  

From the above inequality, we conclude that $\|\mathcal{A}_\pi f\|_{L^2_\pi} \leq \|f\|_{L^2_\pi}$.
Taking the supremum over all $ f $ with $ \|f\|_{L^2_\pi} \leq 1 $, we obtain $
\|\mathcal{A}_\pi\| \leq 1.$ Thus, $\mathcal{A}_\pi$ is a bounded linear operator with norm at most 1, completing the proof. \qed
```

which implies that it is a *bounded linear operator* on $L^2_\pi(\mathcal{X})$. This bound follows from the fact that the transition kernel $p(x, dy)$ defines a Markov operator that preserves probability measures. More precisely, this inequality indicates that $\mathcal{A}_\pi$ is a contraction (or non-expansive) operator in $L^2_\pi(\mathcal{X})$. Moreover, in the special case of reversible Markov chains, where the detailed balance condition holds, the operator $\mathcal{A}_\pi$ becomes self-adjoint in $L^2_\pi(\mathcal{X})$. Self-adjointness is crucial for a clear spectral decomposition of $\mathcal{A}_\pi$, which in turn directly relates to the convergence rate to equilibrium. Specifically, the existence of a spectral gap (i.e., the difference between the leading eigenvalue $1$ and the second-largest eigenvalue in absolute value) quantitatively determines the mixing rate of the Markov chain.

Thus, the invariant measure not only characterizes the equilibrium state but also underpins the analytical properties of the Markov transfer operator, such as stability, convergence, and mixing behavior.

## Ergodicity and Mixing

Once an invariant measure exists, a natural question arises: does the system converge to it? This section introduces ergodicity and mixing, which formalize the idea that a system forgets its initial condition and settles into stable, predictable behavior. *Ergodicity* describes convergence in distribution; *mixing* quantifies the rate at which statistical dependencies vanish.



### Ergodicity

Ergodicity expresses the idea that a system will, over time, explore the entire space in a statistically uniform way. Whether we are modeling physical particles, economic states, or sequences of decisions, we often want to know: *does the system stabilize?* *Does it "average out" in a way that no part of the state space is neglected in the long run?*

Ergodicity can be understood from either the dual perspective, which describes how probability distributions evolve over time, or the primal perspective, which focuses on the evolution of functions (observables). 


---

From the **dual operator’s** view, ergodicity corresponds to **convergence of the law of the process** toward the invariant measure $\pi$, regardless of the initial distribution.

```{admonition} Definition: *Ergodic Markov Chain (dual viewpoint)*  
:class: note

A time-homogeneous Markov chain $\{X_t\}_{t \in \mathbb{N}}$ with transition kernel $p(x, \cdot)$ and invariant measure $\pi$ is said to be **ergodic** if, for any initial distribution $\mu_0$, we have

$$
\mathcal{A}^{*t} \mu_0 \xrightarrow[t \to \infty]{} \pi,
$$

in the sense of weak convergence of measures. Equivalently, for all bounded measurable functions $f : \mathcal{X} \to \mathbb{R}$,

$$
\mathbb{E}_{\mu_0}[f(X_t)] \xrightarrow[t \to \infty]{} \int_{\mathcal{X}} f(x) \, d\pi(x).
$$
```

This means that the system forgets its initial condition in distribution, and time averages of observables converge to expectations under $\pi$.

---

From the **primal operator’s** view, ergodicity means that the **only fixed points** of the Markov transfer operator $\mathcal{A}_F$ in $L^2_\pi(\mathcal{X})$ are constant functions. That is, if $\mathcal{A}_F f = f$, then $f(x) = \int_{\mathcal{X}} f(y) \, d\pi(y)$ $\pi$-almost everywhere.

```{admonition} Definition: *Ergodic Markov Chain (primal viewpoint)*  
:class: note

Let $\pi$ be an invariant probability measure for a time-homogeneous Markov chain with transition kernel $p(x, \cdot)$. Consider the Markov transfer operator $\mathcal{A}_F$ acting on $L^2_\pi(\mathcal{X})$, the space of square-integrable functions with respect to $\pi$.

The chain is said to be **ergodic** if the only functions $f \in L^2_\pi(\mathcal{X})$ satisfying

$$
\mathcal{A}_F f = f
$$

are constant $\pi$-almost everywhere:

$$
f(x) = \int_{\mathcal{X}} f(y) \, d\pi(y).
$$
```

This expresses that repeated application of the transition kernel flattens any observable $f$ toward its mean under $\pi$:

$$
\mathcal{A}_F^t f(x) \xrightarrow[t \to \infty]{} \int_{\mathcal{X}} f(y) \, d\pi(y).
$$

This formulation is particularly useful in spectral analysis, where ergodicity implies that the eigenspace associated to the eigenvalue $1$ consists only of constant functions.

> **Note:** Unlike the dual formulation, this definition requires the operator $\mathcal{A}_F$ to be defined on a Hilbert space such as $L^2_\pi(\mathcal{X})$, which assumes the existence of an invariant measure and square-integrability of observables.



In summary, the **dual operator $\mathcal{A}^*$** pushes probability measures toward $\pi$, while the **primal operator $\mathcal{A}_F$** pushes observables toward constants. Both viewpoints describe the same asymptotic behavior: the system stabilizes to a long-term statistical equilibrium that is independent of its initial condition.


### Mixing and Decorrelation

Ergodicity guarantees convergence of distributions. **Mixing** is a stronger property: it describes how the dependence between past and future vanishes over time.

#### Strong Mixing (α-Mixing)

```{admonition} Definition: *Strong Mixing*  
:class: note

Let $\sigma(X_t)$ denote the sigma-algebra generated by $X_t$. A Markov chain is said to be **strongly mixing** if its $\alpha$-mixing coefficient

$$
\alpha(t) := \sup_{A \in \sigma(X_0), B \in \sigma(X_t)} \left| \mathbb{P}(A \cap B) - \mathbb{P}(A)\mathbb{P}(B) \right|
$$

satisfies $\alpha(t) \to 0$ as $t \to \infty$.
```

This means that events at time $0$ and time $t$ become asymptotically independent. Strong mixing implies ergodicity, but the converse is not necessarily true.

#### Weak Mixing

```{admonition} Definition: *Weak Mixing*  
:class: note

A Markov chain is said to be **weakly mixing** if for all bounded measurable functions $f, g : \mathcal{X} \to \mathbb{R}$,

$$
\text{Cov}(f(X_0), g(X_t)) \to 0 \quad \text{as } t \to \infty.
$$
```

This implies that the correlation between observations decays, even if independence is not achieved.


### Time Averages and Ergodic Theorem

Under ergodicity, long-term averages along trajectories converge to expectations under $\pi$. This is known as the *ergodic theorem*.

```{admonition} Theorem: *Ergodic Theorem*  
:class: important

Let $\{X_t\}_{t \in \mathbb{N}}$ be an ergodic Markov chain with invariant measure $\pi$. Then for any $f \in L^1_\pi(\mathcal{X})$,

$$
\frac{1}{T} \sum_{t=0}^{T-1} f(X_t) \xrightarrow[T \to \infty]{\text{a.s.}} \int_{\mathcal{X}} f(x) d\pi(x).
$$
```

This result justifies the use of $\pi$ as the *empirical average* distribution seen over time.



## Conclusion and Outlook

In this course, we developed a rigorous framework for analyzing Markov chains using the tools of measure theory and operator theory. We explored how the Markov transfer operator governs the evolution of functions (observables), while its dual describes the evolution of distributions. These perspectives are linked through a fundamental duality that illuminates the long-term behavior of stochastic processes. 

We saw how invariant measures characterize statistical equilibrium, how ergodicity captures asymptotic independence from initial conditions, and how spectral properties relate to convergence and mixing.

These theoretical tools are not just abstract constructions—they play a crucial role in practical areas such as statistical physics, machine learning (e.g., MCMC), and dynamical systems. Whether one is interested in the long-term behavior of a physical system or the convergence of a stochastic algorithm, the operator-theoretic view of Markov chains offers both clarity and depth.

This operator-theoretic viewpoint opens the door to more general frameworks where randomness may be replaced—or combined—with deterministic dynamics. In the next section, we’ll see how the **Koopman operator** generalizes this functional viewpoint to deterministic systems—bridging the gap between stochastic processes and nonlinear dynamics.


```{bibliography}
```
