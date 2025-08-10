# D-WAVE implementation of the reduced OneOpto model.

Author: Hao Mack

We will work in parallel to consider different classical and quantum methods to implement the _simplified OneOpto optimization model_. Our goal is to match the value of each characteristic per each risk group to its target as closely as possible while respecting the many guardrails.

The implementation is in the form of a single Jupyter notebook, and a sample JSON of selected bonds and their bounds created by Mackenson.

## Mathematical formulation

We are given three sets, $C$, $L$, $J$, and two global parameters, $N$, and $\vec{r}$ as the input. Our decision variables are the $\vec{y}$, indicating the particular bonds are included in the portfolio or not.

The objective function is a quadratic optimization function in the $\vec{y}$, subject to four linear inequality constraints.

### Securities

$C$ is denoted as the set of securities. Each security $c \in C$ is expressed as a quintuple $(p_c, m_c, M_c, i_c, \delta_c)$, where:

* $p_c$ is the market price
* $m_c$ and $M_c$ are the minimum and maximum trade value
* $i_c$ is the basket inventory
* $\delta_c$ is the minimum increment of a price

### Risk buckets

There is a set $L$ of risk buckets, and for each element of the risk bucket $\ell \in L$, there is a corresponding set of bond securities belonging to that risk bucket. Thus, there is a "set of sets" structure of the form $(\ell, \mathbb{K}_\ell)$.

* $L$ is the set of risk buckets.
* $\mathbb{K}_l$ is the set of bond securities in risk bucket $\ell$.
    * It is noted that it is not necessarily the case that the risk buckets are mutually exclusive.

### Characteristics and Guardrails

There is a set $J$ consisting of the characteristics of the investments. It is a strict subset of the fields that provide bond information. Optimization bonds are based on a combination of an element in $J$, and another element in either $C$ or $L$. $J$ is strictly a preference for the investor.

There are $|J|(2|L| + |C|)$ guardrails with a two-dimensional index, which one of the indices, $j$, is from $J$.

* $\rho_j$ represent the weight for that characteristic when determining the energy of the optimization problem to minimize.

Each element $j \in J$ is a triple $(\vec{K}_{j}, \vec{b}_{j}, \vec{\beta}_{j})$ which the length of the vector depends on $|L|$ for $\vec{K}_j$ and $\vec{b}_j$, and on $|C|$ for $\vec{\beta}_j$. The elements of each vector element of $\vec{K}_j$ and $\vec{b}_j$ are themselves an interval containing a guardrail and a target value.

* $K_{\ell, j} = [K^{\textrm{low}}_{\ell, j}, K^{\textrm{target}}_{\ell, j}, K^{\textrm{up}}_{\ell, j}]$ are the target and guardrails of a characteristic $j$ in risk bucket $\ell$.
* $b_{\ell, j} = [b^{\textrm{up}}_{\ell, j}, b^\textrm{low}_{\ell, j}]$ denote the benchmark guardrails.
* $\beta_{c, j}$ denote the contribution of a unit of bond $c$ to the target of characteristic $j$.

### Global parameters

Two global hyperparameters will be used to find a solution for various use cases.

* $N$ denotes the maximum number of bonds in the portfolio.
* $R = [m_R, M_R]$ denote the range of the residual cash flow of portfolio. 
    * Namely, $m_R$ denote the minimum cash flow, and $M_R$ denote the maximum cash flow. Note that they are signed quantities.

### Decision variable

Let $\vec{y}$ denote the solution set. Each entry of $\vec{y}$ is indexed by an element $c \in C$. The objective function involves the optimzal allocation of the bonds and how many of each bond to allocate. The amount of bonds to allocate, $x_c$, is simplified as a function of $c$ and $y$: $x_c = \frac{m_c + \min{M_c, i_c}}{2\delta_c}$, the average number of bonds from what is feasible for that $c$.

### Optimization problem

We must match the value of each characteristic $j \in J$ in each risk group $\ell \in L$ to its target, $K^{\textrm{target}}_{\ell, j}$, subject to four constraints. 

Note that if a quantity $s \in [I^{\textrm{low}}, I^{\textrm{target}}, I^{\textrm{up}}]$ or $s \in [I^{\textrm{low}}, I^{\textrm{up}}]$, the inequality becomes $I^{\textrm{low}} \leq s \leq I^{\textrm{up}}$.

$$
\begin{align*}
  E &= \min_{\vec{y}} \sum_{\ell \in L} \sum_{j \in J} \rho_j \left(\sum_{c\in{\mathbb{K}_{\ell}}} \beta_{c, j} x_c - K^{\textrm{target}}_{\ell, j}\right)^2 \\
  \sum_{c \in C} y_c &\leq N \\
  \sum_{c \in \mathbb{K}_{\ell}} \beta_{c,j}y_c &\in K_{\ell, j} \\
  \sum_{c \in C} \frac{p_c \delta_c}{100} x_c &\in R \\
  \sum_{c \in \mathbb{K}_{\ell}} \frac{p_c \delta_c}{100} \beta_{c,j} x_c &\in b_{\ell, j}
\end{align*}
$$