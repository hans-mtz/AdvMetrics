<!-- # ---
# document: "article"
# author: "Hans Martinez"
# title: "Advanced Metrics"
# date: "February 19,2021"
# --- -->

# Advanced Metrics
### Problem Set 1

**Hans Martinez**

Feb, 2021

### Question 1

Comparing different quadrature rules to calculate the integral:

$$ \Phi_{\mu,\sigma} = \int_{-\infty}^{a} \phi(x)dx $$
́
1. Gauss-Legendre between $( \mu - 4\sigma, a )$.
2. Gauss-Chebyshev with change of variables (COV).
3. Gauss-Hermite with COV.

For the Gauss-Chevyshev COV, I used $\rho(y)=ln(y+1)-ln(2)+a$, which maps $[-1,1] \to [-\infty,a]$. For Gauss-Hermite, $\rho(y)=-e^{-y}+a$ which maps $[ - \infty , \infty] \to [-\infty,a]$.

The next table summarizes the results.

| $\mu$ | $\sigma$ |a|Gauss-Legendre|Gauss-Chebyshev|Gauss-Hermite|
|----|-----|-----|--------------|---------------|-------------|
|-1  |0.5  |-2.5 |0.01460856    |0.01703302     |0.01628041   |
|-1  |0.5  |-1.75|0.14208332    |0.14488280     |0.14076644   |
|-1  |0.5  |-1   |0.49766113    |0.50080450     |0.49416151   |
|-1  |0.5  |-0.25|0.85323895    |0.85606443     |0.84864305   |
|-1  |0.5  |0.5  |0.98071371    |0.98223965     |0.97563008   |
|-1  |5    |-16  |0.00000000    |0.00000000     |0.00000000   |
|-1  |5    |-8.5 |0.00039812    |0.00039904     |0.00039085   |
|-1  |5    |-1   |0.49999976    |0.50394920     |0.49813474   |
|-1  |5    |6.5  |1.00001041    |0.52275421     |1.06847568   |
|-1  |5    |14   |1.01309803    |0.00007890     |0.21901183   |
|1   |0.5  |-0.5 |0.01460856    |0.01703302     |0.01628041   |
|1   |0.5  |0.25 |0.14208332    |0.14488280     |0.14076644   |
|1   |0.5  |1    |0.49766113    |0.50080450     |0.49416151   |
|1   |0.5  |1.75 |0.85323895    |0.85606443     |0.84864305   |
|1   |0.5  |2.5  |0.98071371    |0.98223965     |0.97563008   |
|1   |5    |-14  |0.00000000    |0.00000000     |0.00000000   |
|1   |5    |-6.5 |0.00039812    |0.00039904     |0.00039085   |
|1   |5    |1    |0.49999976    |0.50394920     |0.49813474   |
|1   |5    |8.5  |1.00001041    |0.52275421     |1.06847568   |
|1   |5    |16   |1.01309803    |0.00007890     |0.21901183   |


We can observe that integral results from the three quadrature rules are very close across the whole distribution when $\sigma$ is low. However, when the $\sigma$ is high, we start getting odd results in the upper tail of the distribution, like probabilities above 1 for example.

### Question 2
