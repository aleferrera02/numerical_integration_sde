# numerical_integration_sde

Stochastic optimization is a fundamental area in machine learning and large-scale optimization problems. A central algorithm in this context is Stochastic Gradient Descent (SGD), widely used to minimize complex objective functions, especially when the data is large and computing the full gradient becomes impractical. This project focuses on the analysis of the SGD algorithm, examining its interpretation as a Stochastic Differential Equation (SDE). This approach provides a deeper understanding of the algorithm's behavior, particularly in scenarios where the full gradient is unavailable, such as in large-scale optimization problems.

The SDE perspective allows for the study of the dynamics of SGD, with particular attention to its convergence properties and the role of noise. Specifically, we explore how the choice of step size $h$ affects the balance between exploration and convergence speed. Another aspect explored is Noisy Gradient Descent (NGD), a variant of SGD that introduces controlled noise through a “temperature” parameter. This noise helps NGD escape local minima, improving exploration in complex, non-convex landscapes.

Finally, the interpretation of NGD as an SDE provides a rigorous understanding of its dynamics, while also ensuring convergence guarantees. The analysis of both SGD and NGD highlights the power of stochastic methods as practical approximations of deterministic optimization algorithms, offering significant advantages in tackling high-dimensional and non-convex optimization tasks.
