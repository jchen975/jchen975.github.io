---
title: 'From geometric series to fixed point iterations and their applications: Part II'
date: 2025-08-22  
permalink: /posts/2025/08/geo-series-fixed point-2/
tags:
  - math
  - algorithm
---

Part II of a three-part note on how the common existence and uniqueness results in ODEs and MDPs boil down to the Banach fixed point theorem, which itself is based on the simple geometric series. 

_Acknowledgement: the content of this post, in particular the Banach fixed point theorem and its proof is based on Prof. Keith Conrad's [note](https://kconrad.math.uconn.edu/blurbs/analysis/contraction.pdf)._

In Part I of this three-part note, we saw how the geometric series is extended beyond a series of real or complex numbers, into operators and matrices. We introduced the Neumann series, and an application of matrix inverse computation using these tools. In this note, we look at general iterative algorithms and the famed Banach fixed point theorem, and understand how convergence to optimum may be connected to fixed points.  

Iterative algorithms and fixed points
-----
We consider a sequence of points, say \\(\\{x_n\\}_{n=0}^{\infty} \subseteq \mathbb{R}^n\\), whose elements are recursively defined by 

$$
\begin{equation}\label{eq:iterations}
    x_{n+1} \define T(x_{n}), \quad n \inZ_{\geq 0},
\end{equation}
$$

where \\(T : \mathbb{R}^n \goesto \mathbb{R}^n\\) is some fixed function. We may think of the elements in this sequence as _iterations_ of applying the function \\(T\\), since unrolling the recursion we would obtain \\(x\_{n} = T^{n}(x\_0)\\) for all \\(n \inZ_{\geq 0}\\). Common examples of this recursion scheme include (1) the Newton's method to find a root of some differentiable function \\(f : \mathbb{R} \goesto \mathbb{R}\\), which defines the iterations as  

$$
\begin{equation}\label{eq:newton}
    x_{n+1} = x_{n} - \frac{f(x_{n})}{f'(x_{n})},  \quad n \inZ_{\geq 0},
\end{equation}
$$

with \\(x_0\\) being the initial guess of the root, 
and (2) the gradient descent algorithm to find a local minimum on some differentiable objective function \\(f : \mathbb{R}^n \goesto \mathbb{R}\\), which defines the iterations as 

$$
\begin{equation}\label{eq:grad-descent}
    x_{n+1} = x_{n} - \eta \nabla f(x_{n}),  \quad n \inZ_{\geq 0},
\end{equation}
$$

where \\(\eta\\) is the step size, and \\(\nabla f\\) is the gradient of \\(f\\). 
All such iterative algorithms aims to stop at some desired point, which depending on the purpose, can be root of a function (as in Newton's method) or local minimum (as in gradient descent), or even global minimum/maximum. One common definition that can tie these points together, is the notion of _fixed point_. For our purposes, given a function \\(T : V \goesto V\\) with \\(V\\) being a set, \\(x \in V\\) is a fixed point of \\(T\\) if \\(T(x) = x\\). Thus, we can interpret the Newton step \eqref{eq:newton} by defining the map 

$$
    T(x) = x - \frac{f(x)}{f'(x)}. 
$$

Since \\(x^{\star}\\) is a root of \\(f\\) if and only if \\(f(x^{\star}) = 0\\), it amounts to the Newton step \eqref{eq:newton} reaching the point \\(x_{n+1} = T(x_n) = x_{n}\\), that is, a fixed point of \\(T\\). Similarly, we can define 

$$
    T(x) = x - \eta \nabla f(x)
$$

for gradient descent. A local minimum is an \\(x^{\star}\\) that satisfies \\(\nabla f(x^{\star}) = \zero\\), so once again it means that the gradient descent step has reached the point \\(x_{n+1} = T(x_n) = x_n\\), which is again a fixed point of \\(T\\).  

To visualize a fixed point of some function (if one exists) on the two-dimensional Euclidean space, we simply plot \\(x\\) against \\(T(x)\\), then its fixed points are all points that intersect with the line \\(y = x\\). Of course, not all functions have a fixed point, e.g., the real-valued function \\(T(x) = x+1\\) does not have a fixed point, the constant function \\(T(x) = 1\\) does not have a fixed point on the interval \\([0,1)\\), and a fixed point is not even defined for any matrixed-valued function \\(T : \mathbb{R} \goesto \mathbb{R}^{m\times n}\\). Generally, there are many results that provide sufficient conditions for a fixed point to exist, though not every one of them will tell us where it is. The most prominent example is the [Brouwer fixed point theorem](https://en.wikipedia.org/wiki/Brouwer_fixed point_theorem), which says a continuous function mapping a compact, convex set to itself always has at least one fixed point; the generalization of this theorem to set-valued maps, called the [Kakutani fixed point theorem](https://en.wikipedia.org/wiki/Kakutani_fixed point_theorem), lays the foundation for the existence of a Nash equilibrium in game theory. 

Banach fixed point theorem and its proof
-----
The fixed point theorem that we focus on in this note is called the _Banach fixed point theorem_, also called the _contractive mapping theorem_. Unlike Brouwer, which provides sufficient conditions for the _existence_ of a fixed point, the Banach fixed point theorem gives us conditions for when a _unique_ fixed point exists, and how to find it. The latter part is particularly useful in applied scenarios, since we wouldn't be able to do much with only the knowledge that a fixed point exists but not where it is. 

Formally, let us consider a complete metric space \\((M,d)\\), and we want to answer the following questions about a function \\(T : M \goesto M\\), in the following order: 

1. When does \\(T\\) have a fixed point? 
2. When is a fixed point \\(x^{\star}\\) of \\(T\\) unique? 
3. How to compute/find the unique fixed point \\(x^{\star}\\)? 

The answers to these questions lie in the _Banach fixed point theorem_ stated below. 

**Theorem** (Banach fixed point theorem). _Suppose that there exists some positive constant \\(c < 1\\) such that, for every \\(\,x, y \in M\\), the map \\(\,T\\) satisfies_

$$
\begin{equation}\label{eq:contraction}
    d(T(x), T(y)) \leq c d(x,y).
\end{equation}
$$

_Then, there exists a unique fixed point \\(\,x^{\star}\\) of \\(\,T\\) in \\(\,M\\), and every sequence \\(\\{x\_n\\}\_{n \inZ\_{\geq 0}}\\) such that \\(\,x\_n \in M\\) and \\(x_{n+1} \define T(x_n)\\) converges to \\(x^{\star}\\)._

_Proof_. First, we show that if \\(T\\) has a fixed point, then it must be unique. Assume towards contradiction that both \\(x, \tilde{x}\in M\\) are both fixed points of \\(T\\), but \\(x \neq \tilde{x}\\). Then, since \\(x = T(x)\\)  and \\(\tilde{x} = T(\tilde{x})\\), we must have 

$$
    d(x, \tilde{x}) = d(T(x), T(\tilde{x})) \leq cd(x, \tilde{x}) < d(x, \tilde{x}),
$$

which is a contradiction. Next, we show that the sequence \\(\\{x\_n\\}_{n \inZ\_{\geq 0}}\\) is Cauchy; in particular, since \\(M\\) is complete, it converges to a point in \\(M\\). By construction, we have 

$$
d(x_n, x_{n+1}) \leq cd(x_{n-1}, x_n) \leq c^2d(x_{n-2}, x_{n-1}) \leq \cdots \leq c^nd(x_0, x_1)
$$

for \\(n \in \\{1,2,...\\}\\). For any \\(m \inZ_{\geq 0}\\) such that \\(m > n\\), since \\(d\\) is a metric and must satisfy the triangle inequality, by the geometric series formula we also have 

$$
\begin{align*}
    d(x_n, x_m) &\leq d(x_n, x_{n+1}) + \cdots + d(x_{m-1}, x_m) \leq c^n d(x_0, x_1) + \cdots + c^{m-1}d(x_0, x_1) \\
                &= \sum_{i=n}^{m-1} c^i d(x_0, x_1)  \leq \left(\sum_{i=0}^{\infty} c^i - \sum_{i=0}^{n-1} c^i\right) d(x_0, x_1) \\
                &= \left(\frac{1}{1-c} - \frac{1-c^n}{1-c}\right) d(x_0, x_1) \\
                &= \frac{c^n}{1-c} d(x_0, x_1).
\end{align*}
$$

Since \\(c < 1\\), for any \\(\epsilon > 0\\), there is an integer \\(N > 0\\) such that 

$$
\frac{c^N}{1-c} d(x_0, x_1) < \epsilon.
$$

Thus, \\(d(x\_n, x\_m) < \epsilon\\) for all \\(n,m > N\\), so the sequence \\(\{x\_n\}\\) is Cauchy. Finally, let \\(\bar{x} \define \lim\_{n \goesto \infty} x\_n\\), and all we need to do is to prove that \\(T(\bar{x}) = \bar{x}\\) is that \\(T\\) is continuous, since by continuity, we can then conclude that \\(T(\bar{x}) = \lim\_{n \goesto \infty} T(x_n) = \lim\_{n\goesto \infty} x_{n+1} = \bar{x}\\). To see this, take any \\(\epsilon > 0\\) and set \\(\delta \define \frac{\epsilon}{c}\\), then for any \\(x, y \in M\\), by \eqref{eq:contraction} we get that 

$$
    d(x, y) < \delta \implies d(T(x), T(y)) \leq cd(x,y) < c\delta = \epsilon,
$$

i.e., \\(T\\) is not only continuous but uniformly continuous, and the proof is complete. \\(\endproof\\)

As we can see, while the Banach fixed point theorem relies on certain "niceness" properties of the metric space and knowledge of \\(\epsilon\\)-\\(\delta\\) definition of continuity, at its core the proof just utilizes the geometric series identity to obtain the final bound on the iterations, and thereby the fixed point itself. 

We make two additional remarks about the Banach fixed point theorem. First, any map \\(T\\) that satisfies \eqref{eq:contraction} is called a _contraction_ (thus the alternative name of the theorem, "contraction mapping theorem"), and the constant \\(c\\) is called the _contraction constant_. A contraction is by definition a Lipschitz continuous function with Lipschitz constant strictly less than 1. Second, the requirement that \\(T\\) maps a complete metric space \\(M\\) to itself can be restricted to a closed invariant subset of \\(M\\) under \\(T\\). That is, let \\(T\\) be a contraction on \\(M\\), and suppose that \\(\mathcal{S} \subseteq M\\) is a closed set with respect to the metric topology such that \\(T(\mathcal{S}) \subseteq \mathcal{S}\\), then the map \\(T\vert_{\mathcal{S}}: \mathcal{S} \goesto \mathcal{S}\\) is a contraction since \\(\mathcal{S}\\) contains all of its limit points, and thus, the fixed point \\(x^{\star}\\) as well. 


Finally, for practical reasons, we focus on the case of (1) Euclidean spaces, in which the \\(L_p\\) norm is a metric, and (2) \\(T\\) has additional smoothness properties such as continuous differentiability. The (induced) \\(L_p\\) norm of \\(T\\)'s Jacobian evaluated at some point \\(x \inR^n\\), denoted by \\(\Vert J(x)\Vert \\), can sometimes provide a more easily checkable sufficient condition for contractivity. The reason why this condition can be more practical than verifying \eqref{eq:contraction} point-wise is that, in most algorithmic convergence analysis from a fixed point perspective, we often know exactly what the iteration mapping \\(T\\) is, and if it is differentiable, then we can easily compute and bound the \\(L_p\\) norm of its Jacobian.  

**Theorem**. _Let \\(\,\mathcal{D} \subseteq \mathbb{R}^n\\) be a \\(\,T\\)-invariant convex set. Suppose that there exists some positive constant \\(\,c < 1\\) such that \\(\Vert J(x)\Vert \leq c \\) for every \\(\,\,x \in \mathcal{D}\\), then \\(\,T\\) is a contraction with contraction constant \\(\,c\\), with respect to the norm \\(\,\Vert \cdot \Vert\\), i.e., for all \\(\,x, y \in \mathcal{D}\\),_

$$
    \Vert T(x) - T(y) \Vert \leq c \Vert x - y \Vert. 
$$

_Proof_. Take any two points \\(x, y \in \mathcal{D}\\), and define the line segment \\(\Phi(t) \define tx + (1-t)y \\) parameterized by \\(t \in [0,1]\\), i.e., \\(\Phi(0) = y\\) and \\(\Phi(1) = x\\). Since \\(\mathcal{D}\\) is convex and \\(T\\)-invariant, every point on the line segment, and itself under the map \\(T\\), must be contained in \\(\mathcal{D}\\). By chain rule, 

$$
    \frac{\d}{\d t} T(\Phi(t)) = J(\Phi(t)) \nabla \Phi(t) = J(tx + (1-t)y) (x - y), 
$$

by the fundamental theorem of calculus, we then get that 

$$
    T(x) - T(y) = T(\Phi(1)) - T(\Phi(0)) = \int_{0}^1 \frac{\d}{\d t} T(\Phi(t)) \d t = \int_{0}^1  J(tx + (1-t)y) (x - y) \d t, 
$$

and taking the \\(L_p\\) norm of both sides, by triangle inequality and submultiplicativity we get 

$$
\begin{align*}
    \Vert T(x) - T(y) \Vert
    &= \left \Vert \int_{0}^1  J(tx + (1-t)y) (x - y) \d t \right \Vert  \\
    &\leq \int_{0}^1 \left \Vert  J(tx + (1-t)y) \right \Vert \d t \Vert x - y \Vert 
    \leq c \Vert x - y \Vert,
\end{align*}
$$

i.e., \\(T\\) is a contraction with contraction constant \\(c\\). \\(\endproof\\)


Conclusion
-----
In this note, we explored how some iterative numerical algorithms can be understood via fixed point iterations, the basic formulation of the Banach fixed point theorem and its proof, along with a practical sufficient condition for contractivity. This note, of course, only scratches the very top surface of the vast realm of fixed point theories, and we will explore a little more about its applications in the next and final part of this series of notes.  


<br> 

Further readings and references
-----

1. Keith Conrad's expository papers on Banach fixed point theorm [Part 1](https://kconrad.math.uconn.edu/blurbs/analysis/contraction.pdf) [Part 2](https://kconrad.math.uconn.edu/blurbs/analysis/contraction2.pdf)
2. [Banach fixed point theorem in \\(\mathbb{R}^n\\)](https://terpconnect.umd.edu/%7Epetersd/666/fixedpoint.pdf)
3. [Metric space](https://en.wikipedia.org/wiki/Metric_space)

_Last updated: 2025-09-01 5:02 EST_