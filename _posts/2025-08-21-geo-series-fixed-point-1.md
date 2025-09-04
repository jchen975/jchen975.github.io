---
title: 'From geometric series to fixed-point iterations and their applications: Part I'
date: 2025-08-21  
permalink: /blurbs/2025/08/geo-series-fixed-point-1/
tags:
  - math
---

Part I of a three-part note on how the common existence and uniqueness results in ODEs and MDPs boil down to the Banach fixed-point theorem, which itself is based on the simple geometric series. 


Most people will probably learn about the geometric series at some point during high school, when its uses during this first encounter are usually elementary and application-oriented. For example, the geometric series is used in calculating compound interest, rate of return, perhaps some physics phenomenon, though it has far more uses than one might think. In this three-part note, we will see (1) how the geometric series can be extended beyond a series of real or complex numbers, and (2) the central role of geometric series in proving the Banach fixed-point theorem, and (3) the two standard applications of the Banach fixed-point theorem: the Picard–Lindelöf theorem, and convergence of value iteration in Markov decision processes. 

<!-- The geometric series and its extension: the Neumann series
----- -->

To begin, we recall that a geometric series takes the following form: for some number \\(k\\) and  \\(a\\), 

$$
\begin{equation}\label{eq:geo-series-real}
    k + ka + ka^2 + \dots = \sum_{i = 0}^{\infty} ka^i, 
\end{equation}
$$

is a geometric series with a common ratio \\(a\\). If \\(a \< 1\\), it is called the "decay rate" and the series converges to \\(\frac{k}{1-a}\\). One important extension of this idea is the power series, where the leading coefficient \\(k\\) is no longer uniform, but instead itself forms a sequence, i.e., 

$$
    k_0 + k_1 a + k_2 a^2 + \dots = \sum_{i = 0}^{\infty} k_i a^i, 
$$

which is important and a powerful tool in its own right, e.g., in the context of analysis. To stay focused on the purpose of this note, we instead look at the extension into geometric series of a _linear operator_ \\(T: V \goesto V\\), where we take \\(V\\) to be a [complete](https://en.wikipedia.org/wiki/Cauchy_sequence#Completeness), [finite-dimensional](https://en.wikipedia.org/wiki/Dimension_(vector_space)), [normed vector space](https://en.wikipedia.org/wiki/Normed_vector_space), also known as a _Banach space_. The standard \\(n\\)-dimensional Euclidean space \\(\mathbb{R}^n\\) is one such example. We begin by taking \\(k = 1\\) and \\(a < 1\\) in \eqref{eq:geo-series-real}, which along with the convergence result stated previously, becomes 

$$
\begin{equation}\label{eq:geo-series-real-1}
    \sum_{i = 0}^{\infty} a^i = (1-a)^{-1}.
\end{equation}
$$

With the linear operator \\(T: V \goesto V\\), we can define the _Neumann series_: 

$$
\begin{equation}\label{eq:neumann-series}
    \sum_{i = 0}^{\infty} T^i,
\end{equation}
$$

which is a series of _repeated application_ of the operator \\(T\\). 
<!-- To make this notion concrete, let us 
see 
how the gradient descent algorithm on an objective function \\(f : \mathbb{R} \goesot \mathbb{R}\\) might fit into this framework. Let \\(x_0 \inR\\) be some fixed initial condition, then we get 
$$
    x_n = x_{n-1} - \eta f'(x_{n-1}) \define T(x_{n-1}). 
$$
That is, the gradient descent iterations can be viewed as the sequence 
$$
    x_1 \define T(x_0), \; x_2 \define T(x_1) = T^2(x_0), \;, \dots, x_n \define T(x_{n-1}) = \cdots = T^{n}(x_0).
$$
With the simple quadratic objective function \\(f(x) = \frac{1}{2}x^2\\), we have \\(f'(x) = x\\), so 
$$
    x_{n} = x_{n-1} - \eta x_{n-1} = (1 - \eta) x_{n-1} = (1-\eta)^2 x_{n-2} = \cdots (1-\eta)^n x_0, 
$$
and the step size \\(\eta\\) to ensure algorithm convergence can be naturally obtained via the geometric series convergence criterion.   -->
<!-- Back to the general Neumann series \eqref{eq:neumann-series}.  -->
Ideally, the series converges, but what's more, it would be nice to arrive at an analogous result to \eqref{eq:geo-series-real-1}. That is, we need to answer the following question: under what condition on \\(T\\) do we have \eqref{eq:neumann-series-conv}? 

$$
\begin{equation}\label{eq:neumann-series-conv}
    \sum_{i = 0}^{\infty} T^i =(I-T)^{-1}
\end{equation}
$$ 

The answer is, when the operator norm of \\(T\\), denoted by \\(\Vert T \Vert\\), satisfies \\(\Vert T \Vert < 1\\). This condition is clearly analogous to the convergence criterion of geometric series, \\(\vert a \vert \< 1 \\). The Neumann series has important applications in the topic of _resolvent analysis_, a topic in complex analysis with deep connections to linear algebra, and by extension, classical control theory topics such as the frequency response of a linear time-invariant system. An excellent treatment on this topic is given in the reference section at the end. This extension from the (scalar) geometric series to the (vector/operator) Neumann series will be examined further in Part II, where we introduce the Banach fixed-point theorem and its proof.  

In addition, since we are only working with finite-dimensional vector spaces, it is natural to consider the square matrix \\(A\\), which is the matrix representation of operator \\(T\\) in the (canonical) basis. The ties to linear algebra now becomes evident: if the induced matrix norm of \\(A\\), denoted by \\(\Vert A \Vert\\), satisfies \\(\Vert A \Vert < 1\\), we can then write 

$$
    \sum_{i = 0}^{\infty} A^i =(I-A)^{-1},
$$ 

which is a series of matrix powers. Thus, when the series converges, a natural application of this result is in matrix inversion computation: let \\(B \define I - A\\), and suppose that \\(\Vert I - B\Vert < 1\\). Then, to compute the inverse \\(B^{-1}\\), we can write 

$$
    B^{-1} = (I - A)^{-1} = \sum_{i = 0}^{\infty} A^i = \sum_{i = 0}^{\infty} (I-B)^i.
$$ 

In practice, we may truncate the series after sufficiently large number of terms. In addition, a common trick used in computing the inverse of a sum of two matrices, \\((A + B)^{-1}\\), when one of them is invertible, is as follows: suppose that \\(A\\) is invertible and \\(\Vert BA \Vert < 1\\), then 

$$
    (A + B)^{-1} = ((I + BA)A)^{-1} = A^{-1}(I + BA)^{-1} = A^{-1} \sum_{i=0}^{\infty} (-BA)^{i},
$$

which is a particularly useful formula in cases of \\(B\\) being a small perturbation to the "nominal" \\(A\\). 


To conclude this very brief note, we can see that the geometric series is much more than a tool to calculate compound interest in high school, and can in fact be generalized to other mathematical objects. In the next part, we will explore how it is used in the study of fixed-point iterations, which is ubiquitous in algorithm development and convergence analysis. We will then prove the Banach fixed-point theorem, which is an extremely simple yet profound result.  

 
<br> 

Further readings and references
-----
1. [Jordan decomposition theorem and Cayley-Hamilton theorem via resolvent analysis](https://www.maths.usyd.edu.au/u/daners/publ/abstracts/linalg/linalg.pdf)
2. [Operator norm](https://en.wikipedia.org/wiki/Operator_norm) and [induced matrix norm](https://en.wikipedia.org/wiki/Matrix_norm#Matrix_norms_induced_by_vector_norms)

_Last updated: 2025-08-22 12:56 EST_