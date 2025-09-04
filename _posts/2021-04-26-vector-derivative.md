---
title: 'Fréchet derivative and common formulas'
date: 2021-04-26  
permalink: /blurbs/2021/04/vec-derivative/
tags:
  - optimization
  - math
---

Some basic Fréchet derivative formulas one often encounters in optimization problems. 

_Acknowledgement: the content of this post originates from Appendix A.4 of Convex Optimization by Boyd & Vandenberghe, and the Introduction to Nonlinear Systems (ECE1647) lecture notes by Prof. Manfredi Maggiore at University of Toronto._


Definitions and properties
-----

In this brief note, we introduce the _Fréchet derivative_ for vector valued functions in finite-dimensional Banach spaces, and derive some formulas commonly seen for practical applications in optimization and control theory. We focus on functions 

$$
    f : X \goesto Y, 
$$

where \\(X \subset \mathbb{R}^n\\) and \\(Y \subset \mathbb{R}^m\\) are [domains](https://en.wikipedia.org/wiki/Domain_(mathematical_analysis)) (open and connected sets), and \\(n, m \geq 1\\) are integers. Note that we restrict the attention to such domains to have better control over boundary and local behaviors in the ensuing definitions.  
In particular, since most objective functions we will use are real-valued, i.e., \\(m = 1\\), most formulas in this note pertains to this case. 

To begin, we recall the basic derivative definition when \\(n = m = 1\\), which we would see in the first calculus course: a function \\(f : \mathbb{R} \goesto \mathbb{R}\\) is differentiable at \\(x \inR\\) if the limit 

$$
\begin{equation}\label{eq:diff-scalar}
    \lim_{h \goesto 0} \frac{f(x + h) - f(x)}{h}
\end{equation}
$$

exists. The limit is called the _derivative of \\(f\\) at \\(x\\)_, and denoted by \\(f'(x)\\). We are taught that the derivative acts like a linearization (or linear approximation of \\(f\\)) at the point \\(x\\), which provides many convenient tools for analysis and computation. Equivalently, we can write \eqref{eq:diff-scalar} as 

$$
\begin{equation}\label{eq:diff-scalar-1}
    \lim_{h \goesto 0} \frac{\vert f(x + h) - f(x) - f'(x)h \vert}{\vert h\vert} = 0.
\end{equation}
$$

To generalize this definition to vector- and matrix-valued functions, we use the _Fréchet derivative_. (Note: we will only discuss the case of _vector-_ valued functions, since matrix-valued functions can be equivalently understood via the linear, bijective _[vecorization](https://en.wikipedia.org/wiki/Vectorization_(mathematics))_ operation, even though it is only an algebraic "convenience" and does not reflect many deeper geometric meanings of matrix-valued functions). We say \\(f: X \goesto Y\\) is (Fréchet) differentiable at \\(x \in X\\) if there exists a linear operator \\(Df : \mathbb{R}^n \goesto \mathbb{R}^m\\), such that 

$$
\begin{equation}\label{eq:diff-frechet}
    \lim_{h \goesto \zero} \frac{\Vert f(x + h) - f(x) - Df(x)h \Vert}{\Vert h\Vert} = \zero. 
\end{equation}
$$

The Fréchet derivative of \\(f\\) at a point \\(x \in X\\) is denoted by \\Df(x)\\, and is unique if it exists. We can observe the similarity between \eqref{eq:diff-frechet} and \eqref{eq:diff-scalar-1}. Note that the norms in \eqref{eq:diff-frechet} are understood to be \\(L_p\\) norms consistent with the dimensions \\(n,m\\). If \\(f\\) is differentiable at _every_ \\(x \in X\\), then we say \\(f\\) is _differentiable_. Furthremore, if \\(f\\) is differentiable, and the map \\(x \mapsto Df\\) from \\(X\\) to \\(\mathbb{R}^{m \times n}\\) is continuous, then we say \\(f\\) is _continuously differentiable_, and write \\(f\\) is \\(\mathcal{C}^1\\). Under the standard bases of \\(\mathbb{R}^n\\) and \\(\mathbb{R}^m\\), the matrix representation of \\(Df\\) is the _Jacobian matrix_. Let \\(\\{\e\_1, \dots, \e\_n\\}\\) be the set of standard basis vectors in \\(\mathbb{R}^n\\), we define 

$$
    \frac{\partial f_{i}}{\partial x_j}(x) = \lim_{t \goesto 0} \frac{f(x + t \e_j) - f(x)}{t}, 
$$

as the \\(j\\)-th _partial derivative_ of \\(f\\) at \\(x\\), if the limit exists. As a simple example, when \\(f : \mathbb{R}^n \goesto \mathbb{R}\\), we get 

$$
    Df(x) = \begin{bmatrix}
        
        \dfrac{\partial f}{\partial x_1} & \dfrac{\partial f}{\partial x_2} \dots \dfrac{\partial f}{\partial x_n}, 
    \end{bmatrix}
$$

and its transpose is the (more) commonly seen _gradient_, denoted by \\(\nabla f(x)\\). 

The classic chain rule can also be extended to vector-valued functions: let \\(X,Y\\) and \\(Z\\) be domains of \\(\mathbb{R}^n\\), \\(\mathbb{R}^m\\), and \\(\mathbb{R}^p\\), respectively, and let \\(f : Y \goesto Z\\), \\(g : X \goesto Y\\) be differentiable. Then, the composite function \\(f \odot g : X \goesto Z\\) is differentiable and 

$$
    D(f \odot g)(x) = Df(g(x))Dg(x). 
$$

We can further define the second Fréchet derivative, \\(D(Df)\\), and even higher order derivatives, provided \\(f\\) is sufficiently smooth. To be rigorous, we would need to properly define the spaces (bounded multilinear maps) that we work with. However, these definitions get quite messy with the introduction of required advanced tools, which are often unnecessary for practical purposes (see Reference 1 for a better treatment). Instead, we focus on the case of real-valued functions \\(f : \mathbb{R}^n \goesto \mathbb{R}), whose second derivative, under the standard basis, is commonly given by the _Hessian matrix_: 

$$
    \nabla^2 f(x) = D(\nabla f)(x) = \begin{bmatrix}
        \dfrac{\partial^2 f(x)}{\partial x_1\partial x_1} & \cdots & \dfrac{\partial f^2(x)}{\partial x_1 \partial x_n} \\
        \vdots & \ddots & \vdots \\
        \dfrac{\partial^2 f(x)}{\partial x_n \partial x_1} & \cdots & \dfrac{\partial^2 f(x)}{\partial x_n \partial x_n}
\end{bmatrix}.
$$

If \\(f\\) is _twice continuously differentiable_, i.e., \\(f\\) is \\(\mathcal{C}^2\\), then its Hessian is symmetric. 


Some common formulas
-----

Next, we go through some (simple but tedious) example derivative calculations. When I took my first course in optimization, I found these steps helpful to convince myself of the results, but again, they are quite tedious. 


1. Take \\(f(x) = \transpose{a} x\\) for some known vector \\(a \inR^{n}). Then, \\(Df(x) = \transpose{a}\\) and \\(\nabla f(x) = a\\).  
    _Proof_. 

    $$
        D \left(\transpose{a}x\right) = 
        \begin{bmatrix}
            \dfrac{\partial}{\partial x_1}\sum_{i=1}^n a_ix_i &
            \cdots &
            \dfrac{\partial}{\partial x_n}\sum_{i=1}^n a_ix_i
        \end{bmatrix} = 
        \begin{bmatrix}
            a_1 & \cdots & a_n
        \end{bmatrix}
        = \transpose{a}.
    $$

2. Take \\(f(x) = Ax\\) for some known matrix \\(A \inR^{m \times n}\\). Then, \\(Df(x) = A\\). 

    _Proof_. 

    $$
    D (Ax) = D \left(\begin{bmatrix}
        a_{11} & \cdots & a_{1n} \\ 
        \vdots & \ddots & \vdots \\
        a_{n1} & \cdots & a_{nn}
    \end{bmatrix}\begin{bmatrix}
        x_1 \\ \vdots \\ x_n
    \end{bmatrix} \right)
    = \begin{bmatrix}
        D\left( \sum_{i=1}^n a_{1i}x_i \right)\\
        \vdots \\
        D\left(\sum_{i=1}^n a_{ni}x_i \right)
    \end{bmatrix}.
    $$

    From the previous identity, \\(D \\left(\\sum\_{i=1}^n a\_{ji}x\_i\\right)\\) for \\(j \in \\{1,...,n\\}\\) results in the row vector 
    
    $$
    \begin{bmatrix}
        a\_{j1} & \\cdots & a\_{jn}
    \end{bmatrix},
    $$

    which is precisely the \\(j\\)-th row of the matrix \\(A\\). Therefore, \\(D(Ax) = A\\).

3. For a matrix \\(A \inS^n\\), i.e., \\(A = \transpose{A}\\), we have that \\(D \left(\transpose{x}Ax\right) = 2Ax\\). 

    _Proof_. 
    If we treat \\(f(x) = \transpose{x}Ax\\), evidently \\(f\\) is real valued, so we want to find 
    
    $$
    \begin{bmatrix}\dfrac{\partial f(x)}{\partial x_1} & \cdots & \dfrac{\partial f(x)}{\partial x_n}\end{bmatrix}. 
    $$

    Expanding \\(f\\), for \\(k \in \\{1,...,n\\}\\), we have

    $$
    \begin{align*}
        \transpose{x}A{x} &= \sum_{i=1}^n x_i \sum_{j=1}^n a_{ij}x_j \\
        &= x_k\sum_{j=1}^na_{kj}x_j + \sum_{i\neq k} x_i\sum_{j=1}^n a_{ij}x_j \\
        &= x_k\left(a_{kk}x_k + \sum_{j\neq k} a_{kj}x_j\right) + \sum_{i\neq k} \left( x_ia_{ik}x_k + x_i \sum_{j\neq k} a_{ij}x_j\right) \\
        &= a_{kk}x_k^2 + \left(\sum_{j\neq k} a_{kj}x_j\right)x_k + \left(\sum_{i\neq k} a_{ik}x_i\right)x_k + \sum_{i\neq k}\sum_{j\neq k} a_{ij}x_jx_i.
    \end{align*}
    $$

    Take the derivative of each term above with respect to \\(x\_k\\), and we get 

    $$
    \begin{align*}
        \frac{\partial}{\partial x_k} \transpose{x}Ax = 2a_{kk}x_k + \sum_{j\neq k} a_{kj}x_j + \sum_{i\neq k} a_{ik}x_i 
        = 2a_{kk}x_k + 2\sum_{j\neq k} a_{kj}x_j = 2\sum_{j=1}^n a_{kj}x_j,
    \end{align*}
    $$

    which is twice the inner product between the \\(k\\)-th row of \\(A\\) and \\(x\\). We can combine the two sums in the second step above since \\(A = \transpose{A}\\), so \\(a_{kj} = a_{ik}\\).  Therefore, the full derivative \\(D\left(\transpose{x}Ax\right)\\) is \\(2Ax\\). Note that, typically we consider the function \\(\frac{1}{2}\transpose{x}Ax\\) instead, where \\(A \inS^n\\) cancels out the constant 2 in the result, since optimizers are not impacted by such constant factors. 
 
4. It follows that \\(\nabla^2 \left( \transpose{x}Ax\right) = 2A\\).


<br> 

Further readings and references
-----
1. Fréchet derivative in [general Banach spaces](https://desvl.xyz/2020/07/31/frechet-derivative/)

_Last updated: 2025-09-03 15:32 EST_