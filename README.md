# Problem Statement

$ T \in M_n(R) $

$ T = \begin{bmatrix}
t1 & t2 & \dots & tn
\\ t(n+1) & t(n+2) & \dots & t2n 
\\ . & . & \dots & . 
\\ . & . &  & . 
\\ . & . & \dots & . 
\\ . & . & \dots & tnn 
\end{bmatrix}$

$ \text{No. of variables in T} = n^2 $

$ y \in \{ \begin{bmatrix}y_1, y_2, \dots y_{n-1}, 1\end{bmatrix} : y_i \in \{-1, 1\} \} $

$ x \in \{ \begin{bmatrix}x_1, x_2, \dots x_n\end{bmatrix}^T : x_i \in \{-1, 1\} \} $

$ b1 =\begin{bmatrix}1, 1, 1, \dots , (\text{n times})\end{bmatrix}^T $

$t = \begin{bmatrix}
t1
\\ t2
\\ .
\\ .
\\ .
\\ tnn
\end{bmatrix}$

$$
yTx = 1, \quad \forall x, \forall y \tag{1}
$$


$ \text{No. of equations satisfying (1)} = 2^{2n-1} . $

$ \text{Choose } n^2  \text{ equations from above equations.} $

$ \text{No. of linear systems of equations} = \binom{2^{2n-1}}{n^2} . $

$ \text{Consider } $
$$ At = b1.  \tag{2}$$

$ \text{A is the coefficient matrix of the chosen } n^2 \text{ equations.} $



$ \large{\text{1. Find out the number systems with unique solutions. i.e. Find how many A's are there such that } A^{-1} \text{  exists.}} $
$ \large{\text{2. Find the unique solutions. i.e. Find all such corresponding T.}} $



## Example:

$ \text{no. of variables (m) = 4} $

$ n = 2 $

$T = \begin{bmatrix}
t1 & t2
\\ t3 & t4
\end{bmatrix}$

$t = \begin{bmatrix}
t1
\\ t2
\\ t3 
\\ t4
\end{bmatrix}$


$ x \in \left\{ [-1, -1]^T, [-1, 1]^T, [1, -1]^T, [1, 1]^T \right\} $

$ y \in \left\{ [-1, 1], [1, 1] \right\} $

$ yTx = 1, \forall x, \forall y $

$ \text{We get } 2^3 = 8 \text{ equations in 4 variables.} $

$ \text{No. of linear systems} = \binom{8}{4} = 70. $

$ \text{No. of linear systems with unique solution = 16 (after calculations).} $


