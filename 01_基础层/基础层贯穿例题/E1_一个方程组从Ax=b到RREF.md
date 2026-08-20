# E1 一个方程组：从 $A\boldsymbol{x}=\boldsymbol{b}$ 到 RREF

本例把 B4～B7 连起来。重点不是记某次行变换，而是观察每一步解决什么问题。

回到开始时的方程组，其增广矩阵为：

$$
\left[
\begin{array}{ccc|c}
1&2&-1&3\\
2&5&1&8\\
-1&-2&2&-1
\end{array}
\right].
$$

首先用第 $1$ 行消去第 $1$ 列下方的元素：

$$
R_2\leftarrow R_2-2R_1,
\qquad
R_3\leftarrow R_3+R_1.
$$

得到：

$$
\left[
\begin{array}{ccc|c}
1&2&-1&3\\
0&1&3&2\\
0&0&1&2
\end{array}
\right].
$$

现在矩阵已经是行阶梯形。它对应：

$$
\begin{cases}
x+2y-z=3,\\
y+3z=2,\\
z=2.
\end{cases}
$$

从最后一行向上回代：

$$
z=2,
$$

$$
y+3\times2=2
\quad\Longrightarrow\quad
y=-4,
$$

$$
x+2(-4)-2=3
\quad\Longrightarrow\quad
x=13.
$$

所以：

$$
\boldsymbol{x}
=
\begin{bmatrix}
13\\-4\\2
\end{bmatrix}.
$$

这个例子中三个变量列都有主元，没有自由变量，也没有矛盾行，所以方程组有唯一解。

如果继续向上消元，可以得到最简行阶梯形：

$$
\left[
\begin{array}{ccc|c}
1&0&0&13\\
0&1&0&-4\\
0&0&1&2
\end{array}
\right],
$$

此时可以直接读出解，不再需要回代。

## 复盘

重新回答：$A$、$\boldsymbol{x}$、$\boldsymbol{b}$ 分别是什么？主元在哪些列？为什么没有自由变量？为什么最终是唯一解？
