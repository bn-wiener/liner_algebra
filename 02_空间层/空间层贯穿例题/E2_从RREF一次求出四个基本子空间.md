# E2 从 RREF 一次求出四个基本子空间

这个例子保留旧资料“四个基本子空间”的完整矩阵计算，目标是把列空间、行空间、零空间和左零空间放在同一个矩阵上观察。

考虑矩阵：

$$
A=
\begin{bmatrix}
1&2&-1&0\\
2&4&1&3\\
-1&-2&-2&-3
\end{bmatrix}.
$$

它是一个 $3\times4$ 矩阵，因此：

$$
A:\mathbb R^4\longrightarrow\mathbb R^3.
$$

行化简得到：

$$
\operatorname{rref}(A)
=
\begin{bmatrix}
1&2&0&1\\
0&0&1&1\\
0&0&0&0
\end{bmatrix}.
$$

主元位于第 $1$、第 $3$ 列，所以：

$$
r=2.
$$

### 列空间

从原矩阵取第 $1$、第 $3$ 列：

$$
\operatorname{Col}(A)
=
\operatorname{span}
\left\{
\begin{bmatrix}1\\2\\-1\end{bmatrix},
\begin{bmatrix}-1\\1\\-2\end{bmatrix}
\right\}.
$$

列空间位于 $\mathbb R^3$ 中，维数为 $2$。

### 行空间

最简行阶梯形中的非零行给出：

$$
\operatorname{Row}(A)
=
\operatorname{span}
\left\{
\begin{bmatrix}1&2&0&1\end{bmatrix},
\begin{bmatrix}0&0&1&1\end{bmatrix}
\right\}.
$$

行空间位于 $\mathbb R^4$ 中，维数为 $2$。

### 零空间

齐次方程为：

$$
\begin{cases}
x_1+2x_2+x_4=0,\\
x_3+x_4=0.
\end{cases}
$$

令：

$$
x_2=s,
\qquad
x_4=t.
$$

则：

$$
x_1=-2s-t,
\qquad
x_3=-t.
$$

所以：

$$
\mathcal N(A)
=
\operatorname{span}
\left\{
\begin{bmatrix}-2\\1\\0\\0\end{bmatrix},
\begin{bmatrix}-1\\0\\-1\\1\end{bmatrix}
\right\}.
$$

零空间位于 $\mathbb R^4$ 中，维数为：

$$
4-2=2.
$$

### 左零空间

观察原矩阵的三行：

$$
r_3=r_1-r_2.
$$

所以：

$$
r_1-r_2-r_3=0.
$$

因此：

$$
\mathcal N(A^{\mathsf T})
=
\operatorname{span}
\left\{
\begin{bmatrix}1\\-1\\-1\end{bmatrix}
\right\}.
$$

左零空间位于 $\mathbb R^3$ 中，维数为：

$$
3-2=1.
$$

### 正交性验证

取一个行空间基向量和一个零空间基向量：

$$
\begin{bmatrix}1&2&0&1\end{bmatrix}
\begin{bmatrix}-2\\1\\0\\0\end{bmatrix}
=0.
$$

对另一零空间基向量：

$$
\begin{bmatrix}1&2&0&1\end{bmatrix}
\begin{bmatrix}-1\\0\\-1\\1\end{bmatrix}
=0.
$$

同样，左零空间基向量与列空间基向量正交：

$$
\begin{bmatrix}1&-1&-1\end{bmatrix}
\begin{bmatrix}1\\2\\-1\end{bmatrix}
=0,
$$

$$
\begin{bmatrix}1&-1&-1\end{bmatrix}
\begin{bmatrix}-1\\1\\-2\end{bmatrix}
=0.
$$

### 维数核对

输入空间：

$$
\dim\operatorname{Row}(A)
+
\dim\mathcal N(A)
=2+2=4.
$$

输出空间：

$$
\dim\operatorname{Col}(A)
+
\dim\mathcal N(A^{\mathsf T})
=2+1=3.
$$

四个空间的维数恰好填满各自的环境空间。

## 用左零空间判断右端是否可达

对于上面的矩阵，左零空间基向量为：

$$
y=
\begin{bmatrix}1\\-1\\-1\end{bmatrix}.
$$

所以可达右端必须满足：

$$
y^{\mathsf T}b=0.
$$

展开得到：

$$
b_1-b_2-b_3=0.
$$

取：

$$
b_{\text{good}}
=
\begin{bmatrix}1\\2\\-1\end{bmatrix}.
$$

因为：

$$
1-2-(-1)=0,
$$

所以它位于列空间中，方程有解。

再取：

$$
b_{\text{bad}}
=
\begin{bmatrix}1\\2\\0\end{bmatrix}.
$$

因为：

$$
1-2-0=-1\neq0,
$$

所以它不在列空间中，方程无解。

## 怎样系统求四个基本子空间

给定矩阵 $A$，可以按以下顺序计算。

### 先化简矩阵

求出行阶梯形或最简行阶梯形，并记录：

- 主元位置；
- 自由变量；
- 秩。

### 求列空间

根据主元列位置，回到原矩阵选择对应列。

### 求行空间

选择阶梯形矩阵中的所有非零行。

### 求零空间

求解：

$$
Ax=0.
$$

提取自由参数对应的方向向量。

### 求左零空间

求解：

$$
A^{\mathsf T}y=0.
$$

### 用维数检查

若秩为 $r$，结果应满足：

$$
\dim\operatorname{Col}(A)=r,
$$

$$
\dim\operatorname{Row}(A)=r,
$$

$$
\dim\mathcal N(A)=n-r,
$$

$$
\dim\mathcal N(A^{\mathsf T})=m-r.
$$

## 复盘

特别检查：列空间基为什么回到原矩阵取主元列，而行空间基为什么可以直接使用阶梯形矩阵的非零行。
