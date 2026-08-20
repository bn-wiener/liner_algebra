# T22 QR 分解

> **学习定位**：从第二层 Gram–Schmidt 继续前进：把矩阵的列结构写成“标准正交方向 $Q$ + 坐标/上三角结构 $R$”。

## 本节主线

对于满列秩矩阵，QR 分解写成：

$$
A=QR.
$$

$Q$ 记录标准正交方向，$R$ 记录原列向量在这些方向上的坐标。

> 本文对应《阶段 4 任务》的第六步。
>
> 这一节要解决的问题是：怎样把矩阵的普通列向量，改用一组标准正交方向表示，并利用这种表示更方便地求解方程和最小二乘问题？

## 1. 先用一句话理解 QR 分解

设矩阵：

$$
A=
\begin{bmatrix}
|&|&&|\\
a_1&a_2&\cdots&a_n\\
|&|&&|
\end{bmatrix}.
$$

第五步使用 Gram–Schmidt，把普通列向量：

$$
a_1,a_2,\ldots,a_n
$$

变成标准正交向量：

$$
q_1,q_2,\ldots,q_n.
$$

将新方向作为矩阵 $Q$ 的列：

$$
Q=
\begin{bmatrix}
|&|&&|\\
q_1&q_2&\cdots&q_n\\
|&|&&|
\end{bmatrix}.
$$

因为两组向量张成相同的列空间，所以每个原列向量 $a_j$ 都可以用 $q_1,\ldots,q_n$ 表示。

把这些表示系数整理成矩阵 $R$，就得到：

$$
\boxed{A=QR}.
$$

最核心的理解是：

```text
Q：保存新的标准正交方向
R：保存原列向量在这些新方向下的坐标
```

所以 QR 分解本质上是在为 $A$ 的列空间更换一组更方便的基。

### 1.1 不要先把 QR 看成两个矩阵相乘

如果只看公式：

$$
A=QR,
$$

很容易把 QR 分解理解成一种机械计算。更重要的问题其实是：

> $A$、$Q$、$R$ 分别把什么坐标变成什么对象？

设：

$$
A\in\mathbb{R}^{m\times n}.
$$

那么 $A$ 接收的是输入坐标：

$$
x\in\mathbb{R}^n,
$$

并产生输出向量：

$$
Ax\in\mathbb{R}^m.
$$

把 $A$ 写成列向量形式：

$$
A=
\begin{bmatrix}
|&|&&|\\
a_1&a_2&\cdots&a_n\\
|&|&&|
\end{bmatrix},
$$

则：

$$
Ax=x_1a_1+x_2a_2+\cdots+x_na_n.
$$

这说明 $x$ 中的数字是相对于原列向量 $a_1,\ldots,a_n$ 使用的组合系数。

经过 Gram–Schmidt 后，得到同一个列空间中的标准正交基：

$$
Q=
\begin{bmatrix}
|&|&&|\\
q_1&q_2&\cdots&q_n\\
|&|&&|
\end{bmatrix}.
$$

同一个输出向量 $Ax$，也一定可以写成：

$$
Ax=c_1q_1+c_2q_2+\cdots+c_nq_n=Qc.
$$

这里的 $c$ 不是新的输出向量，而是 $Ax$ 在标准正交基 $q_1,\ldots,q_n$ 下的坐标。

因此，QR 分解中真正发生的是：

```text
输入坐标 x（相对于 A 的原列）
        │
        │ R：换算成 Q 基下的坐标
        ▼
列空间坐标 c=Rx
        │
        │ Q：用标准正交基重建输出向量
        ▼
输出向量 Ax=Qc=QRx
```

所以：

$$
\boxed{R\text{ 负责坐标换算，}Q\text{ 负责用新基生成实际输出向量。}}
$$

### 1.2 为什么标准正交基能够直接读取坐标

设某个向量 $y$ 属于 $Q$ 的列空间，并且：

$$
y=c_1q_1+\cdots+c_nq_n=Qc.
$$

分别用 $q_i^T$ 与 $y$ 做内积：

$$
\begin{aligned}
q_i^Ty
&=q_i^T(c_1q_1+\cdots+c_nq_n)\\
&=c_1q_i^Tq_1+\cdots+c_nq_i^Tq_n.
\end{aligned}
$$

因为 $q_1,\ldots,q_n$ 标准正交：

$$
q_i^Tq_j=
\begin{cases}
1,&i=j,\\
0,&i\neq j,
\end{cases}
$$

上式只剩下：

$$
q_i^Ty=c_i.
$$

把所有坐标放在一起：

$$
\boxed{c=Q^Ty}.
$$

因此，对于位于列空间中的向量：

$$
\boxed{
y\xrightarrow{\ Q^T\ }c
\xrightarrow{\ Q\ }y
}.
$$

$Q^T$ 读取标准正交坐标，$Q$ 使用这些坐标重建向量。

### 1.3 从坐标变换直接推出 $R=Q^TA$

现在令：

$$
y=Ax.
$$

$y$ 在 $Q$ 基下的坐标为：

$$
c=Q^Ty=Q^TAx.
$$

另一方面，我们把“输入坐标 $x$ 变成 $Q$ 基坐标 $c$”的矩阵定义为 $R$：

$$
c=Rx.
$$

由于对任意 $x$ 都有：

$$
Rx=Q^TAx,
$$

所以两个线性变换的矩阵必须相同：

$$
\boxed{R=Q^TA}.
$$

这个等式的含义不是抽象的矩阵技巧，而是：

> $R$ 的第 $j$ 列，正是原列向量 $a_j$ 在标准正交基 $Q$ 下的坐标。

因为 $R=Q^TA$，它的第 $i$ 行、第 $j$ 列为：

$$
\boxed{r_{ij}=q_i^Ta_j}.
$$

### 1.4 再从 $R=Q^TA$ 推出 $A=QR$

由：

$$
R=Q^TA
$$

可得：

$$
QR=QQ^TA.
$$

这里 $QQ^T$ 是投影到 $Q$ 的列空间的矩阵，而：

$$
\operatorname{Col}(Q)=\operatorname{Col}(A).
$$

$A$ 的每一列本来就已经位于这个空间中，因此投影不会改变它：

$$
QQ^TA=A.
$$

于是：

$$
\boxed{A=QR}.
$$

这条推导可以压缩为：

$$
\boxed{
R=Q^TA
\quad\Longrightarrow\quad
QR=QQ^TA=A
}.
$$

注意，只有当向量位于 $\operatorname{Col}(Q)$ 中时，$QQ^T$ 才保持它不变。对于一般的 $b\in\mathbb{R}^m$：

$$
QQ^Tb
$$

只是 $b$ 在列空间中的正交投影，不一定等于 $b$。

### 1.5 严格推导 $R$ 为什么是上三角矩阵

“第 $j$ 列只使用前 $j$ 个 $q$”是正确结论，但还需要说明为什么。

Gram–Schmidt 具有下面的逐步空间关系：

$$
\operatorname{span}\{q_1,\ldots,q_j\}
=
\operatorname{span}\{a_1,\ldots,a_j\}.
$$

当 $i>j$ 时，$q_i$ 是处理 $a_i$ 时得到的新方向。它与此前已经得到的整个空间正交：

$$
q_i\perp
\operatorname{span}\{q_1,\ldots,q_{i-1}\}.
$$

因为 $j<i$，所以：

$$
a_j\in
\operatorname{span}\{a_1,\ldots,a_j\}
\subseteq
\operatorname{span}\{a_1,\ldots,a_{i-1}\}.
$$

而前 $i-1$ 个 $a$ 与前 $i-1$ 个 $q$ 张成相同空间，因此：

$$
a_j\in
\operatorname{span}\{q_1,\ldots,q_{i-1}\}.
$$

于是 $q_i$ 与 $a_j$ 正交：

$$
q_i^Ta_j=0.
$$

再使用：

$$
r_{ij}=q_i^Ta_j,
$$

得到：

$$
\boxed{i>j\quad\Longrightarrow\quad r_{ij}=0}.
$$

$i>j$ 恰好对应主对角线下方的位置，所以 $R$ 必然是上三角矩阵。

这不是为了计算方便而人为规定的形状，而是 Gram–Schmidt 按列依次建立新方向的自然结果。

### 1.6 对角线元素为什么是正交化剩余量的长度

第 $j$ 步先从 $a_j$ 中减去它在旧方向上的投影：

$$
u_j
=
a_j-
\sum_{i=1}^{j-1}(q_i^Ta_j)q_i.
$$

然后单位化：

$$
q_j=\frac{u_j}{\|u_j\|}.
$$

因此：

$$
u_j=\|u_j\|q_j.
$$

把投影部分移回去：

$$
a_j
=
\sum_{i=1}^{j-1}(q_i^Ta_j)q_i
+\|u_j\|q_j.
$$

与：

$$
a_j=r_{1j}q_1+\cdots+r_{jj}q_j
$$

逐项比较，得到：

$$
r_{ij}=q_i^Ta_j\quad(i<j),
$$

以及：

$$
\boxed{r_{jj}=\|u_j\|}.
$$

若 $A$ 满列秩，则每一步都会产生非零的新方向：

$$
u_j\neq0,
$$

所以：

$$
r_{jj}>0.
$$

如果某个 $r_{jj}=0$，就表示 $a_j$ 完全落在前面列向量张成的空间中，没有增加新方向。

---

## 2. $Q$ 和 $R$ 分别记录什么

### 2.1 $Q$ 记录方向

$Q$ 的列是标准正交向量：

$$
q_1,q_2,\ldots,q_n.
$$

因此：

$$
\boxed{Q^TQ=I_n}.
$$

而且 Gram–Schmidt 保持列空间不变：

$$
\boxed{
\operatorname{Col}(Q)=\operatorname{Col}(A)
}.
$$

### 2.2 $R$ 记录坐标

第 $j$ 个原列向量可以写成：

$$
a_j
=
r_{1j}q_1+
r_{2j}q_2+
\cdots+
r_{nj}q_n.
$$

这些系数构成 $R$ 的第 $j$ 列：

$$
\begin{bmatrix}
r_{1j}\\
r_{2j}\\
\vdots\\
r_{nj}
\end{bmatrix}.
$$

所以矩阵乘法：

$$
QR
$$

就是使用 $R$ 中的坐标，对 $Q$ 的列进行线性组合，从而重建 $A$ 的所有列。

---

## 3. 从 Gram–Schmidt 直接看出 $A=QR$

设 $A$ 有三个线性无关的列：

$$
A=
\begin{bmatrix}
|&|&|\\
a_1&a_2&a_3\\
|&|&|
\end{bmatrix}.
$$

我们先不直接写 $R$，而是完整执行 Gram–Schmidt，并观察每一步自然产生了哪些系数。

### 3.1 处理第一列：得到 $r_{11}$

第一列前面没有旧方向需要减去，所以：

$$
u_1=a_1.
$$

将它单位化：

$$
q_1=\frac{u_1}{\|u_1\|}
=
\frac{a_1}{\|a_1\|}.
$$

把单位化时除掉的长度记作：

$$
\boxed{r_{11}=\|u_1\|=\|a_1\|}.
$$

因为 $q_1=u_1/r_{11}$，所以反过来：

$$
u_1=r_{11}q_1.
$$

又因为 $u_1=a_1$，得到：

$$
\boxed{a_1=r_{11}q_1}.
$$

因此，$R$ 的第一列已经确定：

$$
\begin{bmatrix}
r_{11}\\
0\\
0
\end{bmatrix}.
$$

下面两个零表示 $a_1$ 不需要使用后来才产生的 $q_2$ 和 $q_3$。

### 3.2 处理第二列：得到 $r_{12}$ 和 $r_{22}$

第二列 $a_2$ 通常含有一部分 $q_1$ 方向。它在 $q_1$ 方向上的坐标是：

$$
\boxed{r_{12}=q_1^Ta_2}.
$$

因此，$a_2$ 在 $q_1$ 方向上的投影为：

$$
r_{12}q_1.
$$

从 $a_2$ 中减去这部分：

$$
u_2=a_2-r_{12}q_1.
$$

$u_2$ 是 $a_2$ 中不能被旧方向 $q_1$ 表示的剩余部分。它将产生第二个新方向。

把 $u_2$ 的长度记作：

$$
\boxed{r_{22}=\|u_2\|}.
$$

然后单位化：

$$
q_2=\frac{u_2}{r_{22}}.
$$

反过来便有：

$$
u_2=r_{22}q_2.
$$

代回 $u_2=a_2-r_{12}q_1$：

$$
r_{22}q_2=a_2-r_{12}q_1.
$$

移项得到：

$$
\boxed{a_2=r_{12}q_1+r_{22}q_2}.
$$

因此，$R$ 的第二列为：

$$
\begin{bmatrix}
r_{12}\\
r_{22}\\
0
\end{bmatrix}
=
\begin{bmatrix}
q_1^Ta_2\\
\|a_2-(q_1^Ta_2)q_1\|\\
0
\end{bmatrix}.
$$

### 3.3 处理第三列：得到 $r_{13}$、$r_{23}$ 和 $r_{33}$

第三列 $a_3$ 可能同时含有 $q_1$ 和 $q_2$ 两个旧方向。

先读取它在 $q_1$ 方向上的坐标：

$$
\boxed{r_{13}=q_1^Ta_3}.
$$

再读取它在 $q_2$ 方向上的坐标：

$$
\boxed{r_{23}=q_2^Ta_3}.
$$

所以 $a_3$ 在旧空间中的投影是：

$$
r_{13}q_1+r_{23}q_2.
$$

减去这个投影：

$$
u_3
=
a_3-r_{13}q_1-r_{23}q_2.
$$

把剩余部分的长度记作：

$$
\boxed{r_{33}=\|u_3\|}.
$$

单位化得到：

$$
q_3=\frac{u_3}{r_{33}}.
$$

因此：

$$
u_3=r_{33}q_3.
$$

代回 $u_3$ 的定义并移项：

$$
\boxed{
a_3
=
r_{13}q_1+r_{23}q_2+r_{33}q_3
}.
$$

所以，$R$ 的第三列为：

$$
\begin{bmatrix}
r_{13}\\
r_{23}\\
r_{33}
\end{bmatrix}
=
\begin{bmatrix}
q_1^Ta_3\\
q_2^Ta_3\\
\|a_3-(q_1^Ta_3)q_1-(q_2^Ta_3)q_2\|
\end{bmatrix}.
$$

### 3.4 把三列结果放在一起

前三步分别得到：

$$
a_1=r_{11}q_1,
$$

$$
a_2=r_{12}q_1+r_{22}q_2,
$$

$$
a_3=r_{13}q_1+r_{23}q_2+r_{33}q_3.
$$

将它们并排写成矩阵：

$$
\begin{bmatrix}
|&|&|\\
a_1&a_2&a_3\\
|&|&|
\end{bmatrix}
=
\begin{bmatrix}
|&|&|\\
q_1&q_2&q_3\\
|&|&|
\end{bmatrix}
\begin{bmatrix}
r_{11}&r_{12}&r_{13}\\
0&r_{22}&r_{23}\\
0&0&r_{33}
\end{bmatrix}.
$$

也就是：

$$
\boxed{A=QR}.
$$

其中每个元素都能追溯到 Gram–Schmidt 的某一步：

$$
\boxed{
R=
\begin{bmatrix}
\|a_1\|&q_1^Ta_2&q_1^Ta_3\\
0&\|u_2\|&q_2^Ta_3\\
0&0&\|u_3\|
\end{bmatrix}
}.
$$

由于：

$$
q_1^Ta_1=\|a_1\|,
$$

$$
q_2^Ta_2=\|u_2\|,
$$

$$
q_3^Ta_3=\|u_3\|,
$$

还可以统一写成：

$$
\boxed{
R=
\begin{bmatrix}
q_1^Ta_1&q_1^Ta_2&q_1^Ta_3\\
0&q_2^Ta_2&q_2^Ta_3\\
0&0&q_3^Ta_3
\end{bmatrix}
}.
$$

以一般的第 $j$ 步为例，这个对角线等式可以直接证明。因为：

$$
a_j
=
\sum_{i=1}^{j-1}r_{ij}q_i+u_j,
$$

所以：

$$
\begin{aligned}
q_j^Ta_j
&=
q_j^T
\left(
\sum_{i=1}^{j-1}r_{ij}q_i+u_j
\right)\\
&=
\sum_{i=1}^{j-1}r_{ij}q_j^Tq_i
+q_j^Tu_j.
\end{aligned}
$$

$q_j$ 与所有旧方向 $q_i$ 正交，因此第一项为零。又因为：

$$
u_j=\|u_j\|q_j,
$$

所以：

$$
q_j^Tu_j
=
\|u_j\|q_j^Tq_j
=
\|u_j\|.
$$

最终得到：

$$
\boxed{q_j^Ta_j=\|u_j\|=r_{jj}}.
$$

### 3.5 推广到第 $j$ 列

在处理第 $j$ 列时，先计算它在所有旧方向上的坐标：

$$
\boxed{
r_{ij}=q_i^Ta_j,
\qquad i=1,2,\ldots,j-1
}.
$$

减去所有旧方向的投影：

$$
u_j
=
a_j-
\sum_{i=1}^{j-1}r_{ij}q_i.
$$

再定义：

$$
\boxed{r_{jj}=\|u_j\|},
$$

以及：

$$
q_j=\frac{u_j}{r_{jj}}.
$$

因为 $u_j=r_{jj}q_j$，所以：

$$
a_j
=
\sum_{i=1}^{j-1}r_{ij}q_i+r_{jj}q_j
=
\sum_{i=1}^{j}r_{ij}q_i.
$$

因此，$R$ 的第 $j$ 列是：

$$
\boxed{
\begin{bmatrix}
r_{1j}\\
r_{2j}\\
\vdots\\
r_{jj}\\
0\\
\vdots\\
0
\end{bmatrix}
}.
$$

完整的逐项规则是：

$$
\boxed{
r_{ij}
=
\begin{cases}
q_i^Ta_j,&i<j,\\
\|u_j\|,&i=j,\\
0,&i>j.
\end{cases}
}
$$

由于 $q_j^Ta_j=\|u_j\|$，前两种情况还可以合并：

$$
\boxed{
r_{ij}
=
\begin{cases}
q_i^Ta_j,&i\leq j,\\
0,&i>j.
\end{cases}
}
$$

### 3.6 为什么这些 $r_{ij}$ 正好组成 $Q^TA$

矩阵乘法 $Q^TA$ 的第 $i$ 行、第 $j$ 列就是：

$$
(Q^TA)_{ij}=q_i^Ta_j.
$$

所以：

$$
Q^TA
=
\begin{bmatrix}
q_1^Ta_1&q_1^Ta_2&\cdots&q_1^Ta_n\\
q_2^Ta_1&q_2^Ta_2&\cdots&q_2^Ta_n\\
\vdots&\vdots&\ddots&\vdots\\
q_n^Ta_1&q_n^Ta_2&\cdots&q_n^Ta_n
\end{bmatrix}.
$$

当 $i>j$ 时，$q_i$ 与较早的 $a_j$ 正交，所以：

$$
q_i^Ta_j=0.
$$

因此：

$$
Q^TA
=
\begin{bmatrix}
r_{11}&r_{12}&\cdots&r_{1n}\\
0&r_{22}&\cdots&r_{2n}\\
\vdots&\ddots&\ddots&\vdots\\
0&\cdots&0&r_{nn}
\end{bmatrix}
=R.
$$

最终得到完整关系：

$$
\boxed{R=Q^TA},
$$

$$
\boxed{A=QR}.
$$

这里每个 $r_{ij}$ 都不是事后填入的数字，而是在 Gram–Schmidt 处理第 $j$ 个列向量时自然产生的投影坐标或新方向长度。

### 3.7 三列数值例子：实际算出 $R$ 的六个非零元素

取：

$$
A=
\begin{bmatrix}
1&1&1\\
1&0&1\\
0&1&1
\end{bmatrix}
=
\begin{bmatrix}
|&|&|\\
a_1&a_2&a_3\\
|&|&|
\end{bmatrix},
$$

其中：

$$
a_1=
\begin{bmatrix}
1\\
1\\
0
\end{bmatrix},
\qquad
a_2=
\begin{bmatrix}
1\\
0\\
1
\end{bmatrix},
\qquad
a_3=
\begin{bmatrix}
1\\
1\\
1
\end{bmatrix}.
$$

第一列产生：

$$
r_{11}=\|a_1\|=\sqrt{2},
$$

$$
q_1=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1\\
1\\
0
\end{bmatrix}.
$$

处理第二列，先计算：

$$
r_{12}=q_1^Ta_2=\frac{1}{\sqrt{2}}.
$$

减去 $q_1$ 方向：

$$
\begin{aligned}
u_2
&=a_2-r_{12}q_1\\
&=
\begin{bmatrix}
1\\
0\\
1
\end{bmatrix}
-
\begin{bmatrix}
1/2\\
1/2\\
0
\end{bmatrix}\\
&=
\begin{bmatrix}
1/2\\
-1/2\\
1
\end{bmatrix}.
\end{aligned}
$$

因此：

$$
r_{22}=\|u_2\|=\frac{\sqrt{6}}{2},
$$

$$
q_2=
\frac{1}{\sqrt{6}}
\begin{bmatrix}
1\\
-1\\
2
\end{bmatrix}.
$$

处理第三列，先读取它在两个旧方向上的坐标：

$$
r_{13}=q_1^Ta_3
=
\frac{1+1}{\sqrt{2}}
=\sqrt{2},
$$

$$
r_{23}=q_2^Ta_3
=
\frac{1-1+2}{\sqrt{6}}
=
\frac{2}{\sqrt{6}}.
$$

减去两个旧方向：

$$
\begin{aligned}
u_3
&=a_3-r_{13}q_1-r_{23}q_2\\
&=
\begin{bmatrix}
1\\
1\\
1
\end{bmatrix}
-
\begin{bmatrix}
1\\
1\\
0
\end{bmatrix}
-
\begin{bmatrix}
1/3\\
-1/3\\
2/3
\end{bmatrix}\\
&=
\begin{bmatrix}
-1/3\\
1/3\\
1/3
\end{bmatrix}.
\end{aligned}
$$

所以：

$$
r_{33}=\|u_3\|=\frac{1}{\sqrt{3}},
$$

$$
q_3=
\frac{1}{\sqrt{3}}
\begin{bmatrix}
-1\\
1\\
1
\end{bmatrix}.
$$

现在 $R$ 中的六个非零元素已经全部得到：

$$
\begin{aligned}
r_{11}&=\sqrt{2},&
r_{12}&=\frac{1}{\sqrt{2}},&
r_{13}&=\sqrt{2},\\
r_{22}&=\frac{\sqrt{6}}{2},&
r_{23}&=\frac{2}{\sqrt{6}},&
r_{33}&=\frac{1}{\sqrt{3}}.
\end{aligned}
$$

因此：

$$
Q=
\begin{bmatrix}
1/\sqrt{2}&1/\sqrt{6}&-1/\sqrt{3}\\
1/\sqrt{2}&-1/\sqrt{6}&1/\sqrt{3}\\
0&2/\sqrt{6}&1/\sqrt{3}
\end{bmatrix},
$$

$$
R=
\begin{bmatrix}
\sqrt{2}&1/\sqrt{2}&\sqrt{2}\\
0&\sqrt{6}/2&2/\sqrt{6}\\
0&0&1/\sqrt{3}
\end{bmatrix}.
$$

逐列观察 $QR$：

$$
\begin{aligned}
\text{第一列：}&\quad r_{11}q_1=a_1,\\
\text{第二列：}&\quad r_{12}q_1+r_{22}q_2=a_2,\\
\text{第三列：}&\quad r_{13}q_1+r_{23}q_2+r_{33}q_3=a_3.
\end{aligned}
$$

所以：

$$
QR=
\begin{bmatrix}
|&|&|\\
a_1&a_2&a_3\\
|&|&|
\end{bmatrix}
=A.
$$

---

## 4. 为什么 $R$ 是上三角矩阵

Gram–Schmidt 处理第 $j$ 个向量 $a_j$ 时，只使用已经产生的方向：

$$
q_1,q_2,\ldots,q_j.
$$

因此：

$$
a_j\in\operatorname{span}\{q_1,q_2,\ldots,q_j\}.
$$

它不需要后面才会出现的方向：

$$
q_{j+1},q_{j+2},\ldots,q_n.
$$

所以当 $i>j$ 时：

$$
r_{ij}=0.
$$

在矩阵 $R$ 中，$i>j$ 对应主对角线下方的位置。这些位置全部为零，因此 $R$ 是上三角矩阵：

$$
\boxed{
R=
\begin{bmatrix}
r_{11}&r_{12}&\cdots&r_{1n}\\
0&r_{22}&\cdots&r_{2n}\\
\vdots&\ddots&\ddots&\vdots\\
0&\cdots&0&r_{nn}
\end{bmatrix}
}.
$$

可以简单记忆：

```text
a₁ 只能使用 q₁
a₂ 可以使用 q₁、q₂
a₃ 可以使用 q₁、q₂、q₃
……

每一列只能使用“到当前为止”已经得到的方向
所以对角线下方全部是 0
```

---

## 5. $R$ 中的元素怎样计算

因为 $Q$ 的列标准正交，$a_j$ 在方向 $q_i$ 上的坐标可以直接用内积得到：

$$
\boxed{r_{ij}=q_i^Ta_j}.
$$

对于 $i>j$：

$$
r_{ij}=0.
$$

对于对角线元素：

$$
r_{jj}=\|u_j\|,
$$

其中 $u_j$ 是 Gram–Schmidt 第 $j$ 步减去旧方向投影后剩余的向量。

在通常约定下：

$$
r_{jj}>0.
$$

把所有坐标一次性写成矩阵：

$$
\boxed{R=Q^TA}.
$$

为什么？

$Q^TA$ 的第 $i$ 行、第 $j$ 列正是：

$$
q_i^Ta_j=r_{ij}.
$$

---

## 6. $A$、$Q$、$R$ 的形状

设：

$$
A\in\mathbb{R}^{m\times n},
$$

并假设：

$$
m\geq n,
$$

而且 $A$ 的 $n$ 个列向量线性无关。

经济型 QR 分解中：

$$
Q\in\mathbb{R}^{m\times n},
$$

$$
R\in\mathbb{R}^{n\times n}.
$$

形状检查：

$$
\underbrace{Q}_{m\times n}
\underbrace{R}_{n\times n}
=
\underbrace{A}_{m\times n}.
$$

这里 $Q$ 不一定是方阵。

只要它的列标准正交，就有：

$$
Q^TQ=I_n.
$$

但当 $m>n$ 时：

$$
QQ^T\neq I_m.
$$

$QQ^T$ 是投影到 $\operatorname{Col}(A)$ 的矩阵。

---

## 7. 一个完整的二维 QR 分解例子

设：

$$
A=
\begin{bmatrix}
1&1\\
1&0
\end{bmatrix}.
$$

它的列为：

$$
a_1=
\begin{bmatrix}
1\\
1
\end{bmatrix},
\qquad
a_2=
\begin{bmatrix}
1\\
0
\end{bmatrix}.
$$

第五步已经对这两个向量做过 Gram–Schmidt。

得到：

$$
q_1=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1\\
1
\end{bmatrix},
\qquad
q_2=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1\\
-1
\end{bmatrix}.
$$

所以：

$$
Q=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1&1\\
1&-1
\end{bmatrix}.
$$

### 第一步：计算 $R=Q^TA$

$$
r_{11}=q_1^Ta_1=\sqrt{2},
$$

$$
r_{12}=q_1^Ta_2=\frac{1}{\sqrt{2}},
$$

$$
r_{21}=q_2^Ta_1=0,
$$

$$
r_{22}=q_2^Ta_2=\frac{1}{\sqrt{2}}.
$$

因此：

$$
R=
\begin{bmatrix}
\sqrt{2}&1/\sqrt{2}\\
0&1/\sqrt{2}
\end{bmatrix}.
$$

### 第二步：理解 $R$ 的两列

$R$ 的第一列表示：

$$
a_1=\sqrt{2}\,q_1+0q_2.
$$

$R$ 的第二列表示：

$$
a_2
=
\frac{1}{\sqrt{2}}q_1
+
\frac{1}{\sqrt{2}}q_2.
$$

### 第三步：验证 $A=QR$

$$
\begin{aligned}
QR
&=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1&1\\
1&-1
\end{bmatrix}
\begin{bmatrix}
\sqrt{2}&1/\sqrt{2}\\
0&1/\sqrt{2}
\end{bmatrix}\\
&=
\begin{bmatrix}
1&1\\
1&0
\end{bmatrix}\\
&=A.
\end{aligned}
$$

同时：

$$
Q^TQ=I_2.
$$

所以这个分解满足 QR 分解的全部要求。

### 第四步：跟踪一个具体输入经过 $R$ 和 $Q$ 的过程

取输入：

$$
x=
\begin{bmatrix}
2\\
-1
\end{bmatrix}.
$$

直接使用 $A$ 计算输出：

$$
Ax
=
\begin{bmatrix}
1&1\\
1&0
\end{bmatrix}
\begin{bmatrix}
2\\
-1
\end{bmatrix}
=
\begin{bmatrix}
1\\
2
\end{bmatrix}.
$$

现在改走 QR 分解的两步过程。

先使用 $R$：

$$
\begin{aligned}
c=Rx
&=
\begin{bmatrix}
\sqrt{2}&1/\sqrt{2}\\
0&1/\sqrt{2}
\end{bmatrix}
\begin{bmatrix}
2\\
-1
\end{bmatrix}\\
&=
\begin{bmatrix}
3/\sqrt{2}\\
-1/\sqrt{2}
\end{bmatrix}.
\end{aligned}
$$

这里的 $c$ 表示最终输出在 $q_1,q_2$ 方向上的坐标。因此：

$$
Ax
=
\frac{3}{\sqrt{2}}q_1
-
\frac{1}{\sqrt{2}}q_2.
$$

再使用 $Q$ 重建实际输出：

$$
\begin{aligned}
Qc
&=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1&1\\
1&-1
\end{bmatrix}
\begin{bmatrix}
3/\sqrt{2}\\
-1/\sqrt{2}
\end{bmatrix}\\
&=
\begin{bmatrix}
1\\
2
\end{bmatrix}\\
&=Ax.
\end{aligned}
$$

这个例子具体说明了：

```text
x 中的 2 和 -1
是使用原列 a₁、a₂ 的组合系数

Rx 中的 3/√2 和 -1/√2
是同一个输出在标准正交方向 q₁、q₂ 下的坐标

Q 再根据这些新坐标重建实际输出 Ax
```

因此 $R$ 不是“额外做了一次无关的变换”，而是在原列组合系数与标准正交基坐标之间进行换算。

---

## 8. 一个长方形矩阵的 QR 分解

设：

$$
A=
\begin{bmatrix}
1&1\\
1&0\\
0&1
\end{bmatrix}
\in\mathbb{R}^{3\times2}.
$$

它的列为：

$$
a_1=
\begin{bmatrix}
1\\
1\\
0
\end{bmatrix},
\qquad
a_2=
\begin{bmatrix}
1\\
0\\
1
\end{bmatrix}.
$$

### 第一步：计算 $q_1$

$$
\|a_1\|=\sqrt{2},
$$

$$
q_1
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1\\
1\\
0
\end{bmatrix}.
$$

### 第二步：处理 $a_2$

$$
q_1^Ta_2=\frac{1}{\sqrt{2}}.
$$

减去投影：

$$
\begin{aligned}
u_2
&=a_2-(q_1^Ta_2)q_1\\
&=
\begin{bmatrix}
1\\
0\\
1
\end{bmatrix}
-
\begin{bmatrix}
1/2\\
1/2\\
0
\end{bmatrix}\\
&=
\begin{bmatrix}
1/2\\
-1/2\\
1
\end{bmatrix}.
\end{aligned}
$$

它的长度为：

$$
\|u_2\|
=
\sqrt{\frac{1}{4}+\frac{1}{4}+1}
=
\frac{\sqrt{6}}{2}.
$$

因此：

$$
q_2
=
\frac{1}{\sqrt{6}}
\begin{bmatrix}
1\\
-1\\
2
\end{bmatrix}.
$$

### 第三步：组成 $Q$ 与 $R$

$$
Q=
\begin{bmatrix}
1/\sqrt{2}&1/\sqrt{6}\\
1/\sqrt{2}&-1/\sqrt{6}\\
0&2/\sqrt{6}
\end{bmatrix}
\in\mathbb{R}^{3\times2}.
$$

$$
R=
\begin{bmatrix}
\sqrt{2}&1/\sqrt{2}\\
0&\sqrt{6}/2
\end{bmatrix}
\in\mathbb{R}^{2\times2}.
$$

这里清楚地看到：

- $A$ 是 $3\times2$；
- 经济型 $Q$ 是 $3\times2$，不是方阵；
- $R$ 是 $2\times2$ 上三角矩阵；
- $QR$ 的形状是 $3\times2$，与 $A$ 相同。

---

## 9. 经济型 QR 与完整 QR

设：

$$
A\in\mathbb{R}^{m\times n},
\qquad
m>n,
$$

并且 $A$ 满列秩。

### 9.1 经济型 QR

经济型分解只保留描述 $\operatorname{Col}(A)$ 所需要的 $n$ 个标准正交方向：

$$
A=QR,
$$

其中：

$$
Q\in\mathbb{R}^{m\times n},
\qquad
R\in\mathbb{R}^{n\times n}.
$$

此时：

$$
Q^TQ=I_n.
$$

经济型 $Q$ 通常不是方阵。

### 9.2 完整 QR

把 $Q$ 的列继续补充为整个 $\mathbb{R}^m$ 的标准正交基，得到：

$$
Q_{\text{full}}
\in\mathbb{R}^{m\times m}.
$$

此时：

$$
Q_{\text{full}}^TQ_{\text{full}}
=
Q_{\text{full}}Q_{\text{full}}^T
=I_m.
$$

对应的：

$$
R_{\text{full}}
\in\mathbb{R}^{m\times n}
$$

可以写成：

$$
R_{\text{full}}
=
\begin{bmatrix}
R\\
0
\end{bmatrix}.
$$

于是：

$$
A=Q_{\text{full}}R_{\text{full}}.
$$

两者对比：

| 类型 | $Q$ 的形状 | $R$ 的形状 | 保存的方向 |
|---|---:|---:|---|
| 经济型 QR | $m\times n$ | $n\times n$ | 只保存列空间需要的方向 |
| 完整 QR | $m\times m$ | $m\times n$ | 补全为整个 $\mathbb{R}^m$ 的基 |

对于满列秩最小二乘，经济型 QR 通常已经足够。

---

## 10. 为什么 $Q$ 与 $A$ 的列空间相同

Gram–Schmidt 的每一步都只做两类操作：

1. 从原向量中减去已有方向的线性组合；
2. 对剩余向量进行单位化。

这些操作不会离开原来的列空间，所以：

$$
\operatorname{Col}(Q)
\subseteq
\operatorname{Col}(A).
$$

另一方面，由：

$$
A=QR
$$

可知 $A$ 的每一列都是 $Q$ 的列的线性组合，所以：

$$
\operatorname{Col}(A)
\subseteq
\operatorname{Col}(Q).
$$

因此：

$$
\boxed{
\operatorname{Col}(A)=\operatorname{Col}(Q)
}.
$$

它们描述的是同一个子空间，只是使用的基不同：

- $A$ 使用普通线性无关列；
- $Q$ 使用标准正交列。

---

## 11. 使用 QR 求解方程

先考虑方阵可逆的情况：

$$
A\in\mathbb{R}^{n\times n}.
$$

若：

$$
A=QR,
$$

则方程：

$$
Ax=b
$$

变成：

$$
QRx=b.
$$

左乘 $Q^T$：

$$
Q^TQRx=Q^Tb.
$$

因为：

$$
Q^TQ=I,
$$

所以：

$$
\boxed{Rx=Q^Tb}.
$$

$R$ 是上三角矩阵，因此可以使用回代求解。

假设：

$$
R=
\begin{bmatrix}
r_{11}&r_{12}&\cdots&r_{1n}\\
0&r_{22}&\cdots&r_{2n}\\
\vdots&\ddots&\ddots&\vdots\\
0&\cdots&0&r_{nn}
\end{bmatrix},
\qquad
y=Q^Tb.
$$

方程 $Rx=y$ 的最后一行只有一个未知量：

$$
r_{nn}x_n=y_n,
$$

所以：

$$
x_n=\frac{y_n}{r_{nn}}.
$$

得到 $x_n$ 后，再处理倒数第二行：

$$
r_{n-1,n-1}x_{n-1}+r_{n-1,n}x_n=y_{n-1}.
$$

一般地，从 $i=n$ 到 $i=1$ 依次计算：

$$
\boxed{
x_i
=
\frac{
y_i-\displaystyle\sum_{j=i+1}^{n}r_{ij}x_j
}{r_{ii}}
}.
$$

这就是回代。满列秩保证每个 $r_{ii}\neq0$，因此每一步都可以除以 $r_{ii}$。

原来的问题：

$$
Ax=b
$$

被分成两步：

```text
第一步：计算 y=Qᵀb
第二步：回代求解 Rx=y
```

整个过程中不需要计算 $A^{-1}$。

---

## 12. 使用 QR 求最小二乘解

现在设：

$$
A\in\mathbb{R}^{m\times n},
\qquad
m\geq n,
$$

并且 $A$ 满列秩。

经济型 QR 分解为：

$$
A=QR.
$$

第四步的最小二乘残差是：

$$
r=b-A\widehat{x}.
$$

最佳残差与 $\operatorname{Col}(A)$ 正交。

由于：

$$
\operatorname{Col}(A)=\operatorname{Col}(Q),
$$

所以：

$$
Q^Tr=0.
$$

代入残差：

$$
Q^T(b-QR\widehat{x})=0.
$$

展开：

$$
Q^Tb-Q^TQR\widehat{x}=0.
$$

因为：

$$
Q^TQ=I_n,
$$

所以：

$$
\boxed{
R\widehat{x}=Q^Tb
}.
$$

这就是使用 QR 求最小二乘解的核心公式。

它把一个 $m$ 维输出空间中的近似问题，转化成一个 $n\times n$ 上三角方程。

### 12.1 用正交分解完整推导最小二乘公式

上面的推导从“最佳残差与列空间正交”出发。还可以直接从距离最小化出发，得到一条更完整的推导。

最小二乘问题是：

$$
\min_x\|b-Ax\|^2.
$$

代入 $A=QR$：

$$
\min_x\|b-QRx\|^2.
$$

将 $b$ 分解成列空间内和列空间外两部分：

$$
b=QQ^Tb+(I-QQ^T)b.
$$

其中：

$$
QQ^Tb\in\operatorname{Col}(Q),
$$

而：

$$
(I-QQ^T)b\perp\operatorname{Col}(Q).
$$

于是残差可以写成：

$$
\begin{aligned}
b-QRx
&=QQ^Tb+(I-QQ^T)b-QRx\\
&=Q(Q^Tb-Rx)+(I-QQ^T)b.
\end{aligned}
$$

这两个部分互相正交：

- $Q(Q^Tb-Rx)$ 位于 $\operatorname{Col}(Q)$ 中；
- $(I-QQ^T)b$ 位于 $\operatorname{Col}(Q)$ 的正交补中。

因此根据勾股定理：

$$
\boxed{
\|b-QRx\|^2
=
\|Q(Q^Tb-Rx)\|^2
+
\|(I-QQ^T)b\|^2
}.
$$

由于 $Q^TQ=I$，$Q$ 不改变它所接收向量的长度：

$$
\|Qz\|^2
=z^TQ^TQz
=z^Tz
=\|z\|^2.
$$

所以：

$$
\boxed{
\|b-Ax\|^2
=
\|Q^Tb-Rx\|^2
+
\|(I-QQ^T)b\|^2
}.
$$

现在观察右侧两项：

1. $\|Q^Tb-Rx\|^2$ 会随 $x$ 改变；
2. $\|(I-QQ^T)b\|^2$ 与 $x$ 无关，是目标 $b$ 到列空间的固定距离。

所以最小二乘只能让第一项尽可能小。

当 $A$ 满列秩时，$R$ 是可逆的 $n\times n$ 上三角矩阵，因此可以选择唯一的 $\widehat{x}$ 使：

$$
Q^Tb-R\widehat{x}=0.
$$

由此得到：

$$
\boxed{R\widehat{x}=Q^Tb}.
$$

此时无法消除的残差为：

$$
\boxed{
r=b-A\widehat{x}=(I-QQ^T)b
}.
$$

而最佳近似输出为：

$$
\boxed{
A\widehat{x}=QQ^Tb
}.
$$

因此，QR 最小二乘实际上同时完成了三件事：

```text
Qᵀb       ：读取 b 在列空间中的坐标
R x̂=Qᵀb  ：找到能产生这些坐标的输入参数 x̂
Q R x̂    ：在输出空间中重建 b 的正交投影
```

### 12.2 精确方程与最小二乘方程的区别

如果：

$$
b\in\operatorname{Col}(A),
$$

那么：

$$
(I-QQ^T)b=0.
$$

此时：

$$
A\widehat{x}=b,
$$

QR 解出的就是精确解。

如果：

$$
b\notin\operatorname{Col}(A),
$$

那么：

$$
(I-QQ^T)b\neq0.
$$

这部分位于列空间外，任何 $Ax$ 都无法产生。QR 只能精确匹配列空间内的坐标：

$$
R\widehat{x}=Q^Tb,
$$

并留下不可消除的正交残差：

$$
r=(I-QQ^T)b.
$$

所以“解上三角方程”并不意味着原方程 $Ax=b$ 一定有精确解；它意味着找到了使 $Ax$ 等于 $b$ 的正交投影的参数。

---

## 13. QR 最小二乘的几何意义

经济型 $Q$ 的列是 $\operatorname{Col}(A)$ 的标准正交基。

所以：

$$
Q^Tb
$$

提取了目标 $b$ 在列空间各个标准正交方向上的坐标。

另一方面：

$$
A\widehat{x}=QR\widehat{x}.
$$

$R\widehat{x}$ 是输出 $A\widehat{x}$ 在 $Q$ 这组标准正交基下的坐标。

为了让 $A\widehat{x}$ 等于 $b$ 在列空间中的投影，这两组坐标必须相同：

$$
\boxed{
R\widehat{x}=Q^Tb
}.
$$

可以这样理解：

```text
目标 b
  │
  │ Qᵀ：读取 b 在列空间中的坐标
  ▼
Qᵀb

输入 x̂
  │
  │ R：计算 Ax̂ 在 Q 基下的坐标
  ▼
R x̂

最小二乘要求：R x̂ = Qᵀb
```

---

## 14. 使用 QR 重算第四步的最小二乘例子

第四步使用过：

$$
A=
\begin{bmatrix}
1&0\\
0&1\\
1&1
\end{bmatrix},
\qquad
b=
\begin{bmatrix}
1\\
1\\
3
\end{bmatrix}.
$$

矩阵的两列是：

$$
a_1=
\begin{bmatrix}
1\\
0\\
1
\end{bmatrix},
\qquad
a_2=
\begin{bmatrix}
0\\
1\\
1
\end{bmatrix}.
$$

对它们进行 Gram–Schmidt，得到：

$$
q_1
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1\\
0\\
1
\end{bmatrix},
$$

$$
q_2
=
\frac{1}{\sqrt{6}}
\begin{bmatrix}
-1\\
2\\
1
\end{bmatrix}.
$$

因此：

$$
Q=
\begin{bmatrix}
1/\sqrt{2}&-1/\sqrt{6}\\
0&2/\sqrt{6}\\
1/\sqrt{2}&1/\sqrt{6}
\end{bmatrix},
$$

$$
R=
\begin{bmatrix}
\sqrt{2}&1/\sqrt{2}\\
0&\sqrt{6}/2
\end{bmatrix}.
$$

### 第一步：计算 $Q^Tb$

$$
q_1^Tb
=
\frac{1+3}{\sqrt{2}}
=2\sqrt{2},
$$

$$
q_2^Tb
=
\frac{-1+2+3}{\sqrt{6}}
=
\frac{4}{\sqrt{6}}.
$$

所以：

$$
Q^Tb
=
\begin{bmatrix}
2\sqrt{2}\\
4/\sqrt{6}
\end{bmatrix}.
$$

### 第二步：求解上三角方程

$$
\begin{bmatrix}
\sqrt{2}&1/\sqrt{2}\\
0&\sqrt{6}/2
\end{bmatrix}
\begin{bmatrix}
\widehat{x}_1\\
\widehat{x}_2
\end{bmatrix}
=
\begin{bmatrix}
2\sqrt{2}\\
4/\sqrt{6}
\end{bmatrix}.
$$

先从第二行开始：

$$
\frac{\sqrt{6}}{2}\widehat{x}_2
=
\frac{4}{\sqrt{6}},
$$

得到：

$$
\widehat{x}_2=\frac{4}{3}.
$$

代回第一行：

$$
\sqrt{2}\widehat{x}_1
+
\frac{1}{\sqrt{2}}\cdot\frac{4}{3}
=
2\sqrt{2},
$$

得到：

$$
\widehat{x}_1=\frac{4}{3}.
$$

因此：

$$
\boxed{
\widehat{x}
=
\begin{bmatrix}
4/3\\
4/3
\end{bmatrix}
}.
$$

与第四步使用正规方程得到的结果完全相同。

---

## 15. QR 与正规方程的关系

第四步的正规方程是：

$$
A^TA\widehat{x}=A^Tb.
$$

代入：

$$
A=QR,
$$

得到：

$$
(QR)^T(QR)\widehat{x}=(QR)^Tb.
$$

展开：

$$
R^TQ^TQR\widehat{x}=R^TQ^Tb.
$$

因为：

$$
Q^TQ=I,
$$

所以：

$$
R^TR\widehat{x}=R^TQ^Tb.
$$

满列秩时 $R$ 可逆，因此可以消去 $R^T$，得到：

$$
R\widehat{x}=Q^Tb.
$$

两个方法求的是同一个最小二乘解，但计算路径不同：

```text
正规方程：先构造 AᵀA，再求解
QR 方法：先把 A 分解为 QR，再解上三角方程
```

---

## 16. 为什么数值计算通常更偏向 QR

### 16.1 $Q$ 不会放大向量长度

因为 $Q$ 的列标准正交：

$$
\|Qy\|=\|y\|
$$

对任意 $y\in\mathbb{R}^n$ 成立。

所以 $Q$ 只改变坐标表示，不会额外拉伸或压缩。

### 16.2 避免显式求逆

QR 方法最终求解：

$$
R\widehat{x}=Q^Tb.
$$

使用回代即可，不需要计算：

$$
(A^TA)^{-1}.
$$

### 16.3 避免直接构造 $A^TA$

构造 $A^TA$ 会使条件问题变得更严重。直观上，如果 $A$ 的两个列方向已经非常接近，那么 $A^TA$ 会让这种接近相关的问题更加突出。

从条件数角度，通常有：

$$
\kappa(A^TA)=\kappa(A)^2.
$$

所以正规方程可能损失更多数值精度。

QR 分解直接处理 $A$，通常比构造 $A^TA$ 更稳定。

### 16.4 实际软件不一定使用经典 Gram–Schmidt

Gram–Schmidt 最适合帮助理解 QR 的结构。

实际数值软件常使用：

- Householder 反射；
- Givens 旋转；
- 改进 Gram–Schmidt。

这些方法通常比经典 Gram–Schmidt 更稳定。

本阶段重点是理解 $A=QR$ 的结构，不要求立即掌握所有数值实现。

---

## 17. QR 分解是否唯一

如果：

$$
A=QR,
$$

把某个 $q_i$ 乘以 $-1$，同时把 $R$ 的第 $i$ 行乘以 $-1$，乘积 $QR$ 不会改变。

所以如果不增加约定，QR 分解可能有符号上的不唯一。

通常要求 $R$ 的对角线元素为正：

$$
r_{ii}>0.
$$

在 $A$ 满列秩的情况下，加上这个约定后，经济型 QR 分解唯一。

这不会改变列空间，只是在固定每个标准正交方向的正负号。

---

## 18. 如果 $A$ 的列线性相关，会发生什么

本步主要讨论满列秩矩阵。

如果 $A$ 的列线性相关，Gram–Schmidt 的某一步会得到：

$$
u_j=0.
$$

对应的：

$$
r_{jj}=\|u_j\|=0.
$$

这时 $R$ 的对角线上出现零，$R$ 不可逆，不能直接通过普通回代得到唯一参数。

几何意义是：

> 第 $j$ 列没有提供新的独立方向，输入参数中存在会被 $A$ 消除的自由方向。

处理秩亏矩阵时，通常需要：

- 带列主元的 QR；
- SVD；
- 伪逆。

这些内容将在后续阶段继续学习。

---

## 19. 怎样检查一个 QR 分解

给定候选矩阵 $Q$ 和 $R$，按以下顺序检查。

### 检查一：形状是否正确

经济型分解中：

$$
A\in\mathbb{R}^{m\times n},
$$

$$
Q\in\mathbb{R}^{m\times n},
\qquad
R\in\mathbb{R}^{n\times n}.
$$

### 检查二：$Q$ 的列是否标准正交

$$
Q^TQ=I_n.
$$

### 检查三：$R$ 是否上三角

主对角线下方应全部为零。

### 检查四：是否能够重建 $A$

$$
QR=A.
$$

### 检查五：列空间是否一致

$$
\operatorname{Col}(Q)=\operatorname{Col}(A).
$$

对于由 Gram–Schmidt 得到的分解，还可以检查：

$$
R=Q^TA.
$$

---

## 20. 一套可以直接使用的计算流程

给定满列秩矩阵：

$$
A=
\begin{bmatrix}
|&|&&|\\
a_1&a_2&\cdots&a_n\\
|&|&&|
\end{bmatrix}.
$$

### 第一步：对列向量做 Gram–Schmidt

得到标准正交向量：

$$
q_1,q_2,\ldots,q_n.
$$

### 第二步：组成 $Q$

$$
Q=
\begin{bmatrix}
|&|&&|\\
q_1&q_2&\cdots&q_n\\
|&|&&|
\end{bmatrix}.
$$

### 第三步：计算 $R$

可以逐项计算：

$$
r_{ij}=q_i^Ta_j,
$$

也可以整体计算：

$$
R=Q^TA.
$$

### 第四步：验证分解

$$
Q^TQ=I,
$$

$$
QR=A.
$$

### 第五步：用于最小二乘

计算：

$$
y=Q^Tb,
$$

然后回代求解：

$$
R\widehat{x}=y.
$$

---

## 21. 常见误区

### 误区一：认为经济型 $Q$ 一定是方阵

若：

$$
A\in\mathbb{R}^{m\times n},
\qquad
m>n,
$$

则经济型：

$$
Q\in\mathbb{R}^{m\times n}
$$

不是方阵，但仍满足：

$$
Q^TQ=I_n.
$$

### 误区二：把 $Q^TQ$ 和 $QQ^T$ 都当成单位矩阵

经济型 $Q$ 满足：

$$
Q^TQ=I_n,
$$

而：

$$
QQ^T
$$

是投影到列空间的矩阵，通常不等于 $I_m$。

### 误区三：认为 $R$ 只记录向量长度

$R$ 的对角线元素与正交化残差的长度有关，但非对角元素记录原列向量在旧方向上的坐标。

所以 $R$ 保存的是完整的坐标关系，不只是长度。

### 误区四：忘记 $R$ 为什么上三角

第 $j$ 个原列向量只会使用 $q_1,\ldots,q_j$，不会使用后面才产生的方向。因此 $i>j$ 时 $r_{ij}=0$。

### 误区五：使用 QR 求最小二乘后还去计算 $A^TA$

QR 方法直接求解：

$$
R\widehat{x}=Q^Tb.
$$

不需要再构造正规方程。

### 误区六：忽略列向量的顺序

Gram–Schmidt 按列依次处理，所以交换 $A$ 的列通常会改变 $Q$ 和 $R$。

虽然新的分解仍然描述交换后的矩阵，但不能随意打乱列顺序后继续使用原来的 $R$。

---

## 22. 回答阶段任务中的六个问题

### 1. $Q$ 和 $R$ 分别记录什么信息？

$Q$ 的列记录 $A$ 的列空间中的标准正交方向；$R$ 的每一列记录对应原列向量在这些方向下的坐标。

### 2. 为什么 $Q$ 与 $A$ 的列空间相同？

Gram–Schmidt 只在原列空间内部减去投影并进行单位化，没有增加或删除独立方向。因此：

$$
\operatorname{Col}(Q)=\operatorname{Col}(A).
$$

### 3. 为什么 $R$ 是上三角矩阵？

第 $j$ 个原列向量只需要 $q_1,\ldots,q_j$ 表示，所以 $i>j$ 时 $r_{ij}=0$。

### 4. 经济型 QR 与完整 QR 有什么区别？

经济型 $Q$ 只包含描述列空间所需的 $n$ 个方向；完整 $Q$ 继续补充为整个 $\mathbb{R}^m$ 的 $m$ 个标准正交方向。

### 5. 为什么 QR 能把最小二乘问题化成上三角方程？

因为 $A=QR$、$Q^TQ=I$，残差正交条件可以化为：

$$
R\widehat{x}=Q^Tb.
$$

$R$ 是上三角矩阵，可以通过回代求解。

### 6. 为什么数值计算通常优先使用 QR，而不是显式求逆？

QR 避免显式计算逆矩阵，也避免直接构造可能放大条件问题的 $A^TA$。正交矩阵的计算性质良好，因此 QR 通常更加稳定。

---

## 23. 本步总结

QR 分解把矩阵写成：

$$
\boxed{A=QR}.
$$

其中：

$$
\boxed{Q^TQ=I},
$$

$$
\boxed{
\operatorname{Col}(Q)=\operatorname{Col}(A)
}.
$$

$R$ 是上三角矩阵，并记录原列向量在 $Q$ 的标准正交列下的坐标：

$$
\boxed{R=Q^TA}.
$$

经济型分解的维度为：

$$
\underbrace{A}_{m\times n}
=
\underbrace{Q}_{m\times n}
\underbrace{R}_{n\times n}.
$$

使用 QR 求最小二乘解：

$$
\boxed{
R\widehat{x}=Q^Tb
}.
$$

最终的心智模型是：

```text
A 的普通列向量
       │
       │ Gram–Schmidt
       ▼
Q：相同列空间中的标准正交方向
       │
       │ R：记录原列在这些方向下的坐标
       ▼
A=QR

最小二乘：
b ──Qᵀ──→ 列空间坐标 Qᵀb
x̂ ──R───→ 输出坐标 R x̂

令 R x̂=Qᵀb，再通过回代求解
```

下一步将把最小二乘和 QR 应用于真实数据：从观测点构建设计矩阵，并完成直线拟合与多项式拟合。

## 本节统一

QR 的核心不是机械分解，而是把一般列基换成标准正交基，并把最小二乘问题转化为稳定的上三角方程。它与 SVD 都依赖正交坐标，但回答的问题不同。

## 下一步为什么自然出现

QR 仍然围绕“列空间的一组好基”；如果希望对任意矩形矩阵同时找到最自然的输入方向、输出方向和缩放强度，就进入 SVD。
