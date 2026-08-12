# 第四阶段：正交、投影、最小二乘与 QR

> 核心问题：当目标向量不在矩阵的列空间中，无法精确求解 $A\boldsymbol x=\boldsymbol b$ 时，怎样找到距离目标最近的可实现结果？

第三阶段研究了方程何时有精确解。本阶段处理更常见的现实情况：观测数据存在噪声，约束彼此不完全一致，目标 $\boldsymbol b$ 不属于列空间，因此方程无解。

我们不再要求：

$$
A\boldsymbol x=\boldsymbol b
$$

而是寻找一个 $\hat{\boldsymbol x}$，使：

$$
A\hat{\boldsymbol x}
$$

尽可能接近 $\boldsymbol b$。这个问题把内积、正交、投影、Gram–Schmidt、QR 分解和最小二乘连接成一条主线：

$$
\text{内积与长度}
\rightarrow
\text{正交}
\rightarrow
\text{正交分解}
\rightarrow
\text{投影}
\rightarrow
\text{Gram--Schmidt}
\rightarrow
\text{QR 分解}
\rightarrow
\text{最小二乘}
$$

本阶段最重要的几何图像是：

$$
\boxed{
\boldsymbol b
=
\underbrace{A\hat{\boldsymbol x}}_{\text{列空间中的投影}}
+
\underbrace{\boldsymbol r}_{\text{与列空间正交的残差}}
}
$$

## 学习目标

完成本阶段后，应当能够：

1. 计算向量内积、范数、距离和夹角。
2. 判断向量、向量组和子空间之间的正交关系。
3. 理解正交补，并解释四个基本子空间之间的正交关系。
4. 将向量投影到一条直线或一个子空间。
5. 推导并使用投影矩阵。
6. 判断一个矩阵是否为正交投影矩阵。
7. 使用 Gram–Schmidt 方法将线性无关向量正交标准化。
8. 计算矩阵的 QR 分解，并解释 $Q$ 与 $R$ 的意义。
9. 从几何和代数两个角度推导正规方程。
10. 求解线性最小二乘问题。
11. 使用 QR 分解更稳定地求最小二乘解。
12. 完成直线拟合和简单多项式拟合。
13. 理解残差、拟合误差与数值稳定性。

## 1. 内积：测量方向的一致程度

### 1.1 定义

对于：

$$
\boldsymbol u=
\begin{bmatrix}
u_1\\u_2\\\vdots\\u_n
\end{bmatrix},
\qquad
\boldsymbol v=
\begin{bmatrix}
v_1\\v_2\\\vdots\\v_n
\end{bmatrix}
$$

标准内积定义为：

$$
\boldsymbol u^\mathsf T\boldsymbol v
=
u_1v_1+u_2v_2+\cdots+u_nv_n
$$

内积也称点积。

例如：

$$
\boldsymbol u=
\begin{bmatrix}1\\2\\-1\end{bmatrix},
\qquad
\boldsymbol v=
\begin{bmatrix}2\\0\\2\end{bmatrix}
$$

则：

$$
\boldsymbol u^\mathsf T\boldsymbol v
=1\cdot2+2\cdot0+(-1)\cdot2
=0
$$

### 1.2 内积的基本性质

对于任意向量 $\boldsymbol u,\boldsymbol v,\boldsymbol w$ 和标量 $c$：

$$
\boldsymbol u^\mathsf T\boldsymbol v
=
\boldsymbol v^\mathsf T\boldsymbol u
$$

$$
\boldsymbol u^\mathsf T(\boldsymbol v+\boldsymbol w)
=
\boldsymbol u^\mathsf T\boldsymbol v
+
\boldsymbol u^\mathsf T\boldsymbol w
$$

$$
(c\boldsymbol u)^\mathsf T\boldsymbol v
=
c(\boldsymbol u^\mathsf T\boldsymbol v)
$$

$$
\boldsymbol u^\mathsf T\boldsymbol u\ge0
$$

并且：

$$
\boldsymbol u^\mathsf T\boldsymbol u=0
\iff
\boldsymbol u=\boldsymbol0
$$

### 1.3 几何意义

内积还满足：

$$
\boldsymbol u^\mathsf T\boldsymbol v
=
\|\boldsymbol u\|
\|\boldsymbol v\|
\cos\theta
$$

其中 $\theta$ 是两个非零向量之间的夹角。

因此：

- 内积为正：夹角小于 $90^\circ$，方向总体一致。
- 内积为零：夹角为 $90^\circ$，两个向量正交。
- 内积为负：夹角大于 $90^\circ$，方向总体相反。

内积并不直接等于夹角，它同时受到两个向量长度的影响。

## 2. 范数、距离与夹角

### 2.1 欧几里得范数

向量的长度定义为：

$$
\|\boldsymbol v\|
=
\sqrt{\boldsymbol v^\mathsf T\boldsymbol v}
=
\sqrt{v_1^2+\cdots+v_n^2}
$$

例如：

$$
\left\|
\begin{bmatrix}3\\4\end{bmatrix}
\right\|
=5
$$

### 2.2 单位向量

长度为 1 的向量称为单位向量。对于非零向量 $\boldsymbol v$，其单位化结果为：

$$
\boldsymbol q
=
\frac{\boldsymbol v}{\|\boldsymbol v\|}
$$

单位化保留方向，只去掉长度。

### 2.3 两点之间的距离

向量 $\boldsymbol u$ 与 $\boldsymbol v$ 之间的欧几里得距离为：

$$
d(\boldsymbol u,\boldsymbol v)
=
\|\boldsymbol u-\boldsymbol v\|
$$

最小二乘问题正是要最小化目标向量与可实现输出之间的距离。

### 2.4 夹角

对于非零向量：

$$
\cos\theta
=
\frac{\boldsymbol u^\mathsf T\boldsymbol v}
{\|\boldsymbol u\|\|\boldsymbol v\|}
$$

例如：

$$
\boldsymbol u=
\begin{bmatrix}1\\0\end{bmatrix},
\qquad
\boldsymbol v=
\begin{bmatrix}1\\1\end{bmatrix}
$$

则：

$$
\cos\theta=\frac{1}{\sqrt2}
$$

所以 $\theta=45^\circ$。

### 2.5 柯西—施瓦茨不等式

任意向量满足：

$$
|\boldsymbol u^\mathsf T\boldsymbol v|
\le
\|\boldsymbol u\|\|\boldsymbol v\|
$$

当且仅当两个向量线性相关时取等号。这个不等式保证夹角公式中的余弦值位于 $[-1,1]$。

### 2.6 三角不等式

$$
\|\boldsymbol u+\boldsymbol v\|
\le
\|\boldsymbol u\|+\|\boldsymbol v\|
$$

它说明从原点直接走到终点，不会比先沿 $\boldsymbol u$ 再沿 $\boldsymbol v$ 更远。

## 3. 正交：彼此不干扰的方向

### 3.1 正交向量

若：

$$
\boldsymbol u^\mathsf T\boldsymbol v=0
$$

则称 $\boldsymbol u$ 与 $\boldsymbol v$ 正交，记作：

$$
\boldsymbol u\perp\boldsymbol v
$$

零向量与所有向量都正交，但讨论方向时通常关注非零向量。

### 3.2 正交向量组

若一组非零向量中任意两个不同向量都正交，则称其为正交向量组。

正交非零向量组一定线性无关。

证明：若

$$
c_1\boldsymbol v_1+\cdots+c_k\boldsymbol v_k=\boldsymbol0
$$

两边与 $\boldsymbol v_j$ 做内积，得到：

$$
c_j\|\boldsymbol v_j\|^2=0
$$

因为 $\boldsymbol v_j\neq\boldsymbol0$，所以 $c_j=0$。对每个 $j$ 都成立，因此向量组线性无关。

### 3.3 标准正交向量组

如果正交向量组中的每个向量长度都为 1，则称为标准正交向量组：

$$
\boldsymbol q_i^\mathsf T\boldsymbol q_j
=
\begin{cases}
1,&i=j\\
0,&i\neq j
\end{cases}
$$

若把这些向量作为矩阵 $Q$ 的列，则：

$$
Q^\mathsf TQ=I
$$

### 3.4 勾股定理

若 $\boldsymbol u\perp\boldsymbol v$，则：

$$
\|\boldsymbol u+\boldsymbol v\|^2
=
\|\boldsymbol u\|^2+\|\boldsymbol v\|^2
$$

证明：

$$
\|\boldsymbol u+\boldsymbol v\|^2
=(\boldsymbol u+\boldsymbol v)^\mathsf T
(\boldsymbol u+\boldsymbol v)
$$

展开后交叉项为零。

### 3.5 标准正交基的优势

若 $\boldsymbol q_1,\dots,\boldsymbol q_k$ 是子空间的一组标准正交基，并且：

$$
\boldsymbol v
=c_1\boldsymbol q_1+\cdots+c_k\boldsymbol q_k
$$

则每个坐标无需解方程组，直接由内积得到：

$$
c_i=\boldsymbol q_i^\mathsf T\boldsymbol v
$$

因此标准正交基是最便于计算的坐标系统。

## 4. 正交补与四个基本子空间

### 4.1 正交补

对于子空间 $W\subseteq\mathbb R^n$，其正交补定义为：

$$
W^\perp
=
\{\boldsymbol x\in\mathbb R^n:
\boldsymbol x^\mathsf T\boldsymbol w=0
\text{ 对所有 }\boldsymbol w\in W\}
$$

$W^\perp$ 本身也是子空间。

例如，在 $\mathbb R^3$ 中：

- 一个经过原点的平面的正交补是一条法线方向。
- 一条经过原点的直线的正交补是与它垂直的平面。

### 4.2 行空间与零空间

若 $\boldsymbol x\in N(A)$，则：

$$
A\boldsymbol x=\boldsymbol0
$$

这表示 $A$ 的每一行都与 $\boldsymbol x$ 正交。因此：

$$
N(A)=C(A^\mathsf T)^\perp
$$

输入空间分解为：

$$
\mathbb R^n
=
C(A^\mathsf T)\oplus N(A)
$$

### 4.3 列空间与左零空间

同理：

$$
N(A^\mathsf T)=C(A)^\perp
$$

输出空间分解为：

$$
\mathbb R^m
=
C(A)\oplus N(A^\mathsf T)
$$

这里的 $\oplus$ 表示正交直和：每个向量都能唯一分解成分别属于两个互相正交子空间的向量。

### 4.4 与维数的关系

若 $A$ 是 $m\times n$ 矩阵，秩为 $r$，则：

$$
\dim C(A^\mathsf T)+\dim N(A)
=r+(n-r)=n
$$

$$
\dim C(A)+\dim N(A^\mathsf T)
=r+(m-r)=m
$$

这不仅是维数相加，更说明两个子空间在几何上恰好补全整个空间。

## 5. 投影到一条直线

### 5.1 问题

给定非零向量 $\boldsymbol a$ 和目标 $\boldsymbol b$，希望在直线

$$
\operatorname{span}(\boldsymbol a)
$$

上找到距离 $\boldsymbol b$ 最近的向量 $\boldsymbol p$。

因为 $\boldsymbol p$ 位于该直线上，所以：

$$
\boldsymbol p=c\boldsymbol a
$$

最近点要求残差：

$$
\boldsymbol r=\boldsymbol b-\boldsymbol p
$$

与直线方向 $\boldsymbol a$ 正交：

$$
\boldsymbol a^\mathsf T(\boldsymbol b-c\boldsymbol a)=0
$$

解得：

$$
c
=
\frac{\boldsymbol a^\mathsf T\boldsymbol b}
{\boldsymbol a^\mathsf T\boldsymbol a}
$$

因此：

$$
\boxed{
\operatorname{proj}_{\boldsymbol a}\boldsymbol b
=
\frac{\boldsymbol a^\mathsf T\boldsymbol b}
{\boldsymbol a^\mathsf T\boldsymbol a}
\boldsymbol a
}
$$

### 5.2 例

设：

$$
\boldsymbol a=
\begin{bmatrix}2\\1\end{bmatrix},
\qquad
\boldsymbol b=
\begin{bmatrix}3\\4\end{bmatrix}
$$

投影系数为：

$$
c=\frac{2\cdot3+1\cdot4}{2^2+1^2}=2
$$

所以：

$$
\boldsymbol p=
2\begin{bmatrix}2\\1\end{bmatrix}
=
\begin{bmatrix}4\\2\end{bmatrix}
$$

残差为：

$$
\boldsymbol r
=
\boldsymbol b-\boldsymbol p
=
\begin{bmatrix}-1\\2\end{bmatrix}
$$

验证：

$$
\boldsymbol a^\mathsf T\boldsymbol r
=2(-1)+1(2)=0
$$

### 5.3 单位方向的简化

若 $\boldsymbol q$ 是单位向量，则：

$$
\operatorname{proj}_{\boldsymbol q}\boldsymbol b
=(\boldsymbol q^\mathsf T\boldsymbol b)\boldsymbol q
$$

因为 $\boldsymbol q^\mathsf T\boldsymbol q=1$。

## 6. 投影到子空间

### 6.1 使用标准正交基

若 $W$ 的标准正交基为：

$$
\boldsymbol q_1,\dots,\boldsymbol q_k
$$

则 $\boldsymbol b$ 在 $W$ 上的投影为：

$$
\boldsymbol p
=
\sum_{i=1}^k
(\boldsymbol q_i^\mathsf T\boldsymbol b)\boldsymbol q_i
$$

令：

$$
Q=
\begin{bmatrix}
\boldsymbol q_1&\cdots&\boldsymbol q_k
\end{bmatrix}
$$

则：

$$
\boxed{
\boldsymbol p=QQ^\mathsf T\boldsymbol b
}
$$

投影矩阵为：

$$
P=QQ^\mathsf T
$$

### 6.2 正交分解

任意 $\boldsymbol b\in\mathbb R^m$ 都能唯一写成：

$$
\boldsymbol b=\boldsymbol p+\boldsymbol r
$$

其中：

$$
\boldsymbol p\in W,
\qquad
\boldsymbol r\in W^\perp
$$

由于二者正交：

$$
\|\boldsymbol b\|^2
=
\|\boldsymbol p\|^2+\|\boldsymbol r\|^2
$$

### 6.3 投影为什么给出最近点

对任意 $\boldsymbol w\in W$：

$$
\boldsymbol b-\boldsymbol w
=(\boldsymbol b-\boldsymbol p)
+(\boldsymbol p-\boldsymbol w)
$$

其中：

- $\boldsymbol b-\boldsymbol p\in W^\perp$。
- $\boldsymbol p-\boldsymbol w\in W$。

两者正交，所以：

$$
\|\boldsymbol b-\boldsymbol w\|^2
=
\|\boldsymbol b-\boldsymbol p\|^2
+\|\boldsymbol p-\boldsymbol w\|^2
\ge
\|\boldsymbol b-\boldsymbol p\|^2
$$

当且仅当 $\boldsymbol w=\boldsymbol p$ 时取等号。因此投影是唯一最近点。

## 7. 一般列空间上的投影矩阵

### 7.1 非正交基的情况

设 $A$ 的列线性无关，列空间为 $C(A)$。投影点可以写成：

$$
\boldsymbol p=A\hat{\boldsymbol x}
$$

残差必须与列空间中每一列正交：

$$
A^\mathsf T(\boldsymbol b-A\hat{\boldsymbol x})=\boldsymbol0
$$

因此：

$$
A^\mathsf TA\hat{\boldsymbol x}
=
A^\mathsf T\boldsymbol b
$$

若 $A$ 满列秩，则 $A^\mathsf TA$ 可逆：

$$
\hat{\boldsymbol x}
=(A^\mathsf TA)^{-1}A^\mathsf T\boldsymbol b
$$

所以：

$$
\boxed{
P
=
A(A^\mathsf TA)^{-1}A^\mathsf T
}
$$

### 7.2 投影矩阵的性质

正交投影矩阵满足：

$$
P^2=P
$$

因为已经投影到目标子空间中的向量，再投影一次不会改变。

它还满足：

$$
P^\mathsf T=P
$$

即正交投影矩阵是对称矩阵。

补空间投影矩阵为：

$$
I-P
$$

它把向量投影到 $C(A)^\perp=N(A^\mathsf T)$。

### 7.3 $P$ 与 $I-P$

$$
\boldsymbol p=P\boldsymbol b
$$

$$
\boldsymbol r=(I-P)\boldsymbol b
$$

并且：

$$
P\boldsymbol r=\boldsymbol0
$$

$$
(I-P)\boldsymbol p=\boldsymbol0
$$

这两个矩阵把输出空间分解为列空间和左零空间。

### 7.4 公式的限制

公式

$$
(A^\mathsf TA)^{-1}
$$

要求 $A$ 的列线性无关。若列相关，则 $A^\mathsf TA$ 不可逆，需要删除冗余列、使用 QR 的适当变体，或在第六阶段使用伪逆和 SVD。

## 8. Gram–Schmidt 正交化

### 8.1 目标

给定线性无关向量：

$$
\boldsymbol a_1,\dots,\boldsymbol a_k
$$

希望构造标准正交向量：

$$
\boldsymbol q_1,\dots,\boldsymbol q_k
$$

同时保持张成空间不变：

$$
\operatorname{span}(\boldsymbol a_1,\dots,\boldsymbol a_k)
=
\operatorname{span}(\boldsymbol q_1,\dots,\boldsymbol q_k)
$$

### 8.2 两个向量的过程

首先：

$$
\boldsymbol u_1=\boldsymbol a_1
$$

$$
\boldsymbol q_1
=
\frac{\boldsymbol u_1}{\|\boldsymbol u_1\|}
$$

然后从 $\boldsymbol a_2$ 中减去它在 $\boldsymbol q_1$ 方向上的投影：

$$
\boldsymbol u_2
=
\boldsymbol a_2
-(\boldsymbol q_1^\mathsf T\boldsymbol a_2)\boldsymbol q_1
$$

再单位化：

$$
\boldsymbol q_2
=
\frac{\boldsymbol u_2}{\|\boldsymbol u_2\|}
$$

### 8.3 一般形式

第 $j$ 个正交方向为：

$$
\boldsymbol u_j
=
\boldsymbol a_j
-
\sum_{i=1}^{j-1}
(\boldsymbol q_i^\mathsf T\boldsymbol a_j)\boldsymbol q_i
$$

然后：

$$
\boldsymbol q_j
=
\frac{\boldsymbol u_j}{\|\boldsymbol u_j\|}
$$

### 8.4 完整例子

设：

$$
\boldsymbol a_1=
\begin{bmatrix}1\\1\\0\end{bmatrix},
\qquad
\boldsymbol a_2=
\begin{bmatrix}1\\0\\1\end{bmatrix}
$$

第一步：

$$
\|\boldsymbol a_1\|=\sqrt2
$$

$$
\boldsymbol q_1
=
\frac1{\sqrt2}
\begin{bmatrix}1\\1\\0\end{bmatrix}
$$

第二步，计算投影：

$$
\boldsymbol q_1^\mathsf T\boldsymbol a_2
=
\frac1{\sqrt2}
$$

所以：

$$
\boldsymbol u_2
=
\begin{bmatrix}1\\0\\1\end{bmatrix}
-
\frac12
\begin{bmatrix}1\\1\\0\end{bmatrix}
=
\begin{bmatrix}
\frac12\\-\frac12\\1
\end{bmatrix}
$$

其长度为：

$$
\|\boldsymbol u_2\|
=
\sqrt{\frac32}
$$

单位化得到：

$$
\boldsymbol q_2
=
\frac1{\sqrt6}
\begin{bmatrix}1\\-1\\2\end{bmatrix}
$$

验证：

$$
\boldsymbol q_1^\mathsf T\boldsymbol q_2=0
$$

$$
\|\boldsymbol q_1\|=\|\boldsymbol q_2\|=1
$$

### 8.5 线性相关时会发生什么

如果某个 $\boldsymbol a_j$ 位于前面向量的张成空间中，那么减去全部投影后：

$$
\boldsymbol u_j=\boldsymbol0
$$

此时无法单位化。这恰好检测出原向量组存在冗余。

## 9. QR 分解

### 9.1 从 Gram–Schmidt 到 QR

设满列秩矩阵：

$$
A=
\begin{bmatrix}
\boldsymbol a_1&\cdots&\boldsymbol a_n
\end{bmatrix}
$$

Gram–Schmidt 产生标准正交列：

$$
Q=
\begin{bmatrix}
\boldsymbol q_1&\cdots&\boldsymbol q_n
\end{bmatrix}
$$

每个 $\boldsymbol a_j$ 都可以由前 $j$ 个 $\boldsymbol q_i$ 表示：

$$
\boldsymbol a_j
=
r_{1j}\boldsymbol q_1+\cdots+r_{jj}\boldsymbol q_j
$$

把系数排列成上三角矩阵 $R$，得到：

$$
\boxed{A=QR}
$$

### 9.2 $Q$ 与 $R$ 的意义

- $Q$ 的列是 $C(A)$ 的一组标准正交基。
- $Q^\mathsf TQ=I$。
- $R$ 记录原列向量在正交基 $Q$ 下的坐标。
- $R$ 是上三角矩阵。
- $A$ 与 $Q$ 具有相同列空间。

对于满列秩的高矩阵，$Q$ 通常是 $m\times n$，$R$ 是 $n\times n$。这称为经济型或约化 QR 分解。

### 9.3 继续前面的例子

对于：

$$
A=
\begin{bmatrix}
1&1\\
1&0\\
0&1
\end{bmatrix}
$$

前面已经得到：

$$
Q=
\begin{bmatrix}
\frac1{\sqrt2}&\frac1{\sqrt6}\\
\frac1{\sqrt2}&-\frac1{\sqrt6}\\
0&\frac2{\sqrt6}
\end{bmatrix}
$$

$R$ 的元素为：

$$
r_{ij}=\boldsymbol q_i^\mathsf T\boldsymbol a_j
$$

所以：

$$
R=
\begin{bmatrix}
\sqrt2&\frac1{\sqrt2}\\
0&\sqrt{\frac32}
\end{bmatrix}
$$

并且确有：

$$
A=QR
$$

### 9.4 方阵正交矩阵

若 $Q$ 是方阵并满足：

$$
Q^\mathsf TQ=I
$$

则：

$$
Q^{-1}=Q^\mathsf T
$$

正交矩阵保持内积、长度和角度：

$$
\|Q\boldsymbol x\|=\|\boldsymbol x\|
$$

旋转矩阵和反射矩阵都是正交矩阵。

## 10. 最小二乘问题

### 10.1 从无解到最佳近似

当：

$$
A\boldsymbol x=\boldsymbol b
$$

无精确解时，寻找：

$$
\boxed{
\hat{\boldsymbol x}
=
\arg\min_{\boldsymbol x}
\|A\boldsymbol x-\boldsymbol b\|^2
}
$$

之所以使用平方，是因为：

- 平方保持非负。
- 最小化范数与最小化范数平方得到同一结果。
- 平方形式易于求导和展开。
- 每个残差分量的平方可以累加。

### 10.2 几何解释

所有 $A\boldsymbol x$ 构成列空间 $C(A)$。因此最小二乘是在列空间中寻找距离 $\boldsymbol b$ 最近的向量：

$$
\boldsymbol p=A\hat{\boldsymbol x}
$$

残差：

$$
\boldsymbol r
=
\boldsymbol b-A\hat{\boldsymbol x}
$$

必须与列空间正交：

$$
\boldsymbol r\perp C(A)
$$

### 10.3 正规方程

因为残差与 $A$ 的每一列正交：

$$
A^\mathsf T\boldsymbol r=\boldsymbol0
$$

代入残差：

$$
A^\mathsf T(\boldsymbol b-A\hat{\boldsymbol x})
=
\boldsymbol0
$$

得到正规方程：

$$
\boxed{
A^\mathsf TA\hat{\boldsymbol x}
=
A^\mathsf T\boldsymbol b
}
$$

若 $A$ 满列秩：

$$
\hat{\boldsymbol x}
=(A^\mathsf TA)^{-1}A^\mathsf T\boldsymbol b
$$

### 10.4 为什么叫正规方程

这里的“正规”来自法向或正交条件。正规方程表达的不是原方程已经变得精确，而是误差向量已经垂直于所有可调整方向。

### 10.5 最小二乘解是否唯一

- 投影点 $\boldsymbol p$ 始终唯一。
- 若 $A$ 满列秩，则产生该投影点的 $\hat{\boldsymbol x}$ 唯一。
- 若 $A$ 的列相关，则可能有多个不同的 $\boldsymbol x$ 产生同一个最佳投影点。

第六阶段将使用伪逆从多个最小二乘解中选择范数最小的解。

## 11. 使用 QR 求最小二乘

若：

$$
A=QR
$$

其中 $Q^\mathsf TQ=I$，则：

$$
A^\mathsf TA
=
R^\mathsf TQ^\mathsf TQR
=
R^\mathsf TR
$$

正规方程变为：

$$
R^\mathsf TR\hat{\boldsymbol x}
=
R^\mathsf TQ^\mathsf T\boldsymbol b
$$

若 $R$ 可逆，可约去 $R^\mathsf T$：

$$
\boxed{
R\hat{\boldsymbol x}
=
Q^\mathsf T\boldsymbol b
}
$$

因为 $R$ 是上三角矩阵，可以通过回代快速求解。

QR 方法避免显式形成 $A^\mathsf TA$，通常比直接解正规方程更稳定。

投影也可以直接写成：

$$
\boldsymbol p=QQ^\mathsf T\boldsymbol b
$$

因为 $Q$ 与 $A$ 具有相同列空间。

## 12. 直线拟合：贯穿本阶段的完整例子

给定三个观测点：

$$
(0,1),\qquad(1,2),\qquad(2,2)
$$

希望拟合直线：

$$
y=\beta_0+\beta_1x
$$

### 12.1 建立矩阵模型

每个点给出一个方程：

$$
\begin{cases}
\beta_0=1\\
\beta_0+\beta_1=2\\
\beta_0+2\beta_1=2
\end{cases}
$$

写成：

$$
A\boldsymbol\beta=\boldsymbol b
$$

其中：

$$
A=
\begin{bmatrix}
1&0\\
1&1\\
1&2
\end{bmatrix},
\qquad
\boldsymbol\beta=
\begin{bmatrix}\beta_0\\\beta_1\end{bmatrix},
\qquad
\boldsymbol b=
\begin{bmatrix}1\\2\\2\end{bmatrix}
$$

三个点不在同一条直线上，所以原方程没有精确解。

### 12.2 建立正规方程

$$
A^\mathsf TA
=
\begin{bmatrix}
3&3\\
3&5
\end{bmatrix}
$$

$$
A^\mathsf T\boldsymbol b
=
\begin{bmatrix}
5\\6
\end{bmatrix}
$$

所以：

$$
\begin{bmatrix}
3&3\\
3&5
\end{bmatrix}
\begin{bmatrix}
\hat\beta_0\\\hat\beta_1
\end{bmatrix}
=
\begin{bmatrix}
5\\6
\end{bmatrix}
$$

两式相减得到：

$$
2\hat\beta_1=1
$$

因此：

$$
\hat\beta_1=\frac12
$$

再代回：

$$
\hat\beta_0=\frac76
$$

最佳拟合直线为：

$$
\boxed{
\hat y=\frac76+\frac12x
}
$$

### 12.3 投影与预测值

预测向量为：

$$
\boldsymbol p
=
A\hat{\boldsymbol\beta}
=
\begin{bmatrix}
\frac76\\
\frac53\\
\frac{13}{6}
\end{bmatrix}
$$

它是 $\boldsymbol b$ 在 $C(A)$ 上的正交投影。

### 12.4 残差

定义残差为观测值减预测值：

$$
\boldsymbol r
=
\boldsymbol b-\boldsymbol p
=
\begin{bmatrix}
-\frac16\\
\frac13\\
-\frac16
\end{bmatrix}
$$

验证正规条件：

$$
A^\mathsf T\boldsymbol r
=
\begin{bmatrix}
1&1&1\\
0&1&2
\end{bmatrix}
\begin{bmatrix}
-\frac16\\
\frac13\\
-\frac16
\end{bmatrix}
=
\begin{bmatrix}0\\0\end{bmatrix}
$$

第一行说明残差和为零；第二行说明残差与输入 $x$ 的加权和为零。

### 12.5 误差大小

残差平方和为：

$$
\|\boldsymbol r\|^2
=
\frac1{36}+\frac19+\frac1{36}
=
\frac16
$$

它是所有直线模型中能够达到的最小平方误差。

### 12.6 几何和统计语言的对应

| 线性代数语言 | 拟合语言 |
|---|---|
| $\boldsymbol b$ | 观测值 |
| $C(A)$ | 模型能够产生的全部预测 |
| $A\hat{\boldsymbol\beta}$ | 最佳预测值 |
| $\boldsymbol r$ | 残差 |
| 投影 | 选择最佳拟合模型 |
| $\|\boldsymbol r\|^2$ | 残差平方和 |

## 13. 多项式拟合与一般线性模型

### 13.1 多项式拟合

若拟合二次多项式：

$$
y=\beta_0+\beta_1x+\beta_2x^2
$$

对数据点 $(x_i,y_i)$，设计矩阵为：

$$
A=
\begin{bmatrix}
1&x_1&x_1^2\\
1&x_2&x_2^2\\
\vdots&\vdots&\vdots\\
1&x_m&x_m^2
\end{bmatrix}
$$

然后求：

$$
\min_{\boldsymbol\beta}
\|A\boldsymbol\beta-\boldsymbol y\|^2
$$

虽然模型关于 $x$ 是非线性的二次曲线，但它关于未知参数 $\beta_0,\beta_1,\beta_2$ 是线性的，因此仍是线性最小二乘问题。

### 13.2 多元线性模型

若每个样本有多个特征：

$$
y_i
=
\beta_0+\beta_1x_{i1}+\cdots+\beta_nx_{in}
$$

仍可组织为：

$$
A\boldsymbol\beta\approx\boldsymbol y
$$

矩阵的每一行对应一个样本，每一列对应一个特征或基函数。

### 13.3 模型空间

改变设计矩阵的列，相当于改变允许使用的模型基函数。最小二乘始终是在这些列张成的模型空间中寻找距离数据最近的预测。

## 14. 数值稳定性与计算方法

### 14.1 为什么不宜显式求逆

理论公式：

$$
\hat{\boldsymbol x}
=(A^\mathsf TA)^{-1}A^\mathsf T\boldsymbol b
$$

适合推导，但实际计算不应先求逆矩阵。通常直接求解线性系统：

$$
A^\mathsf TA\hat{\boldsymbol x}
=A^\mathsf T\boldsymbol b
$$

或者优先使用 QR 分解。

### 14.2 正规方程会放大条件问题

矩阵的条件数衡量输入误差或舍入误差被放大的程度。通常有：

$$
\kappa(A^\mathsf TA)\approx\kappa(A)^2
$$

因此形成 $A^\mathsf TA$ 可能让原本已经不稳定的问题变得更敏感。

### 14.3 QR 的优势

QR 分解使用正交变换。正交矩阵保持长度，不会额外放大误差，因此比直接构造正规方程更稳定。

实际数值软件通常使用 Householder 反射或 Givens 旋转计算 QR，而不是直接使用经典 Gram–Schmidt。

### 14.4 经典与改进 Gram–Schmidt

在精确数学中，两种写法等价；在浮点计算中，经典 Gram–Schmidt 可能因舍入误差失去正交性。改进 Gram–Schmidt 会逐步更新剩余向量，通常更稳定。

### 14.5 列接近相关

若 $A$ 的列几乎线性相关，则：

- $A^\mathsf TA$ 接近不可逆。
- 最小二乘系数可能对微小数据变化非常敏感。
- 预测值可能仍较稳定，但参数本身不稳定。
- 应考虑重新选择特征、正则化或使用 SVD。

## 15. 常见误区

### 误区 1：内积为零意味着其中一个向量为零

两个非零向量也可以内积为零，这表示它们正交。

### 误区 2：线性无关就等于正交

正交非零向量一定线性无关，但线性无关向量不一定正交。

### 误区 3：投影系数就是内积

只有投影方向是单位向量时，系数才是 $\boldsymbol q^\mathsf T\boldsymbol b$。一般情况必须除以 $\boldsymbol a^\mathsf T\boldsymbol a$。

### 误区 4：投影后的向量与原向量正交

投影向量 $\boldsymbol p$ 位于目标子空间中。真正与目标子空间正交的是残差：

$$
\boldsymbol r=\boldsymbol b-\boldsymbol p
$$

### 误区 5：$A\hat{\boldsymbol x}=\boldsymbol b$

最小二乘问题中通常没有精确相等。正确关系是：

$$
A\hat{\boldsymbol x}=P\boldsymbol b
$$

### 误区 6：正规方程让原方程重新变得有精确解

正规方程求的是投影点对应的参数。它没有改变 $\boldsymbol b$ 不属于列空间这一事实。

### 误区 7：残差与 $\boldsymbol b$ 正交

残差与 $C(A)$ 正交，因此与 $A$ 的每一列正交；它不一定与 $\boldsymbol b$ 正交。

### 误区 8：所有满足 $P^2=P$ 的矩阵都是正交投影

$P^2=P$ 只表示它是某种投影。正交投影还要求：

$$
P^\mathsf T=P
$$

### 误区 9：$Q^\mathsf TQ=I$ 就一定有 $QQ^\mathsf T=I$

若 $Q$ 是高矩阵且列标准正交，则 $Q^\mathsf TQ=I$，但 $QQ^\mathsf T$ 是投影到 $C(Q)$ 的矩阵，不一定是整个输出空间上的单位矩阵。

### 误区 10：最小残差意味着每个残差分量都最小

最小二乘最小化的是所有残差平方之和。某个单独数据点的误差可能并不是最小。

### 误区 11：QR 分解只是另一种矩阵乘法练习

QR 的核心是用标准正交基描述同一个列空间，使投影和求解变得简单稳定。

## 16. 应用连接

- **线性回归**：用多个特征预测连续目标。
- **曲线拟合**：使用多项式或其他基函数近似实验数据。
- **测量融合**：从带噪声的多次观测中估计位置或系统参数。
- **计算机视觉**：相机标定、位姿估计和三维重建常使用最小二乘。
- **信号处理**：把信号投影到一组基函数上，分离有效成分和噪声。
- **推荐系统**：在低维模型空间中近似用户评分。
- **控制与系统辨识**：根据输入输出数据估计动态模型参数。
- **数值算法**：QR 用于稳定求解方程、最小二乘和特征值计算。
- **傅里叶分析**：正弦和余弦构成正交函数系统，系数由投影得到。
- **机器学习**：特征正交化、线性模型、主成分分析都依赖本阶段思想。

## 17. 阶段练习

### 基础题

1. 设：

   $$
   \boldsymbol u=
   \begin{bmatrix}1\\2\\-1\end{bmatrix},
   \qquad
   \boldsymbol v=
   \begin{bmatrix}2\\0\\2\end{bmatrix}
   $$

   计算内积和两个向量的长度，并判断它们是否正交。

2. 求向量 $(3,4)^\mathsf T$ 的单位化结果。

3. 求向量：

   $$
   \boldsymbol b=
   \begin{bmatrix}3\\4\end{bmatrix}
   $$

   在：

   $$
   \boldsymbol a=
   \begin{bmatrix}2\\1\end{bmatrix}
   $$

   方向上的投影和残差。

4. 判断矩阵：

   $$
   P=
   \begin{bmatrix}
   1&0\\
   0&0
   \end{bmatrix}
   $$

   是否为正交投影矩阵，并说明投影到哪个子空间。

5. 设：

   $$
   \boldsymbol q_1=
   \frac1{\sqrt2}
   \begin{bmatrix}1\\1\\0\end{bmatrix},
   \qquad
   \boldsymbol q_2=
   \frac1{\sqrt6}
   \begin{bmatrix}1\\-1\\2\end{bmatrix}
   $$

   验证它们标准正交。

### 理解题

6. 为什么正交非零向量一定线性无关？

7. 为什么投影残差必须与目标子空间正交？

8. 证明投影矩阵：

   $$
   P=QQ^\mathsf T
   $$

   在 $Q^\mathsf TQ=I$ 时满足 $P^2=P$ 和 $P^\mathsf T=P$。

9. 为什么 $N(A)=C(A^\mathsf T)^\perp$？

10. 为什么最小二乘投影点总是唯一，而最小二乘参数不一定唯一？

11. 为什么使用 QR 求最小二乘通常比直接形成 $A^\mathsf TA$ 更稳定？

### 综合题

12. 使用 Gram–Schmidt 将下列向量正交标准化：

   $$
   \boldsymbol a_1=
   \begin{bmatrix}1\\1\\0\end{bmatrix},
   \qquad
   \boldsymbol a_2=
   \begin{bmatrix}1\\0\\1\end{bmatrix}
   $$

13. 求矩阵的经济型 QR 分解：

   $$
   A=
   \begin{bmatrix}
   1&1\\
   1&0\\
   0&1
   \end{bmatrix}
   $$

14. 将：

   $$
   \boldsymbol b=
   \begin{bmatrix}2\\1\\3\end{bmatrix}
   $$

   投影到第 13 题矩阵 $A$ 的列空间，并求残差。

15. 使用正规方程求下面超定系统的最小二乘解：

   $$
   \begin{bmatrix}
   1&0\\
   1&1\\
   1&2
   \end{bmatrix}
   \begin{bmatrix}x_1\\x_2\end{bmatrix}
   \approx
   \begin{bmatrix}1\\2\\2\end{bmatrix}
   $$

16. 对数据点：

   $$
   (0,2),\qquad(1,3),\qquad(2,5)
   $$

   拟合直线 $y=\beta_0+\beta_1x$，求最佳参数、预测值和残差。

17. 设：

   $$
   A=
   \begin{bmatrix}
   1&1\\
   1&-1\\
   1&1
   \end{bmatrix},
   \qquad
   \boldsymbol b=
   \begin{bmatrix}2\\0\\3\end{bmatrix}
   $$

   1. 判断 $A$ 的两列是否正交。
   2. 根据判断选择合适的方法求最小二乘解。
   3. 求投影与残差。

18. 设 $A$ 满列秩且 $A=QR$。从正规方程出发，推导：

   $$
   R\hat{\boldsymbol x}=Q^\mathsf T\boldsymbol b
   $$

## 18. 参考答案

1. 内积为：

   $$
   \boldsymbol u^\mathsf T\boldsymbol v
   =1\cdot2+2\cdot0+(-1)\cdot2
   =0
   $$

   长度为：

   $$
   \|\boldsymbol u\|=\sqrt6,
   \qquad
   \|\boldsymbol v\|=\sqrt8=2\sqrt2
   $$

   两个非零向量内积为零，所以正交。

2. 因为长度为 5，所以单位向量为：

   $$
   \begin{bmatrix}3/5\\4/5\end{bmatrix}
   $$

3. 投影系数为：

   $$
   c=\frac{\boldsymbol a^\mathsf T\boldsymbol b}
   {\boldsymbol a^\mathsf T\boldsymbol a}
   =\frac{10}{5}=2
   $$

   所以：

   $$
   \boldsymbol p=
   \begin{bmatrix}4\\2\end{bmatrix},
   \qquad
   \boldsymbol r=
   \begin{bmatrix}-1\\2\end{bmatrix}
   $$

4. 有：

   $$
   P^2=P,
   \qquad
   P^\mathsf T=P
   $$

   所以它是正交投影矩阵，将向量投影到 $x$ 轴。

5. 两个向量的长度都为 1，并且：

   $$
   \boldsymbol q_1^\mathsf T\boldsymbol q_2
   =
   \frac{1-1+0}{\sqrt{12}}
   =0
   $$

   所以它们标准正交。

6. 假设正交非零向量的线性组合为零。分别与每个向量做内积后，其他项都消失，只剩对应系数乘该向量长度平方，因此每个系数都必须为零。

7. 若残差还有沿目标子空间的分量，就可以沿该方向移动当前近似点并进一步缩短距离。只有残差完全位于正交补中时，当前点才是最近点。

8. 因为：

   $$
   P^2
   =QQ^\mathsf TQQ^\mathsf T
   =Q(Q^\mathsf TQ)Q^\mathsf T
   =QQ^\mathsf T
   =P
   $$

   并且：

   $$
   P^\mathsf T
   =(QQ^\mathsf T)^\mathsf T
   =QQ^\mathsf T
   =P
   $$

9. $\boldsymbol x\in N(A)$ 意味着 $A\boldsymbol x=\boldsymbol0$，即 $A$ 的每一行与 $\boldsymbol x$ 内积为零。因此 $\boldsymbol x$ 与整个行空间正交，反过来也成立。

10. 闭子空间中的正交投影点唯一。若参数不唯一，则两个参数之差位于 $N(A)$，它们仍然产生同一个投影点。因此输出唯一不代表输入表示唯一。

11. 形成 $A^\mathsf TA$ 会使条件数近似平方，并可能放大舍入误差。QR 主要使用保持长度的正交变换，避免这一额外放大。

12. 结果为：

   $$
   \boldsymbol q_1=
   \frac1{\sqrt2}
   \begin{bmatrix}1\\1\\0\end{bmatrix}
   $$

   $$
   \boldsymbol q_2=
   \frac1{\sqrt6}
   \begin{bmatrix}1\\-1\\2\end{bmatrix}
   $$

13. 经济型 QR 分解为：

   $$
   Q=
   \begin{bmatrix}
   \frac1{\sqrt2}&\frac1{\sqrt6}\\
   \frac1{\sqrt2}&-\frac1{\sqrt6}\\
   0&\frac2{\sqrt6}
   \end{bmatrix}
   $$

   $$
   R=
   \begin{bmatrix}
   \sqrt2&\frac1{\sqrt2}\\
   0&\sqrt{\frac32}
   \end{bmatrix}
   $$

14. 使用 $\boldsymbol p=QQ^\mathsf T\boldsymbol b$。投影坐标为：

   $$
   \boldsymbol q_1^\mathsf T\boldsymbol b
   =\frac3{\sqrt2}
   $$

   $$
   \boldsymbol q_2^\mathsf T\boldsymbol b
   =\frac7{\sqrt6}
   $$

   因此：

   $$
   \boldsymbol p
   =
   \frac3{\sqrt2}\boldsymbol q_1
   +
   \frac7{\sqrt6}\boldsymbol q_2
   =
   \begin{bmatrix}
   \frac83\\
   \frac13\\
   \frac73
   \end{bmatrix}
   $$

   残差为：

   $$
   \boldsymbol r
   =\boldsymbol b-\boldsymbol p
   =
   \begin{bmatrix}
   -\frac23\\
   \frac23\\
   \frac23
   \end{bmatrix}
   $$

   可以验证 $A^\mathsf T\boldsymbol r=\boldsymbol0$。

15. 正规方程为：

   $$
   \begin{bmatrix}
   3&3\\
   3&5
   \end{bmatrix}
   \begin{bmatrix}\hat x_1\\\hat x_2\end{bmatrix}
   =
   \begin{bmatrix}5\\6\end{bmatrix}
   $$

   解得：

   $$
   \hat x_1=\frac76,
   \qquad
   \hat x_2=\frac12
   $$

16. 设计矩阵仍为：

   $$
   A=
   \begin{bmatrix}
   1&0\\
   1&1\\
   1&2
   \end{bmatrix}
   $$

   观测向量为：

   $$
   \boldsymbol b=
   \begin{bmatrix}2\\3\\5\end{bmatrix}
   $$

   有：

   $$
   A^\mathsf T\boldsymbol b=
   \begin{bmatrix}10\\13\end{bmatrix}
   $$

   正规方程：

   $$
   \begin{bmatrix}
   3&3\\
   3&5
   \end{bmatrix}
   \begin{bmatrix}\hat\beta_0\\\hat\beta_1\end{bmatrix}
   =
   \begin{bmatrix}10\\13\end{bmatrix}
   $$

   解得：

   $$
   \hat\beta_0=\frac{11}{6},
   \qquad
   \hat\beta_1=\frac32
   $$

   预测值为：

   $$
   \boldsymbol p=
   \begin{bmatrix}
   \frac{11}{6}\\
   \frac{10}{3}\\
   \frac{29}{6}
   \end{bmatrix}
   $$

   残差为：

   $$
   \boldsymbol r=
   \boldsymbol b-\boldsymbol p
   =
   \begin{bmatrix}
   \frac16\\
   -\frac13\\
   \frac16
   \end{bmatrix}
   $$

17. 两列为：

   $$
   \boldsymbol a_1=
   \begin{bmatrix}1\\1\\1\end{bmatrix},
   \qquad
   \boldsymbol a_2=
   \begin{bmatrix}1\\-1\\1\end{bmatrix}
   $$

   这里需要注意：

   $$
   \boldsymbol a_1^\mathsf T\boldsymbol a_2
   =1-1+1
   =1
   $$

   因此题目中的两列实际上并不正交，不能直接使用正交列公式。这一检查本身就是题目要发现的问题。应改用正规方程：

   $$
   A^\mathsf TA=
   \begin{bmatrix}3&1\\1&3\end{bmatrix},
   \qquad
   A^\mathsf T\boldsymbol b=
   \begin{bmatrix}5\\5\end{bmatrix}
   $$

   解得：

   $$
   \hat{\boldsymbol x}=
   \begin{bmatrix}5/4\\5/4\end{bmatrix}
   $$

   投影为：

   $$
   \boldsymbol p=
   A\hat{\boldsymbol x}
   =
   \begin{bmatrix}5/2\\0\\5/2\end{bmatrix}
   $$

   残差为：

   $$
   \boldsymbol r=
   \boldsymbol b-\boldsymbol p
   =
   \begin{bmatrix}-1/2\\0\\1/2\end{bmatrix}
   $$

18. 代入 $A=QR$：

   $$
   A^\mathsf TA\hat{\boldsymbol x}
   =A^\mathsf T\boldsymbol b
   $$

   得：

   $$
   R^\mathsf TQ^\mathsf TQR\hat{\boldsymbol x}
   =R^\mathsf TQ^\mathsf T\boldsymbol b
   $$

   因为 $Q^\mathsf TQ=I$ 且满列秩时 $R$ 可逆：

   $$
   R^\mathsf TR\hat{\boldsymbol x}
   =R^\mathsf TQ^\mathsf T\boldsymbol b
   $$

   左乘 $(R^\mathsf T)^{-1}$，得到：

   $$
   R\hat{\boldsymbol x}
   =Q^\mathsf T\boldsymbol b
   $$

## 19. 阶段检验

在不查阅资料的情况下，回答以下问题：

1. 内积如何同时描述代数运算和几何方向关系？
2. 范数、距离和夹角如何由内积得到？
3. 为什么正交非零向量组一定线性无关？
4. 什么是正交补？四个基本子空间如何组成两组正交补？
5. 如何推导向量到一条直线的投影公式？
6. 为什么投影点是子空间中的唯一最近点？
7. 正交投影矩阵为什么同时满足 $P^2=P$ 和 $P^\mathsf T=P$？
8. Gram–Schmidt 每一步减去的是什么？为什么张成空间不变？
9. QR 分解中的 $Q$ 和 $R$ 分别表示什么？
10. 为什么最小二乘等价于将 $\boldsymbol b$ 投影到 $C(A)$？
11. 如何从残差正交条件推导正规方程？
12. 最小二乘投影点与最小二乘参数的唯一性有何区别？
13. 为什么 QR 能将最小二乘问题化为上三角方程？
14. 为什么实际计算更偏好 QR，而不是显式形成 $(A^\mathsf TA)^{-1}$？
15. 能否从一组数据建立设计矩阵，并解释每一列对应的模型基函数？

能够清楚回答这些问题，独立完成阶段练习，并从几何上解释投影与残差，就可以进入第五阶段：行列式、特征值、对角化与动态系统。
