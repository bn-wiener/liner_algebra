# 第六阶段：SVD、伪逆、低秩近似与数据分析

> 核心问题：任意矩阵最本质的输入方向和输出方向是什么？矩阵真正保留了多少信息？当矩阵不可逆或数据维度过高时，怎样得到稳定、紧凑而有意义的表示？

前五个阶段分别研究了向量空间、线性变换、方程求解、正交投影和特征结构。本阶段使用奇异值分解将这些主题统一起来。

对于任意实矩阵：

$$
A\in\mathbb R^{m\times n}
$$

都存在分解：

$$
\boxed{
A=U\Sigma V^\mathsf T
}
$$

其几何过程是：

$$
\text{输入空间旋转或反射}
\xrightarrow{V^\mathsf T}
\text{沿正交方向缩放}
\xrightarrow{\Sigma}
\text{输出空间旋转或反射}
\xrightarrow{U}
$$

SVD 不要求矩阵是方阵、对称矩阵或可逆矩阵，因此它比普通特征分解适用范围更广。

本阶段的知识链为：

$$
\text{奇异值与奇异向量}
\rightarrow
\text{SVD}
\rightarrow
\text{四个基本子空间}
\rightarrow
\text{伪逆}
\rightarrow
\text{稳定最小二乘}
\rightarrow
\text{低秩近似}
\rightarrow
\text{PCA 与数据压缩}
$$

## 学习目标

完成本阶段后，应当能够：

1. 解释为什么任意矩阵都存在 SVD。
2. 根据 $A^\mathsf TA$ 或 $AA^\mathsf T$ 求奇异值和奇异向量。
3. 区分完整 SVD、经济型 SVD 和紧致 SVD。
4. 从旋转、缩放和投影角度解释 SVD 的几何意义。
5. 使用 SVD 直接读出矩阵的秩和四个基本子空间。
6. 理解外积展开 $A=\sum_i\sigma_i\boldsymbol u_i\boldsymbol v_i^\mathsf T$。
7. 构造 Moore–Penrose 伪逆。
8. 使用伪逆求精确解、最小二乘解和最小范数解。
9. 理解 $AA^+$ 与 $A^+A$ 的投影意义。
10. 使用奇异值分析矩阵的条件数和数值敏感性。
11. 理解 Eckart–Young 最佳低秩近似定理。
12. 使用截断 SVD 进行压缩、去噪和降维。
13. 从 SVD 推导 PCA 的主成分方向、得分和解释方差。
14. 比较特征分解、QR 和 SVD 的适用范围。

## 1. 为什么需要 SVD

### 1.1 特征分解的限制

普通特征分解研究：

$$
A\boldsymbol v=\lambda\boldsymbol v
$$

它要求输入和输出位于同一个空间，因此 $A$ 必须是方阵。

即使 $A$ 是方阵，也可能：

- 没有足够多的线性无关特征向量。
- 在实数范围内没有实特征值。
- 特征向量不正交。
- 对微小扰动非常敏感。

### 1.2 SVD 的优势

SVD 对任意 $m\times n$ 实矩阵都存在，并且提供：

- 输入空间中的标准正交方向。
- 输出空间中的标准正交方向。
- 每对方向之间的非负缩放强度。
- 秩和四个基本子空间。
- 稳定的最小二乘求解。
- 最佳低秩近似。

### 1.3 两组空间

对于：

$$
A:\mathbb R^n\rightarrow\mathbb R^m
$$

SVD 不要求某个输入方向在变换后仍位于原来的直线上，而是寻找：

$$
A\boldsymbol v_i
=
\sigma_i\boldsymbol u_i
$$

其中：

- $\boldsymbol v_i\in\mathbb R^n$ 是右奇异向量，位于输入空间。
- $\boldsymbol u_i\in\mathbb R^m$ 是左奇异向量，位于输出空间。
- $\sigma_i\ge0$ 是奇异值。

## 2. 从 $A^\mathsf TA$ 出发

### 2.1 对称半正定矩阵

对任意矩阵 $A$：

$$
A^\mathsf TA
$$

都是 $n\times n$ 对称矩阵，并且：

$$
\boldsymbol x^\mathsf TA^\mathsf TA\boldsymbol x
=
\|A\boldsymbol x\|^2
\ge0
$$

所以 $A^\mathsf TA$ 是半正定矩阵。

根据谱定理，它拥有一组标准正交特征向量：

$$
\boldsymbol v_1,\dots,\boldsymbol v_n
$$

以及非负特征值：

$$
\lambda_1,\dots,\lambda_n
$$

### 2.2 奇异值

奇异值定义为：

$$
\boxed{
\sigma_i=\sqrt{\lambda_i(A^\mathsf TA)}
}
$$

通常按降序排列：

$$
\sigma_1\ge\sigma_2\ge\cdots\ge\sigma_p\ge0
$$

其中：

$$
p=\min(m,n)
$$

### 2.3 右奇异向量

$A^\mathsf TA$ 的标准正交特征向量就是右奇异向量：

$$
A^\mathsf TA\boldsymbol v_i
=
\sigma_i^2\boldsymbol v_i
$$

它们给出输入空间中矩阵作用最自然的正交方向。

### 2.4 左奇异向量

当 $\sigma_i>0$ 时，定义：

$$
\boxed{
\boldsymbol u_i
=
\frac{A\boldsymbol v_i}{\sigma_i}
}
$$

因此：

$$
A\boldsymbol v_i
=
\sigma_i\boldsymbol u_i
$$

并且：

$$
A^\mathsf T\boldsymbol u_i
=
\sigma_i\boldsymbol v_i
$$

左奇异向量也是 $AA^\mathsf T$ 的特征向量：

$$
AA^\mathsf T\boldsymbol u_i
=
\sigma_i^2\boldsymbol u_i
$$

### 2.5 两个对称矩阵的共同非零谱

$A^\mathsf TA$ 与 $AA^\mathsf T$ 尺寸可能不同，但具有相同的非零特征值：

$$
\sigma_1^2,\dots,\sigma_r^2
$$

其中：

$$
r=\operatorname{rank}(A)
$$

## 3. SVD 的矩阵形式

### 3.1 完整 SVD

完整 SVD 写成：

$$
A=U\Sigma V^\mathsf T
$$

其中：

$$
U\in\mathbb R^{m\times m},
\qquad
V\in\mathbb R^{n\times n}
$$

都是正交矩阵：

$$
U^\mathsf TU=UU^\mathsf T=I_m
$$

$$
V^\mathsf TV=VV^\mathsf T=I_n
$$

$\Sigma$ 是 $m\times n$ 的矩形对角矩阵：

$$
\Sigma=
\begin{bmatrix}
\sigma_1&&&&\\
&\sigma_2&&&\\
&&\ddots&&\\
&&&\sigma_p&\\
&&&&0
\end{bmatrix}
$$

### 3.2 经济型 SVD

若只保留 $p=\min(m,n)$ 个可能的奇异方向，可写成经济型 SVD。对于 $m\ge n$：

$$
A=U_p\Sigma_pV^\mathsf T
$$

其中：

$$
U_p\in\mathbb R^{m\times n},
\qquad
\Sigma_p\in\mathbb R^{n\times n}
$$

### 3.3 紧致 SVD

若矩阵秩为 $r$，只保留非零奇异值：

$$
\boxed{
A=U_r\Sigma_rV_r^\mathsf T
}
$$

其中：

$$
U_r\in\mathbb R^{m\times r}
$$

$$
\Sigma_r\in\mathbb R^{r\times r}
$$

$$
V_r\in\mathbb R^{n\times r}
$$

紧致 SVD 最直接地暴露矩阵真正有效的 $r$ 个方向。

### 3.4 SVD 不唯一

- 奇异向量可以同时改变符号，不影响分解。
- 重复奇异值对应的子空间内可以选择不同标准正交基。
- 零奇异值对应的左右补空间基也有多种选择。

奇异值本身及相关子空间是确定的，但具体基向量可能不唯一。

## 4. SVD 的几何意义

对于输入 $\boldsymbol x$：

$$
A\boldsymbol x
=
U\Sigma V^\mathsf T\boldsymbol x
$$

可以分成三步。

### 4.1 第一步：$V^\mathsf T$

$V^\mathsf T$ 将输入向量转换到右奇异向量组成的标准正交坐标系。

因为 $V$ 正交，这一步只旋转或反射，不改变长度和角度。

### 4.2 第二步：$\Sigma$

$\Sigma$ 沿各坐标轴分别缩放：

$$
x_i\mapsto\sigma_ix_i
$$

- 大奇异值方向被强烈放大。
- 小奇异值方向被弱化。
- 零奇异值方向被完全消除。

### 4.3 第三步：$U$

$U$ 将缩放后的结果旋转或反射到输出空间中的左奇异向量方向。

### 4.4 单位球的像

输入空间中的单位球：

$$
\|\boldsymbol x\|=1
$$

经过 $A$ 后变成输出空间中的椭球：

- 椭球主轴方向是 $\boldsymbol u_i$。
- 主轴长度是 $\sigma_i$。
- 对应输入方向是 $\boldsymbol v_i$。

因此 SVD 是任意线性变换的主轴描述。

## 5. 手算 SVD 的基本步骤

对于 $A\in\mathbb R^{m\times n}$：

1. 计算 $A^\mathsf TA$。
2. 求 $A^\mathsf TA$ 的特征值和标准正交特征向量。
3. 对特征值开平方得到奇异值。
4. 按奇异值从大到小排列右奇异向量。
5. 对每个非零奇异值计算：

   $$
   \boldsymbol u_i=\frac{A\boldsymbol v_i}{\sigma_i}
   $$

6. 补充左零空间中的标准正交向量，得到完整 $U$。
7. 组织 $U,\Sigma,V$ 并验证：

   $$
   A=U\Sigma V^\mathsf T
   $$

手算中也可以从 $AA^\mathsf T$ 开始，特别是在 $m<n$ 时，选择尺寸较小的对称矩阵通常更方便。

## 6. 贯穿案例：一个秩一矩阵的完整 SVD

考虑：

$$
A=
\begin{bmatrix}
1&1\\
1&1\\
0&0
\end{bmatrix}
\in\mathbb R^{3\times2}
$$

### 6.1 计算 $A^\mathsf TA$

$$
A^\mathsf TA
=
\begin{bmatrix}
2&2\\
2&2
\end{bmatrix}
$$

其特征值为：

$$
\lambda_1=4,
\qquad
\lambda_2=0
$$

所以奇异值为：

$$
\sigma_1=2,
\qquad
\sigma_2=0
$$

### 6.2 右奇异向量

对 $\lambda_1=4$，取：

$$
\boldsymbol v_1
=
\frac1{\sqrt2}
\begin{bmatrix}1\\1\end{bmatrix}
$$

对 $\lambda_2=0$，取：

$$
\boldsymbol v_2
=
\frac1{\sqrt2}
\begin{bmatrix}1\\-1\end{bmatrix}
$$

因此：

$$
V=
\frac1{\sqrt2}
\begin{bmatrix}
1&1\\
1&-1
\end{bmatrix}
$$

### 6.3 左奇异向量

对非零奇异值：

$$
\boldsymbol u_1
=
\frac{A\boldsymbol v_1}{2}
$$

因为：

$$
A\boldsymbol v_1
=
\begin{bmatrix}
\sqrt2\\
\sqrt2\\
0
\end{bmatrix}
$$

所以：

$$
\boldsymbol u_1
=
\frac1{\sqrt2}
\begin{bmatrix}1\\1\\0\end{bmatrix}
$$

再选择与 $\boldsymbol u_1$ 正交的单位向量：

$$
\boldsymbol u_2
=
\frac1{\sqrt2}
\begin{bmatrix}1\\-1\\0\end{bmatrix}
$$

$$
\boldsymbol u_3
=
\begin{bmatrix}0\\0\\1\end{bmatrix}
$$

于是：

$$
U=
\begin{bmatrix}
\frac1{\sqrt2}&\frac1{\sqrt2}&0\\
\frac1{\sqrt2}&-\frac1{\sqrt2}&0\\
0&0&1
\end{bmatrix}
$$

### 6.4 奇异值矩阵

$$
\Sigma=
\begin{bmatrix}
2&0\\
0&0\\
0&0
\end{bmatrix}
$$

因此完整 SVD 为：

$$
\boxed{
A=U\Sigma V^\mathsf T
}
$$

### 6.5 紧致 SVD

因为秩为 1，只保留非零奇异方向：

$$
A
=
\boldsymbol u_1
\begin{bmatrix}2\end{bmatrix}
\boldsymbol v_1^\mathsf T
$$

即：

$$
A
=
2
\left(
\frac1{\sqrt2}
\begin{bmatrix}1\\1\\0\end{bmatrix}
\right)
\left(
\frac1{\sqrt2}
\begin{bmatrix}1&1\end{bmatrix}
\right)
$$

## 7. SVD 与秩

矩阵的秩等于非零奇异值的数量：

$$
\boxed{
\operatorname{rank}(A)
=
\#\{i:\sigma_i>0\}
}
$$

原因是：

- 正交矩阵 $U,V$ 不改变秩。
- $\Sigma$ 的秩就是非零对角元素数量。

对于贯穿案例，只有 $\sigma_1=2$ 非零，所以：

$$
\operatorname{rank}(A)=1
$$

### 7.1 数值秩

在浮点计算中，理论上为零的奇异值可能表现为非常小的正数。因此通常设置阈值：

$$
\sigma_i>\text{tolerance}
$$

才视为有效奇异值。

数值秩依赖于数据尺度、误差水平和应用目标，不应机械地将所有非零浮点数都视为有效方向。

## 8. SVD 与四个基本子空间

设：

$$
\sigma_1\ge\cdots\ge\sigma_r>0
$$

其余奇异值为零。

### 8.1 列空间

$$
\boxed{
C(A)
=
\operatorname{span}
(\boldsymbol u_1,\dots,\boldsymbol u_r)
}
$$

非零奇异值对应的左奇异向量张成所有可达到输出。

### 8.2 左零空间

$$
\boxed{
N(A^\mathsf T)
=
\operatorname{span}
(\boldsymbol u_{r+1},\dots,\boldsymbol u_m)
}
$$

### 8.3 行空间

$$
\boxed{
C(A^\mathsf T)
=
\operatorname{span}
(\boldsymbol v_1,\dots,\boldsymbol v_r)
}
$$

### 8.4 零空间

$$
\boxed{
N(A)
=
\operatorname{span}
(\boldsymbol v_{r+1},\dots,\boldsymbol v_n)
}
$$

### 8.5 贯穿案例

对秩一矩阵：

$$
C(A)
=
\operatorname{span}
\left(
\frac1{\sqrt2}
\begin{bmatrix}1\\1\\0\end{bmatrix}
\right)
$$

$$
N(A^\mathsf T)
=
\operatorname{span}
\left(
\frac1{\sqrt2}
\begin{bmatrix}1\\-1\\0\end{bmatrix},
\begin{bmatrix}0\\0\\1\end{bmatrix}
\right)
$$

$$
C(A^\mathsf T)
=
\operatorname{span}
\left(
\frac1{\sqrt2}
\begin{bmatrix}1\\1\end{bmatrix}
\right)
$$

$$
N(A)
=
\operatorname{span}
\left(
\frac1{\sqrt2}
\begin{bmatrix}1\\-1\end{bmatrix}
\right)
$$

SVD 不需要再次消元，就同时给出了四个基本子空间的标准正交基。

## 9. 外积展开

紧致 SVD 可以写成：

$$
\boxed{
A
=
\sum_{i=1}^r
\sigma_i\boldsymbol u_i\boldsymbol v_i^\mathsf T
}
$$

每一项：

$$
\sigma_i\boldsymbol u_i\boldsymbol v_i^\mathsf T
$$

都是秩一矩阵。

它的作用是：

1. 使用 $\boldsymbol v_i^\mathsf T\boldsymbol x$ 测量输入在 $\boldsymbol v_i$ 方向上的分量。
2. 乘以 $\sigma_i$ 进行缩放。
3. 沿输出方向 $\boldsymbol u_i$ 产生结果。

因此任意秩为 $r$ 的矩阵都是 $r$ 个正交秩一模式之和。

对于贯穿案例：

$$
A=2\boldsymbol u_1\boldsymbol v_1^\mathsf T
$$

只有一个有效模式。

## 10. Moore–Penrose 伪逆

### 10.1 为什么需要伪逆

普通逆矩阵要求矩阵是可逆方阵。但实际问题中的矩阵经常：

- 不是方阵。
- 秩不足。
- 方程无精确解。
- 方程有无穷多个解。

伪逆将逆矩阵的思想推广到任意矩阵。

### 10.2 从 SVD 定义

若：

$$
A=U\Sigma V^\mathsf T
$$

则 Moore–Penrose 伪逆定义为：

$$
\boxed{
A^+
=
V\Sigma^+U^\mathsf T
}
$$

$\Sigma^+$ 的构造方法是：

- 将每个非零奇异值 $\sigma_i$ 替换为 $1/\sigma_i$。
- 零奇异值仍保持为零。
- 将矩形尺寸转置。

紧致形式为：

$$
A^+
=
V_r\Sigma_r^{-1}U_r^\mathsf T
$$

### 10.3 贯穿案例的伪逆

因为：

$$
\sigma_1=2
$$

所以：

$$
A^+
=
\boldsymbol v_1
\frac12
\boldsymbol u_1^\mathsf T
$$

计算得：

$$
\boxed{
A^+
=
\frac14
\begin{bmatrix}
1&1&0\\
1&1&0
\end{bmatrix}
}
$$

### 10.4 Moore–Penrose 四个条件

伪逆是唯一满足下列条件的矩阵：

$$
AA^+A=A
$$

$$
A^+AA^+=A^+
$$

$$
(AA^+)^\mathsf T=AA^+
$$

$$
(A^+A)^\mathsf T=A^+A
$$

前两个条件表达广义逆关系，后两个条件保证相关投影是正交投影。

### 10.5 可逆矩阵的情况

若 $A$ 是可逆方阵，则所有奇异值非零，并且：

$$
A^+=A^{-1}
$$

所以伪逆确实是普通逆矩阵的推广。

## 11. 伪逆与投影

### 11.1 输出空间投影

$$
\boxed{
AA^+
}
$$

是投影到列空间 $C(A)$ 的正交投影矩阵。

使用紧致 SVD：

$$
AA^+
=
U_r\Sigma_rV_r^\mathsf T
V_r\Sigma_r^{-1}U_r^\mathsf T
=
U_rU_r^\mathsf T
$$

### 11.2 输入空间投影

$$
\boxed{
A^+A
}
$$

是投影到行空间 $C(A^\mathsf T)$ 的正交投影矩阵：

$$
A^+A
=
V_rV_r^\mathsf T
$$

### 11.3 贯穿案例

$$
AA^+
=
\frac12
\begin{bmatrix}
1&1&0\\
1&1&0\\
0&0&0
\end{bmatrix}
$$

它把 $\mathbb R^3$ 投影到：

$$
\operatorname{span}
\left(
\begin{bmatrix}1\\1\\0\end{bmatrix}
\right)
$$

而：

$$
A^+A
=
\frac12
\begin{bmatrix}
1&1\\
1&1
\end{bmatrix}
$$

它把 $\mathbb R^2$ 投影到：

$$
\operatorname{span}
\left(
\begin{bmatrix}1\\1\end{bmatrix}
\right)
$$

## 12. 伪逆与最小二乘

对于任意矩阵和目标向量：

$$
\boxed{
\hat{\boldsymbol x}=A^+\boldsymbol b
}
$$

给出一个最小二乘解。

它具有更强的性质：

> $A^+\boldsymbol b$ 是所有最小二乘解中欧几里得范数最小的解。

### 12.1 SVD 形式

$$
A^+\boldsymbol b
=
\sum_{i=1}^r
\frac{\boldsymbol u_i^\mathsf T\boldsymbol b}{\sigma_i}
\boldsymbol v_i
$$

过程是：

1. 将 $\boldsymbol b$ 投影到各左奇异方向。
2. 除以对应非零奇异值。
3. 沿右奇异方向恢复输入系数。

### 12.2 贯穿案例

取：

$$
\boldsymbol b=
\begin{bmatrix}2\\0\\1\end{bmatrix}
$$

则：

$$
\hat{\boldsymbol x}
=
A^+\boldsymbol b
=
\frac14
\begin{bmatrix}
1&1&0\\
1&1&0
\end{bmatrix}
\begin{bmatrix}2\\0\\1\end{bmatrix}
$$

所以：

$$
\boxed{
\hat{\boldsymbol x}
=
\begin{bmatrix}1/2\\1/2\end{bmatrix}
}
$$

对应预测：

$$
\boldsymbol p
=
A\hat{\boldsymbol x}
=
\begin{bmatrix}1\\1\\0\end{bmatrix}
$$

残差：

$$
\boldsymbol r
=
\boldsymbol b-\boldsymbol p
=
\begin{bmatrix}1\\-1\\1\end{bmatrix}
$$

残差属于左零空间：

$$
\boldsymbol r\in N(A^\mathsf T)
$$

### 12.3 所有最小二乘解

零空间为：

$$
N(A)
=
\operatorname{span}
\left(
\begin{bmatrix}1\\-1\end{bmatrix}
\right)
$$

因此所有产生同一最佳预测的参数为：

$$
\boldsymbol x
=
\begin{bmatrix}1/2\\1/2\end{bmatrix}
+t
\begin{bmatrix}1\\-1\end{bmatrix}
$$

其中 $t=0$ 时范数最小，因为伪逆解位于行空间，与零空间正交。

## 13. 精确解、最小二乘解与最小范数解

伪逆统一处理多种情况。

### 13.1 可逆方阵

若 $A$ 可逆：

$$
A^+\boldsymbol b=A^{-1}\boldsymbol b
$$

得到唯一精确解。

### 13.2 满列秩高矩阵

若 $m>n$ 且 $A$ 满列秩：

$$
A^+
=(A^\mathsf TA)^{-1}A^\mathsf T
$$

得到唯一最小二乘解。

### 13.3 满行秩宽矩阵

若 $m<n$ 且 $A$ 满行秩：

$$
A^+
=
A^\mathsf T(AA^\mathsf T)^{-1}
$$

每个右端都有精确解，伪逆选择其中范数最小的一个。

### 13.4 一般秩亏矩阵

SVD 公式：

$$
A^+=V_r\Sigma_r^{-1}U_r^\mathsf T
$$

仍然有效，同时处理无解、无穷多解和冗余方向。

## 14. 奇异值与矩阵范数

### 14.1 谱范数

矩阵的二范数定义为：

$$
\|A\|_2
=
\max_{\|\boldsymbol x\|=1}
\|A\boldsymbol x\|
$$

SVD 给出：

$$
\boxed{
\|A\|_2=\sigma_1
}
$$

最大奇异值是矩阵对单位向量能够产生的最大放大倍数。

### 14.2 最小放大倍数

若 $A$ 满列秩，则：

$$
\min_{\|\boldsymbol x\|=1}
\|A\boldsymbol x\|
=
\sigma_{\min}
$$

最小奇异值越接近零，说明某个输入方向几乎被压扁。

### 14.3 Frobenius 范数

$$
\|A\|_F
=
\sqrt{\sum_{i,j}a_{ij}^2}
$$

SVD 给出：

$$
\boxed{
\|A\|_F^2
=
\sum_i\sigma_i^2
}
$$

因此矩阵总平方能量等于所有奇异值平方之和。

## 15. 条件数与数值敏感性

### 15.1 二范数条件数

对于满秩矩阵：

$$
\boxed{
\kappa_2(A)
=
\frac{\sigma_{\max}}{\sigma_{\min}}
}
$$

若矩阵秩亏，则：

$$
\kappa_2(A)=\infty
$$

### 15.2 几何意义

单位球经过 $A$ 后变成椭球。条件数等于最长主轴与最短非零主轴之比。

- $\kappa$ 接近 1：各方向缩放接近，问题较稳定。
- $\kappa$ 很大：某些方向几乎被压扁，逆问题敏感。
- $\kappa=\infty$：存在完全丢失方向，无法唯一求逆。

### 15.3 小奇异值为什么危险

伪逆中包含：

$$
\frac1{\sigma_i}
$$

当 $\sigma_i$ 很小时，目标向量在 $\boldsymbol u_i$ 方向上的微小噪声会被巨大放大到输入解中。

SVD 不仅告诉我们答案，还揭示哪些方向使答案不可靠。

### 15.4 例

设：

$$
A=
\begin{bmatrix}
1&0\\
0&0.001
\end{bmatrix}
$$

奇异值为：

$$
1,\qquad0.001
$$

所以：

$$
\kappa_2(A)=1000
$$

第二方向上的输出误差在求逆时可能被放大约 1000 倍。

## 16. 截断 SVD 与稳定化

### 16.1 基本思想

若非常小的奇异值主要代表噪声，可以只保留较大的奇异值：

$$
A_k
=
\sum_{i=1}^k
\sigma_i\boldsymbol u_i\boldsymbol v_i^\mathsf T
$$

对应截断伪逆：

$$
A_k^+
=
\sum_{i=1}^k
\frac1{\sigma_i}
\boldsymbol v_i\boldsymbol u_i^\mathsf T
$$

这会牺牲一部分精确性，换取更好的稳定性。

### 16.2 信号与噪声

常见解释是：

- 大奇异值对应稳定、主要的结构。
- 小奇异值对应弱方向，可能包含噪声。

但不能仅凭大小自动判断语义。阈值应结合噪声水平、验证误差和具体应用。

### 16.3 正则化

另一种方法是 Tikhonov 或岭正则化：

$$
\min_{\boldsymbol x}
\left(
\|A\boldsymbol x-\boldsymbol b\|^2
+\lambda\|\boldsymbol x\|^2
\right)
$$

在奇异方向上，其缩放因子从：

$$
\frac1{\sigma_i}
$$

变为：

$$
\frac{\sigma_i}{\sigma_i^2+\lambda}
$$

因此小奇异值方向不会被无限放大。

## 17. 最佳低秩近似

### 17.1 问题

希望使用秩不超过 $k$ 的矩阵近似 $A$：

$$
\operatorname{rank}(B)\le k
$$

并使误差尽可能小。

### 17.2 截断 SVD

保留前 $k$ 个最大奇异值：

$$
\boxed{
A_k
=
\sum_{i=1}^k
\sigma_i\boldsymbol u_i\boldsymbol v_i^\mathsf T
}
$$

### 17.3 Eckart–Young 定理

$A_k$ 是所有秩不超过 $k$ 的矩阵中最佳近似。

在谱范数下：

$$
\boxed{
\|A-A_k\|_2
=
\sigma_{k+1}
}
$$

在 Frobenius 范数下：

$$
\boxed{
\|A-A_k\|_F
=
\sqrt{
\sigma_{k+1}^2+\cdots+\sigma_r^2
}
}
$$

### 17.4 例

若矩阵奇异值为：

$$
5,\qquad2,\qquad0.5
$$

最佳秩一近似的误差为：

$$
\|A-A_1\|_2=2
$$

$$
\|A-A_1\|_F
=
\sqrt{2^2+0.5^2}
=
\frac{\sqrt{17}}2
$$

最佳秩二近似的误差为：

$$
\|A-A_2\|_2
=
\|A-A_2\|_F
=0.5
$$

## 18. 低秩表示与压缩

### 18.1 存储量

一个 $m\times n$ 矩阵需要存储：

$$
mn
$$

个数。

秩 $k$ 截断 SVD 只需存储：

- $U_k$：$mk$ 个数。
- $\Sigma_k$：$k$ 个数。
- $V_k$：$nk$ 个数。

总计：

$$
k(m+n+1)
$$

当：

$$
k\ll\min(m,n)
$$

时可显著压缩。

### 18.2 图像压缩

灰度图像可以表示为像素矩阵 $A$。截断 SVD 使用少数秩一模式近似图像：

$$
A_k
=
\sigma_1\boldsymbol u_1\boldsymbol v_1^\mathsf T
+\cdots+
\sigma_k\boldsymbol u_k\boldsymbol v_k^\mathsf T
$$

较大的奇异值通常保留主要亮度、边缘和大尺度结构；较小模式补充细节和噪声。

### 18.3 压缩率与质量

$k$ 越小：

- 存储更少。
- 计算更快。
- 细节损失更多。

$k$ 越大：

- 近似更精确。
- 存储成本更高。

选择 $k$ 是压缩率与重建质量之间的权衡。

## 19. PCA：寻找数据的主要变化方向

### 19.1 数据矩阵

设有 $m$ 个样本、$n$ 个特征，将数据组织为：

$$
X\in\mathbb R^{m\times n}
$$

每一行是一个样本，每一列是一个特征。

### 19.2 中心化

PCA 前必须减去每个特征的均值：

$$
X_c
=
X-\boldsymbol1\boldsymbol\mu^\mathsf T
$$

中心化使数据围绕原点变化。若不中心化，第一主成分可能主要描述均值位置，而不是样本变化方向。

### 19.3 协方差矩阵

样本协方差矩阵为：

$$
C
=
\frac1{m-1}
X_c^\mathsf TX_c
$$

它是对称半正定矩阵。

若：

$$
X_c=U\Sigma V^\mathsf T
$$

则：

$$
C
=
V
\frac{\Sigma^\mathsf T\Sigma}{m-1}
V^\mathsf T
$$

所以：

- $V$ 的列是主成分方向。
- 第 $i$ 个主成分方差为：

  $$
  \frac{\sigma_i^2}{m-1}
  $$

### 19.4 主成分得分

将数据投影到前 $k$ 个主成分：

$$
Z=X_cV_k
$$

利用 SVD：

$$
Z=U_k\Sigma_k
$$

$Z$ 是降维后的样本表示。

### 19.5 重建

从低维表示近似重建：

$$
\hat X_c
=
ZV_k^\mathsf T
=
U_k\Sigma_kV_k^\mathsf T
$$

再加回均值：

$$
\hat X
=
\hat X_c
+\boldsymbol1\boldsymbol\mu^\mathsf T
$$

### 19.6 解释方差比

前 $k$ 个主成分保留的方差比例为：

$$
\frac{
\sigma_1^2+\cdots+\sigma_k^2
}{
\sigma_1^2+\cdots+\sigma_r^2
}
$$

该比例常用于选择降维维数，但还应结合任务效果和可解释性。

## 20. PCA 完整例子

考虑已经中心化的数据：

$$
X=
\begin{bmatrix}
-2&-4\\
-1&-2\\
0&0\\
1&2\\
2&4
\end{bmatrix}
$$

所有样本都位于直线：

$$
y=2x
$$

### 20.1 秩与主方向

矩阵可以写成外积：

$$
X=
\begin{bmatrix}
-2\\-1\\0\\1\\2
\end{bmatrix}
\begin{bmatrix}1&2\end{bmatrix}
$$

所以秩为 1。

第一右奇异向量为：

$$
\boldsymbol v_1
=
\frac1{\sqrt5}
\begin{bmatrix}1\\2\end{bmatrix}
$$

它正是数据直线方向。

第二方向可取：

$$
\boldsymbol v_2
=
\frac1{\sqrt5}
\begin{bmatrix}-2\\1\end{bmatrix}
$$

该方向与数据直线正交。

### 20.2 奇异值

两个外积向量长度分别为：

$$
\sqrt{10},
\qquad
\sqrt5
$$

所以唯一非零奇异值为：

$$
\sigma_1=\sqrt{50}=5\sqrt2
$$

第二奇异值为零。

### 20.3 方差

样本数量 $m=5$，第一主成分方差为：

$$
\frac{\sigma_1^2}{m-1}
=
\frac{50}{4}
=12.5
$$

第二主成分方差为零。

### 20.4 降维

只保留第一主成分不会损失任何信息，因为数据本来就位于一条直线上。

每个样本由一个标量得分表示：

$$
Z=X\boldsymbol v_1
$$

这将二维数据无损压缩为一维数据。

## 21. 去噪与潜在结构

现实数据不会严格位于低维子空间中。通常：

$$
X
=
X_{\text{signal}}
+X_{\text{noise}}
$$

若信号具有近似低秩结构，而噪声分散在许多方向，则：

- 大奇异值捕获主要结构。
- 小奇异值捕获细节和部分噪声。

截断 SVD：

$$
X_k=U_k\Sigma_kV_k^\mathsf T
$$

可以作为去噪结果。

但必须注意：

- 有意义的弱信号也可能对应小奇异值。
- 结构化噪声也可能产生大奇异值。
- 选择 $k$ 需要结合领域知识和验证。

## 22. SVD、特征分解与 QR 的比较

| 方法 | 适用矩阵 | 核心结构 | 典型用途 |
|---|---|---|---|
| 特征分解 | 方阵；最好是对称矩阵 | 不变方向 | 动态系统、谱分析 |
| 对称特征分解 | 实对称矩阵 | 正交特征方向 | 二次型、PCA、优化 |
| QR 分解 | 一般矩阵 | 正交列空间基与三角坐标 | 最小二乘、线性求解 |
| SVD | 任意矩阵 | 输入输出主方向与缩放 | 伪逆、压缩、降维 |

### 22.1 特征分解与 SVD

若 $A$ 是实对称半正定矩阵：

$$
A=Q\Lambda Q^\mathsf T
$$

此时特征值非负，SVD 与特征分解高度一致。

若 $A$ 是一般矩阵：

- 特征值可以为负数或复数。
- 奇异值始终是非负实数。
- 特征向量描述同一空间中的不变方向。
- 奇异向量描述输入空间与输出空间之间的正交方向对应。

### 22.2 QR 与 SVD

QR 比 SVD 计算更便宜，适合满秩最小二乘。

SVD 计算成本更高，但能：

- 处理秩亏。
- 揭示条件数。
- 给出最小范数解。
- 提供最佳低秩近似。

## 23. 极分解

SVD 还可以导出极分解。以 $m\ge n$ 的经济型 SVD 为例：

$$
A=U_p\Sigma_pV^\mathsf T
$$

定义：

$$
H=(A^\mathsf TA)^{1/2}
=V\Sigma_pV^\mathsf T
$$

是对称半正定矩阵，描述沿正交方向的纯拉伸。

而：

$$
Q=U_pV^\mathsf T
$$

描述旋转、反射或更一般的部分等距映射，并且：

$$
A=QH
$$

当 $A$ 满列秩时：

$$
Q^\mathsf TQ=I_n
$$

若 $A$ 秩亏，$Q$ 应理解为相应有效子空间上的部分等距映射。对于 $m<n$ 的矩阵，可以使用对应的左极分解。

极分解可以理解为矩阵版本的复数极坐标：

$$
\text{线性变换}
=
\text{旋转或反射}
\times
\text{纯拉伸}
$$

## 24. 常见误区

### 误区 1：奇异值就是特征值

奇异值是 $A^\mathsf TA$ 特征值的非负平方根。只有在特殊矩阵中，它们才与 $A$ 的特征值直接一致。

### 误区 2：SVD 只适用于方阵

SVD 最重要的优势之一就是适用于任意矩形矩阵。

### 误区 3：左右奇异向量位于同一个空间

右奇异向量属于输入空间 $\mathbb R^n$，左奇异向量属于输出空间 $\mathbb R^m$。

### 误区 4：零奇异值没有信息

零奇异值准确揭示被完全消除的输入方向和无法达到的输出方向，是理解零空间和左零空间的关键。

### 误区 5：SVD 完全唯一

奇异值唯一，但奇异向量存在符号、重复子空间和零空间基选择的不唯一性。

### 误区 6：伪逆会恢复所有丢失信息

伪逆无法恢复零空间中已经丢失的分量。它只选择与观测一致或最接近观测的最小范数输入。

### 误区 7：$A^+A=I$ 对所有矩阵成立

$A^+A$ 是投影到行空间。只有满列秩时，行空间等于整个输入空间，才有 $A^+A=I_n$。

### 误区 8：$AA^+=I$ 对所有矩阵成立

$AA^+$ 是投影到列空间。只有满行秩时，列空间等于整个输出空间，才有 $AA^+=I_m$。

### 误区 9：最大奇异值越大，矩阵整体越不稳定

不稳定主要由最大与最小非零奇异值的比值决定，而不是仅看最大奇异值。

### 误区 10：小奇异值一定是噪声

小奇异值表示弱方向，但弱方向可能包含重要信息。是否舍弃必须结合任务判断。

### 误区 11：最佳低秩近似对所有误差度量都最佳

Eckart–Young 定理保证截断 SVD 在谱范数和 Frobenius 范数等特定酉不变范数下最佳，不应无限推广到所有指标。

### 误区 12：PCA 不需要中心化

未中心化 PCA 往往优先描述数据相对原点的整体位置，而不是样本围绕均值的变化。

### 误区 13：PCA 的主成分具有天然语义

主成分是最大方差方向，未必自动对应容易命名或具有因果意义的现实因素。

### 误区 14：保留方差越多，下游任务一定越好

高方差方向不一定最有利于分类、预测或决策。PCA 是无监督几何方法，不知道下游目标。

## 25. 应用连接

- **最小二乘**：秩亏或病态系统的稳定求解。
- **图像压缩**：使用少量奇异模式近似像素矩阵。
- **图像去噪**：截断弱奇异方向。
- **推荐系统**：用户—物品矩阵的低秩潜在因子模型。
- **自然语言处理**：潜在语义分析使用截断 SVD 提取语义方向。
- **PCA**：寻找数据主要变化方向并降维。
- **系统辨识**：识别动态系统的有效阶数。
- **控制与逆问题**：条件数分析和正则化。
- **信号处理**：子空间方法、滤波和频谱估计。
- **计算机视觉**：基础矩阵估计、形状分析和低秩几何模型。
- **生物信息学**：基因表达矩阵降维和噪声抑制。
- **科学计算**：大型矩阵的随机化低秩近似。
- **模型压缩**：用低秩因子近似神经网络权重矩阵。

## 26. 阶段练习

### 基础题

1. 求矩阵：

   $$
   A=
   \begin{bmatrix}
   3&0\\
   0&2\\
   0&0
   \end{bmatrix}
   $$

   的奇异值、秩和条件数。

2. 对贯穿矩阵：

   $$
   A=
   \begin{bmatrix}
   1&1\\
   1&1\\
   0&0
   \end{bmatrix}
   $$

   计算 $A^\mathsf TA$，并求其特征值和标准正交特征向量。

3. 写出第 2 题矩阵的完整 SVD 和紧致 SVD。

4. 根据第 3 题 SVD，直接写出四个基本子空间的一组标准正交基。

5. 求第 2 题矩阵的伪逆。

6. 使用伪逆求：

   $$
   A\boldsymbol x\approx
   \begin{bmatrix}2\\0\\1\end{bmatrix}
   $$

   的最小范数最小二乘解、投影和残差。

### 理解题

7. 为什么 $A^\mathsf TA$ 的特征值一定非负？

8. 为什么非零奇异值的数量等于矩阵的秩？

9. 证明 $AA^+$ 是列空间上的正交投影。

10. 证明 $A^+A$ 是行空间上的正交投影。

11. 为什么伪逆解 $A^+\boldsymbol b$ 与零空间正交？

12. 为什么小奇异值会使逆问题对噪声敏感？

13. SVD 与普通特征分解在几何问题上分别回答什么问题？

### 综合题

14. 设矩阵的奇异值为：

   $$
   8,\qquad3,\qquad1,\qquad0
   $$

   1. 求矩阵秩。
   2. 求最佳秩一近似的谱范数和 Frobenius 范数误差。
   3. 求最佳秩二近似的两种误差。
   4. 求矩阵的标准二范数条件数；若只考虑非零奇异子空间，再求有效条件数。

15. 设：

   $$
   A=
   \begin{bmatrix}
   1&0\\
   0&0.01
   \end{bmatrix},
   \qquad
   \boldsymbol b=
   \begin{bmatrix}1\\0.01\end{bmatrix}
   $$

   1. 求 $A^{-1}\boldsymbol b$。
   2. 若第二个观测变为 $0.011$，求新解。
   3. 解释变化被放大的原因。

16. 证明截断 SVD：

   $$
   A_k=
   \sum_{i=1}^k
   \sigma_i\boldsymbol u_i\boldsymbol v_i^\mathsf T
   $$

   的秩不超过 $k$。

17. 一个 $1000\times800$ 矩阵使用秩 20 的截断 SVD 存储：

   1. 原矩阵需要存储多少个数？
   2. 截断 SVD 需要存储多少个数？
   3. 存储量约为原来的百分之多少？

18. 对已经中心化的数据：

   $$
   X=
   \begin{bmatrix}
   -2&-4\\
   -1&-2\\
   0&0\\
   1&2\\
   2&4
   \end{bmatrix}
   $$

   1. 求秩。
   2. 求第一主成分方向。
   3. 求唯一非零奇异值。
   4. 求第一主成分方差。
   5. 判断降到一维是否损失信息。

19. 设：

   $$
   A=U\Sigma V^\mathsf T
   $$

   推导：

   $$
   \|A\|_F^2
   =
   \sum_i\sigma_i^2
   $$

20. 比较下列场景应优先使用特征分解、QR 还是 SVD：

   1. 求实对称矩阵的主方向。
   2. 求满列秩高矩阵的最小二乘解。
   3. 求秩亏矩阵的最小范数解。
   4. 对图像矩阵进行低秩压缩。

## 27. 参考答案

1. 因为矩阵已经是矩形对角形式，奇异值为：

   $$
   \sigma_1=3,
   \qquad
   \sigma_2=2
   $$

   秩为 2，条件数为：

   $$
   \kappa_2(A)=\frac32
   $$

2. 

   $$
   A^\mathsf TA=
   \begin{bmatrix}
   2&2\\
   2&2
   \end{bmatrix}
   $$

   特征值为 4 和 0。对应标准正交特征向量可取：

   $$
   \boldsymbol v_1=
   \frac1{\sqrt2}
   \begin{bmatrix}1\\1\end{bmatrix},
   \qquad
   \boldsymbol v_2=
   \frac1{\sqrt2}
   \begin{bmatrix}1\\-1\end{bmatrix}
   $$

3. 完整 SVD 为：

   $$
   U=
   \begin{bmatrix}
   \frac1{\sqrt2}&\frac1{\sqrt2}&0\\
   \frac1{\sqrt2}&-\frac1{\sqrt2}&0\\
   0&0&1
   \end{bmatrix}
   $$

   $$
   \Sigma=
   \begin{bmatrix}
   2&0\\
   0&0\\
   0&0
   \end{bmatrix}
   $$

   $$
   V=
   \frac1{\sqrt2}
   \begin{bmatrix}
   1&1\\
   1&-1
   \end{bmatrix}
   $$

   紧致 SVD 为：

   $$
   A
   =
   \left(
   \frac1{\sqrt2}
   \begin{bmatrix}1\\1\\0\end{bmatrix}
   \right)
   [2]
   \left(
   \frac1{\sqrt2}
   \begin{bmatrix}1&1\end{bmatrix}
   \right)
   $$

4. 

   $$
   C(A)=
   \operatorname{span}
   \left(
   \frac1{\sqrt2}
   \begin{bmatrix}1\\1\\0\end{bmatrix}
   \right)
   $$

   $$
   N(A^\mathsf T)=
   \operatorname{span}
   \left(
   \frac1{\sqrt2}
   \begin{bmatrix}1\\-1\\0\end{bmatrix},
   \begin{bmatrix}0\\0\\1\end{bmatrix}
   \right)
   $$

   $$
   C(A^\mathsf T)=
   \operatorname{span}
   \left(
   \frac1{\sqrt2}
   \begin{bmatrix}1\\1\end{bmatrix}
   \right)
   $$

   $$
   N(A)=
   \operatorname{span}
   \left(
   \frac1{\sqrt2}
   \begin{bmatrix}1\\-1\end{bmatrix}
   \right)
   $$

5. 

   $$
   A^+
   =
   \frac14
   \begin{bmatrix}
   1&1&0\\
   1&1&0
   \end{bmatrix}
   $$

6. 

   $$
   \hat{\boldsymbol x}
   =
   A^+\boldsymbol b
   =
   \begin{bmatrix}1/2\\1/2\end{bmatrix}
   $$

   投影为：

   $$
   \boldsymbol p
   =
   A\hat{\boldsymbol x}
   =
   \begin{bmatrix}1\\1\\0\end{bmatrix}
   $$

   残差为：

   $$
   \boldsymbol r
   =
   \boldsymbol b-\boldsymbol p
   =
   \begin{bmatrix}1\\-1\\1\end{bmatrix}
   $$

7. 对任意 $\boldsymbol x$：

   $$
   \boldsymbol x^\mathsf TA^\mathsf TA\boldsymbol x
   =
   \|A\boldsymbol x\|^2
   \ge0
   $$

   所以 $A^\mathsf TA$ 半正定，所有特征值非负。

8. 正交矩阵不改变秩，而 $\Sigma$ 的秩等于其非零对角元素数量，因此：

   $$
   \operatorname{rank}(A)
   =
   \operatorname{rank}(\Sigma)
   =
   \#\{i:\sigma_i>0\}
   $$

9. 使用紧致 SVD：

   $$
   AA^+
   =
   U_r\Sigma_rV_r^\mathsf T
   V_r\Sigma_r^{-1}U_r^\mathsf T
   =
   U_rU_r^\mathsf T
   $$

   $U_r$ 的列是列空间的标准正交基，所以这是列空间上的正交投影。

10. 同理：

   $$
   A^+A
   =
   V_rV_r^\mathsf T
   $$

   $V_r$ 的列是行空间的标准正交基，所以这是行空间上的正交投影。

11. $A^+\boldsymbol b$ 是 $V_r$ 各列的线性组合，因此属于行空间。第三和第四阶段已经知道：

   $$
   C(A^\mathsf T)\perp N(A)
   $$

   所以伪逆解与零空间正交。

12. 伪逆沿第 $i$ 个奇异方向除以 $\sigma_i$。若 $\sigma_i$ 很小，$\boldsymbol b$ 在对应左奇异方向上的微小扰动会被乘以很大的 $1/\sigma_i$。

13. 特征分解寻找方阵中经过变换后保持在原直线上的方向。SVD 寻找输入空间的右奇异方向如何映射为输出空间的左奇异方向，并给出非负缩放强度。

14. 非零奇异值有 3 个，所以秩为 3。

   最佳秩一近似：

   $$
   \|A-A_1\|_2=3
   $$

   $$
   \|A-A_1\|_F
   =
   \sqrt{3^2+1^2}
   =
   \sqrt{10}
   $$

   最佳秩二近似：

   $$
   \|A-A_2\|_2=1
   $$

   $$
   \|A-A_2\|_F=1
   $$

   因为存在零奇异值，矩阵秩亏，所以标准二范数条件数为：

   $$
   \kappa_2(A)=\infty
   $$

   若只考虑非零奇异子空间，有效条件数为：

   $$
   \frac81=8
   $$

15. 原解为：

   $$
   \boldsymbol x=
   \begin{bmatrix}1\\1\end{bmatrix}
   $$

   新观测对应：

   $$
   \boldsymbol x_{\text{new}}
   =
   \begin{bmatrix}1\\1.1\end{bmatrix}
   $$

   第二个观测只增加了 $0.001$，但第二个解分量增加了 $0.1$。原因是第二奇异值只有 $0.01$，求逆时该方向被放大 100 倍。

16. 每一项：

   $$
   \sigma_i\boldsymbol u_i\boldsymbol v_i^\mathsf T
   $$

   的秩至多为 1。$k$ 个秩一矩阵之和的秩至多为 $k$，因此：

   $$
   \operatorname{rank}(A_k)\le k
   $$

17. 原矩阵存储量：

   $$
   1000\cdot800=800000
   $$

   截断 SVD 存储量：

   $$
   20(1000+800+1)
   =36020
   $$

   比例约为：

   $$
   \frac{36020}{800000}\times100\%
   \approx4.50\%
   $$

18. 数据矩阵可以写成两个向量的外积，所以秩为 1。

   第一主成分方向为：

   $$
   \boldsymbol v_1
   =
   \frac1{\sqrt5}
   \begin{bmatrix}1\\2\end{bmatrix}
   $$

   唯一非零奇异值为：

   $$
   \sigma_1=5\sqrt2
   $$

   第一主成分方差：

   $$
   \frac{\sigma_1^2}{5-1}
   =
   \frac{50}{4}
   =
   12.5
   $$

   第二奇异值为零，因此降到一维不会损失信息。

19. 因为 Frobenius 范数在正交变换下不变：

   $$
   \|A\|_F
   =
   \|U\Sigma V^\mathsf T\|_F
   =
   \|\Sigma\|_F
   $$

   而 $\Sigma$ 只有奇异值位于对角线上，所以：

   $$
   \|A\|_F^2
   =
   \sum_i\sigma_i^2
   $$

20. 实对称矩阵主方向优先使用对称特征分解；满列秩高矩阵最小二乘优先使用 QR；秩亏矩阵最小范数解使用 SVD；图像低秩压缩使用截断 SVD。

## 28. 阶段检验

在不查阅资料的情况下，回答以下问题：

1. SVD 为什么比普通特征分解适用范围更广？
2. 奇异值与 $A^\mathsf TA$ 的特征值有什么关系？
3. 左奇异向量和右奇异向量分别位于哪个空间？
4. $A\boldsymbol v_i=\sigma_i\boldsymbol u_i$ 有什么几何意义？
5. 完整 SVD、经济型 SVD 和紧致 SVD 有什么区别？
6. 如何从 SVD 直接读出矩阵的秩？
7. 如何从 $U$ 和 $V$ 直接读出四个基本子空间？
8. 外积展开怎样把矩阵拆成多个秩一模式？
9. 伪逆如何通过翻转非零奇异值构造？
10. $AA^+$ 和 $A^+A$ 分别投影到哪个子空间？
11. 为什么 $A^+\boldsymbol b$ 是最小范数最小二乘解？
12. 最大和最小奇异值分别描述什么？
13. 条件数为什么能衡量逆问题敏感性？
14. 小奇异值为什么既代表弱信息方向，也可能放大噪声？
15. 截断 SVD 为什么是最佳低秩近似？
16. 图像压缩中，秩 $k$ 如何控制存储量和重建质量？
17. PCA 为什么必须先中心化数据？
18. PCA 的主成分方向、得分和解释方差如何从 SVD 得到？
19. 何时应优先选择 QR，何时应选择 SVD？
20. SVD 如何将本课程中的几何变换、子空间、投影、方程求解和数据降维统一起来？

能够清楚回答这些问题，独立完成阶段练习，并能从一个矩阵的 SVD 解释其有效方向、信息强度和数值稳定性，就完成了线性代数主干课程的六个阶段。
