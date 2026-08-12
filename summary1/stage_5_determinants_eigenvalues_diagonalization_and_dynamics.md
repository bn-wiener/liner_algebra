# 第五阶段：行列式、特征值、对角化与动态系统

> 核心问题：一个线性变换怎样改变空间的面积、体积和方向？哪些方向经过变换后保持不变？怎样利用这些特殊方向理解反复作用的矩阵？

前四个阶段研究了向量空间、线性变换、方程求解和最佳近似。本阶段开始研究矩阵内部更深的结构。

本阶段有两条相互连接的主线：

$$
\text{矩阵}
\rightarrow
\text{行列式}
\rightarrow
\text{体积缩放与可逆性}
$$

以及：

$$
\text{矩阵}
\rightarrow
\text{特征值与特征向量}
\rightarrow
\text{对角化}
\rightarrow
\text{矩阵幂与动态系统}
$$

它们最终在对称矩阵中汇合：

$$
A=A^\mathsf T
\quad\Longrightarrow\quad
A=Q\Lambda Q^\mathsf T
$$

对称矩阵拥有实特征值和一组标准正交特征向量，因此具有最清晰的几何结构。

## 学习目标

完成本阶段后，应当能够：

1. 从面积、体积和空间定向解释行列式。
2. 使用消元、三角结构和代数余子式计算行列式。
3. 熟练使用行列式的主要性质。
4. 解释行列式为零与不可逆、秩下降和空间压缩之间的关系。
5. 求矩阵的特征值、特征向量和特征空间。
6. 区分代数重数与几何重数。
7. 判断矩阵能否对角化，并构造 $A=P\Lambda P^{-1}$。
8. 使用对角化计算矩阵幂和矩阵函数。
9. 理解相似矩阵表示同一个线性变换的不同坐标形式。
10. 掌握实对称矩阵的谱定理。
11. 使用特征值判断二次型的正定性。
12. 理解 Rayleigh 商与主方向。
13. 使用特征分解分析离散动态系统的长期行为。
14. 理解 Markov 链稳定状态和线性微分方程中的特征结构。

## 1. 行列式：空间尺度的缩放因子

### 1.1 二维几何意义

设：

$$
A=
\begin{bmatrix}
a&b\\
c&d
\end{bmatrix}
$$

两个标准基向量经过 $A$ 变换后成为矩阵的两列：

$$
\boldsymbol a_1=
\begin{bmatrix}a\\c\end{bmatrix},
\qquad
\boldsymbol a_2=
\begin{bmatrix}b\\d\end{bmatrix}
$$

单位正方形被变成由 $\boldsymbol a_1,\boldsymbol a_2$ 张成的平行四边形。

行列式：

$$
\det(A)=ad-bc
$$

其绝对值等于该平行四边形的面积：

$$
|\det(A)|
=
\text{面积缩放倍数}
$$

### 1.2 符号与空间定向

行列式的符号描述定向：

- $\det(A)>0$：空间定向保持。
- $\det(A)<0$：空间定向翻转。
- $\det(A)=0$：面积被压缩为零。

例如关于 $y$ 轴反射：

$$
F=
\begin{bmatrix}
-1&0\\
0&1
\end{bmatrix}
$$

有：

$$
\det(F)=-1
$$

面积大小不变，但空间定向发生翻转。

### 1.3 三维与高维

对于 $3\times3$ 矩阵：

$$
|\det(A)|
$$

表示单位立方体经过变换后的平行六面体体积。

对于 $n\times n$ 矩阵，行列式表示 $n$ 维有向体积的缩放因子。虽然高维图形难以直接画出，缩放规律仍完全成立。

### 1.4 例

设：

$$
A=
\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}
$$

则：

$$
\det(A)=2\cdot2-1\cdot1=3
$$

所以该变换把面积放大为原来的 3 倍，并保持空间定向。

## 2. 行列式的基本性质

### 2.1 单位矩阵

$$
\det(I)=1
$$

因为单位变换不改变体积和定向。

### 2.2 矩阵乘积

$$
\boxed{
\det(AB)=\det(A)\det(B)
}
$$

几何上，先由 $B$ 缩放体积，再由 $A$ 缩放，整体缩放倍数就是两者乘积。

### 2.3 转置

$$
\det(A^\mathsf T)=\det(A)
$$

行列式对行和列具有对称地位。

### 2.4 逆矩阵

若 $A$ 可逆，则：

$$
\det(A^{-1})
=
\frac1{\det(A)}
$$

因为：

$$
1=\det(I)=\det(AA^{-1})
=\det(A)\det(A^{-1})
$$

### 2.5 行列式与标量

对于 $n\times n$ 矩阵：

$$
\det(cA)=c^n\det(A)
$$

因为矩阵的每一列都被放大 $c$ 倍，总体积被放大 $c^n$ 倍。

注意它一般不是 $c\det(A)$。

### 2.6 三角矩阵

若 $A$ 是上三角或下三角矩阵，则：

$$
\det(A)=a_{11}a_{22}\cdots a_{nn}
$$

即行列式等于对角元素乘积。

## 3. 行变换与行列式

三种初等行变换对行列式有不同影响。

### 3.1 交换两行

交换两行会改变空间定向：

$$
\det(B)=-\det(A)
$$

### 3.2 某行乘常数

若将某一行乘 $c$：

$$
\det(B)=c\det(A)
$$

只有一个方向被缩放了 $c$ 倍。

### 3.3 某行加上另一行的倍数

$$
R_i\leftarrow R_i+cR_j
$$

不会改变行列式。

几何上，这种剪切操作改变形状，但不改变体积。

### 3.4 使用消元计算行列式

设：

$$
A=
\begin{bmatrix}
1&2&1\\
2&5&3\\
1&0&2
\end{bmatrix}
$$

执行：

$$
R_2\leftarrow R_2-2R_1,
\qquad
R_3\leftarrow R_3-R_1
$$

得到：

$$
\begin{bmatrix}
1&2&1\\
0&1&1\\
0&-2&1
\end{bmatrix}
$$

再执行：

$$
R_3\leftarrow R_3+2R_2
$$

得到：

$$
\begin{bmatrix}
1&2&1\\
0&1&1\\
0&0&3
\end{bmatrix}
$$

以上操作均为行倍加，不改变行列式。因此：

$$
\det(A)=1\cdot1\cdot3=3
$$

若消元中交换行或缩放行，必须同步记录其影响。

## 4. 代数余子式展开

### 4.1 子式与代数余子式

删除第 $i$ 行和第 $j$ 列后得到的行列式称为子式，记为：

$$
M_{ij}
$$

代数余子式定义为：

$$
C_{ij}=(-1)^{i+j}M_{ij}
$$

符号排列为：

$$
\begin{bmatrix}
+&-&+&\cdots\\
-&+&-&\cdots\\
+&-&+&\cdots\\
\vdots&\vdots&\vdots&\ddots
\end{bmatrix}
$$

### 4.2 沿一行展开

沿第 $i$ 行展开：

$$
\det(A)
=
\sum_{j=1}^n a_{ij}C_{ij}
$$

也可以沿任意一列展开。

### 4.3 三阶例子

设：

$$
A=
\begin{bmatrix}
1&2&0\\
3&-1&2\\
0&4&1
\end{bmatrix}
$$

沿第一行展开：

$$
\det(A)
=
1
\begin{vmatrix}
-1&2\\
4&1
\end{vmatrix}
-
2
\begin{vmatrix}
3&2\\
0&1
\end{vmatrix}
$$

所以：

$$
\det(A)
=1(-1-8)-2(3-0)
=-15
$$

### 4.4 何时使用余子式

- 矩阵中某行或某列有许多零时，余子式展开很方便。
- 对一般大矩阵，消元比递归展开高效得多。
- 余子式更适合理论证明，而不是大规模数值计算。

## 5. 行列式与可逆性

对于 $n\times n$ 方阵，下列条件等价：

1. $A$ 可逆。
2. $\det(A)\neq0$。
3. $\operatorname{rank}(A)=n$。
4. $A$ 的列线性无关。
5. $A$ 的列张成 $\mathbb R^n$。
6. $N(A)=\{\boldsymbol0\}$。
7. 对每个 $\boldsymbol b$，方程 $A\boldsymbol x=\boldsymbol b$ 有唯一解。
8. $0$ 不是 $A$ 的特征值。

### 5.1 为什么行列式为零意味着不可逆

若：

$$
\det(A)=0
$$

则变换后的单位体积为零，说明整个空间被压缩到更低维子空间。至少有一个输入方向被消除，因此不同输入可能产生同一输出，无法唯一恢复。

### 5.2 可逆矩阵定理的统一认识

这些等价条件分别来自不同视角：

- 行列式视角：体积是否被压扁。
- 消元视角：是否每列都有主元。
- 子空间视角：零空间是否只有零向量。
- 方程视角：是否每个目标都有唯一原像。
- 特征视角：是否存在被压缩为零的非零方向。

## 6. 特征值与特征向量

### 6.1 定义

若存在非零向量 $\boldsymbol v$ 和标量 $\lambda$，满足：

$$
\boxed{
A\boldsymbol v=\lambda\boldsymbol v
}
$$

则：

- $\boldsymbol v$ 是 $A$ 的特征向量。
- $\lambda$ 是对应的特征值。

特征向量经过变换后仍位于原来的直线上，只发生缩放或反向。

### 6.2 不同特征值的几何含义

- $\lambda>1$：沿该方向放大。
- $0<\lambda<1$：沿该方向缩小。
- $\lambda<0$：方向翻转并按 $|\lambda|$ 缩放。
- $\lambda=1$：该方向保持不变。
- $\lambda=0$：该方向被压缩为零，因此矩阵不可逆。

### 6.3 为什么零向量不是特征向量

零向量满足：

$$
A\boldsymbol0=\lambda\boldsymbol0
$$

对任意 $\lambda$ 都成立，无法提供任何特殊方向信息，因此定义中必须要求 $\boldsymbol v\neq\boldsymbol0$。

### 6.4 例

设：

$$
A=
\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}
$$

对向量：

$$
\boldsymbol v_1=
\begin{bmatrix}1\\1\end{bmatrix}
$$

有：

$$
A\boldsymbol v_1
=
\begin{bmatrix}3\\3\end{bmatrix}
=3\boldsymbol v_1
$$

所以 $\boldsymbol v_1$ 是特征向量，对应 $\lambda_1=3$。

对：

$$
\boldsymbol v_2=
\begin{bmatrix}1\\-1\end{bmatrix}
$$

有：

$$
A\boldsymbol v_2
=
\begin{bmatrix}1\\-1\end{bmatrix}
=1\boldsymbol v_2
$$

所以对应 $\lambda_2=1$。

## 7. 特征多项式与特征空间

### 7.1 求特征值

由：

$$
A\boldsymbol v=\lambda\boldsymbol v
$$

得到：

$$
(A-\lambda I)\boldsymbol v=\boldsymbol0
$$

要存在非零解，矩阵 $A-\lambda I$ 必须不可逆：

$$
\boxed{
\det(A-\lambda I)=0
}
$$

这称为特征方程。

### 7.2 特征多项式

$$
p_A(\lambda)=\det(A-\lambda I)
$$

称为矩阵 $A$ 的特征多项式。它是 $n$ 次多项式，其根是矩阵的特征值。

有些教材使用：

$$
\det(\lambda I-A)
$$

两种定义的根完全相同，只可能相差整体符号。

### 7.3 二阶例子

对：

$$
A=
\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}
$$

有：

$$
\det(A-\lambda I)
=
\begin{vmatrix}
2-\lambda&1\\
1&2-\lambda
\end{vmatrix}
$$

$$
=(2-\lambda)^2-1
$$

$$
=(\lambda-3)(\lambda-1)
$$

所以特征值为：

$$
\lambda_1=3,
\qquad
\lambda_2=1
$$

### 7.4 求特征向量

对每个特征值，求解：

$$
(A-\lambda I)\boldsymbol v=\boldsymbol0
$$

特征值 $\lambda$ 对应的特征空间为：

$$
E_\lambda
=
N(A-\lambda I)
$$

特征空间包括零向量，但其中所有非零向量才是特征向量。

对于 $\lambda=3$：

$$
A-3I=
\begin{bmatrix}
-1&1\\
1&-1
\end{bmatrix}
$$

所以：

$$
E_3
=
\operatorname{span}
\left(
\begin{bmatrix}1\\1\end{bmatrix}
\right)
$$

对于 $\lambda=1$：

$$
E_1
=
\operatorname{span}
\left(
\begin{bmatrix}1\\-1\end{bmatrix}
\right)
$$

### 7.5 迹与行列式

对于 $n\times n$ 矩阵，按代数重数计：

$$
\operatorname{tr}(A)
=
\lambda_1+\cdots+\lambda_n
$$

$$
\det(A)
=
\lambda_1\cdots\lambda_n
$$

其中：

$$
\operatorname{tr}(A)
=a_{11}+\cdots+a_{nn}
$$

称为矩阵的迹。

对贯穿矩阵：

$$
\operatorname{tr}(A)=4=3+1
$$

$$
\det(A)=3=3\cdot1
$$

## 8. 代数重数与几何重数

### 8.1 代数重数

特征值 $\lambda$ 作为特征多项式根出现的次数，称为代数重数。

例如：

$$
p_A(\lambda)=(\lambda-2)^3(\lambda+1)
$$

则特征值 2 的代数重数为 3。

### 8.2 几何重数

特征空间的维数：

$$
\dim E_\lambda
=
\dim N(A-\lambda I)
$$

称为特征值 $\lambda$ 的几何重数。

总有：

$$
1
\le
\text{几何重数}
\le
\text{代数重数}
$$

### 8.3 缺陷矩阵

考虑：

$$
A=
\begin{bmatrix}
1&1\\
0&1
\end{bmatrix}
$$

特征多项式为：

$$
(1-\lambda)^2
$$

所以 $\lambda=1$ 的代数重数为 2。

但：

$$
A-I=
\begin{bmatrix}
0&1\\
0&0
\end{bmatrix}
$$

其零空间只有一个独立方向：

$$
E_1
=
\operatorname{span}
\left(
\begin{bmatrix}1\\0\end{bmatrix}
\right)
$$

几何重数为 1。矩阵没有足够多的独立特征向量，因此不能对角化。

## 9. 对角化

### 9.1 为什么对角矩阵简单

若：

$$
\Lambda=
\begin{bmatrix}
\lambda_1&&0\\
&\ddots&\\
0&&\lambda_n
\end{bmatrix}
$$

则：

$$
\Lambda\boldsymbol x
=
\begin{bmatrix}
\lambda_1x_1\\
\vdots\\
\lambda_nx_n
\end{bmatrix}
$$

每个坐标方向彼此独立，只进行缩放。

### 9.2 构造对角化

若 $A$ 有 $n$ 个线性无关特征向量：

$$
\boldsymbol v_1,\dots,\boldsymbol v_n
$$

令：

$$
P=
\begin{bmatrix}
\boldsymbol v_1&\cdots&\boldsymbol v_n
\end{bmatrix}
$$

以及：

$$
\Lambda=
\begin{bmatrix}
\lambda_1&&0\\
&\ddots&\\
0&&\lambda_n
\end{bmatrix}
$$

因为：

$$
AP=P\Lambda
$$

所以：

$$
\boxed{
A=P\Lambda P^{-1}
}
$$

### 9.3 换基解释

$P$ 的列是一组特征向量基。执行顺序为：

1. $P^{-1}$：把标准坐标转换成特征向量基下的坐标。
2. $\Lambda$：沿各特征方向独立缩放。
3. $P$：把结果转换回标准坐标。

因此对角化是在寻找一种让线性变换最简单的坐标系统。

### 9.4 可对角化条件

矩阵可对角化，当且仅当它有 $n$ 个线性无关特征向量。

以下条件足以保证可对角化：

- $A$ 有 $n$ 个互不相同的特征值。
- $A$ 是实对称矩阵。
- 对每个特征值，几何重数等于代数重数。

特征值重复并不必然导致不可对角化。例如单位矩阵只有一个特征值 1，但每个非零向量都是特征向量，因此完全可以对角化。

## 10. 相似矩阵

若存在可逆矩阵 $P$，使：

$$
B=P^{-1}AP
$$

则称 $A$ 与 $B$ 相似。

相似矩阵表示同一个线性变换在不同基下的矩阵形式，因此共享许多结构量：

- 相同特征值。
- 相同特征多项式。
- 相同行列式。
- 相同迹。
- 相同秩。

但它们的具体元素、列空间和特征向量坐标一般不同。

对角化就是寻找一个与 $A$ 相似的对角矩阵 $\Lambda$。

## 11. 使用对角化计算矩阵幂

若：

$$
A=P\Lambda P^{-1}
$$

则：

$$
A^2
=
P\Lambda P^{-1}P\Lambda P^{-1}
=
P\Lambda^2P^{-1}
$$

一般地：

$$
\boxed{
A^k=P\Lambda^kP^{-1}
}
$$

而：

$$
\Lambda^k
=
\operatorname{diag}
(\lambda_1^k,\dots,\lambda_n^k)
$$

### 11.1 贯穿矩阵的幂

对：

$$
A=
\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}
$$

取标准正交特征向量：

$$
\boldsymbol q_1
=
\frac1{\sqrt2}
\begin{bmatrix}1\\1\end{bmatrix},
\qquad
\boldsymbol q_2
=
\frac1{\sqrt2}
\begin{bmatrix}1\\-1\end{bmatrix}
$$

于是：

$$
Q=
\frac1{\sqrt2}
\begin{bmatrix}
1&1\\
1&-1
\end{bmatrix}
$$

$$
\Lambda=
\begin{bmatrix}
3&0\\
0&1
\end{bmatrix}
$$

并且：

$$
A=Q\Lambda Q^\mathsf T
$$

所以：

$$
A^k
=
\frac12
\begin{bmatrix}
3^k+1&3^k-1\\
3^k-1&3^k+1
\end{bmatrix}
$$

直接逐次矩阵相乘很复杂，而特征分解把问题变成标量幂。

### 11.2 矩阵函数

若函数 $f$ 可以作用于特征值，则可定义：

$$
f(A)
=
P
\begin{bmatrix}
f(\lambda_1)&&0\\
&\ddots&\\
0&&f(\lambda_n)
\end{bmatrix}
P^{-1}
$$

例如矩阵指数、矩阵平方根和某些微分方程解都可这样构造。

## 12. 复特征值与旋转

实矩阵的特征值不一定全是实数。

例如旋转 $90^\circ$：

$$
R=
\begin{bmatrix}
0&-1\\
1&0
\end{bmatrix}
$$

特征方程为：

$$
\lambda^2+1=0
$$

所以：

$$
\lambda=\pm i
$$

在实平面中，没有任何非零方向经过 $90^\circ$ 旋转后仍留在原直线上，因此没有实特征向量。

复特征值的模描述缩放，辐角描述旋转。共轭复特征值通常对应实空间中的旋转—缩放二维子空间。

本课程后续主要关注实对称矩阵，因为它们保证拥有实特征值和实正交特征向量。

## 13. 对称矩阵与谱定理

### 13.1 对称矩阵

若：

$$
A^\mathsf T=A
$$

则称 $A$ 为实对称矩阵。

例如：

$$
\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}
$$

就是对称矩阵。

### 13.2 对称矩阵的特殊性质

实对称矩阵满足：

1. 所有特征值都是实数。
2. 不同特征值对应的特征向量正交。
3. 可以选择一组标准正交特征向量作为 $\mathbb R^n$ 的基。
4. 一定可以正交对角化。

### 13.3 谱定理

对任意实对称矩阵，都存在正交矩阵 $Q$ 和实对角矩阵 $\Lambda$，使：

$$
\boxed{
A=Q\Lambda Q^\mathsf T
}
$$

因为 $Q$ 正交：

$$
Q^{-1}=Q^\mathsf T
$$

相比一般对角化 $P\Lambda P^{-1}$，正交对角化具有更稳定、更清晰的几何意义：先旋转或反射坐标系，沿正交方向缩放，再恢复坐标。

### 13.4 为什么不同特征值的特征向量正交

设：

$$
A\boldsymbol u=\lambda\boldsymbol u,
\qquad
A\boldsymbol v=\mu\boldsymbol v
$$

并且 $\lambda\neq\mu$。因为 $A=A^\mathsf T$：

$$
\boldsymbol u^\mathsf TA\boldsymbol v
=(A\boldsymbol u)^\mathsf T\boldsymbol v
$$

左侧为：

$$
\mu\boldsymbol u^\mathsf T\boldsymbol v
$$

右侧为：

$$
\lambda\boldsymbol u^\mathsf T\boldsymbol v
$$

因此：

$$
(\mu-\lambda)
\boldsymbol u^\mathsf T\boldsymbol v=0
$$

由于 $\mu\neq\lambda$，所以：

$$
\boldsymbol u^\mathsf T\boldsymbol v=0
$$

## 14. 二次型

### 14.1 定义

形如：

$$
q(\boldsymbol x)
=
\boldsymbol x^\mathsf TA\boldsymbol x
$$

的标量表达式称为二次型。

对于：

$$
A=
\begin{bmatrix}
a&b\\
b&d
\end{bmatrix},
\qquad
\boldsymbol x=
\begin{bmatrix}x\\y\end{bmatrix}
$$

有：

$$
\boldsymbol x^\mathsf TA\boldsymbol x
=
ax^2+2bxy+dy^2
$$

### 14.2 为什么只需研究对称矩阵

任意方阵都可分解为对称部分和反对称部分：

$$
A
=
\frac{A+A^\mathsf T}{2}
+
\frac{A-A^\mathsf T}{2}
$$

反对称部分满足：

$$
\boldsymbol x^\mathsf T
\left(
\frac{A-A^\mathsf T}{2}
\right)
\boldsymbol x=0
$$

因此二次型只由对称部分决定。

### 14.3 主轴变换

若对称矩阵：

$$
A=Q\Lambda Q^\mathsf T
$$

令：

$$
\boldsymbol y=Q^\mathsf T\boldsymbol x
$$

则：

$$
\boldsymbol x^\mathsf TA\boldsymbol x
=
\boldsymbol y^\mathsf T\Lambda\boldsymbol y
=
\lambda_1y_1^2+\cdots+\lambda_ny_n^2
$$

交叉项消失。特征向量给出二次曲面的主轴方向，特征值给出各主轴方向上的曲率或缩放强度。

### 14.4 贯穿矩阵的二次型

对：

$$
A=
\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}
$$

有：

$$
q(x,y)
=
2x^2+2xy+2y^2
$$

在特征向量坐标中：

$$
q
=
3y_1^2+y_2^2
$$

所以两个主方向上的系数分别为 3 和 1。

## 15. 正定、半正定与不定矩阵

对实对称矩阵 $A$：

### 15.1 正定

若对所有非零 $\boldsymbol x$：

$$
\boldsymbol x^\mathsf TA\boldsymbol x>0
$$

则称 $A$ 正定。

等价地：

$$
\lambda_i>0
\quad\text{对所有 }i
$$

### 15.2 半正定

若：

$$
\boldsymbol x^\mathsf TA\boldsymbol x\ge0
$$

则称 $A$ 半正定。等价地，所有特征值非负。

### 15.3 负定与不定

- 所有特征值小于零：负定。
- 所有特征值不大于零：半负定。
- 既有正特征值又有负特征值：不定。

### 15.4 几何与优化意义

- 正定二次型的等值面是椭球。
- 不定二次型具有鞍形结构。
- 优化中，正定 Hessian 表示严格局部极小。
- 半正定矩阵常表示能量、方差和协方差结构。

### 15.5 二阶正定判据

对于对称矩阵：

$$
A=
\begin{bmatrix}
a&b\\
b&d
\end{bmatrix}
$$

正定当且仅当：

$$
a>0,
\qquad
\det(A)>0
$$

更高维中，Sylvester 判据要求所有顺序主子式都为正。

贯穿矩阵的特征值为 3 和 1，均为正，所以它是正定矩阵。

## 16. Rayleigh 商

对于实对称矩阵 $A$ 和非零向量 $\boldsymbol x$，定义：

$$
R_A(\boldsymbol x)
=
\frac{\boldsymbol x^\mathsf TA\boldsymbol x}
{\boldsymbol x^\mathsf T\boldsymbol x}
$$

它表示矩阵沿方向 $\boldsymbol x$ 的有效缩放或二次型强度。

若 $\boldsymbol x$ 是特征向量：

$$
R_A(\boldsymbol x)=\lambda
$$

谱定理给出：

$$
\lambda_{\min}
\le
R_A(\boldsymbol x)
\le
\lambda_{\max}
$$

并且：

$$
\lambda_{\max}
=
\max_{\|\boldsymbol x\|=1}
\boldsymbol x^\mathsf TA\boldsymbol x
$$

$$
\lambda_{\min}
=
\min_{\|\boldsymbol x\|=1}
\boldsymbol x^\mathsf TA\boldsymbol x
$$

所以最大和最小特征向量分别给出二次型最强和最弱的方向。

## 17. 离散动态系统

### 17.1 状态演化

离散线性动态系统写成：

$$
\boldsymbol x_{k+1}=A\boldsymbol x_k
$$

因此：

$$
\boldsymbol x_k=A^k\boldsymbol x_0
$$

若 $A=P\Lambda P^{-1}$，将初始状态分解为特征向量：

$$
\boldsymbol x_0
=
c_1\boldsymbol v_1+\cdots+c_n\boldsymbol v_n
$$

则：

$$
\boxed{
\boldsymbol x_k
=
c_1\lambda_1^k\boldsymbol v_1
+\cdots+
c_n\lambda_n^k\boldsymbol v_n
}
$$

每个特征方向独立演化。

### 17.2 长期行为

- $|\lambda|<1$：对应分量逐渐衰减。
- $|\lambda|>1$：对应分量指数增长。
- $\lambda=1$：对应分量保持不变。
- $\lambda=-1$：对应分量正负交替。
- 复特征值：通常产生旋转或振荡。

### 17.3 稳定性

谱半径定义为：

$$
\rho(A)
=
\max_i|\lambda_i|
$$

若：

$$
\rho(A)<1
$$

则对任意初始状态：

$$
\boldsymbol x_k\rightarrow\boldsymbol0
$$

若存在 $|\lambda|>1$ 的特征值，沿对应方向的微小分量会增长，系统通常不稳定。

当 $\rho(A)=1$ 时，需要进一步考虑对应特征值、重数和 Jordan 结构，不能只凭模等于 1 下结论。

### 17.4 贯穿矩阵的动态行为

取：

$$
\boldsymbol x_0=
\begin{bmatrix}1\\0\end{bmatrix}
$$

它可以写成：

$$
\boldsymbol x_0
=
\frac12
\begin{bmatrix}1\\1\end{bmatrix}
+
\frac12
\begin{bmatrix}1\\-1\end{bmatrix}
$$

所以：

$$
\boldsymbol x_k
=
\frac{3^k}{2}
\begin{bmatrix}1\\1\end{bmatrix}
+
\frac12
\begin{bmatrix}1\\-1\end{bmatrix}
$$

当 $k$ 很大时，第一项占主导，状态方向越来越接近 $(1,1)^\mathsf T$。

这揭示一个普遍规律：若存在唯一模最大的特征值，反复作用后系统通常趋向其主特征向量方向。

## 18. Markov 链

### 18.1 状态转移

设 $\boldsymbol p_k$ 表示时刻 $k$ 的概率分布，使用列向量约定：

$$
\boldsymbol p_{k+1}=M\boldsymbol p_k
$$

若 $M$ 的每一列元素非负且列和为 1，则 $M$ 是列随机矩阵。

### 18.2 稳态分布

稳态分布 $\boldsymbol\pi$ 满足：

$$
M\boldsymbol\pi=\boldsymbol\pi
$$

所以它是特征值 1 对应的特征向量，并且还要满足：

$$
\pi_i\ge0,
\qquad
\sum_i\pi_i=1
$$

### 18.3 例

设：

$$
M=
\begin{bmatrix}
0.8&0.3\\
0.2&0.7
\end{bmatrix}
$$

求解：

$$
(M-I)\boldsymbol\pi=\boldsymbol0
$$

得到比例：

$$
\pi_1:\pi_2=3:2
$$

归一化后：

$$
\boldsymbol\pi=
\begin{bmatrix}0.6\\0.4\end{bmatrix}
$$

另一个特征值为 0.5，其分量会按 $0.5^k$ 衰减。因此长期状态趋向稳态分布。

## 19. 线性微分方程与矩阵指数

连续时间线性系统写成：

$$
\frac{d\boldsymbol x}{dt}
=
A\boldsymbol x
$$

其解为：

$$
\boldsymbol x(t)
=
e^{At}\boldsymbol x(0)
$$

矩阵指数定义为：

$$
e^{At}
=
I+At+\frac{A^2t^2}{2!}
+\frac{A^3t^3}{3!}+\cdots
$$

若：

$$
A=P\Lambda P^{-1}
$$

则：

$$
e^{At}
=
Pe^{\Lambda t}P^{-1}
$$

其中：

$$
e^{\Lambda t}
=
\operatorname{diag}
(e^{\lambda_1t},\dots,e^{\lambda_nt})
$$

因此：

- 实部为负的特征值对应衰减模式。
- 实部为正的特征值对应增长模式。
- 非零虚部对应振荡。

## 20. 贯穿本阶段的完整总结

再次考虑：

$$
A=
\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}
$$

### 20.1 行列式

$$
\det(A)=3
$$

面积放大 3 倍，矩阵可逆。

### 20.2 特征结构

$$
\lambda_1=3,
\qquad
\boldsymbol v_1=
\begin{bmatrix}1\\1\end{bmatrix}
$$

$$
\lambda_2=1,
\qquad
\boldsymbol v_2=
\begin{bmatrix}1\\-1\end{bmatrix}
$$

两个特征方向正交。

### 20.3 谱分解

$$
A=Q
\begin{bmatrix}
3&0\\
0&1
\end{bmatrix}
Q^\mathsf T
$$

其中：

$$
Q=
\frac1{\sqrt2}
\begin{bmatrix}
1&1\\
1&-1
\end{bmatrix}
$$

### 20.4 矩阵幂

$$
A^k
=
\frac12
\begin{bmatrix}
3^k+1&3^k-1\\
3^k-1&3^k+1
\end{bmatrix}
$$

### 20.5 二次型

$$
\boldsymbol x^\mathsf TA\boldsymbol x
=
2x^2+2xy+2y^2
$$

在特征坐标中：

$$
=3y_1^2+y_2^2
$$

两个特征值均为正，因此矩阵正定。

### 20.6 动态行为

反复应用 $A$ 时，特征值 3 对应分量快速增长，特征值 1 对应分量保持不变。除非初始状态恰好完全位于第二特征方向，否则长期方向趋向 $(1,1)^\mathsf T$。

同一个矩阵由此被统一理解为：

$$
\boxed{
\text{体积缩放}
\leftrightarrow
\text{不变方向}
\leftrightarrow
\text{正交主轴}
\leftrightarrow
\text{长期动态}
}
$$

## 21. 常见误区

### 误区 1：行列式就是矩阵所有元素的乘积

只有三角矩阵的行列式等于对角元素乘积。一般矩阵必须使用行列式规则计算。

### 误区 2：$\det(A+B)=\det(A)+\det(B)$

行列式不是普通线性函数，一般不满足该等式。

### 误区 3：行列式为负表示面积为负

实际面积是 $|\det(A)|$。负号表示空间定向翻转。

### 误区 4：特征向量经过变换后完全不变

特征向量的方向保持在同一条直线上，但长度可能改变，方向也可能在负特征值时翻转。

### 误区 5：零向量是所有矩阵的特征向量

特征向量定义中明确要求非零。

### 误区 6：矩阵的每一列都是特征向量

矩阵的列是标准基向量的像，不一定是矩阵自己的特征向量。

### 误区 7：存在特征值就一定能对角化

每个复方阵都有特征值，但只有拥有足够多线性无关特征向量时才能对角化。

### 误区 8：重复特征值一定不能对角化

关键是对应特征空间是否提供足够多独立特征向量，而不是是否重复。

### 误区 9：相似矩阵是相同矩阵

相似矩阵通常元素不同，只是表示同一个线性变换的不同坐标形式。

### 误区 10：对称矩阵只是外观对称

对称性带来实特征值、正交特征向量和稳定的谱分解，是深刻的结构条件。

### 误区 11：所有正定矩阵都必须是对角矩阵

正定矩阵可以有交叉项；通过正交换基后才变成正对角形式。

### 误区 12：最大特征值总能决定全部长期行为

还要考虑初始状态是否含有对应特征向量分量、最大模特征值是否唯一，以及矩阵是否可对角化。

### 误区 13：$\rho(A)=1$ 就一定稳定

模为 1 的重复特征值和 Jordan 块可能导致多项式增长，需要进一步分析。

## 22. 应用连接

- **几何变换**：行列式测量面积、体积和定向变化。
- **积分换元**：Jacobian 行列式修正局部体积尺度。
- **微分方程**：特征值决定系统的增长、衰减与振荡模式。
- **控制系统**：特征值位置决定稳定性与响应速度。
- **Markov 链**：特征值 1 对应稳态分布，其他特征值决定收敛速度。
- **振动分析**：特征向量是固有振型，特征值与固有频率相关。
- **结构工程**：刚度矩阵和质量矩阵的特征问题描述结构响应。
- **量子力学**：可观测量由对称或 Hermitian 算子表示，特征值是可能测量结果。
- **优化**：Hessian 的特征值描述不同方向的曲率和极值类型。
- **数据分析**：协方差矩阵的特征向量给出主成分方向。
- **图论**：图拉普拉斯矩阵的特征值反映连通性和聚类结构。
- **搜索排序**：PageRank 可理解为随机矩阵的主特征向量问题。

## 23. 阶段练习

### 基础题

1. 计算：

   $$
   \det
   \begin{bmatrix}
   3&1\\
   2&4
   \end{bmatrix}
   $$

   并解释其几何意义。

2. 已知 $\det(A)=-2$，其中 $A$ 是 $3\times3$ 矩阵。分别求：

   $$
   \det(A^\mathsf T),
   \qquad
   \det(A^{-1}),
   \qquad
   \det(3A)
   $$

3. 使用消元计算：

   $$
   \det
   \begin{bmatrix}
   1&2&1\\
   2&5&3\\
   1&0&2
   \end{bmatrix}
   $$

4. 求矩阵：

   $$
   A=
   \begin{bmatrix}
   4&1\\
   2&3
   \end{bmatrix}
   $$

   的特征值。

5. 求第 4 题中每个特征值对应的特征空间。

6. 判断：

   $$
   B=
   \begin{bmatrix}
   1&1\\
   0&1
   \end{bmatrix}
   $$

   是否可对角化。

### 理解题

7. 为什么 $\det(A)=0$ 意味着至少有一个非零方向被矩阵压缩为零？

8. 为什么不同特征值对应的特征向量一定线性无关？

9. 特征值的代数重数和几何重数有什么区别？它们与对角化有什么关系？

10. 为什么 $A=P\Lambda P^{-1}$ 可以理解为一次换基？

11. 为什么实对称矩阵一定可以正交对角化？

12. 为什么正定对称矩阵的所有特征值都为正？

13. 对动态系统 $\boldsymbol x_{k+1}=A\boldsymbol x_k$，为什么特征值的模控制长期增长或衰减？

### 综合题

14. 对矩阵：

   $$
   A=
   \begin{bmatrix}
   4&1\\
   2&3
   \end{bmatrix}
   $$

   1. 构造 $A=P\Lambda P^{-1}$。
   2. 使用对角化计算 $A^k$。
   3. 求 $A^5(1,0)^\mathsf T$。

15. 对称矩阵：

   $$
   S=
   \begin{bmatrix}
   5&2\\
   2&2
   \end{bmatrix}
   $$

   1. 求特征值和标准正交特征向量。
   2. 写出谱分解 $S=Q\Lambda Q^\mathsf T$。
   3. 判断 $S$ 是否正定。
   4. 写出二次型。

16. 判断下列二次型的类型：

   $$
   q_1(x,y)=3x^2+2xy+3y^2
   $$

   $$
   q_2(x,y)=x^2-4y^2
   $$

   $$
   q_3(x,y)=(x-y)^2
   $$

17. 设动态系统：

   $$
   \boldsymbol x_{k+1}
   =
   \begin{bmatrix}
   0.8&0\\
   0&1.2
   \end{bmatrix}
   \boldsymbol x_k
   $$

   对初始状态 $(1,1)^\mathsf T$，求 $\boldsymbol x_k$，并描述长期行为。

18. 对 Markov 矩阵：

   $$
   M=
   \begin{bmatrix}
   0.8&0.3\\
   0.2&0.7
   \end{bmatrix}
   $$

   求稳态分布，并验证概率和为 1。

19. 对贯穿矩阵：

   $$
   A=
   \begin{bmatrix}
   2&1\\
   1&2
   \end{bmatrix}
   $$

   1. 求 $A^{10}$。
   2. 求二次型在单位圆上的最大值和最小值。
   3. 指出达到最大值和最小值的方向。

## 24. 参考答案

1. 

   $$
   \det
   \begin{bmatrix}
   3&1\\
   2&4
   \end{bmatrix}
   =12-2=10
   $$

   该变换把面积放大 10 倍，并保持定向。

2. 

   $$
   \det(A^\mathsf T)=-2
   $$

   $$
   \det(A^{-1})=-\frac12
   $$

   因为 $A$ 是 $3\times3$：

   $$
   \det(3A)=3^3\det(A)=27(-2)=-54
   $$

3. 使用讲义中的行倍加消元可得上三角矩阵，对角线为 $1,1,3$，所以：

   $$
   \det(A)=3
   $$

4. 特征方程为：

   $$
   \det(A-\lambda I)
   =
   \begin{vmatrix}
   4-\lambda&1\\
   2&3-\lambda
   \end{vmatrix}
   $$

   $$
   =(4-\lambda)(3-\lambda)-2
   $$

   $$
   =\lambda^2-7\lambda+10
   =(\lambda-5)(\lambda-2)
   $$

   所以特征值为 5 和 2。

5. 对 $\lambda=5$：

   $$
   A-5I=
   \begin{bmatrix}
   -1&1\\
   2&-2
   \end{bmatrix}
   $$

   所以：

   $$
   E_5
   =
   \operatorname{span}
   \left(
   \begin{bmatrix}1\\1\end{bmatrix}
   \right)
   $$

   对 $\lambda=2$：

   $$
   A-2I=
   \begin{bmatrix}
   2&1\\
   2&1
   \end{bmatrix}
   $$

   所以：

   $$
   E_2
   =
   \operatorname{span}
   \left(
   \begin{bmatrix}1\\-2\end{bmatrix}
   \right)
   $$

6. $B$ 只有特征值 1，代数重数为 2，但特征空间：

   $$
   E_1
   =
   \operatorname{span}
   \left(
   \begin{bmatrix}1\\0\end{bmatrix}
   \right)
   $$

   只有一个独立方向，因此不可对角化。

7. 行列式绝对值是体积缩放倍数。若为零，单位体积被压成零体积，列向量线性相关，所以存在非零系数 $\boldsymbol x$ 使 $A\boldsymbol x=\boldsymbol0$。

8. 对不同特征值 $\lambda_1,\dots,\lambda_k$ 的特征向量建立线性组合为零。应用 $A$ 并与原关系消去一个特征值后，可以逐步证明所有系数为零，因此它们线性无关。

9. 代数重数是特征值作为特征多项式根的次数；几何重数是对应特征空间的维数。矩阵可对角化当且仅当所有特征空间总共提供 $n$ 个独立特征向量，即每个特征值的几何重数达到其代数重数。

10. $P^{-1}$ 把标准坐标转换到特征向量基，$\Lambda$ 在该基中沿各坐标方向独立缩放，$P$ 再转换回标准坐标。

11. 谱定理保证实对称矩阵有实特征值，不同特征值的特征空间彼此正交；每个重复特征空间内部还可选择标准正交基，最终得到完整标准正交特征基。

12. 若 $A\boldsymbol v=\lambda\boldsymbol v$ 且 $\boldsymbol v\neq0$，则：

   $$
   \boldsymbol v^\mathsf TA\boldsymbol v
   =
   \lambda\|\boldsymbol v\|^2
   $$

   正定要求左侧为正，而 $\|\boldsymbol v\|^2>0$，所以 $\lambda>0$。

13. 在特征向量分解中，每次作用 $A$ 都把对应分量乘 $\lambda$。经过 $k$ 次后变为 $\lambda^k$，所以其模决定增长、衰减或保持。

14. 取：

   $$
   P=
   \begin{bmatrix}
   1&1\\
   1&-2
   \end{bmatrix},
   \qquad
   \Lambda=
   \begin{bmatrix}
   5&0\\
   0&2
   \end{bmatrix}
   $$

   有：

   $$
   P^{-1}
   =
   \begin{bmatrix}
   \frac23&\frac13\\
   \frac13&-\frac13
   \end{bmatrix}
   $$

   因此：

   $$
   A^k
   =
   P
   \begin{bmatrix}
   5^k&0\\
   0&2^k
   \end{bmatrix}
   P^{-1}
   $$

   展开得：

   $$
   A^k
   =
   \frac13
   \begin{bmatrix}
   2\cdot5^k+2^k&5^k-2^k\\
   2\cdot5^k-2\cdot2^k&5^k+2\cdot2^k
   \end{bmatrix}
   $$

   对 $k=5$：

   $$
   A^5
   \begin{bmatrix}1\\0\end{bmatrix}
   =
   \frac13
   \begin{bmatrix}
   2\cdot3125+32\\
   2\cdot3125-64
   \end{bmatrix}
   =
   \begin{bmatrix}
   2094\\
   2062
   \end{bmatrix}
   $$

15. 特征方程：

   $$
   \det(S-\lambda I)
   =(5-\lambda)(2-\lambda)-4
   $$

   $$
   =\lambda^2-7\lambda+6
   =(\lambda-6)(\lambda-1)
   $$

   对 $\lambda=6$，可取特征向量 $(2,1)^\mathsf T$；对 $\lambda=1$，可取 $(1,-2)^\mathsf T$。单位化：

   $$
   \boldsymbol q_1=
   \frac1{\sqrt5}
   \begin{bmatrix}2\\1\end{bmatrix},
   \qquad
   \boldsymbol q_2=
   \frac1{\sqrt5}
   \begin{bmatrix}1\\-2\end{bmatrix}
   $$

   所以：

   $$
   Q=
   \frac1{\sqrt5}
   \begin{bmatrix}
   2&1\\
   1&-2
   \end{bmatrix},
   \qquad
   \Lambda=
   \begin{bmatrix}
   6&0\\
   0&1
   \end{bmatrix}
   $$

   $$
   S=Q\Lambda Q^\mathsf T
   $$

   两个特征值均为正，所以 $S$ 正定。二次型为：

   $$
   \boldsymbol x^\mathsf TS\boldsymbol x
   =5x^2+4xy+2y^2
   $$

16. 

   $$
   q_1(x,y)
   =
   \begin{bmatrix}x&y\end{bmatrix}
   \begin{bmatrix}3&1\\1&3\end{bmatrix}
   \begin{bmatrix}x\\y\end{bmatrix}
   $$

   对应特征值为 4 和 2，均为正，所以正定。

   $$
   q_2(x,y)=x^2-4y^2
   $$

   对应特征值 1 和 $-4$，所以不定。

   $$
   q_3(x,y)=(x-y)^2
   $$

   始终非负，但在 $x=y\neq0$ 时等于零，所以半正定但不正定。

17. 因为矩阵已经对角化：

   $$
   \boldsymbol x_k
   =
   \begin{bmatrix}
   0.8^k&0\\
   0&1.2^k
   \end{bmatrix}
   \begin{bmatrix}1\\1\end{bmatrix}
   =
   \begin{bmatrix}
   0.8^k\\
   1.2^k
   \end{bmatrix}
   $$

   第一分量趋于零，第二分量增长，因此系统不稳定，长期方向趋近第二坐标轴。

18. 解：

   $$
   M\boldsymbol\pi=\boldsymbol\pi
   $$

   得：

   $$
   0.2\pi_1=0.3\pi_2
   $$

   所以：

   $$
   \pi_1:\pi_2=3:2
   $$

   归一化：

   $$
   \boldsymbol\pi=
   \begin{bmatrix}0.6\\0.4\end{bmatrix}
   $$

   两个分量非负且总和为 1。

19. 

   $$
   A^{10}
   =
   \frac12
   \begin{bmatrix}
   3^{10}+1&3^{10}-1\\
   3^{10}-1&3^{10}+1
   \end{bmatrix}
   $$

   因为：

   $$
   3^{10}=59049
   $$

   所以：

   $$
   A^{10}
   =
   \begin{bmatrix}
   29525&29524\\
   29524&29525
   \end{bmatrix}
   $$

   单位圆上二次型的最大值是最大特征值 3，在方向 $(1,1)^\mathsf T$ 上取得；最小值是最小特征值 1，在方向 $(1,-1)^\mathsf T$ 上取得。方向应理解为对应的单位化向量及其反方向。

## 25. 阶段检验

在不查阅资料的情况下，回答以下问题：

1. 行列式的绝对值和符号分别表示什么？
2. 三种初等行变换怎样影响行列式？
3. 为什么 $\det(AB)=\det(A)\det(B)$ 具有自然的几何解释？
4. 行列式为零与不可逆、秩下降和零空间非平凡有什么关系？
5. 特征向量为什么代表线性变换的不变方向？
6. 如何从 $\det(A-\lambda I)=0$ 得到特征值？
7. 特征空间、代数重数和几何重数分别是什么？
8. 矩阵可对角化的充要条件是什么？
9. $A=P\Lambda P^{-1}$ 的三个因子分别完成什么坐标操作？
10. 为什么对角化能够简化矩阵幂和矩阵函数？
11. 实对称矩阵为什么格外重要？
12. 谱定理如何解释对称矩阵的几何作用？
13. 如何使用特征值判断二次型的正定性？
14. Rayleigh 商的最大值和最小值与特征值有什么关系？
15. 特征值怎样决定离散动态系统的长期行为？
16. Markov 链的稳态为什么对应特征值 1？
17. 连续系统中，特征值实部怎样决定增长、衰减和振荡？

能够清楚回答这些问题，独立完成阶段练习，并能把一个矩阵解释为体积缩放、特征方向和动态模式，就可以进入第六阶段：SVD、伪逆、低秩近似与数据分析。
