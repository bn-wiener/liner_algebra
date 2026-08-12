# 第一阶段：向量与空间

> 核心问题：一个空间由哪些独立方向构成？怎样用这些方向表示空间中的所有对象？

本阶段建立线性代数最基础的语言。重点不是记忆符号，而是理解下面这条认知链：

$$
\text{向量}
\rightarrow
\text{线性组合}
\rightarrow
\text{张成空间}
\rightarrow
\text{线性无关}
\rightarrow
\text{基}
\rightarrow
\text{维数与坐标}
$$

## 学习目标

完成本阶段后，应当能够：

1. 从几何、数据和抽象对象三个角度解释向量。
2. 计算向量的加法、数乘和线性组合，并说明其几何意义。
3. 判断一个向量是否属于给定向量组的张成空间。
4. 判断一组向量是否线性相关，并找出其中的冗余方向。
5. 判断一组向量是否构成某个空间的基。
6. 解释基、维数和坐标之间的关系。
7. 在简单的不同基之间表示同一个向量。

## 1. 向量：空间中的对象

### 1.1 三种理解方式

向量可以同时理解为：

- **几何对象**：具有方向和大小的箭头，例如位移、速度和力。
- **数据对象**：一组有顺序的数，例如一个人的身高、体重和年龄。
- **抽象对象**：只要满足向量加法和标量乘法规则，函数、多项式和信号也可以视为向量。

例如：

$$
\boldsymbol v=
\begin{bmatrix}
2\\
3
\end{bmatrix}
$$

在平面中，它可以表示从原点到 $(2,3)$ 的位移；在数据问题中，它也可以表示一个具有两个特征的样本。

向量通常写成列形式：

$$
\boldsymbol v=
\begin{bmatrix}
v_1\\
v_2\\
\vdots\\
v_n
\end{bmatrix}
\in\mathbb R^n
$$

$\mathbb R^n$ 表示由 $n$ 个实数组成的所有向量的集合。

### 1.2 向量的基本运算

设：

$$
\boldsymbol u=
\begin{bmatrix}
u_1\\u_2
\end{bmatrix},
\qquad
\boldsymbol v=
\begin{bmatrix}
v_1\\v_2
\end{bmatrix}
$$

向量加法按对应分量进行：

$$
\boldsymbol u+\boldsymbol v=
\begin{bmatrix}
u_1+v_1\\
u_2+v_2
\end{bmatrix}
$$

标量 $c$ 与向量相乘：

$$
c\boldsymbol v=
\begin{bmatrix}
cv_1\\
cv_2
\end{bmatrix}
$$

几何上：

- 向量加法表示位移的合成。
- 当 $c>1$ 时，向量被拉长。
- 当 $0<c<1$ 时，向量被缩短。
- 当 $c<0$ 时，向量还会反向。
- 当 $c=0$ 时，结果是零向量 $\boldsymbol 0$。

### 例 1：合成位移

一个人先向右 3 米、向上 1 米，再向左 1 米、向上 2 米：

$$
\begin{bmatrix}3\\1\end{bmatrix}
+
\begin{bmatrix}-1\\2\end{bmatrix}
=
\begin{bmatrix}2\\3\end{bmatrix}
$$

最终位移是向右 2 米、向上 3 米。

## 2. 线性组合：用已有方向构造新向量

给定向量 $\boldsymbol v_1,\dots,\boldsymbol v_k$，表达式

$$
c_1\boldsymbol v_1+cdots+c_k\boldsymbol v_k
$$

称为这些向量的**线性组合**，其中 $c_1,dots,c_k$ 是标量。

线性组合的核心问题是：

> 只使用现有的方向，通过缩放和相加，能够构造出哪些向量？

### 例 2：用两个方向表示目标向量

设：

$$
\boldsymbol v_1=
\begin{bmatrix}1\\1\end{bmatrix},
\qquad
\boldsymbol v_2=
\begin{bmatrix}1\\-1\end{bmatrix}
$$

尝试表示：

$$
\boldsymbol b=
\begin{bmatrix}4\\2\end{bmatrix}
$$

令 $c_1\boldsymbol v_1+c_2\boldsymbol v_2=\boldsymbol b$，得到：

$$
\begin{cases}
c_1+c_2=4\\
c_1-c_2=2
\end{cases}
$$

解得 $c_1=3$、$c_2=1$，所以：

$$
\begin{bmatrix}4\\2\end{bmatrix}
=3
\begin{bmatrix}1\\1\end{bmatrix}
+
\begin{bmatrix}1\\-1\end{bmatrix}
$$

这说明 $\boldsymbol b$ 可以由 $\boldsymbol v_1$ 和 $\boldsymbol v_2$ 构造出来。

## 3. 张成空间：所有可到达的向量

向量 $\boldsymbol v_1,\dots,\boldsymbol v_k$ 的所有线性组合构成它们的**张成空间**：

$$
\operatorname{span}(\boldsymbol v_1,\dots,\boldsymbol v_k)
=
\left\{
c_1\boldsymbol v_1+\cdots+c_k\boldsymbol v_k
\right\}
$$

在二维空间中：

- 一个非零向量张成一条经过原点的直线。
- 两个共线向量仍然只张成一条直线。
- 两个不共线向量张成整个 $\mathbb R^2$。

在三维空间中：

- 一个非零向量张成一条直线。
- 两个不共线向量张成一个经过原点的平面。
- 三个不共面的向量张成整个 $\mathbb R^3$。

注意：由向量张成的直线或平面一定经过原点，因为所有系数都取 0 时必然得到零向量。

### 例 3：判断目标是否属于张成空间

设：

$$
\boldsymbol v_1=
\begin{bmatrix}1\\2\end{bmatrix},
\qquad
\boldsymbol v_2=
\begin{bmatrix}2\\4\end{bmatrix}
$$

因为 $\boldsymbol v_2=2\boldsymbol v_1$，所以它们只张成直线：

$$
\operatorname{span}(\boldsymbol v_1,\boldsymbol v_2)
=
\left\{
t\begin{bmatrix}1\\2\end{bmatrix}:t\in\mathbb R
\right\}
$$

$\begin{bmatrix}3\\6\end{bmatrix}$ 位于这条直线上，而 $\begin{bmatrix}3\\5\end{bmatrix}$ 不在其中。

## 4. 线性相关：识别冗余方向

若存在一组**不全为零**的系数，使得：

$$
c_1\boldsymbol v_1+cdots+c_k\boldsymbol v_k=\boldsymbol 0
$$

则称这组向量**线性相关**。如果只有

$$
c_1=\cdots=c_k=0
$$

才能得到零向量，则称它们**线性无关**。

直觉上：

- 线性相关意味着至少有一个向量可以由其他向量表示，存在冗余。
- 线性无关意味着每个向量都提供了一个新的、不可替代的方向。

### 例 4：三个向量中的冗余

设：

$$
\boldsymbol v_1=
\begin{bmatrix}1\\0\end{bmatrix},
\quad
\boldsymbol v_2=
\begin{bmatrix}0\\1\end{bmatrix},
\quad
\boldsymbol v_3=
\begin{bmatrix}1\\1\end{bmatrix}
$$

因为：

$$
\boldsymbol v_1+\boldsymbol v_2-\boldsymbol v_3=\boldsymbol 0
$$

所以三个向量线性相关。$\boldsymbol v_3$ 没有增加新的可到达方向。

### 快速判断

- 向量组中包含零向量，一定线性相关。
- 两个向量互为倍数，一定线性相关。
- $\mathbb R^2$ 中超过两个向量，一定线性相关。
- $\mathbb R^3$ 中超过三个向量，一定线性相关。
- 向量数量少于空间维数时，可能线性无关，但不能张成整个空间。

最后两条将在理解“维数”后变得自然。

## 5. 基：无冗余的完整坐标系统

向量空间 $V$ 的一组向量要成为一组**基**，必须同时满足：

1. **能够张成 $V$**：没有缺少方向。
2. **线性无关**：没有多余方向。

因此，基可以理解为：

> 描述一个空间所需的一组最精简、完整的基本方向。

$\mathbb R^2$ 的标准基是：

$$
\boldsymbol e_1=
\begin{bmatrix}1\\0\end{bmatrix},
\qquad
\boldsymbol e_2=
\begin{bmatrix}0\\1\end{bmatrix}
$$

但基并不唯一。例如：

$$
\boldsymbol b_1=
\begin{bmatrix}1\\1\end{bmatrix},
\qquad
\boldsymbol b_2=
\begin{bmatrix}1\\-1\end{bmatrix}
$$

同样构成 $\mathbb R^2$ 的一组基，因为它们不共线，并且能够张成整个平面。

### 为什么基表示是唯一的

若一组基可以用两种方式表示同一个向量：

$$
\boldsymbol v
=c_1\boldsymbol b_1+\cdots+c_n\boldsymbol b_n
=d_1\boldsymbol b_1+\cdots+d_n\boldsymbol b_n
$$

两式相减得到：

$$
(c_1-d_1)\boldsymbol b_1+cdots+(c_n-d_n)\boldsymbol b_n=\boldsymbol 0
$$

由于基向量线性无关，只能有 $c_i-d_i=0$。因此 $c_i=d_i$，基下的坐标表示是唯一的。

## 6. 维数：独立方向的数量

一个向量空间中任意一组基所包含的向量数量都相同，这个数量称为空间的**维数**。

$$
\dim(\mathbb R^2)=2,
\qquad
\dim(\mathbb R^3)=3
$$

维数的本质不是“写了几个数字”，而是：

> 完整描述这个空间需要多少个相互独立的方向？

### 抽象空间的例子

次数不超过 2 的多项式组成空间：

$$
P_2=\{a+bx+cx^2:a,b,c\in\mathbb R\}
$$

它的一组基是：

$$
1,\quad x,\quad x^2
$$

所以：

$$
\dim(P_2)=3
$$

多项式 $2-3x+5x^2$ 在这组基下，可以用坐标表示为：

$$
\begin{bmatrix}2\\-3\\5\end{bmatrix}
$$

这说明线性代数研究的不只是空间中的箭头，还研究所有具有线性结构的对象。

## 7. 向量空间与子空间

一个集合若对向量加法和标量乘法保持封闭，就具有线性结构。直观地说：

- 集合中的两个对象相加，结果仍在集合中。
- 集合中的对象乘任意标量，结果仍在集合中。
- 因而集合必然包含零向量。

在 $\mathbb R^3$ 中，下列集合是子空间：

$$
W=
\left\{
\begin{bmatrix}x\\y\\0\end{bmatrix}:x,y\in\mathbb R
\right\}
$$

它是经过原点的 $xy$ 平面。

下列集合不是子空间：

$$
S=
\left\{
\begin{bmatrix}x\\y\end{bmatrix}:x+y=1
\right\}
$$

因为零向量不满足 $x+y=1$。几何上，它是一条不经过原点的直线。

阶段 3 将进一步研究矩阵产生的列空间、行空间和零空间；此处只需建立“子空间必须保留线性运算”的直觉。

## 8. 坐标：向量在一组基下的描述

设 $B=(\boldsymbol b_1,\dots,\boldsymbol b_n)$ 是空间的一组基。若：

$$
\boldsymbol v
=c_1\boldsymbol b_1+cdots+c_n\boldsymbol b_n
$$

则称：

$$
[\boldsymbol v]_B=
\begin{bmatrix}
c_1\\
\vdots\\
c_n
\end{bmatrix}
$$

为 $\boldsymbol v$ 在基 $B$ 下的坐标。

继续使用：

$$
\boldsymbol b_1=
\begin{bmatrix}1\\1\end{bmatrix},
\qquad
\boldsymbol b_2=
\begin{bmatrix}1\\-1\end{bmatrix}
$$

由于：

$$
\begin{bmatrix}4\\2\end{bmatrix}
=3\boldsymbol b_1+\boldsymbol b_2
$$

所以：

$$
\left[
\begin{bmatrix}4\\2\end{bmatrix}
\right]_B
=
\begin{bmatrix}3\\1\end{bmatrix}
$$

同一个向量没有改变，改变的是描述它所使用的坐标语言：

$$
\boxed{\text{向量不变，坐标随基改变}}
$$

阶段 2 将进一步把“基向量组成的矩阵”理解为坐标变换工具。

## 9. 一个贯穿本阶段的例子

设：

$$
\boldsymbol v_1=
\begin{bmatrix}1\\1\\0\end{bmatrix},
\quad
\boldsymbol v_2=
\begin{bmatrix}0\\1\\1\end{bmatrix},
\quad
\boldsymbol v_3=
\begin{bmatrix}1\\2\\1\end{bmatrix}
$$

观察到：

$$
\boldsymbol v_3=\boldsymbol v_1+\boldsymbol v_2
$$

因此：

- 三个向量线性相关。
- $\boldsymbol v_3$ 是冗余向量。
- 三个向量的张成空间与 $\boldsymbol v_1,\boldsymbol v_2$ 的张成空间相同。
- $\boldsymbol v_1,\boldsymbol v_2$ 不互为倍数，因此线性无关。
- $\boldsymbol v_1,\boldsymbol v_2$ 构成该张成空间的一组基。
- 该张成空间的维数为 2，是 $\mathbb R^3$ 中一个经过原点的平面。

这个例子把“线性组合、张成、相关性、基和维数”连接成了一个整体。

## 10. 常见误区

### 误区 1：向量就是坐标

坐标只是向量在某组基下的表示。更换基以后坐标会改变，但向量本身没有改变。

### 误区 2：向量数量就是张成空间的维数

向量之间可能存在冗余。三个向量也可能只张成一条直线或一个平面。

### 误区 3：能够张成空间就一定是基

还必须满足线性无关。例如 $\boldsymbol e_1,\boldsymbol e_2,\boldsymbol e_1+\boldsymbol e_2$ 能张成 $\mathbb R^2$，但存在冗余，不是一组基。

### 误区 4：线性无关就一定能张成整个空间

$\mathbb R^3$ 中两个不共线向量可以线性无关，但只能张成一个平面，不能张成整个 $\mathbb R^3$。

### 误区 5：任意直线和平面都是子空间

只有经过原点且对线性运算封闭的直线和平面才是子空间。

## 11. 应用连接

- **数据表示**：一个样本可以表示成特征向量。
- **坐标系统**：选择不同的基，可以让同一个问题更容易描述。
- **数据冗余**：线性相关表示某些特征可以由其他特征构造。
- **信号与图像**：复杂信号可以由一组基本信号线性组合得到。
- **降维**：寻找更少的独立方向表示主要信息，本质上是在寻找合适的低维基。
- **后续矩阵理论**：矩阵的列就是一组向量，矩阵的秩就是其中独立方向的数量。

## 12. 阶段练习

### 基础题

1. 计算：

   $$
   2\begin{bmatrix}1\\-1\end{bmatrix}
   -3\begin{bmatrix}0\\2\end{bmatrix}
   $$

2. 判断 $\begin{bmatrix}6\\9\end{bmatrix}$ 是否属于 $\operatorname{span}\left(\begin{bmatrix}2\\3\end{bmatrix}\right)$。

3. 描述 $\operatorname{span}\left(\begin{bmatrix}1\\0\\0\end{bmatrix},\begin{bmatrix}0\\1\\0\end{bmatrix}\right)$ 的几何形状。

4. 判断下列向量是否线性相关：

   $$
   \begin{bmatrix}1\\0\end{bmatrix},
   \quad
   \begin{bmatrix}0\\1\end{bmatrix},
   \quad
   \begin{bmatrix}2\\3\end{bmatrix}
   $$

5. 判断下列向量能否构成 $\mathbb R^2$ 的基：

   $$
   \begin{bmatrix}1\\2\end{bmatrix},
   \qquad
   \begin{bmatrix}2\\5\end{bmatrix}
   $$

### 理解题

6. 为什么包含零向量的向量组一定线性相关？

7. 为什么 $\mathbb R^3$ 中两个线性无关的向量不能张成整个 $\mathbb R^3$？

8. 给出一组能够张成 $\mathbb R^2$、但不是 $\mathbb R^2$ 的基的向量。

9. 判断集合 $\{(x,y):y=2x\}$ 是否为 $\mathbb R^2$ 的子空间，并说明原因。

10. 判断集合 $\{(x,y):y=2x+1\}$ 是否为 $\mathbb R^2$ 的子空间，并说明原因。

### 综合题

11. 设基：

    $$
    B=\left(
    \begin{bmatrix}1\\1\end{bmatrix},
    \begin{bmatrix}2\\-1\end{bmatrix}
    \right)
    $$

    求 $\boldsymbol v=\begin{bmatrix}5\\1\end{bmatrix}$ 在基 $B$ 下的坐标。

12. 设：

    $$
    \boldsymbol v_1=\begin{bmatrix}1\\0\\1\end{bmatrix},
    \quad
    \boldsymbol v_2=\begin{bmatrix}0\\1\\1\end{bmatrix},
    \quad
    \boldsymbol v_3=\begin{bmatrix}1\\1\\2\end{bmatrix}
    $$

    判断三个向量是否线性相关，并找出它们张成空间的一组基和维数。

## 13. 参考答案

1. $\begin{bmatrix}2\\-8\end{bmatrix}$。
2. 属于，因为 $\begin{bmatrix}6\\9\end{bmatrix}=3\begin{bmatrix}2\\3\end{bmatrix}$。
3. $\mathbb R^3$ 中经过原点的 $xy$ 平面。
4. 线性相关，因为第三个向量等于前两个向量的线性组合：

   $$
   \begin{bmatrix}2\\3\end{bmatrix}
   =2\begin{bmatrix}1\\0\end{bmatrix}
   +3\begin{bmatrix}0\\1\end{bmatrix}
   $$

5. 可以。两个向量不互为倍数，因此线性无关并能张成 $\mathbb R^2$。
6. 只需令零向量对应的系数非零、其他系数为零，就能得到一个非平凡的零线性组合。
7. $\mathbb R^3$ 需要三个独立方向；两个线性无关向量最多张成一个平面。
8. 例如 $\boldsymbol e_1,\boldsymbol e_2,\boldsymbol e_1+\boldsymbol e_2$。它们能张成 $\mathbb R^2$，但线性相关。
9. 是。它是一条经过原点的直线，并且对向量加法和标量乘法封闭。
10. 不是。它不包含零向量，也不对标量乘法封闭。
11. 设：

    $$
    c_1\begin{bmatrix}1\\1\end{bmatrix}
    +c_2\begin{bmatrix}2\\-1\end{bmatrix}
    =\begin{bmatrix}5\\1\end{bmatrix}
    $$

    解得 $c_1=\frac{7}{3}$、$c_2=\frac{4}{3}$，因此：

    $$
    [\boldsymbol v]_B=
    \begin{bmatrix}7/3\\4/3\end{bmatrix}
    $$

12. 因为 $\boldsymbol v_3=\boldsymbol v_1+\boldsymbol v_2$，所以三个向量线性相关。$\boldsymbol v_1,\boldsymbol v_2$ 线性无关，因此它们构成张成空间的一组基，该空间维数为 2。

## 14. 阶段检验

在不查阅资料的情况下，回答以下问题：

1. 向量和向量坐标有什么区别？
2. 线性组合描述了什么操作？
3. 张成空间描述了什么范围？
4. 线性相关为什么意味着信息冗余？
5. 一组向量成为基需要满足哪两个条件？
6. 维数为什么等于独立方向的数量？
7. 为什么子空间必须经过原点？

能够清楚回答这些问题，并独立完成阶段练习，就可以进入第二阶段：矩阵与线性变换。
