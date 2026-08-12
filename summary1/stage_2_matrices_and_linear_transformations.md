# 第二阶段：矩阵与线性变换

> 核心问题：一个空间中的对象，怎样按照统一规则变成另一个对象？

第一阶段研究了空间中的对象、方向和坐标。本阶段开始研究作用在这些对象上的规则，并建立下面这条认知链：

$$
\text{矩阵}
\rightarrow
\text{线性变换}
\rightarrow
\text{矩阵乘向量}
\rightarrow
\text{基本几何变换}
\rightarrow
\text{变换复合}
\rightarrow
\text{逆变换}
\rightarrow
\text{换基}
$$

学习时应始终保留两个同步视角：

- **代数视角**：矩阵是一张按行列排列的数表，可以执行运算。
- **几何视角**：矩阵描述一个线性变换，表示整个空间怎样被拉伸、旋转、反射、剪切或压缩。

本阶段最重要的目标，是看到矩阵时不再只想到“行乘列”，而是首先想到：

> 这个矩阵把输入空间中的基向量送到了哪里？

## 学习目标

完成本阶段后，应当能够：

1. 识别矩阵的尺寸，并说明一个 $m\times n$ 矩阵为什么表示从 $\mathbb R^n$ 到 $\mathbb R^m$ 的变换。
2. 判断一个给定规则是否为线性变换。
3. 根据基向量的像写出线性变换的矩阵。
4. 从列组合、分量计算、方程组和几何变换四个角度解释 $A\boldsymbol x$。
5. 写出二维缩放、旋转、反射、剪切和投影矩阵。
6. 解释矩阵乘法为什么表示线性变换的复合，并正确判断执行顺序。
7. 判断简单变换是否可逆，并求简单矩阵的逆。
8. 区分线性变换、仿射变换和坐标变换。
9. 使用换基矩阵在不同坐标系统之间转换表示。
10. 将矩阵计算结果翻译回空间、数据或实际问题中的含义。

## 1. 矩阵：组织多个输入与输出关系

### 1.1 矩阵的形式与尺寸

一个 $m\times n$ 矩阵写成：

$$
A=
\begin{bmatrix}
a_{11}&a_{12}&\cdots&a_{1n}\\
a_{21}&a_{22}&\cdots&a_{2n}\\
\vdots&\vdots&\ddots&\vdots\\
a_{m1}&a_{m2}&\cdots&a_{mn}
\end{bmatrix}
$$

其中：

- $m$ 是行数，对应输出向量的分量数量。
- $n$ 是列数，对应输入向量的分量数量。
- $a_{ij}$ 表示第 $i$ 行、第 $j$ 列的元素。

若 $A$ 是 $m\times n$ 矩阵，则它可以与 $n$ 维列向量相乘，并产生 $m$ 维输出：

$$
A:\mathbb R^n\rightarrow\mathbb R^m
$$

例如：

$$
A=
\begin{bmatrix}
1&2&0\\
-1&3&4
\end{bmatrix}
$$

是一个 $2\times3$ 矩阵。它接收三维输入，产生二维输出：

$$
A
\begin{bmatrix}x_1\\x_2\\x_3\end{bmatrix}
=
\begin{bmatrix}
x_1+2x_2\\
-x_1+3x_2+4x_3
\end{bmatrix}
$$

### 1.2 行视角与列视角

同一个矩阵可以按行或按列观察：

$$
A=
\begin{bmatrix}
-&\boldsymbol r_1^\mathsf T&-\\
-&\boldsymbol r_2^\mathsf T&-\\
&\vdots&\\
-&\boldsymbol r_m^\mathsf T&-
\end{bmatrix}
=
\begin{bmatrix}
|&|&&|\\
\boldsymbol a_1&\boldsymbol a_2&\cdots&\boldsymbol a_n\\
|&|&&|
\end{bmatrix}
$$

- 行视角强调每个输出分量怎样由输入分量计算出来。
- 列视角强调输出怎样由矩阵的列向量线性组合得到。

在本阶段，列视角尤其重要，因为第 $j$ 列就是第 $j$ 个标准基向量经过变换后的结果。

### 1.3 矩阵的基本运算

尺寸相同的矩阵可以相加：

$$
(A+B)_{ij}=a_{ij}+b_{ij}
$$

矩阵可以与标量相乘：

$$
(cA)_{ij}=ca_{ij}
$$

矩阵的转置将行列互换：

$$
(A^\mathsf T)_{ij}=a_{ji}
$$

例如：

$$
\begin{bmatrix}
1&2&3\\
4&5&6
\end{bmatrix}^{\mathsf T}
=
\begin{bmatrix}
1&4\\
2&5\\
3&6
\end{bmatrix}
$$

转置会在后续阶段用于内积、正交投影、最小二乘和对称矩阵。

## 2. 线性变换：保持线性组合的规则

### 2.1 定义

设：

$$
T:\mathbb R^n\rightarrow\mathbb R^m
$$

若对任意向量 $\boldsymbol u,\boldsymbol v$ 和任意标量 $c$，都有：

$$
T(\boldsymbol u+\boldsymbol v)
=T(\boldsymbol u)+T(\boldsymbol v)
$$

以及：

$$
T(c\boldsymbol u)=cT(\boldsymbol u)
$$

则称 $T$ 为线性变换。

这两个条件也可以合并为：

$$
T(a\boldsymbol u+b\boldsymbol v)
=aT(\boldsymbol u)+bT(\boldsymbol v)
$$

也就是说，线性变换会保留线性组合的结构。

### 2.2 线性变换必然保持的结构

任何线性变换都满足：

$$
T(\boldsymbol0)=\boldsymbol0
$$

证明很直接：

$$
T(\boldsymbol0)
=T(0\boldsymbol v)
=0T(\boldsymbol v)
=\boldsymbol0
$$

此外，线性变换还具有以下几何特征：

- 经过原点的直线变换后仍为经过原点的直线，或者被压缩成原点。
- 平行关系保持不变。
- 原点不会被移动。
- 向量之间的线性组合关系保持不变。

长度和角度不一定保持。缩放会改变长度，剪切会改变角度。

### 2.3 如何判断一个规则是否线性

例 1：

$$
T(x,y)=(2x+y,x-y)
$$

每个输出分量都是输入分量的齐次一次组合，没有常数项，因此它是线性变换。其矩阵为：

$$
A=
\begin{bmatrix}
2&1\\
1&-1
\end{bmatrix}
$$

例 2：

$$
T(x,y)=(2x+y+1,x-y)
$$

因为：

$$
T(0,0)=(1,0)\neq(0,0)
$$

所以它不是线性变换。

例 3：

$$
T(x,y)=(x^2,y)
$$

因为一般有：

$$
T(c\boldsymbol x)\neq cT(\boldsymbol x)
$$

所以它不是线性变换。

快速排除规则：若变换包含常数平移项、变量乘积、平方、绝对值等非线性结构，通常不是线性变换。不过最终判断仍应回到线性条件。

## 3. 基向量的像决定整个线性变换

### 3.1 核心推导

二维空间中的任意向量都可以写成：

$$
\boldsymbol x=x_1\boldsymbol e_1+x_2\boldsymbol e_2
$$

其中：

$$
\boldsymbol e_1=
\begin{bmatrix}1\\0\end{bmatrix},
\qquad
\boldsymbol e_2=
\begin{bmatrix}0\\1\end{bmatrix}
$$

利用线性性质：

$$
T(\boldsymbol x)
=x_1T(\boldsymbol e_1)+x_2T(\boldsymbol e_2)
$$

因此，只要知道 $T(\boldsymbol e_1)$ 和 $T(\boldsymbol e_2)$，就知道 $T$ 对所有二维向量的作用。

将两个基向量的像作为列向量排列：

$$
A=
\begin{bmatrix}
|&|\\
T(\boldsymbol e_1)&T(\boldsymbol e_2)\\
|&|
\end{bmatrix}
$$

便得到线性变换在标准基下的矩阵。

在 $n$ 维空间中同理：

$$
A=
\begin{bmatrix}
T(\boldsymbol e_1)&T(\boldsymbol e_2)&\cdots&T(\boldsymbol e_n)
\end{bmatrix}
$$

### 3.2 由几何规则构造矩阵

设变换 $T$ 将：

$$
\boldsymbol e_1\mapsto
\begin{bmatrix}2\\1\end{bmatrix},
\qquad
\boldsymbol e_2\mapsto
\begin{bmatrix}-1\\3\end{bmatrix}
$$

那么：

$$
A=
\begin{bmatrix}
2&-1\\
1&3
\end{bmatrix}
$$

对于：

$$
\boldsymbol x=
\begin{bmatrix}4\\2\end{bmatrix}
=4\boldsymbol e_1+2\boldsymbol e_2
$$

有：

$$
T(\boldsymbol x)
=4T(\boldsymbol e_1)+2T(\boldsymbol e_2)
=
\begin{bmatrix}6\\10\end{bmatrix}
$$

这体现了最重要的原则：输入使用什么系数组合原来的基，输出就使用同样的系数组合变换后的基。

## 4. 矩阵乘向量：同一运算的四个视角

设：

$$
A=
\begin{bmatrix}
2&1\\
0&3
\end{bmatrix},
\qquad
\boldsymbol x=
\begin{bmatrix}4\\-1\end{bmatrix}
$$

那么：

$$
A\boldsymbol x=
\begin{bmatrix}7\\-3\end{bmatrix}
$$

### 4.1 分量计算视角

每个输出分量由矩阵的一行与输入向量共同计算：

$$
A\boldsymbol x=
\begin{bmatrix}
2\cdot4+1\cdot(-1)\\
0\cdot4+3\cdot(-1)
\end{bmatrix}
=
\begin{bmatrix}7\\-3\end{bmatrix}
$$

### 4.2 列向量线性组合视角

$$
A\boldsymbol x
=4
\begin{bmatrix}2\\0\end{bmatrix}
-
\begin{bmatrix}1\\3\end{bmatrix}
=
\begin{bmatrix}7\\-3\end{bmatrix}
$$

输入向量的分量，就是组合矩阵各列的系数。

### 4.3 线性变换视角

矩阵 $A$ 将整个平面中的每一点 $\boldsymbol x$ 送到新位置 $A\boldsymbol x$。它不仅改变一个向量，而是同时规定整个坐标网格如何变化。

### 4.4 方程组视角

若已知输出 $\boldsymbol b$，方程：

$$
A\boldsymbol x=\boldsymbol b
$$

是在寻找一个输入，使其经过变换后到达 $\boldsymbol b$。这将成为第三阶段研究线性方程组的主线。

### 4.5 输入输出尺寸

若 $A$ 是 $m\times n$ 矩阵，则：

$$
\underbrace{A}_{m\times n}
\underbrace{\boldsymbol x}_{n\times1}
=
\underbrace{\boldsymbol b}_{m\times1}
$$

内部尺寸 $n$ 必须一致，输出尺寸由外部的 $m$ 决定。

## 5. 二维空间中的基本线性变换

学习每一种变换时，都应完成三件事：

1. 计算 $A\boldsymbol e_1$ 和 $A\boldsymbol e_2$。
2. 画出两个变换后的基向量。
3. 画出单位正方形变换后的图形。

### 5.1 恒等变换

$$
I=
\begin{bmatrix}
1&0\\
0&1
\end{bmatrix}
$$

它保持所有向量不变：

$$
I\boldsymbol x=\boldsymbol x
$$

### 5.2 缩放

$$
S=
\begin{bmatrix}
s_x&0\\
0&s_y
\end{bmatrix}
$$

- $s_x$ 控制水平方向缩放。
- $s_y$ 控制竖直方向缩放。
- 负缩放还会使对应方向发生反射。
- 某个缩放系数为零时，对应方向会被完全压扁。

例如：

$$
\begin{bmatrix}
2&0\\
0&\frac12
\end{bmatrix}
$$

将水平方向放大 2 倍，将竖直方向缩小到原来的一半。

### 5.3 旋转

逆时针旋转 $\theta$ 的矩阵为：

$$
R(\theta)=
\begin{bmatrix}
\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta
\end{bmatrix}
$$

它来自两个标准基向量旋转后的结果：

$$
R(\theta)\boldsymbol e_1=
\begin{bmatrix}\cos\theta\\\sin\theta\end{bmatrix},
\qquad
R(\theta)\boldsymbol e_2=
\begin{bmatrix}-\sin\theta\\\cos\theta\end{bmatrix}
$$

例如逆时针旋转 $90^\circ$：

$$
R\left(\frac\pi2\right)=
\begin{bmatrix}
0&-1\\
1&0
\end{bmatrix}
$$

### 5.4 反射

关于 $x$ 轴反射：

$$
F_x=
\begin{bmatrix}
1&0\\
0&-1
\end{bmatrix}
$$

关于 $y$ 轴反射：

$$
F_y=
\begin{bmatrix}
-1&0\\
0&1
\end{bmatrix}
$$

关于直线 $y=x$ 反射：

$$
F_{y=x}=
\begin{bmatrix}
0&1\\
1&0
\end{bmatrix}
$$

最后一个矩阵会交换向量的两个坐标。

### 5.5 剪切

水平剪切：

$$
H_x=
\begin{bmatrix}
1&k\\
0&1
\end{bmatrix}
$$

它将：

$$
(x,y)\mapsto(x+ky,y)
$$

竖直高度不变，水平方向的移动量与高度成正比。单位正方形会变成平行四边形。

竖直剪切：

$$
H_y=
\begin{bmatrix}
1&0\\
k&1
\end{bmatrix}
$$

### 5.6 投影

投影到 $x$ 轴：

$$
P_x=
\begin{bmatrix}
1&0\\
0&0
\end{bmatrix}
$$

它将：

$$
(x,y)\mapsto(x,0)
$$

投影到 $y$ 轴：

$$
P_y=
\begin{bmatrix}
0&0\\
0&1
\end{bmatrix}
$$

投影会丢失一个方向的信息，因此通常无法撤销。第四阶段将研究投影到任意直线和子空间的统一公式。

## 6. 矩阵乘法：线性变换的复合

### 6.1 为什么出现矩阵乘法

设：

$$
B:\mathbb R^n\rightarrow\mathbb R^p,
\qquad
A:\mathbb R^p\rightarrow\mathbb R^m
$$

先执行 $B$，再执行 $A$：

$$
\boldsymbol x
\xrightarrow{B}
B\boldsymbol x
\xrightarrow{A}
A(B\boldsymbol x)
$$

复合变换由矩阵 $AB$ 表示：

$$
(AB)\boldsymbol x=A(B\boldsymbol x)
$$

因此，$AB$ 的实际执行顺序是：

$$
\boxed{\text{先 }B\text{，后 }A}
$$

### 6.2 尺寸规则

若 $A$ 是 $m\times p$ 矩阵，$B$ 是 $p\times n$ 矩阵，则：

$$
\underbrace{A}_{m\times p}
\underbrace{B}_{p\times n}
=
\underbrace{AB}_{m\times n}
$$

只有内部尺寸相同，乘法才有定义。

### 6.3 列视角

若：

$$
B=
\begin{bmatrix}
\boldsymbol b_1&\boldsymbol b_2&\cdots&\boldsymbol b_n
\end{bmatrix}
$$

那么：

$$
AB=
\begin{bmatrix}
A\boldsymbol b_1&A\boldsymbol b_2&\cdots&A\boldsymbol b_n
\end{bmatrix}
$$

也就是说，$AB$ 的第 $j$ 列，是 $A$ 对 $B$ 的第 $j$ 列进行变换的结果。

### 6.4 元素计算视角

乘积矩阵的元素为：

$$
(AB)_{ij}=\sum_{k=1}^{p}a_{ik}b_{kj}
$$

这就是通常所说的“$A$ 的第 $i$ 行乘 $B$ 的第 $j$ 列”。它是复合变换在坐标层面的计算规则，而不是矩阵乘法存在的根本原因。

### 6.5 例：先缩放，再旋转

先沿 $x$ 方向放大 2 倍：

$$
S=
\begin{bmatrix}
2&0\\
0&1
\end{bmatrix}
$$

再逆时针旋转 $90^\circ$：

$$
R=
\begin{bmatrix}
0&-1\\
1&0
\end{bmatrix}
$$

整体矩阵为：

$$
RS=
\begin{bmatrix}
0&-1\\
2&0
\end{bmatrix}
$$

对于 $\boldsymbol x=(1,1)^\mathsf T$：

$$
S\boldsymbol x=
\begin{bmatrix}2\\1\end{bmatrix},
\qquad
R(S\boldsymbol x)=
\begin{bmatrix}-1\\2\end{bmatrix}
$$

### 6.6 为什么矩阵乘法一般不交换

对于上面的矩阵：

$$
SR=
\begin{bmatrix}
0&-2\\
1&0
\end{bmatrix}
\neq
RS
$$

先旋转再缩放，与先缩放再旋转，会产生不同结果。因此一般有：

$$
AB\neq BA
$$

不过矩阵乘法满足结合律：

$$
(AB)C=A(BC)
$$

因为连续三个变换的整体结果不依赖于先把哪两个组合成一个整体。

矩阵乘法还满足分配律：

$$
A(B+C)=AB+AC
$$

## 7. 逆矩阵：撤销一个线性变换

### 7.1 逆矩阵的定义

对于方阵 $A$，若存在矩阵 $A^{-1}$，使：

$$
A^{-1}A=AA^{-1}=I
$$

则称 $A$ 可逆，$A^{-1}$ 是 $A$ 的逆矩阵。

几何上，$A^{-1}$ 会撤销 $A$：

$$
\boldsymbol x
\xrightarrow{A}
A\boldsymbol x
\xrightarrow{A^{-1}}
\boldsymbol x
$$

### 7.2 哪些变换可逆

- 非零缩放可逆，逆变换是按倒数缩放。
- 旋转可逆，逆变换是向相反方向旋转。
- 反射可逆，再反射一次即可恢复。
- 剪切可逆，使用相反参数剪切即可恢复。
- 投影不可逆，因为被压掉的方向信息已经丢失。

可逆性的本质不是“能否套公式”，而是：

> 不同输入是否始终产生不同输出，并且每个目标输出是否都能由某个输入产生？

### 7.3 二阶矩阵的逆

设：

$$
A=
\begin{bmatrix}
a&b\\
c&d
\end{bmatrix}
$$

当 $ad-bc\neq0$ 时：

$$
A^{-1}
=
\frac{1}{ad-bc}
\begin{bmatrix}
d&-b\\
-c&a
\end{bmatrix}
$$

其中 $ad-bc$ 将在第五阶段作为行列式深入研究。它为零时，两个变换后的基方向线性相关，空间被压缩到较低维，因此无法恢复原输入。

### 7.4 例：求逆并验证

设：

$$
A=
\begin{bmatrix}
2&1\\
1&1
\end{bmatrix}
$$

因为 $2\cdot1-1\cdot1=1$，所以：

$$
A^{-1}=
\begin{bmatrix}
1&-1\\
-1&2
\end{bmatrix}
$$

验证：

$$
AA^{-1}
=
\begin{bmatrix}
2&1\\
1&1
\end{bmatrix}
\begin{bmatrix}
1&-1\\
-1&2
\end{bmatrix}
=I
$$

### 7.5 逆矩阵与方程求解

若：

$$
A\boldsymbol x=\boldsymbol b
$$

并且 $A$ 可逆，则：

$$
\boldsymbol x=A^{-1}\boldsymbol b
$$

这给出了理论上的解。实际数值计算通常使用消元或矩阵分解，而不是先显式求逆，因为后者往往计算量更大且数值稳定性更差。

### 7.6 逆的顺序

若 $A$ 和 $B$ 都可逆，则：

$$
(AB)^{-1}=B^{-1}A^{-1}
$$

要撤销“先 $B$ 后 $A$”，必须先撤销 $A$，再撤销 $B$，所以顺序反转。

## 8. 换基：改变坐标语言，不改变对象

### 8.1 基矩阵

设二维空间的一组基为：

$$
\mathcal B=(\boldsymbol b_1,\boldsymbol b_2)
$$

把基向量作为列组成基矩阵：

$$
P=
\begin{bmatrix}
|&|\\
\boldsymbol b_1&\boldsymbol b_2\\
|&|
\end{bmatrix}
$$

若向量 $\boldsymbol x$ 在基 $\mathcal B$ 下的坐标为 $[\boldsymbol x]_{\mathcal B}$，则：

$$
\boldsymbol x=P[\boldsymbol x]_{\mathcal B}
$$

这表示使用新基坐标中的系数，对新基向量进行线性组合，得到标准坐标中的向量。

反过来：

$$
[\boldsymbol x]_{\mathcal B}=P^{-1}\boldsymbol x
$$

### 8.2 例：在两套坐标之间转换

设：

$$
\boldsymbol b_1=
\begin{bmatrix}1\\1\end{bmatrix},
\qquad
\boldsymbol b_2=
\begin{bmatrix}1\\-1\end{bmatrix}
$$

于是：

$$
P=
\begin{bmatrix}
1&1\\
1&-1
\end{bmatrix}
$$

若：

$$
[\boldsymbol x]_{\mathcal B}=
\begin{bmatrix}3\\1\end{bmatrix}
$$

则标准坐标为：

$$
\boldsymbol x
=P[\boldsymbol x]_{\mathcal B}
=3\boldsymbol b_1+\boldsymbol b_2
=
\begin{bmatrix}4\\2\end{bmatrix}
$$

### 8.3 线性变换在新基下的矩阵

设线性变换在标准基下由 $A$ 表示。若输入和输出都改用基 $\mathcal B$，则：

$$
[T]_{\mathcal B}=P^{-1}AP
$$

理解这个公式要按执行顺序阅读：

1. $P$：把新基坐标转换为标准坐标。
2. $A$：在标准坐标中执行线性变换。
3. $P^{-1}$：把结果转换回新基坐标。

因此：

$$
[\boldsymbol x]_{\mathcal B}
\xrightarrow{P}
\boldsymbol x
\xrightarrow{A}
A\boldsymbol x
\xrightarrow{P^{-1}}
[A\boldsymbol x]_{\mathcal B}
$$

$A$ 与 $P^{-1}AP$ 描述的是同一个线性变换，只是使用了不同的坐标语言。第五阶段的对角化，就是寻找一种让变换矩阵尽可能简单的基。

### 8.4 主动变换与被动换基

这两个概念很容易混淆：

- **主动变换**：向量真的发生了变化，例如旋转向量。
- **被动换基**：向量本身不变，只改变描述它的坐标系统。

若相机转动后重新描述同一个物体的位置，通常涉及被动坐标变化；若物体本身在固定坐标系中转动，则是主动变换。

## 9. 线性变换与仿射变换

平移：

$$
T(\boldsymbol x)=A\boldsymbol x+\boldsymbol t
$$

当 $\boldsymbol t\neq\boldsymbol0$ 时，它不是普通线性变换，因为：

$$
T(\boldsymbol0)=\boldsymbol t\neq\boldsymbol0
$$

这种“线性变换加平移”称为仿射变换。计算机图形学常使用齐次坐标，把仿射变换统一写成更高维矩阵：

$$
\begin{bmatrix}
x'\\y'\\1
\end{bmatrix}
=
\begin{bmatrix}
a&b&t_x\\
c&d&t_y\\
0&0&1
\end{bmatrix}
\begin{bmatrix}
x\\y\\1
\end{bmatrix}
$$

这并没有让二维平移突然成为二维线性变换，而是把问题嵌入三维齐次坐标空间中统一处理。

## 10. 一个贯穿本阶段的例子

考虑矩阵：

$$
A=
\begin{bmatrix}
2&1\\
1&1
\end{bmatrix}
$$

### 10.1 尺寸与输入输出

$A$ 是 $2\times2$ 矩阵，因此表示：

$$
A:\mathbb R^2\rightarrow\mathbb R^2
$$

### 10.2 基向量的像

$$
A\boldsymbol e_1=
\begin{bmatrix}2\\1\end{bmatrix},
\qquad
A\boldsymbol e_2=
\begin{bmatrix}1\\1\end{bmatrix}
$$

所以单位正方形会变成由 $(2,1)$ 和 $(1,1)$ 张成的平行四边形。

### 10.3 任意输入

对于：

$$
\boldsymbol x=
\begin{bmatrix}3\\-1\end{bmatrix}
$$

有：

$$
A\boldsymbol x
=3
\begin{bmatrix}2\\1\end{bmatrix}
-
\begin{bmatrix}1\\1\end{bmatrix}
=
\begin{bmatrix}5\\2\end{bmatrix}
$$

### 10.4 可逆性

两个列向量不共线，变换没有把整个平面压到一条直线上，因此可以恢复输入：

$$
A^{-1}=
\begin{bmatrix}
1&-1\\
-1&2
\end{bmatrix}
$$

验证：

$$
A^{-1}
\begin{bmatrix}5\\2\end{bmatrix}
=
\begin{bmatrix}3\\-1\end{bmatrix}
$$

### 10.5 不同坐标系统

如果找到一组更适合该变换的基，那么同一个 $A$ 可以表示成更简单的矩阵。这个想法将在特征向量与对角化中发展为：寻找变换保持方向不变的特殊基。

这个例子把本阶段的主要概念连接起来：

$$
\boxed{
\text{矩阵尺寸}
\rightarrow
\text{基向量的像}
\rightarrow
\text{任意向量的像}
\rightarrow
\text{可逆性}
\rightarrow
\text{换基}
}
$$

## 11. 常见误区

### 误区 1：矩阵只是数表

矩阵中的数字是线性变换在特定输入基和输出基下的坐标表示。矩阵背后还有一个与坐标选择无关的线性变换。

### 误区 2：矩阵的行表示基向量的像

在列向量约定下，矩阵的**列**表示输入基向量的像。行用于计算各个输出分量。

### 误区 3：只要公式中有矩阵就是线性变换

$A\boldsymbol x+\boldsymbol b$ 在 $\boldsymbol b\neq\boldsymbol0$ 时是仿射变换，不是普通线性变换。

### 误区 4：矩阵乘法可以交换

矩阵乘法描述有顺序的变换复合。改变顺序，通常会改变结果。

### 误区 5：$AB$ 表示先 $A$ 后 $B$

列向量写在右侧时：

$$
AB\boldsymbol x=A(B\boldsymbol x)
$$

因此是先 $B$ 后 $A$。

### 误区 6：不可逆只是因为不会求逆公式

不可逆是结构问题。只要变换压缩掉某个方向，使不同输入产生相同输出，任何公式都无法恢复已丢失的信息。

### 误区 7：换基改变了向量

换基通常只改变坐标表示，向量本身没有改变。主动变换才会改变向量。

### 误区 8：矩阵的尺寸只影响计算格式

矩阵尺寸明确规定输入和输出所在的空间。$m\times n$ 矩阵接收 $n$ 维输入，产生 $m$ 维输出。

## 12. 应用连接

- **计算机图形学**：旋转、缩放、投影和齐次坐标用于图像与三维模型变换。
- **机器人与导航**：不同传感器、关节和世界坐标系之间需要连续坐标变换。
- **数据处理**：矩阵将一个特征向量转换成新的特征表示。
- **神经网络**：线性层执行 $A\boldsymbol x+\boldsymbol b$，再与非线性激活组合。
- **物理系统**：状态向量经过系统矩阵演化到下一状态。
- **图像处理**：一幅灰度图像可以表示成矩阵，颜色变换和几何变换也可由矩阵描述。
- **经济模型**：投入产出矩阵描述行业之间的线性依赖关系。
- **编码与压缩**：换基可以把信息转换到更适合存储或分析的坐标系统。

## 13. 阶段练习

### 基础题

1. 设 $A$ 是 $3\times2$ 矩阵。它接收几维输入，产生几维输出？写出对应的映射形式。

2. 计算：

   $$
   \begin{bmatrix}
   1&2\\
   -1&3
   \end{bmatrix}
   \begin{bmatrix}4\\-2\end{bmatrix}
   $$

3. 将第 2 题的乘法写成矩阵列向量的线性组合。

4. 已知线性变换满足：

   $$
   T(\boldsymbol e_1)=
   \begin{bmatrix}1\\2\end{bmatrix},
   \qquad
   T(\boldsymbol e_2)=
   \begin{bmatrix}-3\\1\end{bmatrix}
   $$

   写出 $T$ 的矩阵，并计算 $T(2,-1)$。

5. 写出关于 $x$ 轴反射的矩阵，并计算向量 $(3,-2)$ 反射后的结果。

6. 写出逆时针旋转 $90^\circ$ 的矩阵，并计算向量 $(2,1)$ 旋转后的结果。

7. 判断下面哪些矩阵乘法有定义，并写出结果尺寸：

   - $(2\times3)(3\times4)$
   - $(3\times2)(3\times1)$
   - $(4\times2)(2\times5)$

### 理解题

8. 判断下列规则是否为线性变换，并说明原因：

   $$
   T_1(x,y)=(x+2y,3x-y)
   $$

   $$
   T_2(x,y)=(x+1,y)
   $$

   $$
   T_3(x,y)=(|x|,y)
   $$

9. 为什么知道 $T(\boldsymbol e_1),\dots,T(\boldsymbol e_n)$ 就能知道 $T$ 对所有向量的作用？

10. 为什么投影到 $x$ 轴的变换不可逆？请给出两个不同输入产生相同输出的例子。

11. 用几何语言解释为什么矩阵乘法一般不满足交换律。

12. 为什么 $(AB)^{-1}=B^{-1}A^{-1}$，而不是 $A^{-1}B^{-1}$？

### 综合题

13. 先沿 $x$ 方向放大 3 倍，再关于直线 $y=x$ 反射。

   1. 分别写出两个变换矩阵。
   2. 写出整体变换矩阵。
   3. 计算向量 $(1,2)$ 的最终结果。
   4. 比较改变执行顺序后的结果。

14. 设：

   $$
   A=
   \begin{bmatrix}
   1&2\\
   2&4
   \end{bmatrix}
   $$

   1. 画出或描述两个标准基向量的像。
   2. 判断该变换是否可逆。
   3. 找出两个不同输入，使它们经过 $A$ 后产生相同输出。
   4. 解释不可逆的几何原因。

15. 设基：

   $$
   \mathcal B=\left(
   \begin{bmatrix}1\\1\end{bmatrix},
   \begin{bmatrix}1\\-1\end{bmatrix}
   \right)
   $$

   1. 写出基矩阵 $P$。
   2. 已知 $[\boldsymbol x]_{\mathcal B}=(2,3)^\mathsf T$，求 $\boldsymbol x$ 的标准坐标。
   3. 已知 $\boldsymbol y=(6,2)^\mathsf T$，求 $[\boldsymbol y]_{\mathcal B}$。

16. 设线性变换：

   $$
   T(x,y)=(x+y,2x-y)
   $$

   1. 写出标准矩阵 $A$。
   2. 求 $T(3,-1)$。
   3. 判断 $T$ 是否可逆，并求 $A^{-1}$。
   4. 从输出 $(2,7)$ 恢复输入。

## 14. 参考答案

1. 接收二维输入，产生三维输出：

   $$
   A:\mathbb R^2\rightarrow\mathbb R^3
   $$

2. 结果为：

   $$
   \begin{bmatrix}
   1&2\\
   -1&3
   \end{bmatrix}
   \begin{bmatrix}4\\-2\end{bmatrix}
   =
   \begin{bmatrix}0\\-10\end{bmatrix}
   $$

3. 列组合为：

   $$
   4\begin{bmatrix}1\\-1\end{bmatrix}
   -2\begin{bmatrix}2\\3\end{bmatrix}
   =
   \begin{bmatrix}0\\-10\end{bmatrix}
   $$

4. 矩阵为：

   $$
   A=
   \begin{bmatrix}
   1&-3\\
   2&1
   \end{bmatrix}
   $$

   因此：

   $$
   T(2,-1)=
   2\begin{bmatrix}1\\2\end{bmatrix}
   -\begin{bmatrix}-3\\1\end{bmatrix}
   =
   \begin{bmatrix}5\\3\end{bmatrix}
   $$

5. 关于 $x$ 轴反射矩阵为：

   $$
   \begin{bmatrix}1&0\\0&-1\end{bmatrix}
   $$

   所以 $(3,-2)$ 变为 $(3,2)$。

6. 旋转矩阵为：

   $$
   \begin{bmatrix}0&-1\\1&0\end{bmatrix}
   $$

   所以 $(2,1)$ 变为 $(-1,2)$。

7. 第一项有定义，结果为 $2\times4$；第二项无定义，因为内部尺寸 $2$ 与 $3$ 不同；第三项有定义，结果为 $4\times5$。

8. $T_1$ 是线性变换，因为每个输出分量都是输入分量的齐次一次组合。$T_2$ 不是，因为 $T_2(0,0)=(1,0)$。$T_3$ 不是，例如 $T_3(-1,0)=(1,0)$，但 $-T_3(1,0)=(-1,0)$，不满足齐次性。

9. 任意向量都可以表示为：

   $$
   \boldsymbol x=x_1\boldsymbol e_1+\cdots+x_n\boldsymbol e_n
   $$

   线性性质保证：

   $$
   T(\boldsymbol x)
   =x_1T(\boldsymbol e_1)+\cdots+x_nT(\boldsymbol e_n)
   $$

10. 投影会丢失 $y$ 分量。例如 $(1,2)$ 和 $(1,5)$ 都被映射到 $(1,0)$，所以无法根据输出唯一恢复输入。

11. 矩阵乘法表示有顺序的几何操作。先旋转后水平缩放，与先水平缩放后旋转，会沿不同方向拉伸图形，因此结果通常不同。

12. $AB$ 表示先 $B$ 后 $A$。撤销时必须先撤销最后执行的 $A$，再撤销 $B$，所以逆矩阵顺序反转。

13. 沿 $x$ 方向放大 3 倍：

   $$
   S=
   \begin{bmatrix}3&0\\0&1\end{bmatrix}
   $$

   关于 $y=x$ 反射：

   $$
   F=
   \begin{bmatrix}0&1\\1&0\end{bmatrix}
   $$

   先缩放再反射，整体矩阵为：

   $$
   FS=
   \begin{bmatrix}0&1\\3&0\end{bmatrix}
   $$

   因此：

   $$
   FS\begin{bmatrix}1\\2\end{bmatrix}
   =
   \begin{bmatrix}2\\3\end{bmatrix}
   $$

   改变顺序后：

   $$
   SF=
   \begin{bmatrix}0&3\\1&0\end{bmatrix},
   \qquad
   SF\begin{bmatrix}1\\2\end{bmatrix}
   =
   \begin{bmatrix}6\\1\end{bmatrix}
   $$

14. 两个标准基向量的像分别是 $(1,2)$ 和 $(2,4)$，二者共线。该变换不可逆，因为整个平面被压到由 $(1,2)$ 张成的直线上。例如：

   $$
   A\begin{bmatrix}2\\0\end{bmatrix}
   =
   \begin{bmatrix}2\\4\end{bmatrix}
   =
   A\begin{bmatrix}0\\1\end{bmatrix}
   $$

15. 基矩阵为：

   $$
   P=
   \begin{bmatrix}
   1&1\\
   1&-1
   \end{bmatrix}
   $$

   标准坐标为：

   $$
   \boldsymbol x
   =P\begin{bmatrix}2\\3\end{bmatrix}
   =
   \begin{bmatrix}5\\-1\end{bmatrix}
   $$

   对 $\boldsymbol y=(6,2)^\mathsf T$，解：

   $$
   c_1+c_2=6,
   \qquad
   c_1-c_2=2
   $$

   得 $c_1=4,c_2=2$，所以：

   $$
   [\boldsymbol y]_{\mathcal B}=
   \begin{bmatrix}4\\2\end{bmatrix}
   $$

16. 标准矩阵为：

   $$
   A=
   \begin{bmatrix}
   1&1\\
   2&-1
   \end{bmatrix}
   $$

   有：

   $$
   T(3,-1)=
   \begin{bmatrix}2\\7\end{bmatrix}
   $$

   因为 $1\cdot(-1)-1\cdot2=-3\neq0$，所以可逆：

   $$
   A^{-1}
   =
   \begin{bmatrix}
   \frac13&\frac13\\
   \frac23&-\frac13
   \end{bmatrix}
   $$

   因此从输出 $(2,7)$ 恢复输入：

   $$
   A^{-1}
   \begin{bmatrix}2\\7\end{bmatrix}
   =
   \begin{bmatrix}3\\-1\end{bmatrix}
   $$

## 15. 阶段检验

在不查阅资料的情况下，回答以下问题：

1. 为什么 $m\times n$ 矩阵表示从 $\mathbb R^n$ 到 $\mathbb R^m$ 的变换？
2. 判断线性变换需要检查哪两个条件？为什么线性变换必须保持原点？
3. 为什么矩阵的第 $j$ 列等于第 $j$ 个标准基向量的像？
4. 如何从分量、列组合、几何变换和方程组四个角度解释 $A\boldsymbol x$？
5. 怎样根据基向量的像写出旋转、反射、剪切或投影矩阵？
6. 为什么 $AB$ 表示先 $B$ 后 $A$？
7. 为什么矩阵乘法通常不可交换，却满足结合律？
8. 一个线性变换不可逆时，几何上发生了什么？
9. 为什么 $(AB)^{-1}=B^{-1}A^{-1}$？
10. $P$ 与 $P^{-1}$ 在换基过程中分别完成什么工作？
11. 主动变换与被动换基有什么区别？
12. 为什么二维平移不是二维线性变换，却可以用三维齐次矩阵表示？

能够清楚回答这些问题，独立完成阶段练习，并能根据一个二维矩阵画出基向量和单位正方形的变化，就可以进入第三阶段：方程组、消元、秩与四个基本子空间。
