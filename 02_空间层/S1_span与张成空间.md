# S1 span 与张成空间

> **学习定位**：本节属于第二层“空间层——有什么”。正文重点回答“哪些向量可以互相表示、存在多少个独立方向、这些方向形成什么空间”。涉及矩阵作用时只给必要的变换直觉，正式的线性变换理论留到第三层。

## 本节从什么问题产生

第一层已经会把矩阵乘向量看成列向量的线性组合。现在要把问题从“某一个线性组合怎么算”提升为：“允许系数任意变化时，**全部能得到的向量**是什么？”这就自然产生 span。

## 从线性组合到“所有可达到的向量”

给定向量 $\boldsymbol v_1,\dots,\boldsymbol v_k$，表达式

$$
c_1\boldsymbol v_1+\cdots+c_k\boldsymbol v_k
$$

称为这些向量的**线性组合**，其中 $c_1,\dots,c_k$ 是标量。

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

## 张成空间：所有可到达的向量

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

## 补充：span 与 $A\boldsymbol{x}=\boldsymbol{b}$ 的连接

设：

$$
A\boldsymbol{x}=\boldsymbol{b}.
$$

理解“列空间决定方程是否有解”的关键，是把矩阵乘法写成列向量的线性组合。

---

## 矩阵乘法就是列向量的线性组合

设：

$$
A=
\begin{bmatrix}
\boldsymbol{a}_1&
\boldsymbol{a}_2&
\cdots&
\boldsymbol{a}_n
\end{bmatrix},
$$

并且

$$
\boldsymbol{x}
=
\begin{bmatrix}
x_1\\
x_2\\
\vdots\\
x_n
\end{bmatrix}.
$$

那么：

$$
A\boldsymbol{x}
=
x_1\boldsymbol{a}_1
+x_2\boldsymbol{a}_2
+\cdots
+x_n\boldsymbol{a}_n.
$$

因此，不管怎样选择输入 $\boldsymbol{x}$，矩阵能够产生的输出，都只能是列向量

$$
\boldsymbol{a}_1,\ldots,\boldsymbol{a}_n
$$

的线性组合。

---

## 列空间就是所有可能输出的集合

列空间定义为：

$$
\operatorname{Col}(A)
=
\operatorname{span}
\left\{
\boldsymbol{a}_1,\ldots,\boldsymbol{a}_n
\right\}.
$$

也可以写成：

$$
\operatorname{Col}(A)
=
\left\{
A\boldsymbol{x}:
\boldsymbol{x}\in\mathbb R^n
\right\}.
$$

所以：

$$
\boxed{
\operatorname{Col}(A)
=
\text{矩阵 }A\text{ 能产生的全部输出}
}
$$

因此，方程

$$
A\boldsymbol{x}=\boldsymbol{b}
$$

实际上是在问：

> 有没有某个输入 $\boldsymbol{x}$，能够让矩阵 $A$ 产生目标输出 $\boldsymbol{b}$？

这等价于问：

> $\boldsymbol{b}$ 是否属于矩阵所有可能输出的集合？

也就是：

$$
\boxed{
A\boldsymbol{x}=\boldsymbol{b}\text{ 有解}
\Longleftrightarrow
\boldsymbol{b}\in\operatorname{Col}(A)
}
$$

---

## “有解”也可以理解为“能不能用列向量拼出 $\boldsymbol{b}$”

由于：

$$
A\boldsymbol{x}
=
x_1\boldsymbol{a}_1+\cdots+x_n\boldsymbol{a}_n,
$$

所以：

$$
A\boldsymbol{x}=\boldsymbol{b}
$$

等价于：

$$
x_1\boldsymbol{a}_1+\cdots+x_n\boldsymbol{a}_n
=
\boldsymbol{b}.
$$

因此求解 $\boldsymbol{x}$，其实就是在寻找一组系数：

$$
x_1,\ldots,x_n,
$$

使得矩阵的列向量能够线性组合成 $\boldsymbol{b}$。

如果能够拼出来，那么有解。

如果不能拼出来，那么无解。

---

## 一个几何例子

设：

$$
A=
\begin{bmatrix}
1&0\\
0&1\\
0&0
\end{bmatrix}.
$$

它定义：

$$
A:\mathbb R^2\rightarrow\mathbb R^3.
$$

两列分别为：

$$
\boldsymbol{a}_1=
\begin{bmatrix}
1\\
0\\
0
\end{bmatrix},
\qquad
\boldsymbol{a}_2=
\begin{bmatrix}
0\\
1\\
0
\end{bmatrix}.
$$

因此：

$$
A
\begin{bmatrix}
x_1\\
x_2
\end{bmatrix}
=
\begin{bmatrix}
x_1\\
x_2\\
0
\end{bmatrix}.
$$

所以所有输出的第三个分量都只能为 0。

也就是说：

$$
\operatorname{Col}(A)
$$

就是 $\mathbb R^3$ 中的 $xy$ 平面。

如果：

$$
\boldsymbol{b}
=
\begin{bmatrix}
3\\
2\\
0
\end{bmatrix},
$$

那么：

$$
\boldsymbol{b}\in\operatorname{Col}(A),
$$

所以方程有解。

如果：

$$
\boldsymbol{b}
=
\begin{bmatrix}
3\\
2\\
1
\end{bmatrix},
$$

由于第三个分量不为 0，因此：

$$
\boldsymbol{b}\notin\operatorname{Col}(A),
$$

所以无论怎样选择输入，都无法得到这个输出。

因此方程无解。

---

## 从线性变换角度理解

对于：

$$
A:\mathbb R^n\rightarrow\mathbb R^m,
$$

输入空间的标准基为：

$$
\boldsymbol{e}_1,\ldots,\boldsymbol{e}_n.
$$

矩阵的第 $i$ 列就是：

$$
A\boldsymbol{e}_i.
$$

因此矩阵的列向量，就是输入空间各标准基经过线性变换后的结果。

任意输入：

$$
\boldsymbol{x}
=
x_1\boldsymbol{e}_1+\cdots+x_n\boldsymbol{e}_n
$$

经过变换后：

$$
A\boldsymbol{x}
=
x_1A\boldsymbol{e}_1+\cdots+x_nA\boldsymbol{e}_n.
$$

也就是：

$$
A\boldsymbol{x}
=
x_1\boldsymbol{a}_1+\cdots+x_n\boldsymbol{a}_n.
$$

所以整个输入空间经过 $A$ 映射后能够到达的区域，就是这些列向量张成的空间，也就是列空间。

因此可以把列空间理解为：

$$
\boxed{
\text{整个输入空间经过 }A\text{ 映射后的可达输出范围}
}
$$

---

## 本节统一认识

span 不是“把几个向量放在一起”，而是这些向量的**全部线性组合所形成的集合**。因此判断 $\boldsymbol b$ 是否属于某个 span，本质上就是判断能否找到一组系数把 $\boldsymbol b$ 拼出来；当向量作为矩阵列时，这立即变成 $A\boldsymbol x=\boldsymbol b$ 是否有解。

## 下一步

有了“能够生成什么”，下一个问题就是：生成这些向量时是否存在冗余方向？这引出线性相关与线性无关。
