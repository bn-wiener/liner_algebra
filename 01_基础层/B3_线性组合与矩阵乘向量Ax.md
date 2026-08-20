# B3 线性组合与矩阵乘向量 $A\boldsymbol{x}$

> **基础层的任务：会表示、会计算。** 本层始终围绕 $A\boldsymbol{x}=\boldsymbol{b}$ 展开。允许提前看到“空间”和“变换”的影子，但正式的空间定义与变换理论留到后两层。当前目标是：看到一个线性问题，能够把它表示出来、算出来，并读懂计算结果。

## 1. 从“怎样用已有向量构造新向量”开始

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

## 2. 矩阵乘向量的核心公式

若 $A=[\boldsymbol{a}_1\ \boldsymbol{a}_2\ \cdots\ \boldsymbol{a}_n]$，则

$$
A\boldsymbol{x}=x_1\boldsymbol{a}_1+x_2\boldsymbol{a}_2+\cdots+x_n\boldsymbol{a}_n.
$$

“行乘列”告诉你怎样逐项算；“列向量线性组合”告诉你这个运算究竟是什么。

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

### 分量计算视角

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

### 列向量线性组合视角

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

### 线性变换视角

矩阵 $A$ 将整个平面中的每一点 $\boldsymbol x$ 送到新位置 $A\boldsymbol x$。它不仅改变一个向量，而是同时规定整个坐标网格如何变化。

### 方程组视角

若已知输出 $\boldsymbol b$，方程：

$$
A\boldsymbol x=\boldsymbol b
$$

是在寻找一个输入，使其经过变换后到达 $\boldsymbol b$。这将成为后续研究线性方程组的主线。

### 输入输出尺寸

若 $A$ 是 $m\times n$ 矩阵，则：

$$
\underbrace{A}_{m\times n}
\underbrace{\boldsymbol x}_{n\times1}
=
\underbrace{\boldsymbol b}_{m\times1}
$$

内部尺寸 $n$ 必须一致，输出尺寸由外部的 $m$ 决定。

## 3. 从 $A\boldsymbol{x}$ 自然产生三个问题

1. 已知 $\boldsymbol{x}$，怎样计算 $A\boldsymbol{x}$？
2. 已知目标 $\boldsymbol{b}$，能否找到 $\boldsymbol{x}$ 使 $A\boldsymbol{x}=\boldsymbol{b}$？
3. 哪些不同输入会得到相同输出？

第一层后半部分主要研究后两个问题。

当输入向量的坐标为：

```text
x =
(x₁,x₂,…,xₙ)ᵀ
```

矩阵写成列形式：

```text
A = [a₁ a₂ … aₙ]
```

那么：

```text
Ax = x₁a₁ + x₂a₂ + … + xₙaₙ
```

这个式子有两种同时成立的理解。

从线性组合角度看：

> 用 `x` 中的各个分量作为权重，组合矩阵的各列。

从线性变换角度看：

> 输入向量 `x` 经过矩阵 `A` 所表示的线性变换，产生输出向量。

所以，矩阵乘向量并不是一种孤立的计算规则，而是把下面三个层次连接起来：

```text
输入坐标
→ 组合输入基方向
→ 线性变换
→ 组合基向量的像
→ 输出坐标
```

更完整的表达是：

```text
[T(v)]C = A[v]B
```

在标准基下，通常简写成：

```text
b = Ax
```

---

## 4. 本章必须掌握

- 会用行乘列计算 $A\boldsymbol{x}$。
- 会把 $A\boldsymbol{x}$ 写成列向量线性组合。
- 会做输入/输出尺寸检查。
- 能从 $A\boldsymbol{x}=\boldsymbol{b}$ 看出未知量是组合系数 $\boldsymbol{x}$。
