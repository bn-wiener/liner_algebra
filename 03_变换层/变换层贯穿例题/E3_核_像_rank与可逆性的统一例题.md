# E3 核、像、rank 与可逆性的统一例题

## 十五、具体例子

取：

$$
A=
\begin{bmatrix}
1&2&3\\
2&4&6
\end{bmatrix}
$$

这是一个：

$$
A:\mathbb{R}^3\rightarrow\mathbb{R}^2
$$

的线性变换。

因为第二行是第一行的 2 倍，所以：

$$
\operatorname{rank}(A)=1
$$

---

### 15.1 列空间

矩阵的列向量为：

$$
\begin{bmatrix}
1\\
2
\end{bmatrix},
\quad
\begin{bmatrix}
2\\
4
\end{bmatrix},
\quad
\begin{bmatrix}
3\\
6
\end{bmatrix}
$$

后两个都是第一个的倍数。

所以：

$$
\operatorname{Col}(A)
=
\operatorname{span}
\left\{
\begin{bmatrix}
1\\
2
\end{bmatrix}
\right\}
$$

因此列空间是 $\mathbb{R}^2$ 中的一条直线。

矩阵所有输出都只能落在这个方向上。

---

### 15.2 行空间

矩阵两行为：

$$
(1,2,3)
$$

和：

$$
(2,4,6)
$$

第二行是第一行的 2 倍。

所以：

$$
\operatorname{Row}(A)
=
\operatorname{span}\{(1,2,3)\}
$$

它是 $\mathbb{R}^3$ 中的一条直线。

这意味着：

> 虽然输入空间有 3 个维度，但矩阵实际上只能感知 1 个独立方向。

---

### 15.3 零空间

求：

$$
Ax=0
$$

即：

$$
x_1+2x_2+3x_3=0
$$

因此：

$$
x_1=-2x_2-3x_3
$$

所以：

$$
x=
x_2
\begin{bmatrix}
-2\\
1\\
0
\end{bmatrix}
+
x_3
\begin{bmatrix}
-3\\
0\\
1
\end{bmatrix}
$$

于是：

$$
N(A)
=
\operatorname{span}
\left\{
\begin{bmatrix}
-2\\
1\\
0
\end{bmatrix},
\begin{bmatrix}
-3\\
0\\
1
\end{bmatrix}
\right\}
$$

因此零空间是 $\mathbb{R}^3$ 中一个二维平面。

检查正交关系：

$$
(1,2,3)
\cdot
(-2,1,0)
=
-2+2
=
0
$$

以及：

$$
(1,2,3)
\cdot
(-3,0,1)
=
-3+3
=
0
$$

所以：

$$
\operatorname{Row}(A)
\perp
N(A)
$$

---

### 15.4 左零空间

求：

$$
A^Ty=0
$$

即：

$$
\begin{bmatrix}
1&2\\
2&4\\
3&6
\end{bmatrix}
\begin{bmatrix}
y_1\\
y_2
\end{bmatrix}
=
0
$$

得到：

$$
y_1+2y_2=0
$$

所以：

$$
N(A^T)
=
\operatorname{span}
\left\{
\begin{bmatrix}
-2\\
1
\end{bmatrix}
\right\}
$$

而列空间的方向为：

$$
\begin{bmatrix}
1\\
2
\end{bmatrix}
$$

检查：

$$
\begin{bmatrix}
1\\
2
\end{bmatrix}^T
\begin{bmatrix}
-2\\
1
\end{bmatrix}
=
-2+2
=
0
$$

所以：

$$
\operatorname{Col}(A)
\perp
N(A^T)
$$

---


## 十六、这个例子的真正含义

矩阵：

$$
A:\mathbb{R}^3\rightarrow\mathbb{R}^2
$$

输入空间有 3 个维度，但是：

$$
\operatorname{rank}(A)=1
$$

说明实际上只有一个有效输入方向。

所以：

$$
\mathbb{R}^3
=
\underbrace{\operatorname{Row}(A)}_{1\text{维}}
\oplus
\underbrace{N(A)}_{2\text{维}}
$$

可以理解为：

```text
3维输入空间
   │
   ├── 1维：A 能够感知
   │
   └── 2维：A 完全看不见
```

而输出空间：

$$
\mathbb{R}^2
$$

被分解为：

$$
\mathbb{R}^2
=
\underbrace{\operatorname{Col}(A)}_{1\text{维}}
\oplus
\underbrace{N(A^T)}_{1\text{维}}
$$

也就是：

```text
2维输出空间
   │
   ├── 1维：A 能够产生
   │
   └── 1维：A 永远无法产生
```

这就是 $\operatorname{rank}(A)=1$ 的真正含义。

---


## 十八、四个基本子空间与 Ax=b

以后看到：

$$
Ax=b
$$

可以直接问三个问题。

### 18.1 有没有解

判断：

$$
b\in\operatorname{Col}(A)?
$$

如果是，则有解。

如果不是，则无解。

---

### 18.2 解是否唯一

判断：

$$
N(A)=\{0\}?
$$

如果：

$$
N(A)=\{0\}
$$

那么解至多唯一。

如果：

$$
N(A)\neq\{0\}
$$

那么只要存在一个解，就存在无穷多个解。

---

### 18.3 为什么有些 b 无法产生

输出空间可以分解为：

$$
\mathbb{R}^m
=
\operatorname{Col}(A)
\oplus
N(A^T)
$$

因此任意：

$$
b\in\mathbb{R}^m
$$

可以写成：

$$
b=b_c+b_l
$$

其中：

$$
b_c\in\operatorname{Col}(A)
$$

$$
b_l\in N(A^T)
$$

如果：

$$
b_l\neq0
$$

那么 $A$ 无法产生这一部分，所以：

$$
Ax=b
$$

无解。

---
