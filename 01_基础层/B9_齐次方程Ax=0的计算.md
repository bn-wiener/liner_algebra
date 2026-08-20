# B9 齐次方程 $A\boldsymbol{x}=\boldsymbol{0}$ 的计算

> **基础层的任务：会表示、会计算。** 本层始终围绕 $A\boldsymbol{x}=\boldsymbol{b}$ 展开。允许提前看到“空间”和“变换”的影子，但正式的空间定义与变换理论留到后两层。当前目标是：看到一个线性问题，能够把它表示出来、算出来，并读懂计算结果。

## 本章核心

非齐次方程问“怎样得到目标 $\boldsymbol{b}$”；齐次方程把右端设为零，问哪些输入会相互抵消而得到零输出。第一层先完整掌握怎样求这些解；它们的全体在第二层会正式命名为零空间 $\operatorname{Null}(A)$。

## 为什么齐次方程是核心

齐次方程是右端向量为零的矩阵方程：

$$
A\boldsymbol{x}=\boldsymbol{0}.
$$

它看起来比一般的非齐次方程：

$$
A\boldsymbol{x}=\boldsymbol{b}
$$

简单，但它揭示了矩阵最重要的一种信息：哪些输入经过矩阵后完全消失。

如果一个非零向量 $\boldsymbol{z}$ 满足：

$$
A\boldsymbol{z}=\boldsymbol{0},
$$

那么输入沿着 $\boldsymbol{z}$ 方向发生变化时，输出完全没有变化。这个方向就是矩阵无法识别的输入方向。

本章要建立的主线是：

$$
\text{齐次方程}
\longrightarrow
\text{自由变量}
\longrightarrow
\text{零空间}
\longrightarrow
\text{被消除的输入方向}.
$$

## 齐次方程的基本性质

### 齐次方程始终有零解

令：

$$
\boldsymbol{x}=\boldsymbol{0},
$$

则：

$$
A\boldsymbol{0}=\boldsymbol{0}.
$$

因此齐次方程至少有一个解，称为平凡解：

$$
\boldsymbol{x}=\boldsymbol{0}.
$$

如果还有非零解，则称这些解为非平凡解。

所以齐次方程不可能无解，它只有两种情况：

- 只有零解；
- 有零解和无穷多个非零解。

### 齐次解集对线性运算封闭

如果 $\boldsymbol{x}_1$ 和 $\boldsymbol{x}_2$ 都是齐次方程的解，那么：

$$
A\boldsymbol{x}_1=\boldsymbol{0},
\qquad
A\boldsymbol{x}_2=\boldsymbol{0}.
$$

对任意实数 $c,d$，有：

$$
A(c\boldsymbol{x}_1+d\boldsymbol{x}_2)
=cA\boldsymbol{x}_1+dA\boldsymbol{x}_2
=\boldsymbol{0}.
$$

因此 $c\boldsymbol{x}_1+d\boldsymbol{x}_2$ 仍是解。齐次解集是一个子空间，这个子空间就是零空间。

## 如何通过消元求零空间

求 $\mathcal N(A)$ 的任务就是求解：

$$
A\boldsymbol{x}=\boldsymbol{0}.
$$

基本步骤如下：

1. 写出齐次方程 $A\boldsymbol{x}=\boldsymbol{0}$。
2. 对系数矩阵 $A$ 行化简。
3. 找出主元变量和自由变量。
4. 给每个自由变量设置参数。
5. 用参数表示主元变量。
6. 将通解拆成参数与方向向量的线性组合。
7. 这些方向向量构成零空间的一组基。

由于右端是零，行变换不会产生矛盾行。真正需要判断的是：是否有自由变量。

## 例一：两个自由方向

考虑：

$$
A=
\begin{bmatrix}
1&2&-1\\
2&4&-2
\end{bmatrix}.
$$

求解：

$$
A\boldsymbol{x}=\boldsymbol{0}.
$$

设：

$$
\boldsymbol{x}=\begin{bmatrix}x_1\\x_2\\x_3\end{bmatrix}.
$$

对应方程：

$$
\begin{cases}
x_1+2x_2-x_3=0,\\
2x_1+4x_2-2x_3=0.
\end{cases}
$$

第二行是第一行的 $2$ 倍，因此消元得到：

$$
\begin{bmatrix}
1&2&-1\\
0&0&0
\end{bmatrix}.
$$

主元变量是 $x_1$，自由变量是 $x_2,x_3$。令：

$$
x_2=s,
\qquad
x_3=t.
$$

第一行给出：

$$
x_1=-2s+t.
$$

因此：

$$
\boldsymbol{x}
=
\begin{bmatrix}-2s+t\\s\\t\end{bmatrix}
=s\begin{bmatrix}-2\\1\\0\end{bmatrix}
+t\begin{bmatrix}1\\0\\1\end{bmatrix}.
$$

所以：

$$
\mathcal N(A)
=
\operatorname{span}
\left\{
\begin{bmatrix}-2\\1\\0\end{bmatrix},
\begin{bmatrix}1\\0\\1\end{bmatrix}
\right\}.
$$

零空间有两个独立方向，维数为 $2$。

验证其中一个方向：

$$
A\begin{bmatrix}-2\\1\\0\end{bmatrix}
=
\begin{bmatrix}0\\0\end{bmatrix}.
$$

另一个方向也满足同样性质：

$$
A\begin{bmatrix}1\\0\\1\end{bmatrix}
=
\begin{bmatrix}0\\0\end{bmatrix}.
$$

## 参数为什么对应零空间基向量

上例的通解是：

$$
\boldsymbol{x}=s\boldsymbol{v}_1+t\boldsymbol{v}_2,
$$

其中：

$$
\boldsymbol{v}_1=
\begin{bmatrix}-2\\1\\0\end{bmatrix},
\qquad
\boldsymbol{v}_2=
\begin{bmatrix}1\\0\\1\end{bmatrix}.
$$

参数 $s$ 只控制 $\boldsymbol{v}_1$ 方向，参数 $t$ 只控制 $\boldsymbol{v}_2$ 方向。因为 $s,t$ 可以独立取值，所以这两个方向彼此独立。

更一般地，如果齐次方程的通解为：

$$
\boldsymbol{x}
=s_1\boldsymbol{v}_1+\cdots+s_k\boldsymbol{v}_k,
$$

并且 $s_1,\ldots,s_k$ 是自由参数，那么：

$$
\mathcal N(A)=\operatorname{span}\{\boldsymbol{v}_1,\ldots,\boldsymbol{v}_k\}.
$$

在标准的参数化过程中，每个自由变量对应一个基础解向量。它们构成零空间的一组基，前提是每个参数被独立设置为 $1$，其余参数设置为 $0$。

## 例二：只有零解

考虑：

$$
A=
\begin{bmatrix}
1&2\\
3&4
\end{bmatrix}.
$$

齐次方程为：

$$
\begin{cases}
x_1+2x_2=0,\\
3x_1+4x_2=0.
\end{cases}
$$

消元：

$$
\begin{bmatrix}
1&2\\
3&4
\end{bmatrix}
\xrightarrow{R_2\leftarrow R_2-3R_1}
\begin{bmatrix}
1&2\\
0&-2
\end{bmatrix}.
$$

两列都有主元，因此没有自由变量。由第二行得到 $x_2=0$，再由第一行得到 $x_1=0$。

所以：

$$
\mathcal N(A)=\{\boldsymbol{0}\}.
$$

这说明该矩阵不会把任何非零输入方向压缩成零，也就是说它是单射。

## 进一步例子：求一组齐次解方向

考虑：

$$
A=
\begin{bmatrix}
1&2&0&1\\
0&0&1&1
\end{bmatrix}.
$$

求 $\mathcal N(A)$。

齐次方程为：

$$
\begin{cases}
x_1+2x_2+x_4=0,\\
x_3+x_4=0.
\end{cases}
$$

主元变量为 $x_1,x_3$，自由变量为 $x_2,x_4$。令：

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

通解为：

$$
\boldsymbol{x}
=
\begin{bmatrix}
-2s-t\\s\\-t\\t
\end{bmatrix}
=
s\begin{bmatrix}-2\\1\\0\\0\end{bmatrix}
+t\begin{bmatrix}-1\\0\\-1\\1\end{bmatrix}.
$$

因此：

$$
\mathcal N(A)
=
\operatorname{span}
\left\{
\begin{bmatrix}-2\\1\\0\\0\end{bmatrix},
\begin{bmatrix}-1\\0\\-1\\1\end{bmatrix}
\right\}.
$$

这两个向量线性无关，因为它们分别对应两个独立自由参数。零空间维数为 $2$，也符合：

$$
4-\operatorname{rank}(A)=4-2=2.
$$

## 如何验证一组齐次解方向

得到候选向量后，应进行三项检查。

### 每个向量确实在零空间中

对每个候选向量 $\boldsymbol{v}_i$ 计算：

$$
A\boldsymbol{v}_i.
$$

必须得到零向量。

### 候选向量彼此线性无关

如果参数来自不同自由变量，通常可以直接确认线性无关。也可以把向量作为列组成矩阵，再检查其秩。

### 向量数量等于零空间维数

若 $A$ 有 $n$ 列、秩为 $r$，应有：

$$
\text{零空间基向量数量}=n-r.
$$

如果数量少了，说明漏掉了自由变量；如果数量多了，说明方向之间存在依赖。

## 常见误区

### 误区 1：齐次方程可能无解

错误。齐次方程始终有零解。只有非齐次方程才可能因为矛盾而无解。

### 误区 2：零空间中的向量是输出向量

错误。零空间位于输入空间 $\mathbb R^n$。它的向量经过矩阵作用后才变成输出零向量。

### 误区 3：零空间只有零向量就表示矩阵没有输出

错误。它表示没有非零输入会被压缩成零，即输入不会发生混淆。矩阵仍然可以产生很多非零输出。

### 误区 4：有自由变量时只写一个解

自由变量可以任意取值，必须写出参数形式，才能表示全部零空间。

### 误区 5：参数向量不需要验证

建议始终计算 $A\boldsymbol{v}_i$。这是发现符号错误和抄写错误最快的方法。

### 误区 6：零空间基可以直接取化简矩阵的非零列

零空间基不是主元列。它来自齐次方程的参数化通解；列空间基才需要回到原矩阵取主元列。

## 本节练习

### 练习 1：判断零空间位置

若 $A$ 是 $4\times7$ 矩阵，零空间位于哪个空间？一个零空间向量有多少个分量？

### 练习 2：求零空间

求下面矩阵的零空间一组基：

$$
A=
\begin{bmatrix}
1&2&1\\
2&4&2
\end{bmatrix}.
$$

### 练习 3：判断是否有非零零空间

若 $A$ 是 $3\times4$ 矩阵，且 $\operatorname{rank}(A)=3$，判断 $\mathcal N(A)$ 是否包含非零向量，并求其维数。

### 练习 4：参数形式

已知：

$$
\operatorname{rref}(A)=
\begin{bmatrix}
1&0&2&0\\
0&1&-1&0
\end{bmatrix}.
$$

求 $\mathcal N(A)$ 的一组基。

### 练习 5：验证方向

验证：

$$
A=
\begin{bmatrix}
1&2&0&1\\
0&0&1&1
\end{bmatrix}
$$

的两个向量：

$$
\boldsymbol{v}_1=
\begin{bmatrix}-2\\1\\0\\0\end{bmatrix},
\qquad
\boldsymbol{v}_2=
\begin{bmatrix}-1\\0\\-1\\1\end{bmatrix}
$$

是否属于零空间。

### 练习 6：非齐次解的结构

设某方程 $A\boldsymbol{x}=\boldsymbol b$ 的一个特解为：

$$
\boldsymbol{x}_p=
\begin{bmatrix}1\\0\\2\end{bmatrix},
$$

且：

$$
\mathcal N(A)=
\operatorname{span}
\left\{
\begin{bmatrix}-1\\1\\0\end{bmatrix}
\right\}.
$$

写出全部解，并说明解是否唯一。

## 练习答案

### 练习 1

零空间位于输入空间：

$$
\mathcal N(A)\subseteq\mathbb R^7.
$$

每个零空间向量有 $7$ 个分量。

### 练习 2

齐次方程为：

$$
x_1+2x_2+x_3=0.
$$

令 $x_2=s,x_3=t$，则：

$$
x_1=-2s-t.
$$

所以：

$$
\mathcal N(A)
=
\operatorname{span}
\left\{
\begin{bmatrix}-2\\1\\0\end{bmatrix},
\begin{bmatrix}-1\\0\\1\end{bmatrix}
\right\}.
$$

### 练习 3

由秩—零度定理：

$$
\dim\mathcal N(A)=4-3=1.
$$

所以零空间包含非零向量。

### 练习 4

齐次方程为：

$$
x_1+2x_3=0,
\qquad
x_2-x_3=0.
$$

令 $x_3=t,x_4=s$，得到：

$$
x_1=-2t,
\qquad
x_2=t.
$$

因此：

$$
\mathcal N(A)
=
\operatorname{span}
\left\{
\begin{bmatrix}-2\\1\\1\\0\end{bmatrix},
\begin{bmatrix}0\\0\\0\\1\end{bmatrix}
\right\}.
$$

### 练习 5

分别计算：

$$
A\boldsymbol{v}_1
=
\begin{bmatrix}
-2+2\\
0
\end{bmatrix}
=
\begin{bmatrix}0\\0\end{bmatrix},
$$

$$
A\boldsymbol{v}_2
=
\begin{bmatrix}
-1+1\\
-1+1
\end{bmatrix}
=
\begin{bmatrix}0\\0\end{bmatrix}.
$$

所以二者都属于 $\mathcal N(A)$。

### 练习 6

全部解为：

$$
\boldsymbol{x}
=
\boldsymbol{x}_p+t
\begin{bmatrix}-1\\1\\0\end{bmatrix}
=
\begin{bmatrix}1\\0\\2\end{bmatrix}
+t\begin{bmatrix}-1\\1\\0\end{bmatrix},
\qquad t\in\mathbb R.
$$

零空间包含非零向量，因此解不唯一；实际上有无穷多解。

## 第一层与第二层的边界

当前必须会消元、识别自由变量、写参数向量形式、验证参数方向、判断是否只有零解。至于“齐次解全体为什么形成子空间、维数为什么等于自由变量数、怎样与其他子空间配对”，第二层再系统学习。
