# 第五步：Gram–Schmidt 正交化

> 本文对应《阶段 4 任务》的第五步。
>
> 这一节要解决的问题是：怎样把一组普通的线性无关向量，变成一组标准正交向量，同时保持它们张成的子空间不变？

## 1. 为什么需要正交化

设两个线性无关向量：

$$
a_1,
\qquad
a_2
$$

张成一个二维子空间。

只要它们线性无关，就可以作为这个子空间的一组基。但是如果它们不垂直，计算坐标和投影就不方便。

```text
普通基                         标准正交基

       a₂                            q₂
      ↗                              ↑
     /                               │
    /                                │
───●────────→ a₁              ──────●──────→ q₁

方向不重复，但彼此倾斜          长度为 1，并且互相垂直
```

对于普通基，要计算一个向量的坐标，通常需要解方程组。

对于标准正交基 $q_1,q_2$，坐标可以直接通过内积得到：

$$
c_1=q_1^Tb,
\qquad
c_2=q_2^Tb.
$$

所以我们希望完成下面的转换：

$$
\boxed{
\{a_1,a_2,\ldots,a_n\}
\longrightarrow
\{q_1,q_2,\ldots,q_n\}
}
$$

并且同时满足：

$$
q_i^Tq_j=
\begin{cases}
1, & i=j,\\
0, & i\neq j,
\end{cases}
$$

以及：

$$
\operatorname{span}\{q_1,\ldots,q_n\}
=
\operatorname{span}\{a_1,\ldots,a_n\}.
$$

完成这件事的方法就是 Gram–Schmidt 正交化。

---

## 2. 先抓住最核心的动作：减去影子

假设先保留第一个方向 $a_1$。

第二个向量 $a_2$ 通常包含两部分：

1. 沿着 $a_1$ 的部分；
2. 垂直于 $a_1$ 的新方向。

```text
                 a₂
                ↗
               /|
              / |
             /  |  新方向 u₂
            /   |
───────────●────●────────→ a₁ 的方向
           0   a₂ 在 a₁ 上的投影
```

如果从 $a_2$ 中减去它在第一个方向上的投影，剩下的部分就与第一个方向垂直：

$$
\boxed{
u_2
=
a_2-
\operatorname{proj}_{a_1}a_2
}.
$$

可以把这个过程理解为：

> $a_2$ 中与旧方向重复的部分被删除，剩下的就是它真正贡献的新方向。

Gram–Schmidt 的每一步都在重复这个动作：

```text
取一个新向量
      ↓
减去它在所有已有方向上的投影
      ↓
只保留与已有方向都垂直的新部分
      ↓
把新部分单位化
```

---

## 3. 正交化和单位化是两件不同的事

学习 Gram–Schmidt 时，经常会同时看到两类向量：

$$
u_1,u_2,\ldots,u_n
$$

和：

$$
q_1,q_2,\ldots,q_n.
$$

它们的职责不同。

### 3.1 正交化

正交化产生：

$$
u_1,u_2,\ldots,u_n,
$$

要求不同向量之间互相垂直：

$$
u_i^Tu_j=0,
\qquad
i\neq j.
$$

但它们的长度不一定为 $1$。

### 3.2 单位化

将每个非零 $u_i$ 除以自己的长度：

$$
q_i=\frac{u_i}{\|u_i\|}.
$$

这样得到：

$$
\|q_i\|=1.
$$

因为单位化只改变长度，不改变方向，所以：

$$
\operatorname{span}\{u_i\}
=
\operatorname{span}\{q_i\}.
$$

简单区分：

| 操作 | 做了什么 | 得到什么 |
|---|---|---|
| 正交化 | 减掉与旧方向重复的投影 | 互相垂直的 $u_i$ |
| 单位化 | 除以向量自身长度 | 长度为 $1$ 的 $q_i$ |

Gram–Schmidt 最终要得到的是标准正交向量 $q_i$，所以两个步骤都需要。

---

## 4. 第一个向量怎样处理

第一个向量前面没有已有方向，因此不需要减去任何投影。

直接令：

$$
u_1=a_1.
$$

然后单位化：

$$
\boxed{
q_1
=
\frac{a_1}{\|a_1\|}
}.
$$

这里要求：

$$
a_1\neq0.
$$

单位化没有创造新方向，只是把第一个方向的长度调整为 $1$。

---

## 5. 第二个向量怎样处理

现在已经有单位方向 $q_1$。

$a_2$ 在 $q_1$ 方向上的投影是：

$$
\operatorname{proj}_{q_1}a_2
=(q_1^Ta_2)q_1.
$$

从 $a_2$ 中减去这个投影：

$$
\boxed{
u_2
=
a_2-(q_1^Ta_2)q_1
}.
$$

$u_2$ 就是 $a_2$ 在旧方向之外贡献的新方向。

如果 $u_2\neq0$，再单位化：

$$
\boxed{
q_2
=
\frac{u_2}{\|u_2\|}
}.
$$

现在：

$$
q_1^Tq_2=0,
\qquad
\|q_1\|=\|q_2\|=1.
$$

---

## 6. 为什么减去投影后一定正交

由：

$$
u_2=a_2-(q_1^Ta_2)q_1
$$

计算 $q_1^Tu_2$：

$$
\begin{aligned}
q_1^Tu_2
&=q_1^T\left(a_2-(q_1^Ta_2)q_1\right)\\
&=q_1^Ta_2-(q_1^Ta_2)q_1^Tq_1.
\end{aligned}
$$

因为 $q_1$ 是单位向量：

$$
q_1^Tq_1=1,
$$

所以：

$$
\begin{aligned}
q_1^Tu_2
&=q_1^Ta_2-(q_1^Ta_2)\\
&=0.
\end{aligned}
$$

因此：

$$
\boxed{u_2\perp q_1}.
$$

这不是巧合。

$a_2$ 在 $q_1$ 方向上的全部分量已经被完整减掉，所以剩余部分在 $q_1$ 方向上的分量只能是零。

---

## 7. 一个完整的二维例子

设：

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

这两个向量线性无关，但并不正交。

### 第一步：处理 $a_1$

$$
u_1=a_1=
\begin{bmatrix}
1\\
1
\end{bmatrix}.
$$

它的长度是：

$$
\|u_1\|=\sqrt{2}.
$$

所以：

$$
q_1
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1\\
1
\end{bmatrix}.
$$

### 第二步：计算 $a_2$ 在 $q_1$ 上的投影

投影系数：

$$
q_1^Ta_2
=
\frac{1}{\sqrt{2}}.
$$

投影向量：

$$
\begin{aligned}
(q_1^Ta_2)q_1
&=
\frac{1}{\sqrt{2}}
\cdot
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1\\
1
\end{bmatrix}\\
&=
\begin{bmatrix}
1/2\\
1/2
\end{bmatrix}.
\end{aligned}
$$

### 第三步：减去重复方向

$$
\begin{aligned}
u_2
&=a_2-(q_1^Ta_2)q_1\\
&=
\begin{bmatrix}
1\\
0
\end{bmatrix}
-
\begin{bmatrix}
1/2\\
1/2
\end{bmatrix}\\
&=
\begin{bmatrix}
1/2\\
-1/2
\end{bmatrix}.
\end{aligned}
$$

### 第四步：单位化 $u_2$

$$
\|u_2\|
=
\sqrt{\frac{1}{4}+\frac{1}{4}}
=
\frac{1}{\sqrt{2}}.
$$

所以：

$$
\begin{aligned}
q_2
&=\frac{u_2}{\|u_2\|}\\
&=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1\\
-1
\end{bmatrix}.
\end{aligned}
$$

### 第五步：检查结果

长度：

$$
\|q_1\|=\|q_2\|=1.
$$

内积：

$$
\begin{aligned}
q_1^Tq_2
&=
\frac{1}{2}
\begin{bmatrix}
1&1
\end{bmatrix}
\begin{bmatrix}
1\\
-1
\end{bmatrix}\\
&=0.
\end{aligned}
$$

最终得到：

$$
\boxed{
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
\end{bmatrix}
}.
$$

它们构成 $\mathbb{R}^2$ 的一组标准正交基。

---

## 8. 为什么张成空间没有改变

正交化改变了基向量的样子，但没有改变它们能够表示的空间。

先看二维情况。

因为：

$$
q_1=\frac{a_1}{\|a_1\|},
$$

所以 $q_1$ 与 $a_1$ 方向相同。

又因为：

$$
u_2=a_2-(q_1^Ta_2)q_1,
$$

$u_2$ 是 $a_2$ 与 $q_1$ 的线性组合，因此：

$$
u_2\in\operatorname{span}\{a_1,a_2\}.
$$

反过来，把公式整理为：

$$
a_2=u_2+(q_1^Ta_2)q_1.
$$

所以 $a_2$ 也能由 $q_1$ 和 $u_2$ 表示。

因此：

$$
\operatorname{span}\{a_1,a_2\}
=
\operatorname{span}\{q_1,u_2\}.
$$

单位化不改变方向，所以：

$$
\operatorname{span}\{q_1,u_2\}
=
\operatorname{span}\{q_1,q_2\}.
$$

最终：

$$
\boxed{
\operatorname{span}\{a_1,a_2\}
=
\operatorname{span}\{q_1,q_2\}
}.
$$

可以用一句话理解：

> Gram–Schmidt 只是重新组织已有方向，没有增加新方向，也没有删除真正独立的方向。

---

## 9. 第三个向量怎样处理

假设已经得到两个标准正交方向：

$$
q_1,
\qquad
q_2.
$$

第三个向量 $a_3$ 可能同时包含：

- $q_1$ 方向上的分量；
- $q_2$ 方向上的分量；
- 与前两个方向都垂直的新分量。

因此必须同时减去前两个方向上的投影：

$$
\boxed{
u_3
=
a_3
-(q_1^Ta_3)q_1
-(q_2^Ta_3)q_2
}.
$$

然后单位化：

$$
\boxed{
q_3
=
\frac{u_3}{\|u_3\|}
}.
$$

因为 $q_1$ 和 $q_2$ 互相垂直，所以两个旧方向上的分量可以分别计算并减去。

完成后：

$$
q_1^Tq_3=0,
\qquad
q_2^Tq_3=0.
$$

---

## 10. 三维例子：逐层删除旧方向

设：

$$
a_1=
\begin{bmatrix}
1\\
0\\
0
\end{bmatrix},
\qquad
a_2=
\begin{bmatrix}
1\\
1\\
0
\end{bmatrix},
\qquad
a_3=
\begin{bmatrix}
1\\
1\\
1
\end{bmatrix}.
$$

### 处理第一个向量

$a_1$ 已经是单位向量，所以：

$$
q_1=
\begin{bmatrix}
1\\
0\\
0
\end{bmatrix}.
$$

### 处理第二个向量

$a_2$ 在 $q_1$ 上的投影为：

$$
(q_1^Ta_2)q_1
=
\begin{bmatrix}
1\\
0\\
0
\end{bmatrix}.
$$

减去投影：

$$
u_2
=
\begin{bmatrix}
1\\
1\\
0
\end{bmatrix}
-
\begin{bmatrix}
1\\
0\\
0
\end{bmatrix}
=
\begin{bmatrix}
0\\
1\\
0
\end{bmatrix}.
$$

$u_2$ 已经是单位向量，所以：

$$
q_2=
\begin{bmatrix}
0\\
1\\
0
\end{bmatrix}.
$$

### 处理第三个向量

$a_3$ 在 $q_1$ 上的投影为：

$$
(q_1^Ta_3)q_1
=
\begin{bmatrix}
1\\
0\\
0
\end{bmatrix}.
$$

$a_3$ 在 $q_2$ 上的投影为：

$$
(q_2^Ta_3)q_2
=
\begin{bmatrix}
0\\
1\\
0
\end{bmatrix}.
$$

同时减去：

$$
\begin{aligned}
u_3
&=a_3-(q_1^Ta_3)q_1-(q_2^Ta_3)q_2\\
&=
\begin{bmatrix}
1\\
1\\
1
\end{bmatrix}
-
\begin{bmatrix}
1\\
0\\
0
\end{bmatrix}
-
\begin{bmatrix}
0\\
1\\
0
\end{bmatrix}\\
&=
\begin{bmatrix}
0\\
0\\
1
\end{bmatrix}.
\end{aligned}
$$

所以：

$$
q_3=
\begin{bmatrix}
0\\
0\\
1
\end{bmatrix}.
$$

这个例子清楚地展示了 Gram–Schmidt 的作用：

```text
a₁ 提供第一个方向
a₂ 减去第一个方向后，留下第二个方向
a₃ 减去前两个方向后，留下第三个方向
```

---

## 11. 一般的 Gram–Schmidt 公式

设线性无关向量：

$$
a_1,a_2,\ldots,a_n.
$$

第一步：

$$
u_1=a_1,
\qquad
q_1=\frac{u_1}{\|u_1\|}.
$$

第二步：

$$
u_2=a_2-(q_1^Ta_2)q_1,
\qquad
q_2=\frac{u_2}{\|u_2\|}.
$$

第三步：

$$
u_3
=
a_3-(q_1^Ta_3)q_1-(q_2^Ta_3)q_2,
$$

$$
q_3=\frac{u_3}{\|u_3\|}.
$$

一般地，第 $j$ 步为：

$$
\boxed{
u_j
=
a_j-
\sum_{i=1}^{j-1}(q_i^Ta_j)q_i
}.
$$

如果 $u_j\neq0$，则：

$$
\boxed{
q_j
=
\frac{u_j}{\|u_j\|}
}.
$$

求和中的每一项：

$$
(q_i^Ta_j)q_i
$$

都是 $a_j$ 在已有单位方向 $q_i$ 上的投影。

所以整个求和：

$$
\sum_{i=1}^{j-1}(q_i^Ta_j)q_i
$$

就是 $a_j$ 在已有子空间：

$$
\operatorname{span}\{q_1,\ldots,q_{j-1}\}
$$

上的投影。

因此也可以写成：

$$
\boxed{
u_j
=
a_j-
\operatorname{proj}_{\operatorname{span}\{q_1,\ldots,q_{j-1}\}}a_j
}.
$$

这说明 Gram–Schmidt 就是在反复计算“投影后的残差”。

---

## 12. 为什么新向量与所有旧方向都正交

对于任意 $k<j$，计算：

$$
q_k^Tu_j
=
q_k^T
\left(
a_j-
\sum_{i=1}^{j-1}(q_i^Ta_j)q_i
\right).
$$

展开：

$$
q_k^Tu_j
=
q_k^Ta_j-
\sum_{i=1}^{j-1}(q_i^Ta_j)(q_k^Tq_i).
$$

由于已有向量标准正交：

$$
q_k^Tq_i=
\begin{cases}
1, & i=k,\\
0, & i\neq k,
\end{cases}
$$

求和中只有 $i=k$ 的一项不为零，所以：

$$
\begin{aligned}
q_k^Tu_j
&=q_k^Ta_j-(q_k^Ta_j)\\
&=0.
\end{aligned}
$$

因此：

$$
\boxed{
u_j\perp q_1,q_2,\ldots,q_{j-1}
}.
$$

单位化不会改变方向，所以 $q_j$ 也与所有旧方向正交。

---

## 13. 从矩阵形状理解 $Q$

设：

$$
A=
\begin{bmatrix}
|&|&&|\\
a_1&a_2&\cdots&a_n\\
|&|&&|
\end{bmatrix}
\in\mathbb{R}^{m\times n},
$$

并假设 $A$ 的列线性无关。

Gram–Schmidt 得到：

$$
Q=
\begin{bmatrix}
|&|&&|\\
q_1&q_2&\cdots&q_n\\
|&|&&|
\end{bmatrix}
\in\mathbb{R}^{m\times n}.
$$

这里 $Q$ 不一定是方阵。

- 每个 $q_i$ 属于 $\mathbb{R}^m$，所以 $Q$ 有 $m$ 行；
- 一共有 $n$ 个基向量，所以 $Q$ 有 $n$ 列。

只要列标准正交，就有：

$$
\boxed{Q^TQ=I_n}.
$$

但如果 $m>n$，一般：

$$
QQ^T\neq I_m.
$$

此时 $QQ^T$ 是投影到 $\operatorname{Col}(A)$ 的矩阵。

Gram–Schmidt 保证：

$$
\boxed{
\operatorname{Col}(Q)=\operatorname{Col}(A)
}.
$$

所以 $Q$ 用一组更容易计算的标准正交列，描述了与 $A$ 完全相同的列空间。

---

## 14. 如果原向量线性相关，会发生什么

设：

$$
a_1=
\begin{bmatrix}
1\\
0
\end{bmatrix},
\qquad
a_2=
\begin{bmatrix}
2\\
0
\end{bmatrix}.
$$

$a_2=2a_1$，所以第二个向量没有提供新方向。

第一步得到：

$$
q_1=
\begin{bmatrix}
1\\
0
\end{bmatrix}.
$$

$a_2$ 在 $q_1$ 上的投影是：

$$
(q_1^Ta_2)q_1
=2
\begin{bmatrix}
1\\
0
\end{bmatrix}
=
\begin{bmatrix}
2\\
0
\end{bmatrix}
=a_2.
$$

减去投影：

$$
u_2
=
a_2-(q_1^Ta_2)q_1
=0.
$$

这时无法进行单位化，因为：

$$
\frac{u_2}{\|u_2\|}
=
\frac{0}{0}
$$

没有定义。

几何意义是：

> $a_2$ 完全位于已有子空间中，减去已有方向上的投影后什么也没有剩下，因此它没有贡献新的独立方向。

所以：

$$
\boxed{
u_j=0
\Longleftrightarrow
a_j\in\operatorname{span}\{a_1,\ldots,a_{j-1}\}
}.
$$

Gram–Schmidt 不仅能构造标准正交基，还能暴露线性相关性。

如果目标只是寻找列空间的一组基，可以跳过这个零向量；如果预期所有列都线性无关，出现零向量就表示输入不满足条件。

---

## 15. 经典 Gram–Schmidt 与改进 Gram–Schmidt

在精确数学中，经典方法和改进方法得到相同的结果。

区别主要出现在浮点数计算中。

### 15.1 经典 Gram–Schmidt

经典形式一次写成：

$$
u_j
=
a_j-
\sum_{i=1}^{j-1}(q_i^Ta_j)q_i.
$$

所有投影系数都直接使用原始的 $a_j$ 计算。

### 15.2 改进 Gram–Schmidt

改进方法先令：

$$
v=a_j.
$$

然后逐个删除已有方向：

$$
v\leftarrow v-(q_1^Tv)q_1,
$$

$$
v\leftarrow v-(q_2^Tv)q_2,
$$

$$
\vdots
$$

最后：

$$
u_j=v.
$$

区别可以概括为：

```text
经典方法：计算完所有投影，再统一相减
改进方法：减掉一个方向后，立即用更新后的剩余向量继续计算
```

在有限精度计算中，改进 Gram–Schmidt 通常更不容易丢失正交性，因此比经典形式更稳定。

本阶段需要理解这种差别，但暂时不要求深入浮点误差分析。

---

## 16. 怎样检查 Gram–Schmidt 的结果

得到 $q_1,ldots,q_n$ 后，应检查三件事。

### 检查一：每个向量长度为 $1$

$$
\|q_i\|=1.
$$

### 检查二：不同向量互相正交

$$
q_i^Tq_j=0,
\qquad
i\neq j.
$$

把这些条件写成矩阵形式就是：

$$
Q^TQ=I.
$$

### 检查三：张成空间没有改变

$$
\operatorname{Col}(Q)=\operatorname{Col}(A).
$$

在手算中，可以通过下面的关系检查：

- 每个 $q_j$ 都能由 $a_1,ldots,a_j$ 线性表示；
- 每个 $a_j$ 也能由 $q_1,ldots,q_j$ 线性表示。

---

## 17. Gram–Schmidt 与投影、最小二乘的联系

第三步学习了：如果 $Q$ 的列标准正交，那么投影公式非常简单：

$$
p=QQ^Tb.
$$

第四步学习了：对于普通满列秩矩阵 $A$，投影公式是：

$$
p=A(A^TA)^{-1}A^Tb.
$$

Gram–Schmidt 把 $A$ 的普通列变成标准正交列 $Q$，同时保持：

$$
\operatorname{Col}(Q)=\operatorname{Col}(A).
$$

因此，投影到 $A$ 的列空间也可以改成：

$$
p=QQ^Tb.
$$

知识关系是：

```text
A 的普通线性无关列
          │
          │ Gram–Schmidt
          ▼
Q 的标准正交列
          │
          ├─ QᵀQ=I
          ├─ Col(Q)=Col(A)
          └─ 投影可写成 QQᵀb
```

下一步会把这个过程写成矩阵分解：

$$
A=QR.
$$

---

## 18. 一套可以直接使用的手算流程

给定线性无关向量：

$$
a_1,a_2,\ldots,a_n.
$$

按照下面的顺序进行。

### 第一步：处理第一个向量

$$
u_1=a_1,
\qquad
q_1=\frac{u_1}{\|u_1\|}.
$$

### 第二步：处理第 $j$ 个向量

先减去所有已有方向上的投影：

$$
u_j
=
a_j-
\sum_{i=1}^{j-1}(q_i^Ta_j)q_i.
$$

### 第三步：检查是否得到零向量

如果：

$$
u_j=0,
$$

说明 $a_j$ 没有提供新方向。

### 第四步：单位化

如果 $u_j\neq0$：

$$
q_j=\frac{u_j}{\|u_j\|}.
$$

### 第五步：验证结果

$$
Q^TQ=I,
$$

并检查：

$$
\operatorname{Col}(Q)=\operatorname{Col}(A).
$$

---

## 19. 常见误区

### 误区一：只减去最近的一个旧方向

处理 $a_j$ 时，必须减去它在所有已有方向上的投影：

$$
(q_1^Ta_j)q_1,
\ldots,
(q_{j-1}^Ta_j)q_{j-1}.
$$

只减去其中一个，不能保证新向量与其他旧方向正交。

### 误区二：忘记已有方向必须是单位向量

公式：

$$
(q_i^Ta_j)q_i
$$

成立的前提是：

$$
\|q_i\|=1.
$$

如果使用未单位化的正交向量 $u_i$，投影应写成：

$$
\frac{u_i^Ta_j}{u_i^Tu_i}u_i.
$$

### 误区三：正交后忘记单位化

$u_i$ 互相垂直，但不一定长度为 $1$。必须继续计算：

$$
q_i=\frac{u_i}{\|u_i\|}.
$$

### 误区四：出现零向量后仍然进行单位化

如果 $u_j=0$，就不能除以 $\|u_j\|$。这表示原向量线性相关，没有产生新方向。

### 误区五：认为 $Q$ 必须是方阵

如果 $A\in\mathbb{R}^{m\times n}$ 且 $m>n$，由其列得到的：

$$
Q\in\mathbb{R}^{m\times n}
$$

通常是长方形矩阵。

只要列标准正交，就有：

$$
Q^TQ=I_n.
$$

### 误区六：认为正交化改变了列空间

Gram–Schmidt 只是在同一个子空间中更换基：

$$
\operatorname{Col}(Q)=\operatorname{Col}(A).
$$

---

## 20. 回答阶段任务中的五个问题

### 1. Gram–Schmidt 每一步减去的是什么？

减去新向量在所有已有标准正交方向上的投影：

$$
\sum_{i=1}^{j-1}(q_i^Ta_j)q_i.
$$

也就是减去与已有子空间重复的部分。

### 2. 为什么减去投影后会得到正交方向？

因为已有方向上的全部分量已经被完整删除。对每个 $k<j$ 都有：

$$
q_k^Tu_j=0.
$$

### 3. 为什么正交化不会改变原来的列空间？

每个新向量都是原向量的线性组合；反过来，每个原向量也能由已经得到的新向量重新表示。因此两组向量互相包含在对方的张成空间中。

### 4. 为什么最后还要单位化？

正交只保证方向垂直，不保证长度为 $1$。单位化以后才能得到标准正交基，使坐标和投影可以直接通过内积计算。

### 5. 原向量线性相关时会发生什么？

某一步会得到：

$$
u_j=0.
$$

这表示 $a_j$ 完全位于已有向量张成的空间中，没有提供新的独立方向，因此不能再单位化为新的 $q_j$。

---

## 21. 本步总结

Gram–Schmidt 的核心不是背诵一串公式，而是理解一个动作：

$$
\boxed{
\text{新方向}
=
\text{原向量}
-
\text{已有子空间上的投影}
}.
$$

第一个方向：

$$
u_1=a_1,
\qquad
q_1=\frac{u_1}{\|u_1\|}.
$$

第 $j$ 个方向：

$$
\boxed{
u_j
=
a_j-
\sum_{i=1}^{j-1}(q_i^Ta_j)q_i
},
$$

$$
\boxed{
q_j
=
\frac{u_j}{\|u_j\|}
}.
$$

最终得到：

$$
\boxed{Q^TQ=I},
$$

同时保持：

$$
\boxed{
\operatorname{Col}(Q)=\operatorname{Col}(A)
}.
$$

最终的心智模型是：

```text
原向量 aⱼ
│
├─ 已有方向上的部分：删除
│
└─ 与已有方向垂直的部分：保留为 uⱼ
                              │
                              │ 单位化
                              ▼
                              qⱼ
```

下一步将把 Gram–Schmidt 中的结果组织成矩阵：

$$
A=QR.
$$

其中 $Q$ 保存标准正交方向，$R$ 保存原列向量在这些新方向下的坐标。
