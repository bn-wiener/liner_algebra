# 第三阶段：方程组、消元、秩与子空间

> 核心问题：给定目标 $\boldsymbol b$，是否存在输入 $\boldsymbol x$，使 $A\boldsymbol x=\boldsymbol b$？如果存在，全部解具有什么结构？

第二阶段把矩阵理解为线性变换。本阶段反过来研究：已知变换规则和目标输出，能否恢复输入？这一问题把线性方程组、消元、秩和子空间连接成一条完整主线：

$$
\text{线性方程组}
\rightarrow
\text{增广矩阵}
\rightarrow
\text{高斯消元}
\rightarrow
\text{主元与自由变量}
\rightarrow
\text{解的结构}
\rightarrow
\text{秩}
\rightarrow
\text{四个基本子空间}
$$

本阶段要建立一个重要习惯：每当看到

$$
A\boldsymbol x=\boldsymbol b
$$

都同时提出四个问题：

1. **方程视角**：各个未知量是否满足全部方程？
2. **列空间视角**：$\boldsymbol b$ 能否由 $A$ 的列线性组合得到？
3. **变换视角**：$\boldsymbol b$ 是否属于线性变换 $A$ 的输出范围？
4. **解空间视角**：哪些输入差异不会改变输出？

## 学习目标

完成本阶段后，应当能够：

1. 在线性方程组、矩阵方程和列向量线性组合之间自由转换。
2. 使用高斯消元将矩阵化为行阶梯形和简化行阶梯形。
3. 识别主元位置、主元变量和自由变量。
4. 判断方程组无解、存在唯一解或存在无穷多解。
5. 使用参数向量形式表示方程组的全部解。
6. 解释齐次方程的零空间及其几何意义。
7. 证明并使用一般解结构 $\boldsymbol x=\boldsymbol x_p+\boldsymbol x_n$。
8. 求列空间、行空间、零空间和左零空间的一组基。
9. 理解秩的多种等价含义，以及行秩等于列秩。
10. 使用秩—零度定理计算子空间维数。
11. 根据矩阵尺寸和秩快速判断解的可能结构。
12. 将这些概念应用于约束、网络、数据冗余和工程模型。

## 1. 线性方程组的三种语言

考虑方程组：

$$
\begin{cases}
x+2y=5\\
3x-y=4
\end{cases}
$$

### 1.1 方程语言

我们寻找一组数 $(x,y)$，使它同时满足两条约束。每个方程在二维平面中表示一条直线，方程组的解是两条直线的公共交点。

### 1.2 矩阵语言

方程组可以写成：

$$
\begin{bmatrix}
1&2\\
3&-1
\end{bmatrix}
\begin{bmatrix}x\\y\end{bmatrix}
=
\begin{bmatrix}5\\4\end{bmatrix}
$$

即：

$$
A\boldsymbol x=\boldsymbol b
$$

### 1.3 列向量语言

同一个问题还可以写成：

$$
x\begin{bmatrix}1\\3\end{bmatrix}
+y\begin{bmatrix}2\\-1\end{bmatrix}
=
\begin{bmatrix}5\\4\end{bmatrix}
$$

这时问题变成：是否能用矩阵的列向量线性组合出目标 $\boldsymbol b$？

因此：

$$
\boxed{
A\boldsymbol x=\boldsymbol b\text{ 有解}
\iff
\boldsymbol b\in C(A)
}
$$

其中 $C(A)$ 表示 $A$ 的列空间。

### 1.4 几何上的三种可能

在二维的两个方程中，可能出现：

- 两条直线相交：唯一解。
- 两条直线重合：无穷多解。
- 两条直线平行且不同：无解。

高维空间中仍然只有这三种解的数量类型：无解、唯一解或无穷多解。

## 2. 增广矩阵与初等行变换

### 2.1 增广矩阵

矩阵方程：

$$
A\boldsymbol x=\boldsymbol b
$$

可以写成增广矩阵：

$$
[A\mid\boldsymbol b]
$$

例如：

$$
\begin{cases}
x+2y-z=3\\
2x+5y+z=8\\
-x-2y+2z=-1
\end{cases}
$$

对应：

$$
\left[
\begin{array}{ccc|c}
1&2&-1&3\\
2&5&1&8\\
-1&-2&2&-1
\end{array}
\right]
$$

竖线左侧是系数矩阵，右侧是目标向量。

### 2.2 三种初等行变换

允许使用三类操作：

1. 交换两行：$R_i\leftrightarrow R_j$。
2. 某一行乘非零常数：$R_i\leftarrow cR_i$，其中 $c\neq0$。
3. 某一行加上另一行的倍数：$R_i\leftarrow R_i+cR_j$。

这些操作不会改变方程组的解集，因为它们只是在用等价方式重新组织原有约束。

例如，将第二个方程减去第一个方程的 2 倍，不会增加或删除任何真正的解。

### 2.3 行变换改变了什么

行变换保持：

- 方程组的解集。
- 系数矩阵的零空间。
- 行空间。
- 秩。

行变换一般不保持：

- 矩阵本身的列空间。
- 原矩阵各列之间的具体几何位置。

因此，寻找列空间的基时，要从**原矩阵**中选取主元列，而不能直接使用阶梯形矩阵的主元列。

### 2.4 初等矩阵

每次初等行变换也可以看成左乘一个可逆矩阵 $E$：

$$
EA
$$

$E$ 是对单位矩阵执行同样行操作得到的初等矩阵。连续消元可以写成：

$$
E_k\cdots E_2E_1A=R
$$

其中 $R$ 是消元后的矩阵。这为后续 LU 分解奠定基础。

## 3. 高斯消元

### 3.1 消元的目标

高斯消元通过行变换，将增广矩阵变为容易回代的阶梯结构。基本过程是：

1. 在当前列寻找非零主元。
2. 必要时交换行，将主元移到适当位置。
3. 用主元消去其下方元素。
4. 向右下方继续处理。
5. 得到行阶梯形后进行回代。

### 3.2 完整例子

从增广矩阵开始：

$$
\left[
\begin{array}{ccc|c}
1&2&-1&3\\
2&5&1&8\\
-1&-2&2&-1
\end{array}
\right]
$$

执行：

$$
R_2\leftarrow R_2-2R_1,
\qquad
R_3\leftarrow R_3+R_1
$$

得到：

$$
\left[
\begin{array}{ccc|c}
1&2&-1&3\\
0&1&3&2\\
0&0&1&2
\end{array}
\right]
$$

从最后一行回代：

$$
z=2
$$

第二行给出：

$$
y+3z=2
\quad\Rightarrow\quad
y=-4
$$

第一行给出：

$$
x+2y-z=3
\quad\Rightarrow\quad
x=13
$$

所以唯一解为：

$$
\boldsymbol x=
\begin{bmatrix}13\\-4\\2\end{bmatrix}
$$

### 3.3 为什么需要交换行

如果主元位置为零，不能直接用它消元。例如：

$$
\begin{bmatrix}
0&1\\
2&3
\end{bmatrix}
$$

应先交换两行。数值计算中，即使主元不为零但非常小，也常通过选取绝对值较大的主元减少舍入误差，这称为选主元策略。

## 4. 行阶梯形、简化行阶梯形与主元

### 4.1 行阶梯形

矩阵处于行阶梯形，通常要求：

1. 所有非零行位于零行上方。
2. 每个非零行的首个非零元素位于上一行首个非零元素的右侧。
3. 每个主元下方元素均为零。

例如：

$$
\begin{bmatrix}
1&2&0&3\\
0&1&4&2\\
0&0&0&1\\
0&0&0&0
\end{bmatrix}
$$

### 4.2 简化行阶梯形

简化行阶梯形还要求：

1. 每个主元等于 1。
2. 每个主元是所在列唯一的非零元素。

例如：

$$
\begin{bmatrix}
1&0&-2&0\\
0&1&3&0\\
0&0&0&1\\
0&0&0&0
\end{bmatrix}
$$

一个矩阵的行阶梯形可能不唯一，但简化行阶梯形是唯一的。

### 4.3 主元与自由变量

主元所在列对应主元变量；非主元列对应自由变量。

例如：

$$
\begin{bmatrix}
1&2&0&-1\\
0&0&1&3\\
0&0&0&0
\end{bmatrix}
$$

主元位于第 1、3 列，因此：

- $x_1,x_3$ 是主元变量。
- $x_2,x_4$ 是自由变量。

自由变量可以独立取值，主元变量由它们决定。

### 4.4 秩

主元数量称为矩阵的秩：

$$
\operatorname{rank}(A)=\text{主元数量}
$$

秩也等于：

- 线性无关列的最大数量。
- 线性无关行的最大数量。
- 列空间的维数。
- 行空间的维数。
- 线性变换输出空间的维数。

这些看似不同的定义实际上描述同一件事：矩阵中包含多少个独立方向。

## 5. 方程组解的三种情况

### 5.1 无解

如果消元后出现：

$$
\begin{bmatrix}
0&0&\cdots&0&|&c
\end{bmatrix},
\qquad c\neq0
$$

就对应矛盾方程：

$$
0=c
$$

因此方程组无解。

例如：

$$
\begin{cases}
x+y=1\\
2x+2y=3
\end{cases}
$$

第二行减去第一行的 2 倍得到 $0=1$，所以无解。

### 5.2 唯一解

方程组一致，并且每个未知量都是主元变量时，存在唯一解。

对于 $n$ 个未知量，这要求：

$$
\operatorname{rank}(A)=n
$$

同时还必须满足方程一致。

### 5.3 无穷多解

方程组一致，但至少存在一个自由变量时，存在无穷多解。

每个自由变量对应一个独立的解方向。因此，自由变量数量就是齐次解空间的维数。

### 5.4 秩判据

方程组有解当且仅当：

$$
\operatorname{rank}(A)
=
\operatorname{rank}([A\mid\boldsymbol b])
$$

若增广列引入新主元，则说明 $\boldsymbol b$ 不在原列空间中，方程无解。

有解时：

- 若 $\operatorname{rank}(A)=n$，解唯一。
- 若 $\operatorname{rank}(A)<n$，存在自由变量，因此有无穷多解。

## 6. 参数向量形式

考虑方程：

$$
\begin{bmatrix}
1&2&0&1\\
0&0&1&1
\end{bmatrix}
\begin{bmatrix}
x_1\\x_2\\x_3\\x_4
\end{bmatrix}
=
\begin{bmatrix}3\\2\end{bmatrix}
$$

对应：

$$
\begin{cases}
x_1+2x_2+x_4=3\\
x_3+x_4=2
\end{cases}
$$

令自由变量：

$$
x_2=s,
\qquad
x_4=t
$$

则：

$$
x_1=3-2s-t,
\qquad
x_3=2-t
$$

所以：

$$
\boldsymbol x=
\begin{bmatrix}
3-2s-t\\
s\\
2-t\\
t
\end{bmatrix}
$$

拆成向量形式：

$$
\boldsymbol x=
\begin{bmatrix}3\\0\\2\\0\end{bmatrix}
+s\begin{bmatrix}-2\\1\\0\\0\end{bmatrix}
+t\begin{bmatrix}-1\\0\\-1\\1\end{bmatrix}
$$

其中第一个向量是一个特解，后两个方向来自对应齐次方程的零空间。

## 7. 齐次方程与零空间

### 7.1 定义

齐次方程为：

$$
A\boldsymbol x=\boldsymbol0
$$

所有齐次解构成 $A$ 的零空间：

$$
N(A)={\boldsymbol x\in\mathbb R^n:A\boldsymbol x=\boldsymbol0}
$$

零空间位于输入空间 $\mathbb R^n$ 中。

### 7.2 零解与非零解

齐次方程总有零解：

$$
\boldsymbol x=\boldsymbol0
$$

若存在自由变量，则还存在非零解。对于有 $n$ 列的矩阵：

$$
\operatorname{rank}(A)<n
\quad\Rightarrow\quad
N(A)\text{ 中存在非零向量}
$$

### 7.3 几何意义

零空间描述所有被矩阵完全消除的输入方向：

$$
\boldsymbol x\in N(A)
\quad\Rightarrow\quad
A\boldsymbol x=\boldsymbol0
$$

如果 $\boldsymbol x_1-\boldsymbol x_2\in N(A)$，那么：

$$
A\boldsymbol x_1=A\boldsymbol x_2
$$

因此，零空间还描述哪些不同输入无法被输出区分。

### 7.4 求零空间的一组基

对 $A$ 进行行化简，找出自由变量，再将解写成参数向量形式。每个自由参数对应一个基础解向量，这些向量构成零空间的一组基。

例如：

$$
\operatorname{rref}(A)=
\begin{bmatrix}
1&2&0&1\\
0&0&1&1
\end{bmatrix}
$$

则：

$$
N(A)=\operatorname{span}\left(
\begin{bmatrix}-2\\1\\0\\0\end{bmatrix},
\begin{bmatrix}-1\\0\\-1\\1\end{bmatrix}
\right)
$$

## 8. 非齐次方程的完整解结构

### 8.1 特解加齐次解

若 $\boldsymbol x_p$ 是：

$$
A\boldsymbol x=\boldsymbol b
$$

的一个特解，而 $\boldsymbol x_n\in N(A)$，那么：

$$
A(\boldsymbol x_p+\boldsymbol x_n)
=A\boldsymbol x_p+A\boldsymbol x_n
=\boldsymbol b
$$

所以 $\boldsymbol x_p+\boldsymbol x_n$ 仍是非齐次方程的解。

反过来，若 $\boldsymbol x$ 和 $\boldsymbol x_p$ 都是特解，则：

$$
A(\boldsymbol x-\boldsymbol x_p)
=\boldsymbol b-\boldsymbol b
=\boldsymbol0
$$

因此 $\boldsymbol x-\boldsymbol x_p\in N(A)$。

所以全部解恰好为：

$$
\boxed{
\boldsymbol x=\boldsymbol x_p+\boldsymbol x_n,
\qquad
\boldsymbol x_n\in N(A)
}
$$

### 8.2 几何意义

- 齐次解集 $N(A)$ 是经过原点的子空间。
- 非齐次解集若非空，则是零空间的平移。
- 若零空间是一条直线，非齐次解集也是一条平行直线，但通常不经过原点。
- 若零空间只有零向量，非齐次方程至多有一个解。

## 9. 列空间：矩阵能够产生的所有输出

### 9.1 定义

设：

$$
A=
\begin{bmatrix}
\boldsymbol a_1&\boldsymbol a_2&\cdots&\boldsymbol a_n
\end{bmatrix}
$$

列空间为：

$$
C(A)=\operatorname{span}(\boldsymbol a_1,\dots,\boldsymbol a_n)
$$

它位于输出空间 $\mathbb R^m$ 中，表示矩阵能够产生的全部输出：

$$
C(A)={A\boldsymbol x:\boldsymbol x\in\mathbb R^n}
$$

因此：

$$
A\boldsymbol x=\boldsymbol b\text{ 有解}
\iff
\boldsymbol b\in C(A)
$$

### 9.2 求列空间的一组基

步骤如下：

1. 对矩阵 $A$ 进行行化简。
2. 找出阶梯形中的主元列位置。
3. 回到**原矩阵 $A$**，选择对应位置的列。

例如：

$$
A=
\begin{bmatrix}
1&2&-1\\
2&4&1\\
-1&-2&-2
\end{bmatrix}
$$

若主元位于第 1、3 列，则列空间的一组基为原矩阵的第 1、3 列：

$$
C(A)=\operatorname{span}\left(
\begin{bmatrix}1\\2\\-1\end{bmatrix},
\begin{bmatrix}-1\\1\\-2\end{bmatrix}
\right)
$$

### 9.3 为什么不能使用化简后的列

行变换会对列向量在输出空间中的具体位置进行变换，因此一般会改变列空间。阶梯形只能告诉我们哪些**列位置**是独立的；真正的列空间基必须从原矩阵中取出。

## 10. 行空间：矩阵中的独立约束

### 10.1 定义

矩阵各行的张成空间称为行空间：

$$
C(A^\mathsf T)
$$

对于 $m\times n$ 矩阵，行向量有 $n$ 个分量，因此行空间位于 $\mathbb R^n$ 中。

行空间描述方程组中所有可由原约束线性组合得到的约束方向。

### 10.2 求行空间的一组基

初等行变换不会改变行空间，所以阶梯形矩阵的所有非零行可以作为原矩阵行空间的一组基。

若：

$$
\operatorname{rref}(A)=
\begin{bmatrix}
1&2&0&1\\
0&0&1&1\\
0&0&0&0
\end{bmatrix}
$$

则行空间的一组基可以取为：

$$
\begin{bmatrix}1&2&0&1\end{bmatrix},
\qquad
\begin{bmatrix}0&0&1&1\end{bmatrix}
$$

### 10.3 行秩等于列秩

主元数量既等于独立行的数量，也等于独立列的数量，因此：

$$
\dim C(A)=\dim C(A^\mathsf T)=\operatorname{rank}(A)
$$

这就是行秩等于列秩。它说明矩阵独立输入—输出方向的数量与独立约束的数量一致。

## 11. 左零空间：输出必须满足的兼容条件

### 11.1 定义

左零空间是 $A^\mathsf T$ 的零空间：

$$
N(A^\mathsf T)
=
\{\boldsymbol y\in\mathbb R^m:A^\mathsf T\boldsymbol y=\boldsymbol0\}
$$

等价地：

$$
\boldsymbol y^\mathsf TA=\boldsymbol0^\mathsf T
$$

左零空间中的向量描述矩阵各行之间的线性依赖关系。

### 11.2 与方程一致性的关系

若 $A\boldsymbol x=\boldsymbol b$ 有解，并且 $\boldsymbol y\in N(A^\mathsf T)$，则：

$$
\boldsymbol y^\mathsf T\boldsymbol b
=\boldsymbol y^\mathsf TA\boldsymbol x
=0
$$

因此，每个左零向量都会对可达到的输出 $\boldsymbol b$ 施加一个兼容条件：

$$
\boldsymbol y^\mathsf T\boldsymbol b=0
$$

如果某个 $\boldsymbol b$ 不满足该条件，它就不可能位于列空间中。

### 11.3 例：识别矛盾约束

若矩阵三行满足：

$$
\boldsymbol r_3=\boldsymbol r_1-\boldsymbol r_2
$$

那么：

$$
\boldsymbol r_1-\boldsymbol r_2-\boldsymbol r_3=\boldsymbol0
$$

对应左零向量：

$$
\boldsymbol y=
\begin{bmatrix}1\\-1\\-1\end{bmatrix}
$$

方程 $A\boldsymbol x=\boldsymbol b$ 有解时，右端也必须满足同样关系：

$$
b_1-b_2-b_3=0
$$

## 12. 四个基本子空间

对于 $m\times n$ 矩阵 $A$，四个基本子空间如下：

| 子空间 | 所在空间 | 含义 | 维数 |
|---|---|---|---|
| 列空间 $C(A)$ | $\mathbb R^m$ | 所有可能输出 | $r$ |
| 零空间 $N(A)$ | $\mathbb R^n$ | 被变换消除的输入方向 | $n-r$ |
| 行空间 $C(A^\mathsf T)$ | $\mathbb R^n$ | 独立输入约束方向 | $r$ |
| 左零空间 $N(A^\mathsf T)$ | $\mathbb R^m$ | 输出必须满足的兼容关系 | $m-r$ |

其中：

$$
r=\operatorname{rank}(A)
$$

可以按输入空间和输出空间组织：

$$
\mathbb R^n:
\qquad
C(A^\mathsf T)
\quad\text{与}\quad
N(A)
$$

$$
\mathbb R^m:
\qquad
C(A)
\quad\text{与}\quad
N(A^\mathsf T)
$$

下一阶段会证明它们分别互为正交补：

$$
N(A)=C(A^\mathsf T)^\perp
$$

$$
N(A^\mathsf T)=C(A)^\perp
$$

直观上，零空间中的方向与每个行向量都正交；左零空间中的方向与每个列向量都正交。

## 13. 秩—零度定理

### 13.1 定理

若 $A$ 是 $m\times n$ 矩阵，则：

$$
\boxed{
\operatorname{rank}(A)
+
\operatorname{nullity}(A)
=n
}
$$

其中：

$$
\operatorname{nullity}(A)=\dim N(A)
$$

对转置矩阵还有：

$$
\operatorname{rank}(A)
+
\dim N(A^\mathsf T)
=m
$$

### 13.2 为什么成立

$A$ 有 $n$ 列，对应 $n$ 个未知量。消元后，每个变量恰好属于两类之一：

- 主元变量，共 $r$ 个。
- 自由变量，共 $n-r$ 个。

每个自由变量产生零空间中的一个独立基础方向，所以：

$$
\dim N(A)=n-r
$$

因此：

$$
r+(n-r)=n
$$

### 13.3 例

若 $A$ 是 $3\times5$ 矩阵，且秩为 3，则：

$$
\dim N(A)=5-3=2
$$

因为 $A^\mathsf T$ 的秩也为 3，所以：

$$
\dim N(A^\mathsf T)=3-3=0
$$

这意味着：

- 输入中有两个不会影响输出的自由方向。
- 输出空间没有额外兼容约束。
- 列空间等于整个 $\mathbb R^3$，所以任意 $\boldsymbol b\in\mathbb R^3$ 都可达到。

## 14. 根据矩阵形状预测解的结构

设 $A$ 是 $m\times n$ 矩阵，秩为 $r$。

### 14.1 高矩阵：$m>n$

方程多于未知量，系统可能约束过多。若列满秩 $r=n$，一致时解唯一，但并非每个 $\boldsymbol b\in\mathbb R^m$ 都可达到。

这类系统常出现在多次观测同一组未知参数的问题中。第四阶段将处理通常无精确解的超定系统。

### 14.2 宽矩阵：$m<n$

未知量多于方程。因为：

$$
r\le m<n
$$

所以必有自由变量。只要系统一致，就有无穷多解。

### 14.3 方阵：$m=n$

若 $r=n$，则矩阵可逆，对每个 $\boldsymbol b$ 都有唯一解。若 $r<n$，则矩阵不可逆：某些 $\boldsymbol b$ 无解，其余可达到的 $\boldsymbol b$ 有无穷多解。

### 14.4 满列秩与满行秩

- 满列秩：$r=n$，零空间只有零向量，输入不会发生混淆。
- 满行秩：$r=m$，列空间等于 $\mathbb R^m$，每个目标输出都可达到。
- 方阵同时满行秩和满列秩时可逆。

## 15. 贯穿四个基本子空间的完整例子

考虑：

$$
A=
\begin{bmatrix}
1&2&-1&0\\
2&4&1&3\\
-1&-2&-2&-3
\end{bmatrix}
$$

第三行满足：

$$
\boldsymbol r_3=\boldsymbol r_1-\boldsymbol r_2
$$

所以三行中最多只有两个独立方向。

### 15.1 行化简与秩

执行：

$$
R_2\leftarrow R_2-2R_1,
\qquad
R_3\leftarrow R_3+R_1
$$

再继续化简，可得：

$$
\operatorname{rref}(A)=
\begin{bmatrix}
1&2&0&1\\
0&0&1&1\\
0&0&0&0
\end{bmatrix}
$$

主元位于第 1、3 列，所以：

$$
\operatorname{rank}(A)=2
$$

### 15.2 列空间

从原矩阵选取第 1、3 列：

$$
C(A)=\operatorname{span}\left(
\begin{bmatrix}1\\2\\-1\end{bmatrix},
\begin{bmatrix}-1\\1\\-2\end{bmatrix}
\right)
$$

列空间位于 $\mathbb R^3$ 中，维数为 2，因此是一个经过原点的平面。

### 15.3 行空间

取简化行阶梯形中的非零行：

$$
C(A^\mathsf T)=\operatorname{span}\left(
\begin{bmatrix}1&2&0&1\end{bmatrix},
\begin{bmatrix}0&0&1&1\end{bmatrix}
\right)
$$

行空间位于 $\mathbb R^4$ 中，维数为 2。

### 15.4 零空间

齐次方程给出：

$$
\begin{cases}
x_1+2x_2+x_4=0\\
x_3+x_4=0
\end{cases}
$$

令：

$$
x_2=s,
\qquad
x_4=t
$$

得到：

$$
\boldsymbol x
=s\begin{bmatrix}-2\\1\\0\\0\end{bmatrix}
+t\begin{bmatrix}-1\\0\\-1\\1\end{bmatrix}
$$

所以：

$$
N(A)=\operatorname{span}\left(
\begin{bmatrix}-2\\1\\0\\0\end{bmatrix},
\begin{bmatrix}-1\\0\\-1\\1\end{bmatrix}
\right)
$$

其维数为：

$$
4-2=2
$$

### 15.5 左零空间

由行关系：

$$
\boldsymbol r_1-\boldsymbol r_2-\boldsymbol r_3=\boldsymbol0
$$

得到：

$$
N(A^\mathsf T)=\operatorname{span}\left(
\begin{bmatrix}1\\-1\\-1\end{bmatrix}
\right)
$$

其维数为：

$$
3-2=1
$$

### 15.6 判断目标是否可达到

取：

$$
\boldsymbol b=
\begin{bmatrix}1\\2\\-1\end{bmatrix}
$$

它正是 $A$ 的第一列，所以显然属于列空间。一个特解为：

$$
\boldsymbol x_p=
\begin{bmatrix}1\\0\\0\\0\end{bmatrix}
$$

全部解为：

$$
\boldsymbol x=
\begin{bmatrix}1\\0\\0\\0\end{bmatrix}
+s\begin{bmatrix}-2\\1\\0\\0\end{bmatrix}
+t\begin{bmatrix}-1\\0\\-1\\1\end{bmatrix}
$$

再取：

$$
\boldsymbol b_{\text{bad}}=
\begin{bmatrix}1\\2\\0\end{bmatrix}
$$

左零向量给出的兼容条件为：

$$
b_1-b_2-b_3=0
$$

但：

$$
1-2-0=-1\neq0
$$

因此 $\boldsymbol b_{\text{bad}}$ 不在列空间中，方程无解。

### 15.7 维数核对

输入空间：

$$
\dim C(A^\mathsf T)+\dim N(A)=2+2=4
$$

输出空间：

$$
\dim C(A)+\dim N(A^\mathsf T)=2+1=3
$$

这个例子将消元、主元、秩、四个基本子空间、方程一致性和全部解结构连接成了一个整体。

## 16. 常见误区

### 误区 1：行变换后的矩阵与原矩阵具有相同列空间

行变换通常会改变列空间。寻找列空间基时，应根据主元位置回到原矩阵选择对应列。

### 误区 2：零空间位于输出空间

若 $A$ 是 $m\times n$ 矩阵，零空间中的向量必须能作为 $A$ 的输入，所以：

$$
N(A)\subseteq\mathbb R^n
$$

### 误区 3：有自由变量就一定有无穷多解

还必须先保证系统一致。若出现矛盾行，方程仍然无解。

### 误区 4：齐次方程可能无解

齐次方程始终至少有零解。真正的问题是是否存在非零解。

### 误区 5：秩等于矩阵中非零元素的数量

秩是独立方向或主元的数量，与非零元素总数没有直接等价关系。

### 误区 6：零行没有意义

零行暴露了约束之间的依赖关系，也决定左零空间的维数。增广矩阵中零行对应的右端若非零，还会直接产生矛盾。

### 误区 7：特解就是全部解

只要零空间中存在非零方向，就可以在特解上加任意齐次解，得到无穷多个解。

### 误区 8：主元列可以从简化矩阵中直接作为列空间基

主元**位置**从简化矩阵确定，但基向量必须取自原矩阵。

### 误区 9：行空间位于 $\mathbb R^m$

$m\times n$ 矩阵的每一行有 $n$ 个分量，因此行空间位于 $\mathbb R^n$。

### 误区 10：方阵不可逆时，所有右端都无解

不可逆方阵的列空间虽然不是整个输出空间，但其中仍包含无穷多个向量。位于列空间中的右端有无穷多解，列空间外的右端无解。

## 17. 应用连接

- **工程约束**：结构力学、电路和控制系统中的未知量通常由线性方程约束。
- **电路分析**：基尔霍夫定律产生线性方程，约束之间的依赖会反映在秩中。
- **网络流**：节点流量守恒可以写成矩阵方程，循环流构成零空间方向。
- **化学方程式配平**：原子守恒条件构成齐次方程，非零零空间向量给出化学计量系数。
- **数据冗余**：相关特征使矩阵秩下降，表示某些数据维度没有提供独立信息。
- **可观测性**：零空间描述无法从输出中识别的输入变化。
- **一致性检查**：左零空间给出测量值必须满足的约束，可用于检测矛盾数据。
- **计算机图形学**：不可逆投影会丢失深度信息，零空间准确描述被丢失的方向。
- **后续最小二乘**：当 $\boldsymbol b$ 不在列空间中时，只能寻找列空间中距离它最近的向量。
- **矩阵分解**：消元过程可以组织成 LU 分解，为大规模方程求解提供高效方法。

## 18. 阶段练习

### 基础题

1. 将方程组写成 $A\boldsymbol x=\boldsymbol b$ 和列向量线性组合形式：

   $$
   \begin{cases}
   x-2y=3\\
   3x+y=7
   \end{cases}
   $$

2. 对下面的增广矩阵进行消元并求解：

   $$
   \left[
   \begin{array}{cc|c}
   1&2&5\\
   3&4&11
   \end{array}
   \right]
   $$

3. 判断下面矩阵是否处于行阶梯形或简化行阶梯形：

   $$
   R=
   \begin{bmatrix}
   1&0&2\\
   0&1&-1\\
   0&0&0
   \end{bmatrix}
   $$

4. 对于：

   $$
   R=
   \begin{bmatrix}
   1&2&0&-1\\
   0&0&1&3\\
   0&0&0&0
   \end{bmatrix}
   $$

   指出主元变量、自由变量、秩和零度。

5. 判断方程组的解类型：

   $$
   \left[
   \begin{array}{ccc|c}
   1&0&2&3\\
   0&1&-1&2\\
   0&0&0&1
   \end{array}
   \right]
   $$

6. 若 $A$ 是 $4\times7$ 矩阵，且 $\operatorname{rank}(A)=3$，求 $\dim N(A)$ 和 $\dim N(A^\mathsf T)$。

### 理解题

7. 为什么初等行变换保持方程组解集，却一般不保持列空间？

8. 为什么齐次方程存在自由变量时一定有非零解？

9. 证明：若 $\boldsymbol x_1$ 和 $\boldsymbol x_2$ 都满足 $A\boldsymbol x=\boldsymbol b$，则 $\boldsymbol x_1-\boldsymbol x_2\in N(A)$。

10. 为什么一个宽矩阵 $A\in\mathbb R^{m\times n}$，当 $n>m$ 时，齐次方程必有非零解？

11. 解释为什么左零空间中的向量会给出右端 $\boldsymbol b$ 的兼容条件。

12. 对一个 $5\times3$ 矩阵，秩最大是多少？如果达到最大秩，四个基本子空间的维数分别是多少？

### 综合题

13. 求下面方程组的全部解，并写成特解加齐次解：

   $$
   \begin{cases}
   x_1+2x_2-x_3=1\\
   2x_1+4x_2+x_3=5
   \end{cases}
   $$

14. 设：

   $$
   A=
   \begin{bmatrix}
   1&2&3\\
   2&4&6\\
   1&1&1
   \end{bmatrix}
   $$

   1. 求 $\operatorname{rref}(A)$ 和 $\operatorname{rank}(A)$。
   2. 求列空间的一组基。
   3. 求行空间的一组基。
   4. 求零空间的一组基。
   5. 求左零空间的一组基。

15. 设：

   $$
   A=
   \begin{bmatrix}
   1&0&1\\
   0&1&1\\
   1&1&2
   \end{bmatrix}
   $$

   对下列两个右端分别判断 $A\boldsymbol x=\boldsymbol b$ 是否有解；若有解，求全部解：

   $$
   \boldsymbol b_1=
   \begin{bmatrix}2\\3\\5\end{bmatrix},
   \qquad
   \boldsymbol b_2=
   \begin{bmatrix}2\\3\\6\end{bmatrix}
   $$

16. 已知一个 $6\times9$ 矩阵的秩为 4：

   1. 求四个基本子空间的维数。
   2. 判断齐次方程是否有非零解。
   3. 判断是否每个 $\boldsymbol b\in\mathbb R^6$ 都能由 $A\boldsymbol x$ 产生。
   4. 若某个非齐次方程有解，它是唯一解还是无穷多解？

## 19. 参考答案

1. 矩阵形式为：

   $$
   \begin{bmatrix}
   1&-2\\
   3&1
   \end{bmatrix}
   \begin{bmatrix}x\\y\end{bmatrix}
   =
   \begin{bmatrix}3\\7\end{bmatrix}
   $$

   列组合形式为：

   $$
   x\begin{bmatrix}1\\3\end{bmatrix}
   +y\begin{bmatrix}-2\\1\end{bmatrix}
   =
   \begin{bmatrix}3\\7\end{bmatrix}
   $$

2. 执行 $R_2\leftarrow R_2-3R_1$：

   $$
   \left[
   \begin{array}{cc|c}
   1&2&5\\
   0&-2&-4
   \end{array}
   \right]
   $$

   得 $y=2,x=1$。

3. $R$ 同时是行阶梯形和简化行阶梯形。两个主元均为 1，且是各自主元列中的唯一非零元素。

4. 主元变量为 $x_1,x_3$，自由变量为 $x_2,x_4$，秩为 2，零度为 $4-2=2$。

5. 最后一行表示 $0=1$，所以无解。

6. 根据秩—零度定理：

   $$
   \dim N(A)=7-3=4
   $$

   $$
   \dim N(A^\mathsf T)=4-3=1
   $$

7. 行变换是对方程约束进行可逆线性组合，所以不改变同时满足这些约束的输入。它相当于左乘一个可逆矩阵，会整体变换每个列向量在输出空间中的位置，因此列空间一般改变。

8. 自由变量可以任取非零值，再由主元方程决定主元变量，从而构造出至少一个非零齐次解。

9. 因为：

   $$
   A(\boldsymbol x_1-\boldsymbol x_2)
   =A\boldsymbol x_1-A\boldsymbol x_2
   =\boldsymbol b-\boldsymbol b
   =\boldsymbol0
   $$

   所以 $\boldsymbol x_1-\boldsymbol x_2\in N(A)$。

10. 因为：

   $$
   \operatorname{rank}(A)\le m<n
   $$

   所以：

   $$
   \dim N(A)=n-\operatorname{rank}(A)>0
   $$

11. 若 $\boldsymbol y\in N(A^\mathsf T)$ 且 $A\boldsymbol x=\boldsymbol b$，则：

   $$
   \boldsymbol y^\mathsf T\boldsymbol b
   =\boldsymbol y^\mathsf TA\boldsymbol x
   =0
   $$

   所以任何可达到的 $\boldsymbol b$ 都必须满足该条件。

12. 最大秩为 3。达到最大秩时：

   $$
   \dim C(A)=3,
   \qquad
   \dim N(A)=0
   $$

   $$
   \dim C(A^\mathsf T)=3,
   \qquad
   \dim N(A^\mathsf T)=5-3=2
   $$

13. 增广矩阵为：

   $$
   \left[
   \begin{array}{ccc|c}
   1&2&-1&1\\
   2&4&1&5
   \end{array}
   \right]
   $$

   执行 $R_2\leftarrow R_2-2R_1$，得到：

   $$
   \left[
   \begin{array}{ccc|c}
   1&2&-1&1\\
   0&0&3&3
   \end{array}
   \right]
   $$

   所以 $x_3=1$。令 $x_2=s$，则：

   $$
   x_1=2-2s
   $$

   全部解为：

   $$
   \boldsymbol x=
   \begin{bmatrix}2\\0\\1\end{bmatrix}
   +s\begin{bmatrix}-2\\1\\0\end{bmatrix}
   $$

14. 行化简：

   $$
   \operatorname{rref}(A)=
   \begin{bmatrix}
   1&0&-1\\
   0&1&2\\
   0&0&0
   \end{bmatrix}
   $$

   所以秩为 2，主元列为第 1、2 列。

   列空间的一组基为：

   $$
   \begin{bmatrix}1\\2\\1\end{bmatrix},
   \qquad
   \begin{bmatrix}2\\4\\1\end{bmatrix}
   $$

   行空间的一组基为：

   $$
   \begin{bmatrix}1&0&-1\end{bmatrix},
   \qquad
   \begin{bmatrix}0&1&2\end{bmatrix}
   $$

   齐次方程为：

   $$
   x_1-x_3=0,
   \qquad
   x_2+2x_3=0
   $$

   令 $x_3=t$，得到零空间的一组基：

   $$
   N(A)=\operatorname{span}\left(
   \begin{bmatrix}1\\-2\\1\end{bmatrix}
   \right)
   $$

   因为第二行是第一行的 2 倍，并且：

   $$
   2\boldsymbol r_1-\boldsymbol r_2=\boldsymbol0
   $$

   可取一个左零向量 $(2,-1,0)^\mathsf T$。解 $A^\mathsf T\boldsymbol y=0$ 得：

   $$
   N(A^\mathsf T)=\operatorname{span}\left(
   \begin{bmatrix}2\\-1\\0\end{bmatrix}
   \right)
   $$

   左零空间维数为 $3-2=1$，所以一个基础向量已经完整。

15. 第三行等于前两行之和，所以可达到的右端必须满足：

   $$
   b_3=b_1+b_2
   $$

   对 $\boldsymbol b_1=(2,3,5)^\mathsf T$，条件成立。方程为：

   $$
   x_1+x_3=2,
   \qquad
   x_2+x_3=3
   $$

   令 $x_3=t$，全部解为：

   $$
   \boldsymbol x=
   \begin{bmatrix}2\\3\\0\end{bmatrix}
   +t\begin{bmatrix}-1\\-1\\1\end{bmatrix}
   $$

   对 $\boldsymbol b_2=(2,3,6)^\mathsf T$，因为 $6\neq2+3$，所以无解。

16. 已知 $m=6,n=9,r=4$：

   $$
   \dim C(A)=4,
   \qquad
   \dim C(A^\mathsf T)=4
   $$

   $$
   \dim N(A)=9-4=5,
   \qquad
   \dim N(A^\mathsf T)=6-4=2
   $$

   零空间维数为 5，所以齐次方程存在非零解。列空间只有 4 维，不等于 $\mathbb R^6$，所以并非每个右端都可达到。只要非齐次方程有解，因为存在 5 个自由方向，所以一定有无穷多解。

## 20. 阶段检验

在不查阅资料的情况下，回答以下问题：

1. $A\boldsymbol x=\boldsymbol b$ 如何同时表示方程组、列组合和逆向变换问题？
2. 三种初等行变换为什么不改变方程组解集？
3. 行阶梯形与简化行阶梯形有什么区别？
4. 主元变量和自由变量分别承担什么角色？
5. 如何从增广矩阵判断无解、唯一解或无穷多解？
6. 为什么非齐次方程的全部解等于一个特解加零空间？
7. 列空间和零空间分别描述矩阵的什么能力？它们位于哪个空间？
8. 为什么列空间基要从原矩阵取，而行空间基可以从阶梯形取？
9. 左零空间如何描述行依赖和右端兼容条件？
10. 秩有哪些等价含义？
11. 四个基本子空间的维数分别是多少？
12. 如何从主元变量和自由变量解释秩—零度定理？
13. 满列秩、满行秩和可逆方阵分别意味着什么？
14. 对于一个给定矩阵，能否独立求出四个基本子空间的一组基？

能够清楚回答这些问题，独立完成阶段练习，并把消元结果解释为子空间结构，就可以进入第四阶段：正交、投影、最小二乘与 QR 分解。
