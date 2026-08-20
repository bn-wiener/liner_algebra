# E7 从正交基到 QR 再到 SVD

## 7. 一个完整的二维 QR 分解例子

设：

$$
A=
\begin{bmatrix}
1&1\\
1&0
\end{bmatrix}.
$$

它的列为：

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

第五步已经对这两个向量做过 Gram–Schmidt。

得到：

$$
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
\end{bmatrix}.
$$

所以：

$$
Q=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1&1\\
1&-1
\end{bmatrix}.
$$

### 第一步：计算 $R=Q^TA$

$$
r_{11}=q_1^Ta_1=\sqrt{2},
$$

$$
r_{12}=q_1^Ta_2=\frac{1}{\sqrt{2}},
$$

$$
r_{21}=q_2^Ta_1=0,
$$

$$
r_{22}=q_2^Ta_2=\frac{1}{\sqrt{2}}.
$$

因此：

$$
R=
\begin{bmatrix}
\sqrt{2}&1/\sqrt{2}\\
0&1/\sqrt{2}
\end{bmatrix}.
$$

### 第二步：理解 $R$ 的两列

$R$ 的第一列表示：

$$
a_1=\sqrt{2}\,q_1+0q_2.
$$

$R$ 的第二列表示：

$$
a_2
=
\frac{1}{\sqrt{2}}q_1
+
\frac{1}{\sqrt{2}}q_2.
$$

### 第三步：验证 $A=QR$

$$
\begin{aligned}
QR
&=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1&1\\
1&-1
\end{bmatrix}
\begin{bmatrix}
\sqrt{2}&1/\sqrt{2}\\
0&1/\sqrt{2}
\end{bmatrix}\\
&=
\begin{bmatrix}
1&1\\
1&0
\end{bmatrix}\\
&=A.
\end{aligned}
$$

同时：

$$
Q^TQ=I_2.
$$

所以这个分解满足 QR 分解的全部要求。

### 第四步：跟踪一个具体输入经过 $R$ 和 $Q$ 的过程

取输入：

$$
x=
\begin{bmatrix}
2\\
-1
\end{bmatrix}.
$$

直接使用 $A$ 计算输出：

$$
Ax
=
\begin{bmatrix}
1&1\\
1&0
\end{bmatrix}
\begin{bmatrix}
2\\
-1
\end{bmatrix}
=
\begin{bmatrix}
1\\
2
\end{bmatrix}.
$$

现在改走 QR 分解的两步过程。

先使用 $R$：

$$
\begin{aligned}
c=Rx
&=
\begin{bmatrix}
\sqrt{2}&1/\sqrt{2}\\
0&1/\sqrt{2}
\end{bmatrix}
\begin{bmatrix}
2\\
-1
\end{bmatrix}\\
&=
\begin{bmatrix}
3/\sqrt{2}\\
-1/\sqrt{2}
\end{bmatrix}.
\end{aligned}
$$

这里的 $c$ 表示最终输出在 $q_1,q_2$ 方向上的坐标。因此：

$$
Ax
=
\frac{3}{\sqrt{2}}q_1
-
\frac{1}{\sqrt{2}}q_2.
$$

再使用 $Q$ 重建实际输出：

$$
\begin{aligned}
Qc
&=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1&1\\
1&-1
\end{bmatrix}
\begin{bmatrix}
3/\sqrt{2}\\
-1/\sqrt{2}
\end{bmatrix}\\
&=
\begin{bmatrix}
1\\
2
\end{bmatrix}\\
&=Ax.
\end{aligned}
$$

这个例子具体说明了：

```text
x 中的 2 和 -1
是使用原列 a₁、a₂ 的组合系数

Rx 中的 3/√2 和 -1/√2
是同一个输出在标准正交方向 q₁、q₂ 下的坐标

Q 再根据这些新坐标重建实际输出 Ax
```

因此 $R$ 不是“额外做了一次无关的变换”，而是在原列组合系数与标准正交基坐标之间进行换算。

---


## 十三、一个具体数值例子

先选择一个容易看清本质的矩阵：

$$
A=
\begin{bmatrix}
3&0\\
0&2
\end{bmatrix}
$$

这里已经是 SVD 的最简单形式：

$$
U=I,\qquad V=I
$$

$$
\Sigma=
\begin{bmatrix}
3&0\\
0&2
\end{bmatrix}
$$

所以：

$$
A=U\Sigma V^T
$$

即：

$$
A=
I
\begin{bmatrix}
3&0\\
0&2
\end{bmatrix}
I
$$

对于：

$$
x=
\begin{bmatrix}
x_1\\
x_2
\end{bmatrix}
$$

有：

$$
Ax=
\begin{bmatrix}
3x_1\\
2x_2
\end{bmatrix}
$$

单位圆：

$$
x_1^2+x_2^2=1
$$

经过 $A$ 后：

$$
y_1=3x_1,\qquad y_2=2x_2
$$

因此：

$$
\frac{y_1^2}{9}+\frac{y_2^2}{4}=1
$$

变成一个椭圆。

两个半轴：

$$
3,\quad 2
$$

正好就是：

$$
\boxed{
\sigma_1=3,\qquad \sigma_2=2
}
$$

---

## 对照理解

QR 重点是在同一个列空间中构造标准正交基并记录坐标；SVD 则同时寻找最自然的输入正交方向和输出正交方向，并用奇异值把两边配对。
