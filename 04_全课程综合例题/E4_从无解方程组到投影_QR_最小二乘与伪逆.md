# 综合例题 E4：从无解方程组到投影、QR、最小二乘与伪逆

## 1. 题目

给定：

$$
A=
\begin{bmatrix}
1&0\\
1&1\\
1&2
\end{bmatrix},
\qquad
\boldsymbol{b}
=
\begin{bmatrix}
1\\2\\2
\end{bmatrix}.
$$

这里 $A$ 有 3 行、2 列，因此：

$$
A:\mathbb{R}^2\rightarrow\mathbb{R}^3.
$$

要求：

1. 判断 $A\boldsymbol{x}=\boldsymbol{b}$ 是否有精确解；
2. 求最小二乘解 $\hat{\boldsymbol{x}}$；
3. 求投影 $\boldsymbol{p}=A\hat{\boldsymbol{x}}$；
4. 求残差并验证其正交性；
5. 用 QR 分解重新求同一最小二乘解；
6. 写出伪逆并验证 $A^+\boldsymbol{b}=\hat{\boldsymbol{x}}$；
7. 解释正规方程、QR、伪逆实际上在解决同一个问题。

---

## 2. 是否存在精确解

方程：

$$
A\boldsymbol{x}=\boldsymbol{b}
$$

展开为：

$$
x_1=1,
$$

$$
x_1+x_2=2,
$$

$$
x_1+2x_2=2.
$$

从前两式得到：

$$
x_1=1,
\qquad
x_2=1.
$$

代入第三式：

$$
1+2(1)=3\neq 2.
$$

所以原方程无解。

空间层语言：

$$
\boldsymbol{b}\notin\operatorname{Col}(A).
$$

因此不能让 $A\boldsymbol{x}$ 精确等于 $\boldsymbol{b}$。

---

## 3. 改变问题：寻找最近的可达向量

我们寻找：

$$
\hat{\boldsymbol{x}}
$$

使：

$$
\|A\boldsymbol{x}-\boldsymbol{b}\|
$$

最小。

对应投影：

$$
\boldsymbol{p}=A\hat{\boldsymbol{x}}
\in\operatorname{Col}(A).
$$

残差：

$$
\boldsymbol{e}
=
\boldsymbol{b}-\boldsymbol{p}
$$

必须满足：

$$
\boldsymbol{e}\perp\operatorname{Col}(A).
$$

矩阵形式为：

$$
A^{\mathsf T}\boldsymbol{e}=\boldsymbol{0}.
$$

代入：

$$
\boldsymbol{e}
=
\boldsymbol{b}-A\hat{\boldsymbol{x}},
$$

得到：

$$
A^{\mathsf T}
(\boldsymbol{b}-A\hat{\boldsymbol{x}})
=\boldsymbol{0}.
$$

于是：

$$
\boxed{
A^{\mathsf T}A\hat{\boldsymbol{x}}
=
A^{\mathsf T}\boldsymbol{b}
}.
$$

这就是正规方程。

---

## 4. 用正规方程求最小二乘解

先计算：

$$
A^{\mathsf T}A
=
\begin{bmatrix}
3&3\\
3&5
\end{bmatrix}.
$$

再算：

$$
A^{\mathsf T}\boldsymbol{b}
=
\begin{bmatrix}
5\\6
\end{bmatrix}.
$$

因此：

$$
\begin{bmatrix}
3&3\\
3&5
\end{bmatrix}
\begin{bmatrix}
\hat{x}_1\\
\hat{x}_2
\end{bmatrix}
=
\begin{bmatrix}
5\\6
\end{bmatrix}.
$$

解得：

$$
\boxed{
\hat{\boldsymbol{x}}
=
\begin{bmatrix}
\frac{7}{6}\\
\frac{1}{2}
\end{bmatrix}
}.
$$

---

## 5. 求投影和残差

$$
\boldsymbol{p}
=A\hat{\boldsymbol{x}}
=
\begin{bmatrix}
1&0\\
1&1\\
1&2
\end{bmatrix}
\begin{bmatrix}
\frac{7}{6}\\
\frac{1}{2}
\end{bmatrix}.
$$

得到：

$$
\boxed{
\boldsymbol{p}
=
\begin{bmatrix}
\frac{7}{6}\\
\frac{5}{3}\\
\frac{13}{6}
\end{bmatrix}
}.
$$

残差：

$$
\boldsymbol{e}
=
\boldsymbol{b}-\boldsymbol{p}
=
\begin{bmatrix}
-\frac{1}{6}\\
\frac{1}{3}\\
-\frac{1}{6}
\end{bmatrix}.
$$

验证：

$$
A^{\mathsf T}\boldsymbol{e}
=
\boldsymbol{0}.
$$

因为第一列：

$$
\begin{bmatrix}
1\\1\\1
\end{bmatrix}
$$

与残差点积：

$$
-\frac{1}{6}+\frac{1}{3}-\frac{1}{6}=0,
$$

第二列：

$$
\begin{bmatrix}
0\\1\\2
\end{bmatrix}
$$

与残差点积：

$$
\frac{1}{3}-\frac{2}{6}=0.
$$

所以残差确实属于：

$$
\operatorname{Col}(A)^\perp.
$$

---

## 6. 几何解释

原目标：

$$
\boldsymbol{b}
$$

不在列空间中。

因此分解：

$$
\boldsymbol{b}
=
\boldsymbol{p}+\boldsymbol{e},
$$

其中：

$$
\boldsymbol{p}\in\operatorname{Col}(A),
$$

$$
\boldsymbol{e}\in\operatorname{Col}(A)^\perp.
$$

$\boldsymbol{p}$ 是列空间中最接近 $\boldsymbol{b}$ 的向量。

最小二乘不是“硬把无解方程变成有解”，而是改变问题：

> 不再要求精确到达 $\boldsymbol{b}$，而是到达离它最近的可达输出。

---

## 7. QR 分解

对 $A$ 的两列做 Gram–Schmidt，可取：

$$
Q=
\begin{bmatrix}
\frac{1}{\sqrt{3}}&-\frac{1}{\sqrt{2}}\\
\frac{1}{\sqrt{3}}&0\\
\frac{1}{\sqrt{3}}&\frac{1}{\sqrt{2}}
\end{bmatrix},
$$

$$
R=
\begin{bmatrix}
\sqrt{3}&\sqrt{3}\\
0&\sqrt{2}
\end{bmatrix}.
$$

于是：

$$
A=QR.
$$

由于 $Q$ 的列标准正交：

$$
Q^{\mathsf T}Q=I.
$$

最小二乘问题：

$$
\min_{\boldsymbol{x}}
\|QR\boldsymbol{x}-\boldsymbol{b}\|
$$

转化为：

$$
R\hat{\boldsymbol{x}}
=
Q^{\mathsf T}\boldsymbol{b}.
$$

计算：

$$
Q^{\mathsf T}\boldsymbol{b}
=
\begin{bmatrix}
\frac{5}{\sqrt{3}}\\
\frac{1}{\sqrt{2}}
\end{bmatrix}.
$$

于是：

$$
\begin{bmatrix}
\sqrt{3}&\sqrt{3}\\
0&\sqrt{2}
\end{bmatrix}
\begin{bmatrix}
\hat{x}_1\\
\hat{x}_2
\end{bmatrix}
=
\begin{bmatrix}
\frac{5}{\sqrt{3}}\\
\frac{1}{\sqrt{2}}
\end{bmatrix}.
$$

第二行：

$$
\sqrt{2}\hat{x}_2
=
\frac{1}{\sqrt{2}},
$$

所以：

$$
\hat{x}_2=\frac{1}{2}.
$$

第一行：

$$
\sqrt{3}\hat{x}_1
+
\sqrt{3}\cdot\frac{1}{2}
=
\frac{5}{\sqrt{3}}.
$$

两边除以 $\sqrt{3}$：

$$
\hat{x}_1+\frac{1}{2}=\frac{5}{3}.
$$

因此：

$$
\hat{x}_1=\frac{7}{6}.
$$

与正规方程完全一致。

---

## 8. 伪逆

由于 $A$ 满列秩，伪逆可写为：

$$
A^+
=(A^{\mathsf T}A)^{-1}A^{\mathsf T}.
$$

计算得到：

$$
A^+
=
\begin{bmatrix}
\frac{5}{6}&\frac{1}{3}&-\frac{1}{6}\\
-\frac{1}{2}&0&\frac{1}{2}
\end{bmatrix}.
$$

于是：

$$
A^+\boldsymbol{b}
=
\begin{bmatrix}
\frac{5}{6}&\frac{1}{3}&-\frac{1}{6}\\
-\frac{1}{2}&0&\frac{1}{2}
\end{bmatrix}
\begin{bmatrix}
1\\2\\2
\end{bmatrix}.
$$

第一分量：

$$
\frac{5}{6}+\frac{2}{3}-\frac{2}{6}
=
\frac{7}{6}.
$$

第二分量：

$$
-\frac{1}{2}+1
=
\frac{1}{2}.
$$

所以：

$$
\boxed{
A^+\boldsymbol{b}
=
\begin{bmatrix}
\frac{7}{6}\\
\frac{1}{2}
\end{bmatrix}
=
\hat{\boldsymbol{x}}
}.
$$

---

## 9. 三种方法其实在做同一件事

### 正规方程

从残差正交：

$$
A^{\mathsf T}(\boldsymbol{b}-A\hat{\boldsymbol{x}})=0
$$

得到：

$$
A^{\mathsf T}A\hat{\boldsymbol{x}}
=A^{\mathsf T}\boldsymbol{b}.
$$

### QR

直接把列空间换成标准正交基，用：

$$
A=QR
$$

把问题化成上三角系统。

### 伪逆

把“最佳反求”写成统一算子：

$$
\hat{\boldsymbol{x}}=A^+\boldsymbol{b}.
$$

三者不是三个不相关的算法，而是同一个投影问题的不同表达。

---

## 10. 本例的跨层主线

$$
\boxed{
A\boldsymbol{x}=\boldsymbol{b}\text{ 无解}
\rightarrow
\boldsymbol{b}\notin\operatorname{Col}(A)
\rightarrow
\text{投影到列空间}
\rightarrow
\text{残差正交}
\rightarrow
\text{正规方程}
\rightarrow
QR
\rightarrow
A^+
}
$$

真正要掌握的是：

> 最小二乘不是孤立的数值技巧，而是“无精确可达输出时，寻找最近可达输出”的空间问题；QR 和伪逆只是更进一步的变换与分解语言。
