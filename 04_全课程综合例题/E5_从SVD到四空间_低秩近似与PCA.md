# 综合例题 E5：从 SVD 到四空间、低秩近似与 PCA

## 1. 题目

考虑已经中心化的数据矩阵：

$$
X=
\begin{bmatrix}
2&0\\
0&1\\
-2&0\\
0&-1
\end{bmatrix}.
$$

每一行是一条二维样本。两列均值都是 $0$，所以数据已中心化。

要求：

1. 求 $X^{\mathsf T}X$；
2. 求奇异值与右奇异向量；
3. 写出一个 SVD；
4. 解释 rank 和四个基本子空间；
5. 求最佳 rank-1 近似；
6. 从 PCA 角度解释第一主方向；
7. 说明 SVD、低秩近似和 PCA 为什么是同一组方向的不同用途。

---

## 2. 先观察数据几何

四个点是：

$$
(2,0),
\quad
(0,1),
\quad
(-2,0),
\quad
(0,-1).
$$

数据沿 $x$ 方向的展开更大：范围大致到 $\pm2$；沿 $y$ 方向只到 $\pm1$。

所以在计算前就应预期：

> 第一主方向应当是 $x$ 轴方向。

SVD 将把这种直觉变成精确结果。

---

## 3. 计算 $X^{\mathsf T}X$

$$
X^{\mathsf T}X
=
\begin{bmatrix}
2&0&-2&0\\
0&1&0&-1
\end{bmatrix}
\begin{bmatrix}
2&0\\
0&1\\
-2&0\\
0&-1
\end{bmatrix}.
$$

得到：

$$
\boxed{
X^{\mathsf T}X
=
\begin{bmatrix}
8&0\\
0&2
\end{bmatrix}
}.
$$

它已经是对角矩阵。

因此右奇异向量可以直接取：

$$
\boldsymbol{v}_1
=
\begin{bmatrix}
1\\0
\end{bmatrix},
\qquad
\boldsymbol{v}_2
=
\begin{bmatrix}
0\\1
\end{bmatrix}.
$$

对应特征值：

$$
8,
\qquad
2.
$$

所以奇异值：

$$
\boxed{
\sigma_1=\sqrt{8}=2\sqrt{2},
\qquad
\sigma_2=\sqrt{2}
}.
$$

---

## 4. 求左奇异向量

由：

$$
X\boldsymbol{v}_i
=
\sigma_i\boldsymbol{u}_i.
$$

### 第一个方向

$$
X\boldsymbol{v}_1
=
\begin{bmatrix}
2\\0\\-2\\0
\end{bmatrix}.
$$

其长度：

$$
2\sqrt{2}.
$$

所以：

$$
\boldsymbol{u}_1
=
\frac{1}{2\sqrt{2}}
\begin{bmatrix}
2\\0\\-2\\0
\end{bmatrix}
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1\\0\\-1\\0
\end{bmatrix}.
$$

### 第二个方向

$$
X\boldsymbol{v}_2
=
\begin{bmatrix}
0\\1\\0\\-1
\end{bmatrix}.
$$

其长度：

$$
\sqrt{2}.
$$

因此：

$$
\boldsymbol{u}_2
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
0\\1\\0\\-1
\end{bmatrix}.
$$

于是紧致 SVD 可以写为：

$$
X
=
U_r\Sigma_rV_r^{\mathsf T},
$$

其中：

$$
U_r
=
\begin{bmatrix}
\frac{1}{\sqrt{2}}&0\\
0&\frac{1}{\sqrt{2}}\\
-\frac{1}{\sqrt{2}}&0\\
0&-\frac{1}{\sqrt{2}}
\end{bmatrix},
$$

$$
\Sigma_r
=
\begin{bmatrix}
2\sqrt{2}&0\\
0&\sqrt{2}
\end{bmatrix},
$$

$$
V_r
=I_2.
$$

---

## 5. rank 与四个基本子空间

两个奇异值都非零，因此：

$$
\operatorname{rank}(X)=2.
$$

因为 $X$ 有 2 列：

$$
\operatorname{nullity}(X)=0.
$$

所以：

$$
\operatorname{Null}(X)=\{\boldsymbol{0}\}.
$$

行空间：

$$
\operatorname{Row}(X)=\mathbb{R}^2.
$$

右奇异向量：

$$
\boldsymbol{v}_1,\boldsymbol{v}_2
$$

正好给出行空间的一组标准正交基。

列空间由：

$$
\boldsymbol{u}_1,\boldsymbol{u}_2
$$

张成：

$$
\operatorname{Col}(X)
=
\operatorname{span}\{\boldsymbol{u}_1,\boldsymbol{u}_2\}.
$$

因为输出空间是 $\mathbb{R}^4$、rank 为 2，所以左零空间维数为：

$$
4-2=2.
$$

完整 SVD 中，$U$ 剩余的两个正交列正好可以作为：

$$
\operatorname{Null}(X^{\mathsf T})
$$

的一组标准正交基。

这就是 SVD 与四个基本子空间的直接连接。

---

## 6. 最佳 rank-1 近似

SVD 展开为：

$$
X
=
\sigma_1\boldsymbol{u}_1\boldsymbol{v}_1^{\mathsf T}
+
\sigma_2\boldsymbol{u}_2\boldsymbol{v}_2^{\mathsf T}.
$$

最佳 rank-1 近似只保留最大奇异值项：

$$
X_1
=
\sigma_1\boldsymbol{u}_1\boldsymbol{v}_1^{\mathsf T}.
$$

因为：

$$
\boldsymbol{v}_1
=
\begin{bmatrix}
1\\0
\end{bmatrix},
$$

所以：

$$
X_1
=
\begin{bmatrix}
2&0\\
0&0\\
-2&0\\
0&0
\end{bmatrix}.
$$

这相当于把所有数据只保留在 $x$ 方向上的变化。

原矩阵中第二方向的变化被舍弃。

---

## 7. PCA 解释

因为数据已经中心化，PCA 的主方向可以由：

$$
X^{\mathsf T}X
$$

的特征向量得到。

这里：

$$
X^{\mathsf T}X
=
\begin{bmatrix}
8&0\\
0&2
\end{bmatrix}.
$$

最大特征值为：

$$
8,
$$

对应方向：

$$
\boxed{
\boldsymbol{v}_1
=
\begin{bmatrix}
1\\0
\end{bmatrix}
}.
$$

所以第一主成分方向就是 $x$ 轴。

数据投影到第一主方向的坐标为：

$$
X\boldsymbol{v}_1
=
\begin{bmatrix}
2\\0\\-2\\0
\end{bmatrix}.
$$

这把每个二维样本压缩成一个一维数值。

---

## 8. 解释方差比例

两个方向对应能量：

$$
\sigma_1^2=8,
\qquad
\sigma_2^2=2.
$$

总量：

$$
8+2=10.
$$

第一方向解释比例：

$$
\frac{8}{10}=0.8.
$$

所以只保留第一主方向，就保留了约 $80\%$ 的总平方变化量。

第二方向贡献：

$$
\frac{2}{10}=0.2.
$$

这就是 PCA 中“解释方差比例”的几何来源。

---

## 9. 同一组奇异方向的三种用途

### SVD

描述变换：

$$
X=U\Sigma V^{\mathsf T}.
$$

### 低秩近似

只保留最大奇异值方向：

$$
X_1
=
\sigma_1\boldsymbol{u}_1\boldsymbol{v}_1^{\mathsf T}.
$$

### PCA

把：

$$
\boldsymbol{v}_1
$$

解释为数据最主要变化方向，并把样本投影到该方向。

所以它们不是三套完全不同的方法。

它们共享的是：

> 找到最重要的正交方向，并按奇异值大小排序。

---

## 10. 本例统一主线

$$
\boxed{
X^{\mathsf T}X
\rightarrow
\text{右奇异向量}
\rightarrow
\text{奇异值}
\rightarrow
X=U\Sigma V^{\mathsf T}
\rightarrow
\text{四个基本子空间}
\rightarrow
\text{rank-1 近似}
\rightarrow
\text{PCA}
}
$$

这个例题展示了整个三层体系的一个终点：

- 第一层的矩阵运算仍然存在；
- 第二层的 rank、正交基、四空间全部被回收；
- 第三层把它们统一成一个“正交换基 → 缩放 → 正交换基”的变换结构。
