# E8 SVD、伪逆、低秩近似与 PCA 统一例题

## 十九、一个完整数值例子

考虑：

$$
A=
\begin{bmatrix}
3&0\\
0&1\\
0&0
\end{bmatrix}.
$$

它实际上已经是 SVD 形式：

$$
U=I_3,
$$

$$
V=I_2,
$$

$$
\Sigma=
\begin{bmatrix}
3&0\\
0&1\\
0&0
\end{bmatrix}.
$$

奇异值：

$$
\sigma_1=3,
\qquad
\sigma_2=1.
$$

---

### 输入

设：

$$
x=
\begin{bmatrix}
2\\
4
\end{bmatrix}.
$$

那么：

$$
Ax
=

\begin{bmatrix}
6\\
4\\
0
\end{bmatrix}.
$$

几何意义：

```text
输入 x1 方向 → 放大3倍 → 输出u1
输入 x2 方向 → 放大1倍 → 输出u2
```

---

## 二十四、一个非常直观的例子

假设：

$$
A
=

10u_1v_1^T
+
3u_2v_2^T
+
0.1u_3v_3^T.
$$

也就是说：

$$
\sigma_1=10,
$$

$$
\sigma_2=3,
$$

$$
\sigma_3=0.1.
$$

那么第三种模式非常弱。

如果做 rank-2 近似：

$$
A_2
=

10u_1v_1^T
+
3u_2v_2^T.
$$

舍弃：

$$
0.1u_3v_3^T.
$$

此时：

$$
|A-A_2|_2=0.1.
$$

所以损失很小。

---

## 20. PCA 完整例子

考虑已经中心化的数据：

$$
X=
\begin{bmatrix}
-2&-4\\
-1&-2\\
0&0\\
1&2\\
2&4
\end{bmatrix}
$$

所有样本都位于直线：

$$
y=2x
$$

### 20.1 秩与主方向

矩阵可以写成外积：

$$
X=
\begin{bmatrix}
-2\\-1\\0\\1\\2
\end{bmatrix}
\begin{bmatrix}1&2\end{bmatrix}
$$

所以秩为 1。

第一右奇异向量为：

$$
\boldsymbol v_1
=
\frac1{\sqrt5}
\begin{bmatrix}1\\2\end{bmatrix}
$$

它正是数据直线方向。

第二方向可取：

$$
\boldsymbol v_2
=
\frac1{\sqrt5}
\begin{bmatrix}-2\\1\end{bmatrix}
$$

该方向与数据直线正交。

### 20.2 奇异值

两个外积向量长度分别为：

$$
\sqrt{10},
\qquad
\sqrt5
$$

所以唯一非零奇异值为：

$$
\sigma_1=\sqrt{50}=5\sqrt2
$$

第二奇异值为零。

### 20.3 方差

样本数量 $m=5$，第一主成分方差为：

$$
\frac{\sigma_1^2}{m-1}
=
\frac{50}{4}
=12.5
$$

第二主成分方差为零。

### 20.4 降维

只保留第一主成分不会损失任何信息，因为数据本来就位于一条直线上。

每个样本由一个标量得分表示：

$$
Z=X\boldsymbol v_1
$$

这将二维数据无损压缩为一维数据。
