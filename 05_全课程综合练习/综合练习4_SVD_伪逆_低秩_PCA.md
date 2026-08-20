# 综合练习 4：SVD、伪逆、低秩与 PCA

## 第一部分：SVD 与伪逆

给定：

$$
A=
\begin{bmatrix}
3&0\\
0&1\\
0&0
\end{bmatrix}.
$$

以及：

$$
\boldsymbol{b}
=
\begin{bmatrix}
3\\1\\2
\end{bmatrix}.
$$

回答：

1. 写出一个紧致 SVD。
2. 求 rank。
3. 求 $\operatorname{Null}(A)$。
4. 求 $\operatorname{Col}(A)$。
5. 写出伪逆 $A^+$。
6. 求 $A^+\boldsymbol{b}$。
7. 求 $AA^+\boldsymbol{b}$，并解释它的空间意义。
8. 求非零奇异值条件数：

$$
\kappa(A)=\frac{\sigma_{\max}}{\sigma_{\min}}.
$$

9. 写出最佳 rank-1 近似 $A_1$。

## 第二部分：PCA

考虑中心化数据：

$$
X=
\begin{bmatrix}
3&0\\
1&1\\
-1&-1\\
-3&0
\end{bmatrix}.
$$

1. 计算 $X^{\mathsf T}X$。
2. 求其最大特征值对应方向。
3. 解释这个方向为什么是第一主成分方向。
4. 写出样本在第一主方向上的一维投影坐标。
5. 说明 PCA 与 SVD 的关系。
