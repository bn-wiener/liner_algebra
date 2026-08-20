# 练习5：伪逆、条件数、低秩近似与 PCA

   5. 判断降到一维是否损失信息。

6. 使用伪逆求：

   $$
   A\boldsymbol x\approx
   \begin{bmatrix}2\\0\\1\end{bmatrix}
   $$

   的最小范数最小二乘解、投影和残差。

### 理解题

9. 证明 $AA^+$ 是列空间上的正交投影。

10. 证明 $A^+A$ 是行空间上的正交投影。

11. 为什么伪逆解 $A^+\boldsymbol b$ 与零空间正交？

12. 为什么小奇异值会使逆问题对噪声敏感？

14. 设矩阵的奇异值为：

   $$
   8,\qquad3,\qquad1,\qquad0
   $$

15. 设：

   $$
   A=
   \begin{bmatrix}
   1&0\\
   0&0.01
   \end{bmatrix},
   \qquad
   \boldsymbol b=
   \begin{bmatrix}1\\0.01\end{bmatrix}
   $$

16. 证明截断 SVD：

   $$
   A_k=
   \sum_{i=1}^k
   \sigma_i\boldsymbol u_i\boldsymbol v_i^\mathsf T
   $$

   的秩不超过 $k$。

17. 一个 $1000\times800$ 矩阵使用秩 20 的截断 SVD 存储：

18. 对已经中心化的数据：

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

19. 设：

   $$
   A=U\Sigma V^\mathsf T
   $$

   推导：

   $$
   \|A\|_F^2
   =
   \sum_i\sigma_i^2
   $$

20. 比较下列场景应优先使用特征分解、QR 还是 SVD：
