# E1 从一组向量到 span、相关性、基与维数

本例直接采用旧资料第一阶段的贯穿例子，并把它作为 S1～S6 的集中复习。

设：

$$
\boldsymbol v_1=
\begin{bmatrix}1\\1\\0\end{bmatrix},
\quad
\boldsymbol v_2=
\begin{bmatrix}0\\1\\1\end{bmatrix},
\quad
\boldsymbol v_3=
\begin{bmatrix}1\\2\\1\end{bmatrix}
$$

观察到：

$$
\boldsymbol v_3=\boldsymbol v_1+\boldsymbol v_2
$$

因此：

- 三个向量线性相关。
- $\boldsymbol v_3$ 是冗余向量。
- 三个向量的张成空间与 $\boldsymbol v_1,\boldsymbol v_2$ 的张成空间相同。
- $\boldsymbol v_1,\boldsymbol v_2$ 不互为倍数，因此线性无关。
- $\boldsymbol v_1,\boldsymbol v_2$ 构成该张成空间的一组基。
- 该张成空间的维数为 2，是 $\mathbb R^3$ 中一个经过原点的平面。

这个例子把“线性组合、张成、相关性、基和维数”连接成了一个整体。

## 复盘

完成这个例子时不要只得到一个数值答案，而要按顺序回答：

1. 原向量组张成什么集合？
2. 是否存在冗余方向？
3. 怎样删去冗余得到一组基？
4. 为什么这组基能保证唯一坐标？
5. 基向量数量为什么就是维数？
