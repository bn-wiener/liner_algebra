# SVD（奇异值分解）

SVD（Singular Value Decomposition，奇异值分解）可以理解为：

> 把任意线性变换 $A$ 拆成“旋转 / 换正交坐标系 → 沿互相垂直的方向缩放 → 再旋转 / 换正交坐标系”三个最简单的步骤。

如果已经理解了特征值和对角化，那么 SVD 最值得建立的心智模型是：

$$
\boxed{
\text{SVD 是“对角化思想”对任意矩阵的推广}
}
$$

并且它比普通对角化更强：

$$
\boxed{
\text{任意 }m\times n\text{ 矩阵都可以做 SVD}
}
$$

---

# 一、先从几何上理解 SVD

考虑一个线性变换：

$$
T(x)=Ax
$$

例如二维情况下，一个单位圆：

$$
\|x\|=1
$$

经过矩阵 $A$ 之后，一般会变成一个椭圆。

也就是说：

$$
\boxed{
\text{单位圆}
\xrightarrow{A}
\text{椭圆}
}
$$

SVD 告诉我们，这个复杂的变形其实可以拆成三个非常简单的动作：

$$
\boxed{
x
\xrightarrow{V^T}
\text{旋转 / 换输入正交坐标系}
\xrightarrow{\Sigma}
\text{沿坐标轴缩放}
\xrightarrow{U}
\text{旋转 / 换输出正交坐标系}
}
$$

因此：

$$
\boxed{
A=U\Sigma V^T
}
$$

从右往左作用：

$$
x
\rightarrow
V^Tx
\rightarrow
\Sigma V^Tx
\rightarrow
U\Sigma V^Tx
$$

最终：

$$
Ax=U\Sigma V^Tx
$$

---

# 二、三个矩阵分别干什么？

这是理解 SVD 最关键的部分。

$$
A=U\Sigma V^T
$$

不要先记公式，先记三个角色：

$$
\boxed{
V^T：把输入向量转换到右奇异向量基下
}
$$

$$
\boxed{
\Sigma：对各个独立方向分别缩放
}
$$

$$
\boxed{
U：把输出奇异向量基下的坐标还原到输出空间
}
$$

所以 SVD 本质上是在寻找：

> $A$ 对哪些输入方向的作用最简单？

这些方向就是右奇异向量。

---

# 三、从一个二维几何图开始

假设：

$$
A:\mathbb R^2\rightarrow\mathbb R^2
$$

我们在输入空间中找到两个互相垂直的方向：

$$
v_1,\quad v_2
$$

满足：

$$
v_1^Tv_2=0
$$

经过 $A$ 后：

$$
Av_1=\sigma_1u_1
$$

$$
Av_2=\sigma_2u_2
$$

其中：

$$
u_1\perp u_2
$$

因此有：

$$
\boxed{
Av_i=\sigma_i u_i
}
$$

这其实就是 SVD 最核心的公式。

它说的是：

> 输入空间中的特殊方向 $v_i$，经过 $A$ 后，变成输出空间中的方向 $u_i$，长度被放大 $\sigma_i$ 倍。

---

# 四、这和特征向量有什么区别？

特征向量满足：

$$
Av_i=\lambda_i v_i
$$

意思是：

$$
\boxed{
\text{输入方向 = 输出方向}
}
$$

方向不改变，只缩放。

而 SVD 是：

$$
Av_i=\sigma_i u_i
$$

这里一般：

$$
u_i\neq v_i
$$

也就是说：

$$
\boxed{
\text{输入方向 }v_i
\rightarrow
\text{输出方向 }u_i
}
$$

所以特征值研究的是：

> 哪些方向经过 $A$ 后仍然保持自己的方向？

而 SVD 研究的是：

> 能不能找到一组正交输入方向，使得经过 $A$ 后，它们对应一组正交输出方向？

答案是：

$$
\boxed{
\text{永远可以}
}
$$

这就是 SVD 强大的地方。

---

# 五、为什么会有两个基 $V$ 和 $U$？

这是理解 SVD 和对角化区别的关键。

矩阵：

$$
A\in\mathbb R^{m\times n}
$$

描述：

$$
A:\mathbb R^n\rightarrow\mathbb R^m
$$

注意：

- 输入空间是 $\mathbb R^n$
- 输出空间是 $\mathbb R^m$

甚至维数都可能不一样。

因此不能总要求：

$$
Av=\lambda v
$$

因为：

$$
v\in\mathbb R^n
$$

而：

$$
Av\in\mathbb R^m
$$

它们可能根本不属于同一个空间。

所以 SVD 使用两组基。

输入空间：

$$
\boxed{
v_1,v_2,\cdots,v_n
}
$$

输出空间：

$$
\boxed{
u_1,u_2,\cdots,u_m
}
$$

于是：

$$
Av_i=\sigma_i u_i
$$

就完全合理了。

---

# 六、从“基与坐标”理解 SVD

假设：

$$
V=
\begin{bmatrix}
|&|&&|\\
v_1&v_2&\cdots&v_n\\
|&|&&|
\end{bmatrix}
$$

因为 $v_i$ 是一组标准正交基，所以：

$$
V^{-1}=V^T
$$

任意输入：

$$
x=c_1v_1+c_2v_2+\cdots+c_nv_n
$$

写成矩阵：

$$
x=Vc
$$

所以：

$$
c=V^Tx
$$

这意味着：

$$
\boxed{
V^Tx
}
$$

并不是神秘的“旋转”。

更本质地说，它是在求：

$$
\boxed{
x\text{ 在 }(v_1,\cdots,v_n)\text{ 基下的坐标}
}
$$

这和对角化中的：

$$
P^{-1}x
$$

把 $x$ 转到特征向量基坐标，是完全同一个逻辑。

区别只是 SVD 的 $V$ 是正交矩阵，所以：

$$
V^{-1}=V^T
$$

---

# 七、然后 $\Sigma$ 做什么？

得到：

$$
c=
\begin{bmatrix}
c_1\\
c_2\\
\vdots
\end{bmatrix}
$$

之后：

$$
\Sigma c
$$

假设二维情况：

$$
\Sigma=
\begin{bmatrix}
\sigma_1&0\\
0&\sigma_2
\end{bmatrix}
$$

那么：

$$
\Sigma
\begin{bmatrix}
c_1\\
c_2
\end{bmatrix}
=
\begin{bmatrix}
\sigma_1c_1\\
\sigma_2c_2
\end{bmatrix}
$$

更准确地说：

$$
\boxed{
[Ax]_U=\Sigma[x]_V
}
$$

也就是说，$\Sigma$ 把输入空间中 $V$ 基下的坐标，变成输出空间中 $U$ 基下的坐标。

---

# 八、最后为什么还要乘 $U$？

经过 $\Sigma$ 后：

$$
\begin{bmatrix}
\sigma_1c_1\\
\sigma_2c_2
\end{bmatrix}
$$

这些其实是输出向量 $Ax$ 在：

$$
u_1,u_2
$$

这组基下的坐标。

而：

$$
U=
\begin{bmatrix}
|&|\\
u_1&u_2\\
|&|
\end{bmatrix}
$$

所以：

$$
U
\begin{bmatrix}
\sigma_1c_1\\
\sigma_2c_2
\end{bmatrix}
=
\sigma_1c_1u_1+\sigma_2c_2u_2
$$

最终回到输出空间中的实际向量。

因此完整过程是：

$$
\boxed{
x
\xrightarrow{V^T}
[x]_V
\xrightarrow{\Sigma}
[Ax]_U
\xrightarrow{U}
Ax
}
$$

这比“旋转—缩放—旋转”更加本质。

---

# 九、SVD 最本质的坐标解释

现在可以写出：

$$
A=U\Sigma V^T
$$

真正的含义是：

$$
\boxed{
\text{标准输入坐标}
\xrightarrow{V^T}
\text{输入奇异向量基坐标}
\xrightarrow{\Sigma}
\text{输出奇异向量基坐标}
\xrightarrow{U}
\text{标准输出坐标}
}
$$

这和对角化：

$$
A=PDP^{-1}
$$

几乎完全平行。

对角化：

$$
x
\xrightarrow{P^{-1}}
[x]_P
\xrightarrow{D}
[Ax]_P
\xrightarrow{P}
Ax
$$

SVD：

$$
x
\xrightarrow{V^T}
[x]_V
\xrightarrow{\Sigma}
[Ax]_U
\xrightarrow{U}
Ax
$$

最大的区别是：

$$
\boxed{
\text{对角化只有一组基 }P
}
$$

而：

$$
\boxed{
\text{SVD 使用输入基 }V\text{ 和输出基 }U
}
$$

---

# 十、为什么 SVD 一定存在？

即使 $A$ 不是方阵，或者 $A$ 不能对角化，我们仍然可以构造：

$$
A^TA
$$

注意：

$$
A\in\mathbb R^{m\times n}
$$

那么：

$$
A^TA\in\mathbb R^{n\times n}
$$

而且 $A^TA$ 一定是实对称矩阵，因为：

$$
(A^TA)^T=A^TA
$$

实对称矩阵一定可以正交对角化，所以：

$$
A^TA=V\Lambda V^T
$$

其中：

$$
\Lambda=
\begin{bmatrix}
\lambda_1&&\\
&\lambda_2&\\
&&\ddots
\end{bmatrix}
$$

---

# 十一、奇异值从哪里来？

由于：

$$
A^TAv_i=\lambda_i v_i
$$

两边左乘 $v_i^T$：

$$
v_i^TA^TAv_i=\lambda_i v_i^Tv_i
$$

因为：

$$
\|v_i\|=1
$$

所以：

$$
\lambda_i=v_i^TA^TAv_i
$$

注意：

$$
v_i^TA^TAv_i
=
(Av_i)^T(Av_i)
=
\|Av_i\|^2
$$

因此：

$$
\boxed{
\lambda_i=\|Av_i\|^2
}
$$

定义：

$$
\boxed{
\sigma_i=\sqrt{\lambda_i}
}
$$

于是：

$$
\boxed{
\sigma_i=\|Av_i\|
}
$$

这就是奇异值真正的几何意义：

$$
\boxed{
\sigma_i=\text{$A$ 在方向 }v_i\text{ 上的伸缩倍率}
}
$$

---

# 十二、为什么 $u_i$ 可以这样得到？

因为：

$$
Av_i
$$

长度为：

$$
\sigma_i
$$

因此把它归一化：

$$
\boxed{
u_i=\frac{Av_i}{\sigma_i}
}
$$

于是：

$$
Av_i=\sigma_i u_i
$$

这正是 SVD 的核心关系。

所以整个 SVD 的构造过程是：

$$
A^TA
\rightarrow
v_i
\rightarrow
\sigma_i
\rightarrow
u_i
$$

具体：

$$
\boxed{
A^TAv_i=\sigma_i^2v_i
}
$$

然后：

$$
\boxed{
u_i=\frac{Av_i}{\sigma_i}
}
$$

---

# 十三、一个具体数值例子

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

# 十四、再看一个带方向变化的例子

考虑：

$$
A=
\begin{bmatrix}
1&1\\
1&1
\end{bmatrix}
$$

先计算：

$$
A^TA
=
\begin{bmatrix}
1&1\\
1&1
\end{bmatrix}
\begin{bmatrix}
1&1\\
1&1
\end{bmatrix}
=
\begin{bmatrix}
2&2\\
2&2
\end{bmatrix}
$$

特征值：

$$
\lambda_1=4,\qquad \lambda_2=0
$$

所以奇异值：

$$
\sigma_1=2,\qquad \sigma_2=0
$$

对应单位特征向量：

$$
v_1=
\frac{1}{\sqrt2}
\begin{bmatrix}
1\\
1
\end{bmatrix}
$$

$$
v_2=
\frac{1}{\sqrt2}
\begin{bmatrix}
1\\
-1
\end{bmatrix}
$$

因此：

$$
V=
\frac1{\sqrt2}
\begin{bmatrix}
1&1\\
1&-1
\end{bmatrix}
$$

对于 $v_1$：

$$
Av_1
=
\frac1{\sqrt2}
\begin{bmatrix}
2\\
2
\end{bmatrix}
=
2
\frac1{\sqrt2}
\begin{bmatrix}
1\\
1
\end{bmatrix}
$$

所以：

$$
\sigma_1=2
$$

而：

$$
Av_2
=
A
\frac1{\sqrt2}
\begin{bmatrix}
1\\
-1
\end{bmatrix}
=
\begin{bmatrix}
0\\
0
\end{bmatrix}
$$

所以：

$$
\sigma_2=0
$$

几何上：

$$
\boxed{
v_1\text{ 方向被放大 2 倍}
}
$$

而：

$$
\boxed{
v_2\text{ 方向被完全压扁}
}
$$

所以整个二维平面最终被压缩成一条直线。

这也立刻告诉我们：

$$
\operatorname{rank}(A)=1
$$

---

# 十五、SVD 和秩之间的关系

设：

$$
A=U\Sigma V^T
$$

其中：

$$
\Sigma=
\begin{bmatrix}
\sigma_1&&&\\
&\sigma_2&&\\
&&\ddots&\\
&&&\sigma_r\\
&&&&0
\end{bmatrix}
$$

如果只有 $r$ 个非零奇异值：

$$
\sigma_1,\ldots,\sigma_r>0
$$

那么：

$$
\boxed{
\operatorname{rank}(A)=r
}
$$

因为正交矩阵 $U,V$ 都不会改变维数。

真正让空间“损失维度”的就是 $\Sigma$ 中的：

$$
\sigma_i=0
$$

因此：

$$
\boxed{
\text{秩 = 没有被 }A\text{ 压扁掉的独立方向数量}
}
$$

---

# 十六、SVD 和零空间

由：

$$
Av_i=\sigma_i u_i
$$

如果：

$$
\sigma_i=0
$$

那么：

$$
Av_i=0
$$

因此：

$$
v_i\in N(A)
$$

所以：

$$
\boxed{
N(A)
=
\operatorname{span}\{v_i:\sigma_i=0\}
}
$$

也就是说：

$$
\boxed{
A\text{ 的零空间}
=
\text{零奇异值对应的右奇异向量张成的空间}
}
$$

---

# 十七、SVD 和列空间

对于所有非零奇异值：

$$
Av_i=\sigma_i u_i
$$

因为：

$$
\sigma_i\neq0
$$

所以：

$$
u_i=\frac1{\sigma_i}Av_i
$$

因此：

$$
u_i\in C(A)
$$

所以：

$$
\boxed{
C(A)
=
\operatorname{span}(u_1,\ldots,u_r)
}
$$

也就是说，非零奇异值对应的 $u_i$ 构成 $A$ 列空间的一组标准正交基。

---

# 十八、SVD 和行空间

同理：

$$
v_1,\ldots,v_r
$$

张成：

$$
C(A^T)
$$

也就是 $A$ 的行空间：

$$
\boxed{
C(A^T)
=
\operatorname{span}(v_1,\ldots,v_r)
}
$$

因此 SVD 几乎一次性把四个基本子空间全部告诉我们：

$$
\boxed{
\begin{aligned}
C(A)&=\operatorname{span}(u_1,\ldots,u_r)\\
N(A^T)&=\operatorname{span}(u_{r+1},\ldots,u_m)\\
C(A^T)&=\operatorname{span}(v_1,\ldots,v_r)\\
N(A)&=\operatorname{span}(v_{r+1},\ldots,v_n)
\end{aligned}
}
$$

---

# 十九、SVD 和 $Ax=b$

回到主线：

$$
Ax=b
$$

利用：

$$
A=U\Sigma V^T
$$

得到：

$$
U\Sigma V^Tx=b
$$

左乘：

$$
U^T
$$

因为：

$$
U^TU=I
$$

所以：

$$
\Sigma V^Tx=U^Tb
$$

令：

$$
z=V^Tx
$$

$$
c=U^Tb
$$

那么方程变成：

$$
\boxed{
\Sigma z=c
}
$$

本来复杂的：

$$
Ax=b
$$

经过换基后，变成了一组彼此独立的一维方程。

因为：

$$
\Sigma=
\begin{bmatrix}
\sigma_1&&\\
&\sigma_2&\\
&&\ddots
\end{bmatrix}
$$

所以：

$$
\sigma_i z_i=c_i
$$

即：

$$
z_i=\frac{c_i}{\sigma_i}
$$

最后：

$$
x=Vz
$$

所以 SVD 本质上是在说：

$$
\boxed{
复杂的\ Ax=b
\rightarrow
换一个最合适的输入基和输出基
\rightarrow
变成一组独立的一维方程
}
$$

---

# 二十、这就是伪逆的来源

如果：

$$
\sigma_i\neq0
$$

就求倒数：

$$
\frac1{\sigma_i}
$$

定义：

$$
\Sigma^+
$$

把非零奇异值取倒数：

$$
\sigma_i
\rightarrow
\frac1{\sigma_i}
$$

零仍然保持零。

于是：

$$
\boxed{
A^+=V\Sigma^+U^T
}
$$

这就是 Moore-Penrose 伪逆。

因此：

$$
\boxed{
x=A^+b
}
$$

能够统一处理：

- 可逆方阵
- 奇异方阵
- 超定方程
- 欠定方程
- 最小二乘
- 最小范数解

所以从 $Ax=b$ 的角度：

$$
\boxed{
\text{SVD 几乎是研究线性方程组最完整的坐标系统}
}
$$

---

# 二十一、SVD 和最小二乘

超定方程：

$$
Ax\approx b
$$

最小化：

$$
\min_x\|Ax-b\|^2
$$

普通推导得到正规方程：

$$
A^TAx=A^Tb
$$

而 SVD：

$$
A=U\Sigma V^T
$$

直接给出：

$$
x^*=A^+b
$$

即：

$$
\boxed{
x^*=V\Sigma^+U^Tb
}
$$

数值计算中，SVD 通常比直接求：

$$
(A^TA)^{-1}A^T
$$

更加稳定。

---

# 二十二、SVD 和 PCA

假设数据矩阵：

$$
X
$$

做 SVD：

$$
X=U\Sigma V^T
$$

那么：

$$
X^TX
=
V\Sigma^TU^TU\Sigma V^T
$$

由于：

$$
U^TU=I
$$

所以：

$$
X^TX
=
V\Sigma^T\Sigma V^T
$$

因此 $V$ 就是 $X^TX$ 的特征向量矩阵。

而：

$$
\sigma_i^2
$$

就是相应特征值。

PCA 本质上寻找：

$$
\boxed{
\text{数据方差最大的正交方向}
}
$$

而这些方向正是 SVD 的右奇异向量。

因此：

$$
\boxed{
\text{PCA}\leftrightarrow\text{SVD}
}
$$

---

# 二十三、SVD 和低秩近似

完整 SVD 可以写成：

$$
A
=
\sigma_1u_1v_1^T
+
\sigma_2u_2v_2^T
+\cdots+
\sigma_ru_rv_r^T
$$

也就是：

$$
\boxed{
A=\sum_{i=1}^r\sigma_i u_i v_i^T
}
$$

每一个：

$$
u_i v_i^T
$$

都是一个 rank-1 矩阵。

所以一个复杂矩阵，其实是若干个秩 1 线性变换叠加起来的。

如果：

$$
\sigma_1\ge\sigma_2\ge\cdots
$$

而后面的奇异值非常小，我们可以只保留前 $k$ 个：

$$
A_k
=
\sum_{i=1}^k
\sigma_i u_i v_i^T
$$

于是：

$$
\boxed{
A\approx A_k
}
$$

这就是下面这些方法背后的核心数学原理：

- 图像压缩
- PCA 降维
- 推荐系统
- LSA
- 去噪
- 模型压缩

---

# 二十四、为什么奇异值越大越“重要”？

因为：

$$
Av_i=\sigma_i u_i
$$

如果：

$$
\sigma_1=100
$$

说明 $A$ 对 $v_1$ 方向非常敏感。

而：

$$
\sigma_{10}=0.001
$$

说明 $A$ 对 $v_{10}$ 方向几乎完全压缩。

因此：

$$
\boxed{
\sigma_i
=
A\text{ 对第 }i\text{ 个特殊输入方向的响应强度}
}
$$

大的奇异值：

$$
\rightarrow
\text{主要结构}
$$

小的奇异值：

$$
\rightarrow
\text{弱结构 / 可能的噪声}
$$

零奇异值：

$$
\rightarrow
\text{完全丢失的信息}
$$

---

# 二十五、SVD 和条件数

如果：

$$
\sigma_{\max}
$$

很大，而：

$$
\sigma_{\min}
$$

非常小，那么某些方向被强烈放大，某些方向几乎被压扁。

定义二范数条件数：

$$
\boxed{
\kappa_2(A)
=
\frac{\sigma_{\max}}{\sigma_{\min}}
}
$$

例如：

$$
\sigma_1=1000,\qquad \sigma_2=0.001
$$

则：

$$
\kappa_2(A)=10^6
$$

说明：

$$
Ax=b
$$

可能非常不稳定。

因为恢复 $x$ 时需要：

$$
\frac1{\sigma_2}=1000
$$

原本很小的误差会被巨大放大。

---

# 二十六、SVD 和特征分解到底是什么关系？

## 特征分解

对于可以对角化的方阵：

$$
A=PDP^{-1}
$$

寻找：

$$
Av_i=\lambda_i v_i
$$

强调：

$$
\boxed{
A\text{ 自己不改变的方向}
}
$$

## SVD

对于任意矩阵：

$$
A=U\Sigma V^T
$$

寻找：

$$
Av_i=\sigma_i u_i
$$

强调：

$$
\boxed{
A\text{ 的输入主方向和输出主方向}
}
$$

而 SVD 的 $V$ 来自 $A^TA$ 的特征分解：

$$
\boxed{
A^TA=V\Sigma^T\Sigma V^T
}
$$

如果 $A$ 是方阵且采用对应维度的简化写法，可以写成：

$$
A^TA=V\Sigma^2V^T
$$

同时：

$$
\boxed{
AA^T=U\Sigma\Sigma^TU^T
}
$$

因此：

$$
\boxed{
v_i=A^TA\text{ 的特征向量}
}
$$

$$
\boxed{
u_i=AA^T\text{ 的特征向量}
}
$$

$$
\boxed{
\sigma_i^2=\text{对应特征值}
}
$$

---

# 二十七、SVD 和 QR 的区别

QR：

$$
A=QR
$$

重点是：

$$
\boxed{
\text{给列空间找一组标准正交基}
}
$$

其中 $Q$ 描述列空间的正交方向，$R$ 记录原列向量在这些正交方向下的坐标。

而 SVD：

$$
A=U\Sigma V^T
$$

做得更彻底：

- $V$：输入空间最佳正交基
- $U$：输出空间最佳正交基
- $\Sigma$：各个独立方向上的伸缩量

所以可以粗略理解：

$$
\boxed{
\text{QR：主要整理列空间}
}
$$

而：

$$
\boxed{
\text{SVD：同时整理输入空间和输出空间}
}
$$

---

# 二十八、把 QR、特征分解、SVD 放到一起

对于：

$$
\boxed{
A:\mathbb R^n\rightarrow\mathbb R^m
}
$$

## QR

问：

> $A$ 的列空间有哪些标准正交方向？

得到：

$$
A=QR
$$

## 特征分解

问：

> 哪些方向经过 $A$ 后方向保持不变？

得到：

$$
A=PDP^{-1}
$$

## SVD

问：

> 有没有输入正交方向，使 $A$ 把它们映射成输出正交方向？

得到：

$$
A=U\Sigma V^T
$$

---

# 二十九、最终统一到 $T(x)\leftrightarrow Ax=b$

对于：

$$
T(x)=Ax
$$

SVD 告诉我们：

$$
A=U\Sigma V^T
$$

因此：

$$
T(x)=U\Sigma V^Tx
$$

完整过程是：

先求输入奇异基坐标：

$$
\boxed{
z=V^Tx
}
$$

然后每个方向独立缩放：

$$
\boxed{
w=\Sigma z
}
$$

注意此时更准确地说：

$$
w=[Ax]_U
$$

最后从输出奇异基恢复：

$$
\boxed{
Ax=Uw
}
$$

所以：

$$
\boxed{
x
\xrightarrow{V^T}
[x]_V
\xrightarrow{\Sigma}
[Ax]_U
\xrightarrow{U}
Ax
}
$$

---

# 三十、SVD 的统一心智模型

```text
                         线性变换
                       T(x) = Ax
                           │
             ┌─────────────┴─────────────┐
             │                           │
          输入空间                    输出空间
           R^n                         R^m
             │                           │
      右奇异向量基 V              左奇异向量基 U
             │                           │
             └──────────┐   ┌────────────┘
                        │   │
                        ▼   ▼

                     A = UΣVᵀ

           x ──Vᵀ──> [x]V ──Σ──> [Ax]U ──U──> Ax
                              │
                              │
                    每个方向彼此独立
                              │
                       σ₁, σ₂, ..., σᵣ
                              │
             ┌────────────────┼────────────────┐
             ▼                ▼                ▼
           σᵢ > 0          σᵢ 很小           σᵢ = 0
             │                │                │
         有效方向           弱方向           丢失方向
             │                │                │
          列空间           低秩近似          零空间
          行空间              PCA              rank
```

最后浓缩成一句话：

$$
\boxed{
\text{SVD 就是在输入空间和输出空间各找一组最合适的正交基，}
}
$$

$$
\boxed{
\text{使任意复杂线性变换 }A\text{ 在这两组基之间只剩下独立的轴向缩放。}
}
$$

也就是：

$$
\boxed{
A
=
\underbrace{U}_{\text{输出基}}
\underbrace{\Sigma}_{\text{坐标间的独立伸缩}}
\underbrace{V^T}_{\text{输入换基}}
}
$$

更本质的坐标链是：

$$
\boxed{
x
\xrightarrow{V^T}
[x]_V
\xrightarrow{\Sigma}
[Ax]_U
\xrightarrow{U}
Ax
}
$$

其中最重要的一句话是：

$$
\boxed{
\text{坐标本身没有方向，基才赋予坐标方向。}
}
$$

因此：

$$
\Sigma V^Tx
$$

并不是“仍然沿 $V$ 方向摆放的几何向量”，而是：

$$
\boxed{
Ax\text{ 在 }U\text{ 基下的坐标}
}
$$

最终乘以 $U$，才得到输出空间中的实际向量 $Ax$。
