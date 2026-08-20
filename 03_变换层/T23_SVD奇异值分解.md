# T23 SVD 奇异值分解

> **学习定位**：第三层的最终统一分解：任意实矩阵都可在两组标准正交基之间表现为独立非负缩放。

## 本节主线

SVD 写成：

$$
A=U\Sigma V^{\mathsf T}.
$$

它可以理解为：

$$
\text{输入换基}\xrightarrow{V^{\mathsf T}}
\text{独立缩放}\xrightarrow{\Sigma}
\text{输出换基}\xrightarrow{U}.
$$

## 一、先从几何上理解 SVD

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


## 二、三个矩阵分别干什么？

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


## 三、从一个二维几何图开始

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


## 四、这和特征向量有什么区别？

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


## 五、为什么会有两个基 $V$ 和 $U$？

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


## 六、从“基与坐标”理解 SVD

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


## 七、然后 $\Sigma$ 做什么？

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


## 八、最后为什么还要乘 $U$？

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


## 九、SVD 最本质的坐标解释

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


## 十、为什么 SVD 一定存在？

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


## 十四、再看一个带方向变化的例子

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


## 二十九、最终统一到 $T(x)\leftrightarrow Ax=b$

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

## 本节统一

SVD 不要求 $A$ 为方阵，也不要求 $A$ 可对角化。它寻找的是两套空间中彼此配对的正交方向，因此比普通特征分解更普适。

## 下一步为什么自然出现

下一节深入追问：这些右奇异向量、左奇异向量和奇异值到底从哪里来，它们与 rank 和四个基本子空间有什么关系？
