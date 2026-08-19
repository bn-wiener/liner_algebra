可以。SVD 与秩、四个基本子空间之间的关系，是 SVD 最值得掌握的部分之一。因为一旦理解这里，你会发现：

$$
\boxed{
\text{SVD 几乎把一个矩阵的全部空间结构都直接展示出来了}
}
$$

核心公式先放在这里。设

$$
A\in\mathbb R^{m\times n}
$$

秩为

$$
\operatorname{rank}(A)=r
$$

做完整 SVD：

$$
\boxed{
A=U\Sigma V^T
}
$$

那么：

$$
\boxed{
\begin{aligned}
C(A)&=\operatorname{span}(u_1,\ldots,u_r)\
N(A^T)&=\operatorname{span}(u_{r+1},\ldots,u_m)\
C(A^T)&=\operatorname{span}(v_1,\ldots,v_r)\
N(A)&=\operatorname{span}(v_{r+1},\ldots,v_n)
\end{aligned}
}
$$

而且：

$$
\boxed{
r=\text{非零奇异值的个数}
}
$$

下面把这几件事从几何到代数完整串起来。

---

# 一、先从 SVD 的几何图景开始

设：

$$ {A:\mathbb R^n\rightarrow\mathbb R^m} $$

SVD：

$$
A=U\Sigma V^T
$$

真正应该理解成：

$$
\boxed{x\xrightarrow{V^T}x_V
    \xrightarrow{\Sigma}
    Ax_U
    \xrightarrow{U}
    Ax
}
$$

其中：

* (V)：输入空间 $$(\mathbb R^n)$$ 的一组标准正交基；
* (U)：输出空间 $$(\mathbb R^m)$$ 的一组标准正交基；
* ($$\Sigma$$)：告诉我们每一个输入方向到底被放大多少。

最关键的关系是：

$$
\boxed{
Av_i=\sigma_i u_i
}
$$

所以对于每一个右奇异向量 (v_i)：

$$
v_i
\xrightarrow{A}
\sigma_i u_i
$$

这里马上出现两种完全不同的情况。

如果：

$$
\sigma_i>0
$$

那么：

$$
v_i
\xrightarrow{A}
\sigma_i u_i\neq0
$$

这个方向被 (A) 保留下来了。

如果：

$$
\sigma_i=0
$$

那么：

$$
Av_i=0
$$

这个方向被 (A) **彻底压扁了**。

于是：

$$
\boxed{
\text{秩其实就是没有被 }A\text{ 压扁的独立方向数}
}
$$

这就是 SVD 与秩关系的几何本质。

---

# 二、为什么秩等于非零奇异值个数？

设：

$$
A=U\Sigma V^T
$$

因为 (U,V) 都是正交矩阵，所以：

$$
U^{-1}=U^T,\qquad V^{-1}=V^T
$$

它们都是可逆矩阵。

而左乘、右乘可逆矩阵不会改变矩阵的秩：

$$
\operatorname{rank}(PAQ)=\operatorname{rank}(A)
$$

只要 (P,Q) 可逆。

所以：

$$
\operatorname{rank}(A)
=

\operatorname{rank}(U\Sigma V^T)
$$

得到：

$$
\boxed{
\operatorname{rank}(A)
=

\operatorname{rank}(\Sigma)
}
$$

而：

$$
\Sigma=
\begin{bmatrix}
\sigma_1&&&&\
&\sigma_2&&&\
&&\ddots&&\
&&&\sigma_r&\
&&&&0
\end{bmatrix}
$$

其中：

$$
\sigma_1\geq\sigma_2\geq\cdots\geq\sigma_r>0
$$

后面全部是 (0)。

所以：

$$
\boxed{
\operatorname{rank}(\Sigma)=r
}
$$

因此：

$$
\boxed{
\operatorname{rank}(A)
======================

#{i:\sigma_i>0}
}
$$

也就是：

> **矩阵的秩 = 非零奇异值的个数。**

---

# 三、从空间维数理解这个结论

假设：

$$
A:\mathbb R^n\rightarrow\mathbb R^m
$$

输入空间有 (n) 个正交方向：

$$
v_1,\ldots,v_n
$$

如果只有前 (r) 个奇异值非零：

$$
\sigma_1,\ldots,\sigma_r>0
$$

那么：

$$
Av_1=\sigma_1u_1
$$

$$
\cdots
$$

$$
Av_r=\sigma_ru_r
$$

而后面的：

$$
Av_{r+1}=0
$$

$$
\cdots
$$

$$
Av_n=0
$$

所以真正能够产生输出的只有：

$$
v_1,\ldots,v_r
$$

共 (r) 个独立方向。

输出因此只能落在：

$$
\operatorname{span}(u_1,\ldots,u_r)
$$

这个 (r) 维空间里。

所以：

$$
\dim C(A)=r
$$

也就是：

$$
\boxed{
\operatorname{rank}(A)=r
}
$$

---

# 四、先回顾四个基本子空间

对于：

$$
A\in\mathbb R^{m\times n}
$$

有四个最基本的子空间。

## 1. 列空间

$$
C(A)
$$

它属于：

$$
\mathbb R^m
$$

代表：

$$
\boxed{
A\text{ 能够产生哪些输出}
}
$$

也就是：

$$
C(A)={Ax:x\in\mathbb R^n}
$$

---

## 2. 零空间

$$
N(A)
$$

它属于：

$$
\mathbb R^n
$$

代表：

$$
\boxed{
哪些输入经过 }A\text{ 后完全消失}
}
$$

定义：

$$
N(A)={x:Ax=0}
$$

---

## 3. 行空间

通常写成：

$$
C(A^T)
$$

它属于：

$$
\mathbb R^n
$$

可以理解为：

$$
\boxed{
输入空间中真正能够被 }A\text{ 感知的方向}
}
$$

---

## 4. 左零空间

$$
N(A^T)
$$

它属于：

$$
\mathbb R^m
$$

代表输出空间中：

$$
\boxed{
A\text{ 永远不可能产生的方向}
}
$$

---

# 五、SVD 怎么把四个空间全部找出来？

设：

$$
A=U\Sigma V^T
$$

而：

$$
U=
\begin{bmatrix}
u_1&\cdots&u_r&
u_{r+1}&\cdots&u_m
\end{bmatrix}
$$

写成：

$$
\boxed{
U=
\begin{bmatrix}
U_r&U_0
\end{bmatrix}
}
$$

其中：

$$
U_r=
\begin{bmatrix}
u_1&\cdots&u_r
\end{bmatrix}
$$

而：

$$
U_0=
\begin{bmatrix}
u_{r+1}&\cdots&u_m
\end{bmatrix}
$$

同理：

$$
V=
\begin{bmatrix}
v_1&\cdots&v_r&
v_{r+1}&\cdots&v_n
\end{bmatrix}
$$

写成：

$$
\boxed{
V=
\begin{bmatrix}
V_r&V_0
\end{bmatrix}
}
$$

那么四个子空间就是：

$$
\boxed{
C(A)=C(U_r)
}
$$

$$
\boxed{
N(A^T)=C(U_0)
}
$$

$$
\boxed{
C(A^T)=C(V_r)
}
$$

$$
\boxed{
N(A)=C(V_0)
}
$$

这四条非常重要。

---

# 六、为什么非零奇异值对应的 (u_i) 张成列空间？

从：

$$
Av_i=\sigma_i u_i
$$

开始。

如果：

$$
\sigma_i>0
$$

那么：

$$
u_i
===

\frac1{\sigma_i}Av_i
$$

因为：

$$
Av_i\in C(A)
$$

所以：

$$
u_i\in C(A)
$$

因此：

$$
u_1,\ldots,u_r
$$

全部属于列空间。

另一方面，任意输入：

$$
x
=

\sum_{i=1}^n c_i v_i
$$

那么：

$$
Ax
==

\sum_{i=1}^n c_i Av_i
$$

因为后面的零奇异值对应：

$$
Av_i=0
$$

所以：

$$
Ax
==

\sum_{i=1}^r
c_i\sigma_i u_i
$$

因此任何输出 (Ax) 都一定落在：

$$
\operatorname{span}(u_1,\ldots,u_r)
$$

所以：

$$
C(A)
\subseteq
\operatorname{span}(u_1,\ldots,u_r)
$$

前面我们又知道：

$$
u_1,\ldots,u_r\in C(A)
$$

于是：

$$
\boxed{
C(A)
====

\operatorname{span}(u_1,\ldots,u_r)
}
$$

---

# 七、列空间的几何意义

这一条非常重要：

$$
\boxed{
C(A)=\operatorname{span}(u_1,\ldots,u_r)
}
$$

说明：

> 所有可能的输出 (Ax)，都只能沿着前 (r) 个左奇异向量方向组合。

也就是说：

```text
输入空间
R^n

v1 ----σ1----> u1
v2 ----σ2----> u2
...
vr ----σr----> ur

v(r+1) --> 0
...
vn      --> 0
```

所以最终：

$$
Ax
$$

只可能出现在：

$$
\operatorname{span}(u_1,\ldots,u_r)
$$

里面。

这就是列空间。

---

# 八、为什么零奇异值对应的 (v_i) 张成零空间？

这是最直接的一条。

因为：

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
v_{r+1},\ldots,v_n
$$

全部属于零空间。

而任意：

$$
x\in N(A)
$$

都可以展开成：

$$
x=\sum_{i=1}^n c_i v_i
$$

因为：

$$
Ax=0
$$

所以：

$$
\sum_{i=1}^r
c_i\sigma_i u_i=0
$$

由于：

$$
u_1,\ldots,u_r
$$

线性无关，而且：

$$
\sigma_i>0
$$

只能有：

$$
c_1=\cdots=c_r=0
$$

所以：

$$
x=
\sum_{i=r+1}^n c_i v_i
$$

于是：

$$
\boxed{
N(A)
====

\operatorname{span}(v_{r+1},\ldots,v_n)
}
$$

---

# 九、零空间的几何意义

这是一个非常漂亮的结论：

$$
\boxed{
N(A)
====

\text{被 }A\text{ 完全压扁的输入方向}
}
$$

比如：

$$
\sigma_3=0
$$

那么：

$$
Av_3=0
$$

意味着：

$$
v_3
$$

这个方向上的信息经过 (A) 后彻底消失。

所以：

$$
\dim N(A)=n-r
$$

这马上得到秩-零化度定理：

$$
\boxed{
\operatorname{rank}(A)+\operatorname{nullity}(A)=n
}
$$

即：

$$
\boxed{
r+(n-r)=n
}
$$

SVD 让这个定理几乎一眼就能看出来。

---

# 十、为什么前 (r) 个 (v_i) 张成行空间？

我们知道：

$$
A^T
===

V\Sigma^TU^T
$$

因此 (A^T) 的左奇异向量就是 (V)。

同样根据刚才“列空间由非零奇异值对应的左奇异向量张成”这个结论：

$$
C(A^T)
======

\operatorname{span}(v_1,\ldots,v_r)
$$

因此：

$$
\boxed{
C(A^T)
======

\operatorname{span}(v_1,\ldots,v_r)
}
$$

而：

$$
C(A^T)
$$

就是 (A) 的行空间。

所以：

$$
\boxed{
\text{行空间 = 有效右奇异向量张成的空间}
}
$$

---

# 十一、行空间真正的几何意义

这一点非常适合和零空间一起理解。

输入空间：

$$
\mathbb R^n
$$

被 (V) 分成两部分：

$$
\boxed{
\mathbb R^n
===========

C(A^T)\oplus N(A)
}
$$

也就是：

$$
\boxed{
\mathbb R^n
===========

\operatorname{span}(v_1,\ldots,v_r)
\oplus
\operatorname{span}(v_{r+1},\ldots,v_n)
}
$$

前一部分：

$$
C(A^T)
$$

是真正会影响输出的输入方向。

后一部分：

$$
N(A)
$$

是输入进去也没有任何作用的方向。

所以你可以形成这样一个心智模型：

```text
输入空间 R^n

┌─────────────────────────────┐
│                             │
│ 行空间 C(Aᵀ)                │
│ v1,...,vr                   │
│                             │
│ 真正影响 Ax 的方向          │
│                             │
├─────────────────────────────┤
│                             │
│ 零空间 N(A)                 │
│ v(r+1),...,vn               │
│                             │
│ 被 A 完全忽略的方向         │
│                             │
└─────────────────────────────┘
```

而且这两个空间正交：

$$
\boxed{
C(A^T)\perp N(A)
}
$$

---

# 十二、为什么 (u_{r+1},\ldots,u_m) 张成左零空间？

左零空间是：

$$
N(A^T)
$$

即满足：

$$
A^Ty=0
$$

的所有 (y)。

从：

$$
A^T
===

V\Sigma^TU^T
$$

来看，如果：

$$
y=u_i
$$

其中：

$$
i>r
$$

那么：

$$
U^Tu_i=e_i
$$

而 (\Sigma^T) 对应位置是零，所以：

$$
A^Tu_i=0
$$

因此：

$$
u_i\in N(A^T)
$$

于是：

$$
\boxed{
N(A^T)
======

\operatorname{span}(u_{r+1},\ldots,u_m)
}
$$

---

# 十三、左零空间的几何意义

输出空间：

$$
\mathbb R^m
$$

也被分成两部分：

$$
\boxed{
\mathbb R^m
===========

C(A)\oplus N(A^T)
}
$$

即：

$$
\boxed{
\mathbb R^m
===========

\operatorname{span}(u_1,\ldots,u_r)
\oplus
\operatorname{span}(u_{r+1},\ldots,u_m)
}
$$

其中：

### 列空间

$$
C(A)
$$

是：

> (A) 真正能够产生的输出方向。

而：

### 左零空间

$$
N(A^T)
$$

是：

> (A) 永远无法产生的输出方向。

所以：

```text
输出空间 R^m

┌──────────────────────────────┐
│                              │
│ 列空间 C(A)                  │
│ u1,...,ur                    │
│                              │
│ A 能产生的方向               │
│                              │
├──────────────────────────────┤
│                              │
│ 左零空间 N(Aᵀ)               │
│ u(r+1),...,um                │
│                              │
│ A 永远无法产生的方向         │
│                              │
└──────────────────────────────┘
```

而且：

$$
\boxed{
C(A)\perp N(A^T)
}
$$

---

# 十四、现在把四个空间放到一张图里

这是整部分最重要的图。

设：

$$
A:\mathbb R^n\rightarrow\mathbb R^m
$$

秩：

$$
r
$$

那么：

```text
            输入空间 R^n                        输出空间 R^m

      ┌─────────────────────┐              ┌─────────────────────┐
      │                     │              │                     │
      │   行空间 C(Aᵀ)      │              │    列空间 C(A)      │
      │                     │              │                     │
      │   v1,...,vr         │────── A ───>│   u1,...,ur         │
      │                     │              │                     │
      │   有效输入方向      │              │   可达输出方向      │
      │                     │              │                     │
      └─────────────────────┘              └─────────────────────┘
               │                                    │
               │ 正交                               │ 正交
               ▼                                    ▼
      ┌─────────────────────┐              ┌─────────────────────┐
      │                     │              │                     │
      │    零空间 N(A)      │              │ 左零空间 N(Aᵀ)     │
      │                     │              │                     │
      │ v(r+1),...,vn       │              │ u(r+1),...,um       │
      │                     │              │                     │
      │ 被 A 压成 0         │              │ A 永远到不了        │
      │                     │              │                     │
      └─────────────────────┘              └─────────────────────┘
```

核心对应：

$$
\boxed{
v_i
\xrightarrow{A}
\sigma_i u_i
}
$$

对于：

$$
i\le r
$$

有：

$$
\sigma_i>0
$$

所以：

$$
v_i\rightarrow u_i
$$

对于：

$$
i>r
$$

有：

$$
\sigma_i=0
$$

所以：

$$
v_i\rightarrow0
$$

---

# 十五、完整写成空间分解

SVD 实际上给出了输入空间的正交分解：

$$
\boxed{
\mathbb R^n
===========

C(A^T)\oplus N(A)
}
$$

也就是：

$$
\boxed{
\mathbb R^n
===========

\operatorname{span}(V_r)
\oplus
\operatorname{span}(V_0)
}
$$

输出空间：

$$
\boxed{
\mathbb R^m
===========

C(A)\oplus N(A^T)
}
$$

也就是：

$$
\boxed{
\mathbb R^m
===========

\operatorname{span}(U_r)
\oplus
\operatorname{span}(U_0)
}
$$

因此 SVD 可以写成非常漂亮的结构：

$$
\boxed{
\underbrace{
C(A^T)
}*{\text{有效输入}}
\xrightarrow{A}
\underbrace{
C(A)
}*{\text{有效输出}}
}
$$

而：

$$
\boxed{
N(A)
\xrightarrow{A}
{0}
}
$$

至于：

$$
N(A^T)
$$

则根本不在 (A) 的输出范围内。

---

# 十六、用一个完整数值例子理解四个空间

考虑：

$$
A=
\begin{bmatrix}
1&1\
2&2\
3&3
\end{bmatrix}
$$

这是：

$$
A:\mathbb R^2\rightarrow\mathbb R^3
$$

很明显两列成比例：

$$
\begin{bmatrix}
1\2\3
\end{bmatrix},
\qquad
\begin{bmatrix}
1\2\3
\end{bmatrix}
$$

所以：

$$
\boxed{
\operatorname{rank}(A)=1
}
$$

---

# 十七、求这个矩阵的 SVD

先计算：

$$
A^TA
====

\begin{bmatrix}
14&14\
14&14
\end{bmatrix}
$$

特征值：

$$
\lambda_1=28
$$

$$
\lambda_2=0
$$

所以奇异值：

$$
\boxed{
\sigma_1=\sqrt{28}=2\sqrt7
}
$$

$$
\boxed{
\sigma_2=0
}
$$

因此：

$$
\operatorname{rank}(A)=1
$$

再次验证：

$$
\boxed{
\text{非零奇异值只有一个}
}
$$

---

# 十八、右奇异向量

对于：

$$
\lambda_1=28
$$

对应：

$$
\boxed{
v_1=
\frac1{\sqrt2}
\begin{bmatrix}
1\
1
\end{bmatrix}
}
$$

对于：

$$
\lambda_2=0
$$

对应：

$$
\boxed{
v_2=
\frac1{\sqrt2}
\begin{bmatrix}
1\
-1
\end{bmatrix}
}
$$

---

# 十九、直接看这两个输入方向

对于 (v_1)：

$$
Av_1
====

\frac1{\sqrt2}
\begin{bmatrix}
2\
4\
6
\end{bmatrix}
$$

即：

$$
Av_1
====

\sqrt2
\begin{bmatrix}
1\
2\
3
\end{bmatrix}
$$

其长度：

$$
|Av_1|
======

# \sqrt{28}

2\sqrt7
$$

所以：

$$
u_1
===

\frac1{\sqrt{14}}
\begin{bmatrix}
1\
2\
3
\end{bmatrix}
$$

得到：

$$
\boxed{
Av_1=2\sqrt7,u_1
}
$$

而对于：

$$
v_2=
\frac1{\sqrt2}
\begin{bmatrix}
1\
-1
\end{bmatrix}
$$

有：

$$
Av_2=0
$$

所以：

$$
\boxed{
v_2\in N(A)
}
$$

---

# 二十、这个例子的零空间

因此：

$$
\boxed{
N(A)
====

\operatorname{span}
\left(
\frac1{\sqrt2}
\begin{bmatrix}
1\
-1
\end{bmatrix}
\right)
}
$$

几何上：

输入空间 (\mathbb R^2) 中：

$$
v_1
$$

是有效方向。

而：

$$
v_2
$$

是完全被 (A) 压扁的方向。

所以：

```text
R² 输入空间

             v1
             /
            /
-----------•-----------
          /
         /
        v2

v1：保留
v2：压成 0
```

---

# 二十一、这个例子的行空间

因为：

$$
r=1
$$

所以：

$$
\boxed{
C(A^T)
======

\operatorname{span}(v_1)
}
$$

即：

$$
\boxed{
C(A^T)
======

\operatorname{span}
\left(
\begin{bmatrix}
1\
1
\end{bmatrix}
\right)
}
$$

它代表：

> 输入空间中真正对输出有影响的方向。

---

# 二十二、这个例子的列空间

非零奇异值对应：

$$
u_1=
\frac1{\sqrt{14}}
\begin{bmatrix}
1\
2\
3
\end{bmatrix}
$$

所以：

$$
\boxed{
C(A)
====

\operatorname{span}
\left(
\begin{bmatrix}
1\
2\
3
\end{bmatrix}
\right)
}
$$

因为 (A) 无论输入什么：

$$
x=
\begin{bmatrix}
x_1\
x_2
\end{bmatrix}
$$

输出：

$$
Ax
==

\begin{bmatrix}
x_1+x_2\
2x_1+2x_2\
3x_1+3x_2
\end{bmatrix}
$$

可以写成：

$$
Ax
==

(x_1+x_2)
\begin{bmatrix}
1\
2\
3
\end{bmatrix}
$$

所以所有输出都只能沿着：

$$
\begin{bmatrix}
1\
2\
3
\end{bmatrix}
$$

这个方向。

---

# 二十三、这个例子的左零空间

输出空间是：

$$
\mathbb R^3
$$

而列空间只有 (1) 维。

所以：

$$
\dim N(A^T)=3-1=2
$$

我们需要找两个与：

$$
u_1=
\frac1{\sqrt{14}}
\begin{bmatrix}
1\
2\
3
\end{bmatrix}
$$

正交的方向。

例如可以选：

$$
u_2
===

\frac1{\sqrt5}
\begin{bmatrix}
2\
-1\
0
\end{bmatrix}
$$

因为：

$$
1\cdot2+2\cdot(-1)+3\cdot0=0
$$

再选：

$$
u_3
===

\frac1{\sqrt{70}}
\begin{bmatrix}
3\
6\
-5
\end{bmatrix}
$$

也有：

$$
1\cdot3+2\cdot6+3\cdot(-5)
==========================

# 3+12-15

0
$$

并且 (u_2,u_3) 彼此正交。

所以：

$$
\boxed{
N(A^T)=\operatorname{span}(u_2,u_3)
}
$$

这两个方向就是：

> 输出空间中 (A) 永远到不了的方向。

---

# 二十四、把这个例子的四个空间完整列出来

对于：

$$
A=
\begin{bmatrix}
1&1\
2&2\
3&3
\end{bmatrix}
$$

有：

$$
\operatorname{rank}(A)=1
$$

---

### 输入空间 (\mathbb R^2)

有效部分：

$$
\boxed{
C(A^T)
======

\operatorname{span}
\left(
\begin{bmatrix}
1\
1
\end{bmatrix}
\right)
}
$$

失效部分：

$$
\boxed{
N(A)
====

\operatorname{span}
\left(
\begin{bmatrix}
1\
-1
\end{bmatrix}
\right)
}
$$

所以：

$$
\boxed{
\mathbb R^2=C(A^T)\oplus N(A)
}
$$

---

### 输出空间 (\mathbb R^3)

可达部分：

$$
\boxed{
C(A)
====

\operatorname{span}
\left(
\begin{bmatrix}
1\
2\
3
\end{bmatrix}
\right)
}
$$

不可达部分：

$$
\boxed{
N(A^T)
======

C(A)^\perp
}
$$

维数：

$$
2
$$

所以：

$$
\boxed{
\mathbb R^3=C(A)\oplus N(A^T)
}
$$

---

# 二十五、SVD 与 (Ax=b) 的关系也出来了

现在考虑：

$$
Ax=b
$$

什么时候有解？

必须满足：

$$
b\in C(A)
$$

而：

$$
C(A)
====

\operatorname{span}(u_1,\ldots,u_r)
$$

所以：

$$
b
$$

不能包含任何左零空间方向。

把：

$$
b
$$

展开成 (U) 基：

$$
b
=

\sum_{i=1}^m c_i u_i
$$

要有精确解，就必须：

$$
\boxed{
c_{r+1}=\cdots=c_m=0
}
$$

因为这些：

$$
u_{r+1},\ldots,u_m
$$

属于：

$$
N(A^T)
$$

是 (A) 永远无法产生的方向。

因此：

$$
\boxed{
Ax=b\text{ 有解}
\iff
b\perp N(A^T)
}
$$

等价于：

$$
\boxed{
b\in C(A)
}
$$

---

# 二十六、SVD 与解的唯一性

另外，解是否唯一取决于：

$$
N(A)
$$

如果：

$$
N(A)={0}
$$

那么：

$$
Ax=b
$$

至多一个解。

如果：

$$
N(A)\neq{0}
$$

那么假设：

$$
x_0
$$

是一个解：

$$
Ax_0=b
$$

对于任意：

$$
z\in N(A)
$$

都有：

$$
Az=0
$$

所以：

$$
A(x_0+z)
========

# Ax_0+Az

b
$$

因此：

$$
x_0+z
$$

也是解。

所以：

$$
\boxed{
N(A)\text{ 描述了解的不唯一性}
}
$$

这也是四个子空间和 (Ax=b) 的重要联系。

---

# 二十七、最小二乘为什么和左零空间有关？

如果：

$$
b\notin C(A)
$$

那么：

$$
Ax=b
$$

没有精确解。

但因为：

$$
\mathbb R^m
===========

C(A)\oplus N(A^T)
$$

所以可以把：

$$
b
$$

唯一分解成：

$$
\boxed{
b=b_C+b_N
}
$$

其中：

$$
b_C\in C(A)
$$

$$
b_N\in N(A^T)
$$

由于 (A) 永远产生不了：

$$
b_N
$$

所以最小二乘实际上是在寻找：

$$
Ax=b_C
$$

也就是把：

$$
b
$$

投影到：

$$
C(A)
$$

上。

因此残差：

$$
r=b-A\hat x
$$

满足：

$$
\boxed{
r\in N(A^T)
}
$$

也就是：

$$
\boxed{
A^Tr=0
}
$$

这正是正规方程：

$$
A^T(A\hat x-b)=0
$$

的几何意义。

---

# 二十八、SVD 为什么把这些关系看得特别清楚？

因为如果：

$$
A=U\Sigma V^T
$$

那么：

$$
U=
\begin{bmatrix}
U_r&U_0
\end{bmatrix}
$$

实际上直接给了输出空间：

$$
\boxed{
\mathbb R^m
===========

\underbrace{C(U_r)}*{C(A)}
\oplus
\underbrace{C(U_0)}*{N(A^T)}
}
$$

而：

$$
V=
\begin{bmatrix}
V_r&V_0
\end{bmatrix}
$$

直接给了输入空间：

$$
\boxed{
\mathbb R^n
===========

\underbrace{C(V_r)}*{C(A^T)}
\oplus
\underbrace{C(V_0)}*{N(A)}
}
$$

所以你甚至可以把 SVD 理解成：

$$
\boxed{
\text{SVD 在给输入空间和输出空间做“结构化拆分”}
}
$$

---

# 二十九、四个空间、秩、维数放在一起

设：

$$
A\in\mathbb R^{m\times n}
$$

且：

$$
\operatorname{rank}(A)=r
$$

那么：

| 空间       | 所在空间          | SVD 基                |    维数 | 含义        |
| -------- | ------------- | -------------------- | ----: | --------- |
| (C(A^T)) | (\mathbb R^n) | (v_1,\ldots,v_r)     |   (r) | 有效输入方向    |
| (N(A))   | (\mathbb R^n) | (v_{r+1},\ldots,v_n) | (n-r) | 被压成 0 的输入 |
| (C(A))   | (\mathbb R^m) | (u_1,\ldots,u_r)     |   (r) | 可达输出方向    |
| (N(A^T)) | (\mathbb R^m) | (u_{r+1},\ldots,u_m) | (m-r) | 永远不可达输出方向 |

所以：

$$
\boxed{
\dim C(A)=\dim C(A^T)=r
}
$$

这也顺便解释了为什么：

$$
\boxed{
\text{列秩 = 行秩}
}
$$

因为它们都等于：

$$
\boxed{
\text{非零奇异值的个数}
}
$$

---

# 三十、把所有关系压缩成一张心智模型

建议你以后看到：

$$
A=U\Sigma V^T
$$

脑中直接出现：

```text
                         A : R^n → R^m


       输入空间 R^n                           输出空间 R^m

 ┌──────────────────────┐              ┌──────────────────────┐
 │                      │              │                      │
 │ 行空间 C(Aᵀ)         │              │ 列空间 C(A)          │
 │                      │              │                      │
 │ v₁,...,vᵣ            │──────A─────>│ u₁,...,uᵣ            │
 │                      │              │                      │
 │ σ₁,...,σᵣ > 0        │              │ 可达输出             │
 │                      │              │                      │
 └──────────────────────┘              └──────────────────────┘
          ⊥                                      ⊥
          │                                      │
 ┌──────────────────────┐              ┌──────────────────────┐
 │                      │              │                      │
 │ 零空间 N(A)          │              │ 左零空间 N(Aᵀ)      │
 │                      │              │                      │
 │ vᵣ₊₁,...,vₙ          │              │ uᵣ₊₁,...,uₘ          │
 │                      │              │                      │
 │ σ = 0                │              │ A 永远到不了         │
 │                      │              │                      │
 └──────────────────────┘              └──────────────────────┘
          │
          └────────────── A ──────────────> 0
```

最核心的关系就是：

$$
\boxed{
Av_i=
\begin{cases}
\sigma_i u_i,& i\le r\
0,&i>r
\end{cases}
}
$$

于是：

$$
\boxed{
\begin{aligned}
\text{非零 }\sigma_i
&\Longleftrightarrow
\text{有效方向}\
&\Longleftrightarrow
\text{秩}\
&\Longleftrightarrow
C(A^T)\leftrightarrow C(A)
\end{aligned}
}
$$

而：

$$
\boxed{
\sigma_i=0
\Longleftrightarrow
v_i\in N(A)
}
$$

最终建议记住这一句话：

$$
\boxed{
\text{SVD 把输入空间分成“有效输入 + 无效输入”，}
}
$$

$$
\boxed{
\text{把输出空间分成“可达输出 + 不可达输出”。}
}
$$

其中：

$$
\boxed{
\underbrace{C(A^T)}*{\text{有效输入}}
\xrightarrow{A}
\underbrace{C(A)}*{\text{可达输出}}
}
$$

而：

$$
\boxed{
N(A)\xrightarrow{A}0
}
$$

并且：

$$
\boxed{
N(A^T)=C(A)^\perp
}
$$

这其实就是 **秩、四个基本子空间、(Ax=b)、最小二乘、伪逆** 能被 SVD 统一起来的根本原因。
