# 第三层：变换层——怎么变

## 这一层研究什么

前两层已经分别解决了两个问题：

- 基础层：$A\boldsymbol{x}=\boldsymbol{b}$ **怎么算**；
- 空间层：向量、子空间、基、维数、四个基本子空间 **有什么**。

第三层把矩阵进一步看成线性变换：

$$
T(\boldsymbol{x})=A\boldsymbol{x}.
$$

于是问题从“矩阵里有哪些数”变成：

> 输入向量和整个空间经过 $A$ 后，方向、尺度、维数与信息发生了什么？

## 建议学习顺序

第一组先建立动态视角：

$$
\text{矩阵}\rightarrow\text{线性变换}\rightarrow\text{基向量的像}\rightarrow\text{复合}\rightarrow\text{可逆}\rightarrow\text{换基}\rightarrow\text{核与像}.
$$

第二组寻找变换自身的特殊结构：

$$
\det(A)\rightarrow A\boldsymbol v=\lambda\boldsymbol v\rightarrow\text{特征空间}\rightarrow\text{相似}\rightarrow\text{对角化}\rightarrow\text{谱定理}.
$$

第三组研究重复作用和方向相关的“能量”：

$$
A^k\rightarrow\text{动态系统}\rightarrow\boldsymbol x^{\mathsf T}A\boldsymbol x\rightarrow\text{正定性}.
$$

第四组把整个线性代数重新统一到矩阵分解：

$$
A=QR,
\qquad
A=U\Sigma V^{\mathsf T}.
$$

SVD 最终把空间层的秩、四个基本子空间、正交与最小二乘重新收回到同一个变换图景中。

## 学习方式

每学一个节点都同时问：

1. 它从前面哪个问题产生？
2. 它解决什么变换问题？
3. 怎样计算？
4. 空间上意味着什么？
5. 变换上意味着什么？
6. 它怎样连接前两层？

正文中不放来源编号和重构记录。所有来源、拆分和跨层账目统一放在 `_audit/`。

## 公式格式

所有公式使用 VSCode Markdown 兼容写法：行内公式使用 `$...$`，独立公式使用 `$$...$$`。
