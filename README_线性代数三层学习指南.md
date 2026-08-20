# 线性代数三层学习体系——总学习指南

## 1. 这套资料解决什么问题

这套资料不按传统教材章节机械堆叠概念，而围绕一条主线组织线性代数：

> **基础层解决“怎么算”，空间层解决“有什么”，变换层解决“怎么变”。**

三层不是彼此割裂的三门课。它们观察的是同一个对象，只是问题不同。

给定矩阵 $A$ 和向量 $\boldsymbol{x}$：

- 基础层研究怎样计算 $A\boldsymbol{x}$、怎样解 $A\boldsymbol{x}=\boldsymbol{b}$；
- 空间层研究哪些输入方向存在、哪些输出可达、哪些方向被消除；
- 变换层研究 $\boldsymbol{x}\mapsto A\boldsymbol{x}$ 怎样改变方向、尺度、维数和信息。

最终整套知识在矩阵分解，尤其是 SVD 中重新汇合。

## 2. 推荐学习顺序

### 第一阶段：先获得计算语言

进入：[`01_基础层_会表示会计算/`](01_基础层_会表示会计算/README_基础层学习指南.md)

主线：

$$
\text{向量}
\rightarrow
\text{矩阵}
\rightarrow
A\boldsymbol{x}
\rightarrow
A\boldsymbol{x}=\boldsymbol{b}
\rightarrow
\text{高斯消元}
\rightarrow
\text{主元/自由变量}
\rightarrow
\text{解集结构}.
$$

目标不是背概念，而是看到一个线性问题时能够表示、计算并读懂解。

### 第二阶段：把计算结果看成空间结构

进入：[`02_空间层_有什么/`](02_空间层_有什么/README_空间层学习指南.md)

主线：

$$
\text{线性组合}
\rightarrow
\operatorname{span}
\rightarrow
\text{线性无关}
\rightarrow
\text{基}
\rightarrow
\text{维数}
\rightarrow
\text{四个基本子空间}
\rightarrow
\text{正交}
\rightarrow
\text{投影}
\rightarrow
\text{最小二乘}.
$$

这一层始终回答三个问题：

1. 哪些向量可以互相表示？
2. 有多少个独立方向？
3. 这些方向形成什么空间？

### 第三阶段：把矩阵看成作用于空间的变换

进入：[`03_变换层_怎么变/`](03_变换层_怎么变/README_变换层学习指南.md)

主线：

$$
T(\boldsymbol{x})=A\boldsymbol{x}
\rightarrow
\text{复合/可逆/换基}
\rightarrow
\ker T,\operatorname{Im}T
\rightarrow
\text{特征结构}
\rightarrow
\text{对角化}
\rightarrow
\text{QR/SVD}
\rightarrow
\text{伪逆/低秩/PCA}.
$$

第三层的关键不是增加更多算法，而是回答：

> 怎样选择合适的坐标和方向，把复杂的线性变换拆成简单的独立变化？

## 3. 学完三层后怎样回收

完成三层后，不建议立刻重新从 B1 开始逐页阅读。先进入：

- [`00_全体系导航/01_三层统一知识地图.md`](00_全体系导航/01_三层统一知识地图.md)：看全局链条；
- [`00_全体系导航/02_跨层桥梁总册.md`](00_全体系导航/02_跨层桥梁总册.md)：把同一概念的计算、空间、变换三种语言对齐；
- [`04_全课程综合例题/`](04_全课程综合例题/)：用同一个矩阵跨层分析；
- [`05_全课程综合练习/`](05_全课程综合练习/)：检查是否真的能跨层调用知识；
- [`06_附录/`](06_附录/)：用于查符号、公式、概念差异和矩阵分解。

## 4. 一个统一的观察模板

以后遇到任何矩阵 $A\in\mathbb{R}^{m\times n}$，可以依次问下面的问题。

### 计算问题

- $A\boldsymbol{x}$ 怎样计算？
- $A\boldsymbol{x}=\boldsymbol{b}$ 是否有解？
- 解是否唯一？
- RREF 中有多少主元和自由变量？

### 空间问题

- $\operatorname{Col}(A)$ 是什么？
- $\operatorname{Null}(A)$ 是什么？
- $\operatorname{rank}(A)$ 与 $\operatorname{nullity}(A)$ 是多少？
- 四个基本子空间怎样组织输入空间和输出空间？
- 是否需要正交、投影和最小二乘？

### 变换问题

- $A$ 把基向量送到哪里？
- 哪些方向被压到 $0$？
- 哪些输出能够到达？
- 变换是否可逆？
- 有没有只缩放、不改变方向的特征方向？
- 能否换到更好的基，使矩阵变简单？
- SVD 怎样把变换拆成正交方向上的独立缩放？

## 5. 什么时候算真正掌握

不是会背：

$$
\operatorname{rank}(A)+\operatorname{nullity}(A)=n
$$

就算掌握秩—零度，也不是会算：

$$
A=PDP^{-1}
$$

就算掌握对角化。

真正掌握要求能够解释这些公式为什么自然产生、它们解决什么问题、怎样和前面的概念相连。例如：

- 自由变量为什么自然变成零空间中的自由方向；
- $\boldsymbol{b}\notin\operatorname{Col}(A)$ 为什么导致无精确解；
- 投影为什么会产生最小二乘；
- 特征向量为什么可以构成一种“适合矩阵的基”；
- SVD 为什么能同时回收 rank、四个基本子空间、正交和最小二乘。

## 6. 文件结构说明

正常学习只需要阅读：

```text
README_线性代数三层学习指南.md
00_全体系导航/
01_基础层_会表示会计算/
02_空间层_有什么/
03_变换层_怎么变/
04_全课程综合例题/
05_全课程综合练习/
06_附录/
```

`_audit/` 和 `_global_audit/` 是重构与完整性审计材料，不参与日常学习。

## 7. 最终主线

整套资料最终希望建立下面这一个认知：

$$
\boxed{
\text{线性关系}
\rightarrow
\text{方程与计算}
\rightarrow
\text{空间与方向}
\rightarrow
\text{线性变换}
\rightarrow
\text{寻找合适的基}
\rightarrow
\text{把复杂变化拆成简单方向变化}
}
$$
