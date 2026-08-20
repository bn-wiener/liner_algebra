# $A\boldsymbol{x}=\boldsymbol{b}$ 主线总结

## 1. 正向：给输入，算输出

$$
\boldsymbol{x}\longmapsto A\boldsymbol{x}.
$$

## 2. 逆向：给输出，找输入

$$
A\boldsymbol{x}=\boldsymbol{b}.
$$

通过增广矩阵与高斯消元，我们把“能不能做到、怎样做到、做到的方法是否唯一”统一为计算问题。

## 3. 三种结果

- 出现矛盾行：无解；
- 一致且无自由变量：唯一解；
- 一致且有自由变量：无穷多解。

## 4. 齐次方程为什么不可缺少

如果 $\boldsymbol{x}_1$ 与 $\boldsymbol{x}_2$ 产生同一个 $\boldsymbol{b}$，则

$$
A(\boldsymbol{x}_1-\boldsymbol{x}_2)=\boldsymbol{0}.
$$

因此齐次方程描述了“在不改变输出的前提下，输入还能怎样变化”。

## 5. 全部解

$$
\boldsymbol{x}=\boldsymbol{x}_p+\boldsymbol{x}_h,
\qquad A\boldsymbol{x}_h=\boldsymbol{0}.
$$

第二层将进一步追问：这些自由方向形成什么空间、维数是多少、怎样选一组最简方向。
