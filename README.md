# Cpp-DG
练习使用C++写间断有限元方法(1d, 2d rectangular mesh) 

框架和一些写法基本上是参考G老师给的。预计包含： 

    1. scalar 和 system的基本求解器
    2. 常用Limiter(TVB limiter, WENO limiter(探测器方面：KXRCF和minmod), PP limiter)
    3. OEDG, OFDG
    4. 测试案例:(常见的Accuracy test, Blast wave Problem, Shu-Osher Problem, Double Mach Reflection Problem, Shock-Reflection Problem....)

# 1D 数据结构说明

*注：* cell 代表单元指标；nvar 代表方程个数；ndof 代表基函数的个数；nq 代表所用高斯积分点的个数

    1. 每个单元的多项式系数存储于std::vector<double>(ncells * nvar * ndof), 因此u(cell, var, mode)对应u[(cell * nvar + var)*ndof + mode]. 循环结构即为cell -> var -> dof。
    2. 基函数在高斯点处的值存储于std::vector<double>(nq * ndof). 因此 phi(q, k) 对应phi[q*ndof + k]. 循环结构即为 q -> k.
    3. 每个界面的数值通量存储于std::vector<double>(nfaces * nvar). 因此 flux(face, var) 对应 flux[face*nvar + var]. 循环结构即为 face -> var.
    4. 每个单元的高斯点处的通量函数值存储于std::vector<double>(ncells * nq * nvar). 因此f(cell, q, var)对应 f[(cell*nq + q)*nvar + var]. 循环结构即为 cell -> q -> var.