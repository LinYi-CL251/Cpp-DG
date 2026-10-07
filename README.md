# Cpp-DG
练习使用C++写间断有限元方法(1d, 2d rectangular mesh) 

框架和一些写法基本上是参考G老师给的。预计包含： 

    1. scalar 和 system的基本求解器
    2. 常用Limiter(TVB limiter, WENO limiter(探测器方面：KXRCF和minmod), PP limiter)
    3. OEDG, OFDG
    4. 测试案例:(常见的Accuracy test, Blast wave Problem, Shu-Osher Problem, Double Mach Reflection Problem, Shock-Reflection Problem....)