# Soft Kill Bi-directional Evolution Structural Optimization

## 文件说明

* SKBESO_main.m: 主函数, 包含灵敏度计算、灵敏度滤波、优化循环等功能;
* InpGenerater.py: 输出模型中的节点和单元信息供主程序灵敏度滤波使用,并根据主程序生成的设计变量文件生成下一循环的有限元模型 *.inp文件;
* OdbReader.py: Abaqus-python 2.x 脚本文件, 用于从*.odb文件中提取单元应变能密度到文件, 用于主程序中灵敏度计算;
* PostProcessor.py: Abaqus-python 2.x 脚本文件, 用于循环收敛后迭代过程信息的提取;
* HistoryPlot.m: 用于绘制迭代过程曲线;
* Design.inp: 初始设计模型.

## 几点说明

* 主程序中灵敏度过滤矩阵的生成仅需进行一次;
* 迭代计数器取值通过Iternum.dat文件在主程序和子程序之间传递;
* 首次运行时主程序自动生成DesignVariables文件夹;
* 首次运行主程序自动生成DV_Iter0.dat文件和IterNumw.dat文件;
* 模型全体节点和单元信息仅在第一次运行时生成;