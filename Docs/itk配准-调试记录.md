配准功能的调试记录
=======
### 为什么记录
这是在这个程序的编写中遇到的第二个疑难杂症,比较有记录的价值.<br>

### 问题描述
1.同样的代码,在官方案例(ITK.sln)中打开完全没问题,但是在我们的reader中就是不好使<br>
2.问题出现的点是无法解释的<br>

### 思路
1.尽量在reader中复现所有的官方案例中的环境<br>

### 我都干了什么-结果如何
1.替换工程用的库为最新版的打开Example编译的版本---没有解决,没有变化<br>
2.将所有修改过的代码全部改回原来的样子---没有解决,没有变化<br>
3.考虑观察者的问题,先在官方案例中去掉观察者,运行,没问题,于是在本工程中叶去掉观察者---没有解决,没有变换<br>
4.查看处理器预定义,发现你一条和ITK有关的定义,抄过来.---暴露出问题:发现若干个头文件加载不上,但是并不会爆出编译错误,而当我把宏定义到全局上而不是这个文件的时候,这个错误爆发,报编译错误,通过everything检索,发现找不到的头文件并不在itk的install目录中,而是在src下的Example中,进一步检查发现,工程的ITK环境有许多缺失,例如库目录没写,附加依赖项可能也不完整----策略:先去掉宏定义,补全环境设置,在视情况而定.<br>
ITK_IO_FACTORY_REGISTER_MANAGER<br>
4.补全库目录,没用<br>
5.添加宏定义,编译,出现未定义外部符号和找不到的头文件,检索到头文件的位置,添加上,继续编译<br>
5,检查官方案例的附加依赖项,借助VScode和everything,检索出如下四个以前没有的文件:<br>
itknetlib-4.12-gd.lib<br>
ITKIOMRC-4.12-gd.lib<br>
itkgdcmopenjpeg-4.12-gd.lib<br>
itkgdcmcharls-4.12-gd.lib<br>
推测这是和Example编译选项有关系的,换句话说,这个官方案例里面,使用了一些额外的,如果不编译例子就不会有的特殊东西.<br>


### 结论
1.使用复杂的库的时候,要做到头文件在用处引用,没用的不引用,力争做到每个include的原因都是可追溯的,尽量减少不可控代码.<br>
2.对于官方案例的研究要深入,每行代码都要追溯,它使用了什么lib,加载了什么dll,使用了那些头文件,又没用冗余,使用了那些预编译命令,这些都很重要<br>
