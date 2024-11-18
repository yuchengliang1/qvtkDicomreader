 <div align=center><img width="300" height="150" src="https://s2.ax1x.com/2019/01/06/FbimEd.png"/></div>

#### 已知BUG
请看Docs/BUGs.md<br>

#### 编程参考
编程中参考的软件:<br>
1.[Pmsdview](http://pmsdview-12.updatestar.com/)<br/>
2.[RadiAnt](https://www.radiantviewer.com/)<br/>
Pmsdview是由飞利浦开发的一个程序;RadiAnt是效率极高且功能完善的Dicom浏览器和PACS客户端<br>
这两个程序都是不开源的,本程序的编写中只是使用了这两个软件并借鉴了界面和一些功能.<br>

#### Build:
1. 依赖:DCMTK3.6.2,64位;vtk8.0.0,64位;Qt5.9.1,64位;itk4.12,64位 <br>
2. DCMTK,vtk,itk依赖自带,请不要修改<br>
3. Qt是重量级库,需要自行安装到本地<br>
4. 环境:Visual Studio 2017,v141工具集<br>
5. 若运行时缺少DLL,请检查相关的环境变量是否正确,或者找到这个dll并放在解决方案文件夹的x64/debug或x64/release文件夹中<br>

#### 未来
1. windows版上传依赖,添加编译脚本.
2. 移植到Ubuntu.
3. 删除一些冗余的功能.
4. 使用QtQuick界面.