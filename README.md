作者：杨波，赵耀辉（z184321）

单位：河北工业大学机械工程学院，天津 300401，中国

邮箱：[boyang@hebut.edu.cn](mailto:boyang@hebut.edu.cn)
[18436129578@163.com](mailto:18436129578@163.com)

项目地址：
https://github.com/z184321/KC-D.git

代码简介

为了建立晶界处具有较高有序性的多晶金刚石模型，本程序基于截角正八面体晶粒构建多晶结构。

整个建模过程由三个依次执行的 Python 脚本完成。首先，将具有不同晶体取向的立方体金刚石晶粒截取为正八面体；随后对正八面体的六个顶点进行截角处理，获得截角正八面体晶粒；最后将不同取向的截角正八面体晶粒按照预设空间位置进行排列和拼接，构建包含 16 个晶粒的扩展 Kelvin 多晶金刚石模型。

1. octahedron.py

功能介绍

octahedron.py 用于将原始立方体金刚石晶粒截取为正八面体晶粒。

程序自动读取 POSCAR 文件中的晶格参数并计算模型的几何中心。根据用户指定的截断半径 cutoff，对原子坐标进行筛选，仅保留满足正八面体几何条件的原子。

程序实现的功能

1. 根据晶格向量自动计算模型的几何中心。
2. 支持 Cartesian 笛卡尔坐标和 Direct 分数坐标，并可自动进行坐标转换。
3. 批量处理 O、A、B、C、D、E、F、I 共 8 种不同取向的晶粒。
4. 输出截取后的正八面体金刚石晶粒结构。

运行所需文件

Python 脚本：

octahedron.py

输入文件：

POSCAR_O1
POSCAR_A1
POSCAR_B1
POSCAR_C1
POSCAR_D1
POSCAR_E1
POSCAR_F1
POSCAR_I1

上述文件为原始立方体金刚石晶粒结构。

程序输出文件：

0POSCAR_O
0POSCAR_A
0POSCAR_B
0POSCAR_C
0POSCAR_D
0POSCAR_E
0POSCAR_F
0POSCAR_I

上述文件为截取后的正八面体晶粒结构。

2. Truncated1.py

功能介绍

Truncated1.py 用于对 octahedron.py 生成的正八面体晶粒进行截角处理，从而获得截角正八面体晶粒。

程序通过截断比例参数 truncation_ratio 控制截角程度，默认值为 0.33。程序沿正八面体的 ±x、±y 和 ±z 六个顶点方向切除尖端区域，从而获得截角正八面体结构。

程序实现的功能

1. 计算正八面体晶粒的几何中心，并将原子坐标平移至以晶粒中心为原点的坐标系。
2. 沿 +x、-x、+y、-y、+z 和 -z 六个方向进行等距截断。
3. 批量处理 O、A、B、C、D、E、F、I 共 8 种不同取向的正八面体晶粒。
4. 根据截断后的实际原子数量自动更新 POSCAR 文件中的原子数。

运行所需文件

Python 脚本：

Truncated1.py

输入文件：

0POSCAR_O
0POSCAR_A
0POSCAR_B
0POSCAR_C
0POSCAR_D
0POSCAR_E
0POSCAR_F
0POSCAR_I

上述文件由 octahedron.py 生成。

程序输出文件：

POSCAR_O
POSCAR_A
POSCAR_B
POSCAR_C
POSCAR_D
POSCAR_E
POSCAR_F
POSCAR_I

上述文件为最终得到的截角正八面体金刚石晶粒。

3. Kelvin1.py

功能介绍

Kelvin1.py 用于将 8 种不同取向的截角正八面体晶粒按照指定的空间位置进行排列和拼接，最终生成包含 16 个晶粒的扩展 Kelvin 多晶金刚石模型。

程序能够自动计算晶粒的几何尺寸和排列参数，并对相邻晶粒共享界面附近的原子进行处理，包括删除距离过近的原子以及调整部分晶界原子位置，从而改善晶粒连接区域的原子结构。

程序实现的功能

1. 基于 O、A、B、C、D、E、F、I 共 8 种不同取向的截角正八面体晶粒构建 16 晶粒多晶模型。
2. 自动计算晶粒尺寸以及模型构建所需的几何参数。
3. 检测不同晶粒之间距离过近的原子，并根据最小原子间距阈值进行删除，默认阈值为 1.0 Å。
4. 调整共享界面附近部分原子的位置，使界面区域的 C-C 键长接近金刚石标准键长 1.54 Å。
5. 根据晶粒排列自动计算模拟盒子的尺寸。
6. 对最终模型进行结构质量检查，包括最小原子间距以及 C-C 键长分布等信息。

运行所需文件

Python 脚本：

Kelvin1.py

输入文件：

POSCAR_O
POSCAR_A
POSCAR_B
POSCAR_C
POSCAR_D
POSCAR_E
POSCAR_F
POSCAR_I

上述文件由 Truncated1.py 生成。

程序输出文件：

Extended_Kelvin_Structure_CustomPositions.vasp

该文件为最终生成的 16 晶粒扩展 Kelvin 多晶金刚石模型。

程序运行顺序

三个 Python 脚本需要按照以下顺序依次运行：

octahedron.py
↓
Truncated1.py
↓
Kelvin1.py

对应的结构演化过程为：

立方体金刚石晶粒
↓
正八面体金刚石晶粒
↓
截角正八面体金刚石晶粒
↓
16 晶粒扩展 Kelvin 多晶金刚石模型

运行环境

程序运行需要 Python 环境，并依赖以下 Python 库：

numpy
scipy

可通过以下命令安装所需依赖：

pip install numpy scipy

最终输出

完成上述三个程序后，将得到最终多晶模型文件：

Extended_Kelvin_Structure_CustomPositions.vasp

该文件采用 VASP POSCAR 格式，可进一步用于 VASP、LAMMPS、OVITO 等软件中的结构处理、可视化及原子尺度模拟。
