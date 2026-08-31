作者：杨波，赵耀辉(z184321)
    地址：河北工业大学机械工程学院，天津 300401，中国
    邮箱：<boyang@hebut.edu.cn>，<18436129578@163.com>
    网址：https://github.com/z184321/KC-D.git

## 代码简介

为了建立晶界处有序的多晶金刚石模型，开发了本程序。本程序用于构建基于截角正八面体晶粒的多晶模型。通过三个顺序执行的 Python 脚本，将八种不同取向的晶粒从正方体逐步切割为正八面体、截角八面体，最终拼接为包含 16 个晶粒的多晶模型结构。

---

## octahedron.py

### 功能介绍

将正方体晶粒（原始 POSCAR 文件）从中心截取为正八面体形状。程序自动读取晶格参数计算系统中心，用户指定截断半径（cutoff），保留满足正八面体条件的原子。

### 程序实现的功能

1. 自动从晶格向量计算系统几何中心。
2. 支持笛卡尔坐标和分数坐标输入，自动转换。
3. 批量处理 8 种晶粒类型（O, A, B, C, D, E, F, I）。
4. 输出截取后的正八面体结构。

### 运行必要的文件

- **输入文件**：`POSCAR_O1`、`POSCAR_A1`、`POSCAR_B1`、`POSCAR_C1`、`POSCAR_D1`、`POSCAR_E1`、`POSCAR_F1`、`POSCAR_I1`（原始正方体晶粒）
- **Python 脚本**：`octahedron.py`

### 程序生成的文件

- **输出文件**：`0POSCAR_O`、`0POSCAR_A`、`0POSCAR_B`、`0POSCAR_C`、`0POSCAR_D`、`0POSCAR_E`、`0POSCAR_F`、`0POSCAR_I`（正八面体晶粒）

---

## Truncated1.py

### 功能介绍

对正八面体晶粒进行截角处理，生成截角正八面体。通过截断比例参数（`truncation_ratio`，默认 0.33）控制截断程度，沿 6 个顶点方向切除尖端。

### 程序实现的功能

1. 计算八面体中心，将原子坐标平移到中心。
2. 沿 6 个顶点方向（±x, ±y, ±z）进行等距截断。
3. 批量处理 8 个正八面体晶粒。
4. 自动更新输出文件的原子计数。

### 运行必要的文件

- **输入文件**：`0POSCAR_O`、`0POSCAR_A`、`0POSCAR_B`、`0POSCAR_C`、`0POSCAR_D`、`0POSCAR_E`、`0POSCAR_F`、`0POSCAR_I`（由 `octahedron.py` 生成）
- **Python 脚本**：`Truncated1.py`

### 程序生成的文件

- **输出文件**：`POSCAR_O`、`POSCAR_A`、`POSCAR_B`、`POSCAR_C`、`POSCAR_D`、`POSCAR_E`、`POSCAR_F`、`POSCAR_I`（截角八面体晶粒）

---

## Kelvin1.py

### 功能介绍

将 8 种截角正八面体晶粒按照指定空间位置排列，拼接成包含 16 个晶粒的扩展 Kelvin 多晶模型。程序自动计算晶粒尺寸、调整共享面原子位置、删除过近原子，并输出最终的多晶结构 POSCAR 文件。

### 程序实现的功能

1. 基于 8 种晶粒类型（O, A, B, C, D, E, F, I）生成 16 个晶粒。
2. 自动计算晶粒尺寸和几何参数。
3. 删除晶粒间距离过近的原子（最小间距阈值：1.0 Å）。
4. 调整共享面处原子位置，确保 C-C 键长（1.54 Å）正确。
5. 自动计算盒子参数，确定盒子尺寸。
6. 结构质量检查：输出最小原子间距、C-C 键长分布等信息。

### 运行必要的文件

- **输入文件**：`POSCAR_O`、`POSCAR_A`、`POSCAR_B`、`POSCAR_C`、`POSCAR_D`、`POSCAR_E`、`POSCAR_F`、`POSCAR_I`（由 `Truncated1.py` 生成）
- **Python 脚本**：`Kelvin1.py`

### 程序生成的文件

- **输出文件**：`Extended_Kelvin_Structure_CustomPositions.vasp`（最终多晶模型 POSCAR）

---

## 运行依赖

本软件包需要以下 Python 环境：

```bash
pip install numpy scipy
