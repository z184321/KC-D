import numpy as np
import math

def read_poscar(filename):
    """读取POSCAR文件"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 读取系统名称
    system_name = lines[0].strip()
    
    # 读取缩放因子
    scale = float(lines[1].strip())
    
    # 读取晶格向量
    lattice_vectors = []
    for i in range(2, 5):
        lattice_vectors.append([float(x) for x in lines[i].split()])
    lattice_vectors = np.array(lattice_vectors) * scale
    
    # 读取元素类型和原子数量
    elements = lines[5].split()
    atom_counts = [int(x) for x in lines[6].split()]
    total_atoms = sum(atom_counts)
    
    # 检查坐标类型
    coord_type = lines[7].strip().lower()
    
    # 读取原子坐标
    coords = []
    for i in range(8, 8 + total_atoms):
        if lines[i].strip():  # 跳过空行
            coords.append([float(x) for x in lines[i].split()[:3]])
    
    coords = np.array(coords)
    
    # 如果是分数坐标，转换为笛卡尔坐标
    if coord_type[0] == 'd':
        coords = np.dot(coords, lattice_vectors)
    
    return system_name, lattice_vectors, elements, atom_counts, coords

def write_poscar(filename, system_name, lattice_vectors, elements, atom_counts, coords):
    """写入POSCAR文件"""
    with open(filename, 'w') as f:
        # 写入系统名称
        f.write(f"{system_name}_truncated\n")
        
        # 写入缩放因子
        f.write("1.0\n")
        
        # 写入晶格向量
        for vec in lattice_vectors:
            f.write(f"{vec[0]:20.16f} {vec[1]:20.16f} {vec[2]:20.16f}\n")
        
        # 写入元素类型和原子数量
        f.write(" ".join(elements) + "\n")
        f.write(" ".join(map(str, atom_counts)) + "\n")
        
        # 写入坐标类型
        f.write("Cartesian\n")
        
        # 写入原子坐标
        for coord in coords:
            f.write(f"{coord[0]:20.16f} {coord[1]:20.16f} {coord[2]:20.16f}\n")

def truncate_octahedron(coords, truncation_ratio=0.33):
    """
    截断正八面体
    
    参数:
    coords: 原子坐标数组
    truncation_ratio: 截断比例 (0-1)，值越大截掉的部分越多
    
    返回:
    保留的原子坐标
    """
    # 计算八面体的中心
    center = np.mean(coords, axis=0)
    
    # 将坐标平移到中心
    centered_coords = coords - center
    
    # 八面体的顶点方向（正八面体的6个顶点）
    octahedron_directions = [
        [1, 0, 0], [-1, 0, 0],
        [0, 1, 0], [0, -1, 0], 
        [0, 0, 1], [0, 0, -1]
    ]
    
    # 计算每个顶点方向的截断平面
    # 截断距离 = 最大距离 × (1 - truncation_ratio)
    max_distance = np.max(np.linalg.norm(centered_coords, axis=1))
    truncation_distance = max_distance * (1 - truncation_ratio)
    
    # 筛选保留的原子
    kept_indices = []
    
    for i, coord in enumerate(centered_coords):
        keep_atom = True
        
        # 检查原子是否在所有截断平面之内
        for direction in octahedron_directions:
            # 计算原子在该方向上的投影距离
            projection = np.dot(coord, direction)
            
            # 如果投影距离超过截断距离，则删除该原子
            if projection > truncation_distance:
                keep_atom = False
                break
        
        if keep_atom:
            kept_indices.append(i)
    
    # 返回保留的原子坐标（恢复原始位置）
    kept_coords = coords[kept_indices]
    
    print(f"原始原子数: {len(coords)}")
    print(f"截断后原子数: {len(kept_coords)}")
    print(f"截断比例: {truncation_ratio}")
    
    return kept_coords

def process_file(input_file, output_file, truncation_ratio=0.33):
    """处理单个文件"""
    try:
        system_name, lattice_vectors, elements, atom_counts, coords = read_poscar(input_file)
        print(f"成功读取 {input_file}")
        print(f"系统: {system_name}")
        print(f"总原子数: {len(coords)}")
        
        # 截断正八面体
        truncated_coords = truncate_octahedron(coords, truncation_ratio=truncation_ratio)
        
        # 更新原子计数（假设只有一种元素）
        new_atom_counts = [len(truncated_coords)]
        
        # 写入新的POSCAR文件
        write_poscar(output_file, system_name, lattice_vectors, elements, new_atom_counts, truncated_coords)
        print(f"截角正八面体已保存到 {output_file}\n")
        
    except FileNotFoundError:
        print(f"错误：找不到文件 {input_file}")
    except Exception as e:
        print(f"处理文件 {input_file} 时出错: {e}\n")

def main():
    """主函数，处理所有8个文件"""
    # 定义输入文件和输出文件的映射
    file_pairs = [
        ("0POSCAR_A", "POSCAR_A"),
        ("0POSCAR_B", "POSCAR_B"),
        ("0POSCAR_C", "POSCAR_C"),
        ("0POSCAR_D", "POSCAR_D"),
        ("0POSCAR_E", "POSCAR_E"),
        ("0POSCAR_F", "POSCAR_F"),
        ("0POSCAR_I", "POSCAR_I"),
        ("0POSCAR_O", "POSCAR_O")
    ]
    
    print("开始处理8个正八面体文件...\n")
    
    # 处理每个文件
    for input_file, output_file in file_pairs:
        print(f"处理文件: {input_file} -> {output_file}")
        print("-" * 40)
        process_file(input_file, output_file)
    
    print("所有文件处理完成！")

if __name__ == "__main__":
    main()