import numpy as np
import os

def read_poscar_direct(filename):
    """
    专门读取POSCAR文件，处理Direct坐标并转换为笛卡尔坐标
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 第一行：注释
    comment = lines[0].strip()
    
    # 第二行：缩放因子
    scaling_factor = float(lines[1].split()[0])
    
    # 第三到五行：晶格向量
    lattice_vectors = []
    for i in range(2, 5):
        lattice_vectors.append([float(x) for x in lines[i].split()[:3]])
    lattice_vectors = np.array(lattice_vectors)
    
    # 第六行：元素种类
    elements = lines[5].split()
    
    # 第七行：各元素原子数
    element_counts = [int(x) for x in lines[6].split()]
    
    # 第八行：坐标类型 (Direct 或 Cartesian)
    coord_type = lines[7].split()[0].lower()
    
    # 读取原子坐标
    atoms_frac = []  # 分数坐标
    atom_elements = []  # 记录每个原子对应的元素
    
    idx = 8
    for i, count in enumerate(element_counts):
        for j in range(count):
            if idx < len(lines):
                coords = [float(x) for x in lines[idx].split()[:3]]
                atoms_frac.append(coords)
                atom_elements.append(elements[i])
                idx += 1
    
    atoms_frac = np.array(atoms_frac)
    
    # 将分数坐标转换为笛卡尔坐标
    # 笛卡尔坐标 = 分数坐标 × 晶格向量 × 缩放因子
    atoms_cart = np.dot(atoms_frac, lattice_vectors) * scaling_factor
    
    return atoms_cart, atoms_frac, atom_elements, lattice_vectors, scaling_factor, coord_type, comment, elements, element_counts

def filter_octahedron(atoms_cart, center_x, center_y, center_z, cutoff):
    """
    使用正八面体条件筛选原子
    """
    filtered_indices = []
    
    for i, atom in enumerate(atoms_cart):
        x, y, z = atom
        # 应用正八面体条件
        condition = (abs(x - center_x) + abs(y - center_y) + abs(z - center_z)) <= cutoff
        if condition:
            filtered_indices.append(i)
    
    return filtered_indices

def write_poscar_cartesian(filename, atoms_cart, atom_elements, lattice_vectors, scaling_factor, comment, elements, element_counts):
    """
    将筛选后的结构写入POSCAR文件，使用笛卡尔坐标
    """
    # 计算新的元素计数
    new_element_counts = []
    for elem in elements:
        count = sum(1 for e in atom_elements if e == elem)
        new_element_counts.append(count)
    
    with open(filename, 'w') as f:
        # 写入注释行
        f.write(f"{comment} - Octahedron (Cartesian)\n")
        
        # 写入缩放因子
        f.write(f"  {scaling_factor:.10f}\n")
        
        # 写入晶格向量
        for vec in lattice_vectors:
            f.write(f"  {vec[0]:.10f}  {vec[1]:.10f}  {vec[2]:.10f}\n")
        
        # 写入元素种类
        f.write("  " + "  ".join(elements) + "\n")
        
        # 写入各元素原子数
        f.write("  " + "  ".join(str(x) for x in new_element_counts) + "\n")
        
        # 写入坐标类型 - 使用笛卡尔坐标
        f.write("Cartesian\n")
        
        # 写入原子坐标 (笛卡尔坐标，已除以缩放因子)
        for atom in atoms_cart:
            cart_coords = atom / scaling_factor
            f.write(f"  {cart_coords[0]:.16f}  {cart_coords[1]:.16f}  {cart_coords[2]:.16f}\n")

def write_poscar_direct(filename, atoms_frac, atom_elements, lattice_vectors, scaling_factor, comment, elements, element_counts):
    """
    将筛选后的结构写入POSCAR文件，使用分数坐标
    """
    # 计算新的元素计数
    new_element_counts = []
    for elem in elements:
        count = sum(1 for e in atom_elements if e == elem)
        new_element_counts.append(count)
    
    with open(filename, 'w') as f:
        # 写入注释行
        f.write(f"{comment} - Octahedron (Direct)\n")
        
        # 写入缩放因子
        f.write(f"  {scaling_factor:.10f}\n")
        
        # 写入晶格向量
        for vec in lattice_vectors:
            f.write(f"  {vec[0]:.10f}  {vec[1]:.10f}  {vec[2]:.10f}\n")
        
        # 写入元素种类
        f.write("  " + "  ".join(elements) + "\n")
        
        # 写入各元素原子数
        f.write("  " + "  ".join(str(x) for x in new_element_counts) + "\n")
        
        # 写入坐标类型 - 使用分数坐标
        f.write("Direct\n")
        
        # 写入原子坐标 (分数坐标)
        for atom in atoms_frac:
            f.write(f"  {atom[0]:.16f}  {atom[1]:.16f}  {atom[2]:.16f}\n")

def calculate_system_center_and_size(lattice_vectors, scaling_factor):
    """
    直接从晶格向量计算系统中心和尺寸
    """
    # 计算晶格向量的模（尺寸）
    size_x = np.linalg.norm(lattice_vectors[0]) * scaling_factor
    size_y = np.linalg.norm(lattice_vectors[1]) * scaling_factor
    size_z = np.linalg.norm(lattice_vectors[2]) * scaling_factor
    
    # 计算系统中心（假设盒子是立方体）
    center_x = size_x / 2
    center_y = size_y / 2
    center_z = size_z / 2
    
    return center_x, center_y, center_z, size_x, size_y, size_z

def process_single_file(input_file, output_file, center_x, center_y, center_z, cutoff, coord_type_choice):
    """
    处理单个文件
    """
    print(f"正在处理文件: {input_file}")
    
    try:
        # 读取POSCAR文件
        atoms_cart, atoms_frac, atom_elements, lattice_vectors, scaling_factor, coord_type, comment, elements, element_counts = read_poscar_direct(input_file)
        
        print(f"  找到 {len(atoms_cart)} 个原子")
        print(f"  原始坐标类型: {coord_type}")
        
        # 筛选原子
        filtered_indices = filter_octahedron(atoms_cart, center_x, center_y, center_z, cutoff)
        print(f"  筛选后剩余 {len(filtered_indices)} 个原子")
        
        # 计算筛选比例
        if len(atoms_cart) > 0:
            ratio = len(filtered_indices) / len(atoms_cart) * 100
            print(f"  筛选比例: {ratio:.2f}%")
        
        # 提取筛选后的原子
        filtered_atoms_cart = atoms_cart[filtered_indices]
        filtered_atoms_frac = atoms_frac[filtered_indices]
        filtered_atom_elements = [atom_elements[i] for i in filtered_indices]
        
        # 显示元素分布变化
        for elem in elements:
            original_count = element_counts[elements.index(elem)]
            new_count = sum(1 for e in filtered_atom_elements if e == elem)
            print(f"  {elem}: {original_count} -> {new_count}")
        
        # 写入输出文件
        if coord_type_choice == "1":
            write_poscar_cartesian(output_file, filtered_atoms_cart, filtered_atom_elements, 
                                 lattice_vectors, scaling_factor, comment, elements, element_counts)
        else:
            write_poscar_direct(output_file, filtered_atoms_frac, filtered_atom_elements, 
                              lattice_vectors, scaling_factor, comment, elements, element_counts)
        
        print(f"  结果已保存到: {output_file}")
        print("-" * 60)
        
        return True
        
    except Exception as e:
        print(f"  处理文件 {input_file} 时发生错误: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """
    主函数 - 批量处理8个POSCAR文件，将Direct坐标转换为笛卡尔坐标后截取正八面体
    """
    # 定义输入文件和输出文件列表
    input_files = ["POSCAR_A1", "POSCAR_B1", "POSCAR_C1", "POSCAR_D1", 
                   "POSCAR_E1", "POSCAR_F1", "POSCAR_I1", "POSCAR_O1"]
    
    # 输出文件名改为0POSCAR_A等
    output_files = ["0POSCAR_A", "0POSCAR_B", "0POSCAR_C", "0POSCAR_D",
                    "0POSCAR_E", "0POSCAR_F", "0POSCAR_I", "0POSCAR_O"]
    
    print("POSCAR正八面体批量截取程序 (处理Direct坐标)")
    print("=" * 60)
    print(f"输入文件: {', '.join(input_files)}")
    print(f"输出文件: {', '.join(output_files)}")
    print("=" * 60)
    
    # 检查所有输入文件是否存在
    missing_files = []
    for input_file in input_files:
        if not os.path.exists(input_file):
            missing_files.append(input_file)
    
    if missing_files:
        print(f"错误: 找不到以下文件: {', '.join(missing_files)}")
        return
    
    try:
        # 以POSCAR_O1为基准文件获取系统参数
        print("\n以POSCAR_O1为基准文件获取系统参数...")
        base_file = "POSCAR_O1"
        atoms_cart, atoms_frac, atom_elements, lattice_vectors, scaling_factor, coord_type, comment, elements, element_counts = read_poscar_direct(base_file)
        
        print(f"基准文件信息:")
        print(f"  找到 {len(atoms_cart)} 个原子")
        print(f"  坐标类型: {coord_type}")
        print(f"  晶格缩放因子: {scaling_factor}")
        print(f"  元素组成: {dict(zip(elements, element_counts))}")
        
        # 直接从晶格向量计算系统中心和尺寸
        center_x, center_y, center_z, size_x, size_y, size_z = calculate_system_center_and_size(lattice_vectors, scaling_factor)
        
        print(f"\n系统几何中心: [{center_x:.6f}, {center_y:.6f}, {center_z:.6f}]")
        print(f"系统尺寸: X={size_x:.6f}, Y={size_y:.6f}, Z={size_z:.6f}")
        
        # 获取用户输入参数（以基准文件为准）
        print("\n请输入正八面体的中心坐标和截断值(以POSCAR_O1为基准):")
        print(f"推荐中心坐标: [{center_x:.6f}, {center_y:.6f}, {center_z:.6f}]")
        
        use_default_center = input("使用推荐中心坐标? (y/n): ").strip().lower()
        if use_default_center == 'y':
            user_center_x, user_center_y, user_center_z = center_x, center_y, center_z
            print(f"使用中心坐标: [{user_center_x:.6f}, {user_center_y:.6f}, {user_center_z:.6f}]")
        else:
            user_center_x = float(input("中心点 x 坐标: "))
            user_center_y = float(input("中心点 y 坐标: "))
            user_center_z = float(input("中心点 z 坐标: "))
        
        # 计算推荐的截断值
        min_size = min(size_x, size_y, size_z)
        recommended_cutoff = min_size * 0.4  # 推荐使用系统最小尺寸的40%作为截断值
        
        print(f"推荐截断值: {recommended_cutoff:.6f} (系统最小尺寸的40%)")
        user_cutoff = float(input("请输入截断值 cutoff: "))
        
        # 选择输出坐标类型
        print("\n选择输出坐标类型 (所有文件使用相同的坐标类型):")
        print("1. 笛卡尔坐标 (Cartesian)")
        print("2. 分数坐标 (Direct)")
        coord_choice = input("请输入选择 (1 或 2): ").strip()
        
        # 批量处理所有文件
        print("\n开始批量处理文件...")
        print("=" * 60)
        
        success_count = 0
        for i in range(len(input_files)):
            input_file = input_files[i]
            output_file = output_files[i]
            
            success = process_single_file(input_file, output_file, 
                                        user_center_x, user_center_y, user_center_z, 
                                        user_cutoff, coord_choice)
            
            if success:
                success_count += 1
        
        print("=" * 60)
        print(f"\n批量处理完成!")
        print(f"成功处理: {success_count}/{len(input_files)} 个文件")
        
        if coord_choice == "1":
            print(f"输出坐标类型: 笛卡尔坐标 (Cartesian)")
        else:
            print(f"输出坐标类型: 分数坐标 (Direct)")
        
        print(f"输出文件: {', '.join(output_files)}")
        
    except Exception as e:
        print(f"发生错误: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()