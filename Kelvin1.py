import numpy as np
from scipy.spatial import cKDTree
import os
from typing import Dict, List, Any, Tuple

def read_poscar(filename: str) -> Dict[str, Any]:
    """读取POSCAR文件"""
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    
    # 读取系统名称
    system_name = lines[0]
    
    # 读取缩放因子
    scaling_factor = float(lines[1])
    
    # 读取晶格矢量
    lattice_vectors = []
    for i in range(2, 5):
        lattice_vectors.append([float(x) for x in lines[i].split()])
    lattice_vectors = np.array(lattice_vectors) * scaling_factor
    
    # 读取原子类型和数量
    atom_types = lines[5].split()
    atom_counts = [int(x) for x in lines[6].split()]
    
    # 检查坐标格式
    coord_type = lines[7].strip()[0].upper()
    
    # 读取原子坐标
    atoms = []
    index = 8
    for i, count in enumerate(atom_counts):
        for j in range(count):
            coords = [float(x) for x in lines[index].split()[:3]]
            if coord_type == 'D':  # Direct coordinates
                # 转换为笛卡尔坐标
                cartesian_coords = np.dot(coords, lattice_vectors)
                atoms.append({
                    'type': atom_types[i],
                    'fractional': np.array(coords),
                    'cartesian': cartesian_coords
                })
            else:  # Cartesian coordinates
                atoms.append({
                    'type': atom_types[i],
                    'fractional': np.dot(coords, np.linalg.inv(lattice_vectors)),
                    'cartesian': np.array(coords) * scaling_factor
                })
            index += 1
    
    return {
        'system_name': system_name,
        'scaling_factor': scaling_factor,
        'lattice_vectors': lattice_vectors,
        'atom_types': atom_types,
        'atom_counts': atom_counts,
        'coord_type': coord_type,
        'atoms': atoms
    }

class ExtendedKelvinStructure:
    def __init__(self, poscar_files: Dict[str, str]):
        """
        基于8个POSCAR文件初始化扩展Kelvin结构
        poscar_files: 字典，键为晶粒类型，值为文件名
        """
        self.poscar_files = poscar_files
        self.grain_types = ['O', 'A', 'B', 'C', 'D', 'E', 'F', 'I']
        
        # 读取所有晶粒类型
        self.grain_data = {}
        for grain_type in self.grain_types:
            if grain_type in poscar_files:
                filename = poscar_files[grain_type]
                if os.path.exists(filename):
                    self.grain_data[grain_type] = read_poscar(filename)
                    print(f"成功读取 {grain_type} 类型晶粒: {filename}")
                else:
                    print(f"错误: 找不到文件 {filename}")
                    self.grain_data[grain_type] = None
        
        # 检查所有文件是否都已读取
        if not all(self.grain_data.values()):
            missing = [gtype for gtype in self.grain_types 
                      if self.grain_data.get(gtype) is None]
            print(f"错误: 以下晶粒类型的数据缺失: {missing}")
            return
        
        # 金刚石C-C键长 (0.154 nm = 1.54 Å)
        self.cc_bond_length = 1.54  # Å
        
        # 最小允许原子间距 (0.1 nm = 1.0 Å)
        self.min_distance = 1.0  # Å
        
        # 使用O晶粒计算几何参数
        self.calculate_geometry_parameters()
    
    def calculate_geometry_parameters(self) -> None:
        """计算截角正八面体的几何参数 - 基于O晶粒"""
        o_data = self.grain_data['O']
        atoms = o_data['atoms']
        coords = np.array([atom['cartesian'] for atom in atoms])
        
        # 计算中心点
        center = np.mean(coords, axis=0)
        
        # 计算最大距离（顶点到中心的距离）
        distances = np.linalg.norm(coords - center, axis=1)
        self.radius = np.max(distances)
        
        # 从原子坐标估算正方形面宽度
        x_coords = coords[:, 0]
        y_coords = coords[:, 1]
        z_coords = coords[:, 2]
        
        # 找到在x、y、z方向上形成平面的原子
        x_unique = np.unique(np.round(x_coords, decimals=3))
        y_unique = np.unique(np.round(y_coords, decimals=3))
        z_unique = np.unique(np.round(z_coords, decimals=3))
        
        # 计算可能的正方形面尺寸
        if len(x_unique) > 1:
            self.grain_size = np.max(x_unique) - np.min(x_unique)
        else:
            # 如果无法直接计算，使用经验公式
            self.grain_size = 2 * self.radius / np.sqrt(2)
        
        print(f"\n=== O晶粒几何参数 ===")
        print(f"晶粒尺寸: {self.grain_size:.3f} Å")
        print(f"近似半径: {self.radius:.3f} Å")
    
    def assign_grains_to_lattice_points(self) -> List[Dict[str, Any]]:
        """将晶粒分配到指定的格点上"""
        print("\n=== 分配晶粒到指定位置 ===")
        
        # 定义16个晶粒的位置倍数
        # 注意：位置倍数是-1, 0, 1, 2, 3, 4，对应实际位置为：-1/2, 0, 1/2, 1, 3/2, 2倍的晶粒尺寸
        grain_multipliers = [
            (0, 'A', np.array([0.0, 0.0, 0.0])),    # 0倍 -> 0
            (1, 'B', np.array([2.0, 0.0, 0.0])),    # 2倍 -> 1倍晶粒尺寸
            (2, 'C', np.array([0.0, 2.0, 0.0])),    # 2倍 -> 1倍晶粒尺寸
            (3, 'D', np.array([2.0, 2.0, 0.0])),    # 2倍 -> 1倍晶粒尺寸
            (4, 'O', np.array([1.0, 1.0, 1.0])),    # 1倍 -> 1/2倍晶粒尺寸
            (5, 'F', np.array([3.0, 1.0, 1.0])),    # 3倍 -> 3/2倍晶粒尺寸
            (6, 'E', np.array([1.0, 3.0, 1.0])),    # 3倍 -> 3/2倍晶粒尺寸
            (7, 'I', np.array([3.0, 3.0, 1.0])),    # 3倍 -> 3/2倍晶粒尺寸
            (8, 'A', np.array([2.0, 2.0, 2.0])),    # 2倍 -> 1倍晶粒尺寸
            (9, 'B', np.array([4.0, 2.0, 2.0])),    # 4倍 -> 2倍晶粒尺寸
            (10, 'C', np.array([2.0, 4.0, 2.0])),   # 4倍 -> 2倍晶粒尺寸
            (11, 'D', np.array([4.0, 4.0, 2.0])),   # 4倍 -> 2倍晶粒尺寸
            (12, 'O', np.array([-1.0, -1.0, -1.0])), # -1倍 -> -1/2倍晶粒尺寸
            (13, 'F', np.array([1.0, -1.0, -1.0])),  # 1倍 -> 1/2倍晶粒尺寸
            (14, 'E', np.array([-1.0, 1.0, -1.0])),  # -1倍 -> -1/2倍晶粒尺寸
            (15, 'I', np.array([1.0, 1.0, -1.0]))    # 1倍 -> 1/2倍晶粒尺寸
        ]
        
        # 转换倍数到实际坐标：实际坐标 = 倍数 × (晶粒尺寸/2)
        half_grain_size = self.grain_size / 2
        print(f"晶粒尺寸: {self.grain_size:.3f} Å, 半尺寸: {half_grain_size:.3f} Å")
        
        # 创建晶粒分配列表
        grain_assignments = []
        
        for grain_id, grain_type, multiplier in grain_multipliers:
            # 计算实际中心位置: 位置 = 倍数 × (晶粒尺寸/2)
            center = multiplier * half_grain_size
            
            grain_assignments.append({
                'type': grain_type,
                'center': center,
                'grain_id': grain_id,
                'multiplier': multiplier
            })
            
            print(f"晶粒 {grain_id:2d}: 类型 {grain_type}, "
                  f"倍数({multiplier[0]:.0f},{multiplier[1]:.0f},{multiplier[2]:.0f}), "
                  f"中心({center[0]:.3f},{center[1]:.3f},{center[2]:.3f}) Å")
        
        return grain_assignments
    
    def create_extended_kelvin_structure(self) -> List[Dict[str, Any]]:
        """创建扩展Kelvin结构"""
        # 获取晶粒分配
        grain_assignments = self.assign_grains_to_lattice_points()
        
        grains = []
        
        print(f"\n=== 创建晶粒 ===")
        
        # 预计算所有晶粒的中心偏移量
        for assignment in grain_assignments:
            grain_type = assignment['type']
            grain_center = assignment['center']
            grain_id = assignment['grain_id']
            
            # 获取该类型晶粒的数据
            grain_data = self.grain_data[grain_type]
            base_atoms = grain_data['atoms']
            
            # 计算基础晶粒的中心
            base_coords = np.array([atom['cartesian'] for atom in base_atoms])
            base_center = np.mean(base_coords, axis=0)
            
            # 创建晶粒
            grain_atoms = []
            for atom in base_atoms:
                # 移动原子到目标位置
                new_cartesian = atom['cartesian'] - base_center + grain_center
                grain_atoms.append({
                    'type': atom['type'],
                    'cartesian': new_cartesian,
                    'grain_id': grain_id,
                    'grain_type': grain_type
                })
            
            grains.append({
                'atoms': grain_atoms,
                'center': grain_center,
                'grain_id': grain_id,
                'type': grain_type,
                'num_atoms': len(grain_atoms)
            })
        
        total_atoms = sum(len(grain['atoms']) for grain in grains)
        print(f"总共创建了 {len(grains)} 个晶粒，总原子数 {total_atoms}")
        return grains
    
    def remove_close_inter_grain_atoms(self, grains: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """删除晶粒间距离小于最小距离的原子"""
        print("\n=== 删除距离过近的原子 ===")
        
        # 收集所有原子
        all_atoms = []
        grain_ids = []
        for grain in grains:
            for atom in grain['atoms']:
                all_atoms.append(atom)
                grain_ids.append(grain['grain_id'])
        
        all_coords = np.array([atom['cartesian'] for atom in all_atoms])
        grain_ids = np.array(grain_ids)
        
        # 使用KD树找到距离过近的原子对
        tree = cKDTree(all_coords)
        close_pairs = tree.query_pairs(self.min_distance, output_type='set')
        
        if len(close_pairs) > 0:
            print(f"发现 {len(close_pairs)} 对距离过近的原子")
            
            # 对于每个过近距离对，删除其中一个原子（只处理不同晶粒间的原子）
            atoms_to_remove = set()
            for i, j in close_pairs:
                if grain_ids[i] != grain_ids[j]:
                    # 选择删除距离晶粒中心较远的原子
                    grain_i = grains[grain_ids[i]]
                    grain_j = grains[grain_ids[j]]
                    
                    dist_i = np.linalg.norm(all_coords[i] - grain_i['center'])
                    dist_j = np.linalg.norm(all_coords[j] - grain_j['center'])
                    
                    if dist_i > dist_j:
                        atoms_to_remove.add(i)
                    else:
                        atoms_to_remove.add(j)
            
            print(f"将删除 {len(atoms_to_remove)} 个晶粒间距离过近的原子")
            
            # 从后向前删除原子，避免索引变化
            atoms_to_remove_sorted = sorted(atoms_to_remove, reverse=True)
            for idx in atoms_to_remove_sorted:
                del all_atoms[idx]
            
            # 重新分配原子到晶粒
            grains = self.redistribute_atoms_to_grains(grains, all_atoms)
        else:
            print("未发现距离过近的原子")
        
        return grains
    
    def redistribute_atoms_to_grains(self, grains: List[Dict[str, Any]], 
                                    all_atoms: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """将原子重新分配到晶粒"""
        # 清空所有晶粒的原子列表
        for grain in grains:
            grain['atoms'] = []
            grain['num_atoms'] = 0
        
        # 根据grain_id重新分配原子
        for atom in all_atoms:
            grain_id = atom['grain_id']
            grains[grain_id]['atoms'].append(atom)
            grains[grain_id]['num_atoms'] += 1
        
        return grains
    
    def adjust_shared_faces(self, grains: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """调整共享面处的原子位置，确保C-C键长正确"""
        # 首先删除晶粒间距离过近的原子
        grains = self.remove_close_inter_grain_atoms(grains)
        
        print("\n=== 调整共享面原子位置 ===")
        
        # 收集所有原子
        all_atoms = []
        for grain in grains:
            all_atoms.extend(grain['atoms'])
        
        all_coords = np.array([atom['cartesian'] for atom in all_atoms])
        
        # 使用KD树找到最近的原子对
        tree = cKDTree(all_coords)
        
        # 找到所有距离接近C-C键长的原子对
        pairs = tree.query_pairs(self.cc_bond_length * 1.2, output_type='set')
        
        if len(pairs) > 0:
            print(f"找到 {len(pairs)} 对接近C-C键长的原子")
            
            # 调整这些原子对的位置，使它们的距离等于C-C键长
            adjusted_pairs = 0
            for i, j in pairs:
                pos_i = all_coords[i]
                pos_j = all_coords[j]
                
                # 计算当前距离
                current_distance = np.linalg.norm(pos_i - pos_j)
                
                if abs(current_distance - self.cc_bond_length) > 0.01:  # 如果距离偏差较大
                    # 计算调整向量
                    direction = (pos_j - pos_i) / current_distance
                    adjustment = (current_distance - self.cc_bond_length) / 2
                    
                    # 调整两个原子的位置
                    all_coords[i] += direction * adjustment
                    all_coords[j] -= direction * adjustment
                    adjusted_pairs += 1
            
            print(f"调整了 {adjusted_pairs} 对原子的位置")
            
            # 更新所有原子的坐标
            idx = 0
            for grain in grains:
                for atom in grain['atoms']:
                    atom['cartesian'] = all_coords[idx]
                    idx += 1
        else:
            print("未找到接近C-C键长的原子对")
        
        return grains
    
    def calculate_box_parameters(self, grains: List[Dict[str, Any]]) -> Tuple[np.ndarray, float]:
        """计算盒子参数：以晶粒0作为盒子顶点，盒子大小等于晶粒0与晶粒9在X方向的距离"""
        # 找到晶粒0和晶粒9
        grain0 = next((g for g in grains if g['grain_id'] == 0), None)
        grain9 = next((g for g in grains if g['grain_id'] == 9), None)
        
        if grain0 is None or grain9 is None:
            raise ValueError("未找到晶粒0或晶粒9")
        
        # 计算晶粒0与晶粒9在X方向的距离
        distance_x = grain9['center'][0] - grain0['center'][0]
        
        print(f"\n=== 盒子参数计算 ===")
        print(f"晶粒0与晶粒9在X方向的距离: {distance_x:.3f} Å")
        
        # 以晶粒0的原子最小坐标作为盒子顶点
        grain0_coords = np.array([atom['cartesian'] for atom in grain0['atoms']])
        box_origin = np.min(grain0_coords, axis=0)
        
        print(f"盒子原点(晶粒0最小坐标): ({box_origin[0]:.3f}, {box_origin[1]:.3f}, {box_origin[2]:.3f}) Å")
        print(f"盒子尺寸: {distance_x:.3f} Å 立方体")
        
        return box_origin, distance_x
    
    def write_poscar(self, grains: List[Dict[str, Any]], 
                    filename: str = "Extended_Kelvin_Structure.vasp") -> None:
        """将扩展Kelvin结构写入POSCAR文件"""
        # 计算盒子参数
        box_origin, box_size = self.calculate_box_parameters(grains)
        
        # 收集所有原子并平移，使盒子原点在(0,0,0)
        all_atoms = []
        for grain in grains:
            for atom in grain['atoms']:
                # 平移原子，使盒子原点在(0,0,0)
                translated_coords = atom['cartesian'] - box_origin
                all_atoms.append({
                    'type': atom['type'],
                    'cartesian': translated_coords
                })
        
        # 检查所有原子是否都在盒子内
        all_coords = np.array([atom['cartesian'] for atom in all_atoms])
        min_coords = np.min(all_coords, axis=0)
        max_coords = np.max(all_coords, axis=0)
        
        print(f"原子坐标范围:")
        print(f"  X: [{min_coords[0]:.3f}, {max_coords[0]:.3f}] (盒子尺寸: {box_size:.3f})")
        print(f"  Y: [{min_coords[1]:.3f}, {max_coords[1]:.3f}] (盒子尺寸: {box_size:.3f})")
        print(f"  Z: [{min_coords[2]:.3f}, {max_coords[2]:.3f}] (盒子尺寸: {box_size:.3f})")
        
        # 检查是否有原子超出盒子边界
        if (max_coords[0] > box_size + 0.01 or max_coords[1] > box_size + 0.01 or 
            max_coords[2] > box_size + 0.01 or min_coords[0] < -0.01 or 
            min_coords[1] < -0.01 or min_coords[2] < -0.01):
            print("警告: 部分原子可能超出盒子边界")
        
        # 统计原子类型
        atom_type_count = {}
        for atom in all_atoms:
            atom_type = atom['type']
            atom_type_count[atom_type] = atom_type_count.get(atom_type, 0) + 1
        
        # 按原子类型排序
        atom_types = sorted(atom_type_count.keys())
        atom_counts = [atom_type_count[t] for t in atom_types]
        
        with open(filename, 'w') as f:
            # 系统名称
            f.write("Extended Kelvin Structure with 16 grains (Custom Positions)\n")
            
            # 缩放因子
            f.write("1.0\n")
            
            # 超胞晶格矢量（立方体）
            f.write(f"  {box_size:.10f}  0.0000000000  0.0000000000\n")
            f.write(f"  0.0000000000  {box_size:.10f}  0.0000000000\n")
            f.write(f"  0.0000000000  0.0000000000  {box_size:.10f}\n")
            
            # 原子类型
            f.write("  " + "  ".join(atom_types) + "\n")
            
            # 每个类型的原子数量
            f.write("  " + "  ".join(str(count) for count in atom_counts) + "\n")
            
            # 坐标类型
            f.write("Cartesian\n")
            
            # 所有原子的笛卡尔坐标
            for atom in all_atoms:
                coords = atom['cartesian']
                f.write(f"  {coords[0]:.10f}  {coords[1]:.10f}  {coords[2]:.10f}\n")
        
        print(f"\n=== 结构输出 ===")
        print(f"结构已写入文件: {filename}")
        print(f"盒子尺寸: {box_size:.3f} Å 立方体")
        print(f"总原子数: {len(all_atoms)}")
        print(f"原子类型和数量: {dict(zip(atom_types, atom_counts))}")
    
    def check_structure_quality(self, grains: List[Dict[str, Any]]) -> Tuple[float, int]:
        """检查扩展Kelvin结构的质量"""
        # 收集所有原子
        all_atoms = []
        for grain in grains:
            all_atoms.extend(grain['atoms'])
        
        all_coords = np.array([atom['cartesian'] for atom in all_atoms])
        
        # 使用KD树快速计算最小原子间距
        tree = cKDTree(all_coords)
        distances, _ = tree.query(all_coords, k=2)
        min_distance = np.min(distances[:, 1])
        
        # 计算原子间距分布
        pairs = tree.query_pairs(self.cc_bond_length * 2.0, output_type='set')
        pair_distances = []
        for i, j in pairs:
            pair_distances.append(np.linalg.norm(all_coords[i] - all_coords[j]))
        
        pair_distances = np.array(pair_distances) if pair_distances else np.array([])
        
        # 统计接近C-C键长的原子对数量
        cc_bond_pairs = 0
        too_close_pairs = 0
        if len(pair_distances) > 0:
            cc_bond_pairs = np.sum((pair_distances > self.cc_bond_length * 0.9) & 
                                  (pair_distances < self.cc_bond_length * 1.1))
            too_close_pairs = np.sum(pair_distances < self.min_distance)
        
        print("\n=== 结构质量检查 ===")
        print(f"总原子数: {len(all_atoms)}")
        print(f"最小原子间距: {min_distance:.3f} Å")
        print(f"目标C-C键长: {self.cc_bond_length} Å")
        print(f"接近C-C键长的原子对数量: {cc_bond_pairs}")
        print(f"过近的原子对数量 (<{self.min_distance} Å): {too_close_pairs}")
        
        if len(pair_distances) > 0:
            print(f"原子间距分布: [{np.min(pair_distances):.3f}, {np.max(pair_distances):.3f}] Å")
        
        return min_distance, cc_bond_pairs

def main():
    # 定义8个POSCAR文件
    poscar_files = {
        'O': 'POSCAR_O',
        'A': 'POSCAR_A', 
        'B': 'POSCAR_B',
        'C': 'POSCAR_C',
        'D': 'POSCAR_D',
        'E': 'POSCAR_E',
        'F': 'POSCAR_F',
        'I': 'POSCAR_I'
    }
    
    # 检查所有文件是否存在
    missing_files = []
    for grain_type, filename in poscar_files.items():
        if not os.path.exists(filename):
            missing_files.append(filename)
    
    if missing_files:
        print("错误: 以下文件不存在:")
        for filename in missing_files:
            print(f"  {filename}")
        print("请确保所有POSCAR文件在当前目录下")
        return
    
    print("=" * 60)
    print("扩展Kelvin结构生成器")
    print("基于文献: Constrained minimal-interface structures in polycrystalline copper")
    print("使用8种晶粒类型: O, A, B, C, D, E, F, I")
    print("排列方式: 自定义16个晶粒位置（基于晶粒尺寸倍数）")
    print("=" * 60)
    
    try:
        # 创建扩展Kelvin结构
        kelvin = ExtendedKelvinStructure(poscar_files)
        
        # 生成16个晶粒
        grains = kelvin.create_extended_kelvin_structure()
        
        # 调整共享面处的原子位置
        grains = kelvin.adjust_shared_faces(grains)
        
        # 检查结构质量
        min_distance, cc_bond_pairs = kelvin.check_structure_quality(grains)
        
        # 写入新的POSCAR文件
        kelvin.write_poscar(grains, "Extended_Kelvin_Structure_CustomPositions.vasp")
        
        # 显示总结信息
        print("\n" + "=" * 60)
        print("总结")
        print("=" * 60)
        print(f"晶粒尺寸: {kelvin.grain_size:.3f} Å")
        print(f"晶粒近似半径: {kelvin.radius:.3f} Å")
        print(f"晶粒总数: 16")
        
        # 统计晶粒类型分布
        grain_type_count = {}
        for grain in grains:
            grain_type = grain['type']
            grain_type_count[grain_type] = grain_type_count.get(grain_type, 0) + 1
        
        print("\n晶粒类型分布:")
        for gtype in sorted(grain_type_count.keys()):
            print(f"  {gtype}: {grain_type_count[gtype]}个晶粒")
        
        # 计算晶粒0和晶粒9的距离
        grain0 = next((g for g in grains if g['grain_id'] == 0), None)
        grain9 = next((g for g in grains if g['grain_id'] == 9), None)
        
        if grain0 is not None and grain9 is not None:
            distance_x = grain9['center'][0] - grain0['center'][0]
            print(f"\n晶粒0与晶粒9在X方向的距离: {distance_x:.3f} Å")
            print(f"晶粒0中心: ({grain0['center'][0]:.3f}, {grain0['center'][1]:.3f}, {grain0['center'][2]:.3f}) Å")
            print(f"晶粒9中心: ({grain9['center'][0]:.3f}, {grain9['center'][1]:.3f}, {grain9['center'][2]:.3f}) Å")
        
        print("\n结构生成完成!")
        
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()