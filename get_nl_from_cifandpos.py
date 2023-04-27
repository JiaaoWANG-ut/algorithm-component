##INPUT: cif,pos,cutoff
##OUTPUT: neighborlist element distance

from ase.io import read
from ase.neighborlist import NeighborList
import numpy as np

# 读取cif文件



def get_nlist_from_pos(structure,cutoff_radius,target_position):
    # 初始化NeighborList对象，用于确定原子的近邻
    cutoffs = [cutoff_radius] * len(structure)
    nl = NeighborList(cutoffs, skin=0.3, self_interaction=False)
    nl.update(structure)
    
    # 遍历结构中的所有原子
    for i, site in enumerate(structure):
        # 计算原子到目标位置的距离
        distance = np.linalg.norm(site.position - target_position)
        # 如果距离小于等于截断半径，说明该原子在截断半径内
        if distance <= cutoff_radius:
            # 确定原子的近邻
            indices, offsets = nl.get_neighbors(i)
            neighbors = structure[indices]
            # 计算每个近邻原子到中心点的距离，并将其存储在列表中
            distances = [np.linalg.norm(neighbor.position - target_position) for neighbor in neighbors]
            # 输出该原子及其近邻的元素种类和坐标信息以及到中心点的距离
            print(f"{site.symbol} at {site.position} (distance = {distance:.2f} Å)")
            



structure = read("2.cif")

# 定义截断半径和目标位置
cutoff_radius = 5.0  # 单位为Å


target_position = np.array([0, 0, 0])  # 用目标位置的实际数值代替x、y、z

get_nlist_from_pos(structure,cutoff_radius,target_position)
