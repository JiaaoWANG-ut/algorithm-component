from ase.io import read
from ase.neighborlist import NeighborList
import numpy as np

def get_neighbors(filename, target_position, cutoff_radius):
    # 读取结构文件
    structure = read(filename)

    # 初始化NeighborList对象，用于确定原子的近邻
    cutoffs = [cutoff_radius] * len(structure)
    nl = NeighborList(cutoffs, skin=0.3, self_interaction=False)
    nl.update(structure)

    # 存储每个原子的近邻及其到目标位置的距离
    result = []

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
            # 将该原子及其近邻的元素种类、坐标信息以及到目标位置的距离添加到结果列表中
            for neighbor, distance, offset in zip(neighbors, distances, offsets):
                result.append((site.symbol, site.position, neighbor.symbol, neighbor.position + np.dot(offset, structure.cell), distance))
    
    return result

result = get_neighbors('2.cif', np.array([0, 0, 0]), 10.0)

for site_symbol, site_position, neighbor_symbol, neighbor_position, distance in result:
    print(f'{site_symbol} at {site_position} has neighbor {neighbor_symbol} at {neighbor_position} (distance = {distance:.2f} Å)')


