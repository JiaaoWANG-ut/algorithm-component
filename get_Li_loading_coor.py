from ase.io import read
from ase.neighborlist import NeighborList
import numpy as np
from sadian import *
def find_possible_lithium_sites(target_position):
    """
    在CIF文件中查找可能的Li离子负载点
    :param cif_file: CIF文件名
    :param target_position: 目标位置的三维坐标
    :param cutoff_radius: 截断半径，单位为Å
    :return: 如果有且只有两个氧原子，并且这两个氧原子的距离在1.5-3A之间，返回True，否则返回False
    """
    # 读取CIF文件
    cif_file="example.cif"
    cutoff_radius=10
    structure = read(cif_file)

    # 初始化NeighborList对象，用于确定原子的近邻
    cutoffs = [cutoff_radius] * len(structure)
    nl = NeighborList(cutoffs, skin=0.3, self_interaction=False)
    nl.update(structure)

    o_count = 0
    distance_count = 0
    # 遍历目标位置周围的原子
    for i, site in enumerate(structure):
        distance = np.linalg.norm(site.position - target_position)
        if distance > 1.5 and distance <= 3.0:
            distance_count += 1
            # 如果该原子是氧原子
            if site.symbol == "O":
                o_count += 1
                # 计算该氧原子到目标位置的距离
                
                #print(distance)
                # 如果距离在1.5-3A之间，距离计数器加1


    # 如果有且只有两个氧原子，并且这两个氧原子的距离在1.5-3A之间，返回True
    if o_count == 2 and distance_count == 2:
        print(target_position)
        return True
    else:
        return False



box_size = read_cif("example.cif")
density = 10

# 生成均匀分布的点
target_positions = generate_uniform_points(box_size, density)


#遍历目标位置列表，依次调用find_possible_lithium_sites函数，并打印输出结果

#for pos in target_positions:
#    if find_possible_lithium_sites(pos):
#        #print("yes")
#        print(pos)
#    #else:
#        #print("no")




from tqdm import tqdm
import multiprocessing as mul
import sys

CORES = 4
pool = mul.Pool(CORES)    
arg = target_positions

# 使用tqdm创建进度条，并设置更新速率为每秒一次
with tqdm(total=len(arg), miniters=1) as pbar:
    # 在pool.map()函数中传递update函数作为回调函数
    def update(*a):
        pbar.update()
        sys.stdout.flush()   # 强制刷新终端

    rel = pool.map_async(find_possible_lithium_sites,arg,callback=update).get()
