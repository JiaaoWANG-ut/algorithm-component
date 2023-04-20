from ase.io import read
import numpy as np

def read_cif(filename):
    """
    从CIF文件中读取结构并返回晶格信息
    """
    atoms = read(filename)
    lattice = atoms.get_cell()
    return np.array(lattice)

def calculate_box_volume(box_size):
    """
    计算盒子的体积
    """
    volume = np.abs(np.linalg.det(box_size))
    return volume

def generate_uniform_points(box_size, density):
    """
    在给定的盒子内生成均匀分布的点
    """
    # 计算盒子的体积和需要生成的点的数量
    volume = calculate_box_volume(box_size)
    num_points = int(volume * density)

    # 生成均匀分布的点
    points = np.random.rand(num_points, 3)
    points = np.dot(points, np.transpose(box_size))
    print("生成点数：",len(points))

    return points


# # 设置盒子的大小和密度
# box_size = read_cif("example.cif")
# density = 1

# # 生成均匀分布的点
# points = generate_uniform_points(box_size, density)

# # 打印输出点的数量和坐标
# print("生成了%d个点" % len(points))
# print("点的坐标为：")
# print(points)
