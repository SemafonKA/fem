import numpy as np
from matplotlib import pyplot as plt


class Point:
    x: float
    y: float

    def __init__(self, x: float = 0.0, y: float = 0.0):
        self.x = x
        self.y = y

    def __str__(self):
        return f"[{self.x}, {self.y}]"


class Mesh:
    material_num: int
    bl_point: int
    br_point: int
    tl_point: int
    tr_point: int

    def __init__(
        self, material: int = 0, bl: int = 0, br: int = 0, tl: int = 0, tr: int = 0
    ):
        self.material_num = material
        self.bl_point = bl
        self.br_point = br
        self.tl_point = tl
        self.tr_point = tr

    def __str__(self):
        return f"[{self.material_num}, {self.bl_point}, {self.br_point}, {self.tl_point}, {self.tr_point}]"


def plot_meshes(points: list[Point], meshes: list[Mesh]):
    _, axes = plt.subplots()
    colors = ["green", "blue", "red", "yellow", "black"]

    for mesh in meshes:
        # bottom edge
        axes.plot(
            [points[mesh.bl_point].x, points[mesh.br_point].x],
            [points[mesh.bl_point].y, points[mesh.br_point].y],
            color=colors[mesh.material_num - 1],
        )
        # top edge
        axes.plot(
            [points[mesh.tl_point].x, points[mesh.tr_point].x],
            [points[mesh.tl_point].y, points[mesh.tr_point].y],
            color=colors[mesh.material_num - 1],
        )
        # left edge
        axes.plot(
            [points[mesh.bl_point].x, points[mesh.tl_point].x],
            [points[mesh.bl_point].y, points[mesh.tl_point].y],
            color=colors[mesh.material_num - 1],
        )
        # right edge
        axes.plot(
            [points[mesh.br_point].x, points[mesh.tr_point].x],
            [points[mesh.br_point].y, points[mesh.tr_point].y],
            color=colors[mesh.material_num - 1],
        )
    plt.show()


def main():
    lines: list[str]
    with open("grid.txt", "r") as file:
        lines = file.readlines()

    [kx, ky] = [int(x) for x in lines[0].split()]
    points: list[Point] = []
    for j in range(ky):
        subr: list[Point] = []
        splt = lines[j + 2].split()
        for i in range(kx):
            p = Point(
                x=float(splt[2 * i]),
                y=float(splt[2 * i + 1]),
            )
            subr.append(p)
        points.extend(subr)

    km = int(lines[ky + 3].split()[0])
    meshes: list[Mesh] = []
    for j in range(km):
        splt = lines[j + ky + 5].split()
        m = Mesh(
            int(splt[0]),
            int(splt[1]),
            int(splt[2]),
            int(splt[3]),
            int(splt[4]),
        )
        meshes.append(m)

    plot_meshes(points, meshes)


if __name__ == "__main__":
    main()
