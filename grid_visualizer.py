import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon


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
    colors = ["#ffc107", "#1976d2", "#FF6347", "#FFD700", "#F5F5DC"]

    blpoint = points[0]
    trpoint = points[-1]
    axes.set_xlim(blpoint.x, trpoint.x)
    axes.set_ylim(blpoint.y, trpoint.y)

    for mesh in meshes:
        bl_point = points[mesh.bl_point]
        br_point = points[mesh.br_point]
        tl_point = points[mesh.tl_point]
        tr_point = points[mesh.tr_point]

        axes.add_patch(
            Polygon(
                [
                    (bl_point.x, bl_point.y),
                    (br_point.x, br_point.y),
                    (tr_point.x, tr_point.y),
                    (tl_point.x, tl_point.y),
                ],
                facecolor=colors[mesh.material_num - 1],
                edgecolor="black",
            )
        )

    plt.show()


def main():
    lines: list[str]
    with open("lab1/grid.txt", "r") as file:
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
