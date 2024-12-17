import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


class Solution:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray

    def __init__(self, x: np.ndarray, y: np.ndarray, z: np.ndarray):
        self.x = x
        self.y = y
        self.z = z


def plot_graphics(solution: Solution):
    xi = np.linspace(min(solution.x), max(solution.x))
    yi = np.linspace(min(solution.y), max(solution.y))
    zi = griddata(
        (solution.x, solution.y),
        solution.z,
        (xi[None, :], yi[:, None]),
        method="linear",
    )
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(xi, yi, zi, shading="auto")
    plt.colorbar(label="f(x,y)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Heatmap of f(x,y)")

    plt.contour(xi, yi, zi, colors="black", linewidths=0.5, levels=14)


def main():
    lines: list[str]
    with open("fem_lib/io/solution.txt", "r") as file:
        lines = file.readlines()

    x = np.zeros((0,))
    y = np.zeros((0,))
    z = np.zeros((0,))  # Посчитанная функция

    for line in lines:
        splt = line.split()
        x = np.append(x, [float(splt[0])])
        y = np.append(y, [float(splt[1])])
        z = np.append(z, [float(splt[2])])

    # Дальше работает автоматика

    spline = Solution(x, y, z)
    plot_graphics(spline)

    plt.show()


if __name__ == "__main__":
    main()
