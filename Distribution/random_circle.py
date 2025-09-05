import numpy as np
import matplotlib.pyplot as plt

def random_positions_circle(num_devices, R):
    """
    円内に面積均等分布で端末を配置する
    """
    r = R * np.sqrt(np.random.rand(num_devices))  # 面積均等
    #r = R * np.random.rand(num_devices)
    theta = 2 * np.pi * np.random.rand(num_devices)  # 角度一様
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

# パラメータ
num_devices = 1000
R = 100000000  # 半径50m

# 端末位置生成
x, y = random_positions_circle(num_devices, R)

# 可視化
plt.figure(figsize=(6,6))
plt.scatter(x, y, s=10, color='blue', alpha=0.6)
plt.scatter(0, 0, color='red', marker='x', s=100, label='GW')
circle = plt.Circle((0,0), R, color='gray', fill=False, linestyle='--')
plt.gca().add_artist(circle)
plt.gca().set_aspect('equal', 'box')
plt.xlim(-R-5, R+5)
plt.ylim(-R-5, R+5)
plt.title(f'{num_devices} devices in circle of radius {R}')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.legend()
plt.show()
