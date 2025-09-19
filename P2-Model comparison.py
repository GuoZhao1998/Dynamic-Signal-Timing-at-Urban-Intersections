import numpy as np
from scipy.optimize import minimize_scalar, root_scalar

f1 = 10/60
f2 = 14/60
s1 = 44/60
s2 = 44/60
C = 1.5*60

def dw(λi,xi,fi,C):
    term1_1 = C * (1 - λi)**2 / (2 * (1 - λi*xi))
    term2_1 = xi**2 / (2*fi*(1 - xi))
    term3_1 = 0.65 * (C / fi**2)**(1/3) * xi**(2 + 5*λi)
    return (term1_1 + term2_1 - term3_1)

def webster(λ1,f1,f2,C):
    x1 = f1/s1
    x2 = f2/s2
    λ2 = 1 - λ1
    return dw(λ1,x1,f1,C) + dw(λ2,x2,f2,C)

def ROSCA_delay(λ1,f1,f2,C):
    λ2 = 1 - λ1
    x1 = f1/s1
    x2 = f2/s2
    num = s1*λ1*x1 * dw(λ1,x1,f1,C) + s2*λ2*x2 * dw(λ2,x2,f2,C)
    den = s1*λ1*x1 + s2*λ2*x2
    return num / den

Q0_1, Q0_2 = 0, 0

def dC1(λ1):  
    return ( C*(f1+s1)*(λ1-1)**2 )/ ( 2*s1 ) + (Q0_1**2)/(2*C*s1*(s1 - f1))

def dC2(λ1):  
    return (C*λ1*s2 + Q0_2)**2/(2*C*s2*(s2 - f2))



x_axis = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]# [i/10 for i in range(0,11)]
y1_axis = []
y2_axis = []
y3_axis = []
y4_axis = [ 34.22, 32.57, 31.08, 30.85, 30.18, 31.09, 31.75, 32.27, 38.73]

for λ1 in x_axis:
    x1 = f1/s1
    x2 = f2/s2
    λ2 = 1 - λ1

    y1_w = dw(λ1,x1,f1,C) + dw(λ2,x2,f2,C)
    y2_R = ( s1*λ1*x1*dw(λ1,x1,f1,C) + s2*λ2*x2*dw(λ2,x2,f2,C) )/ (s1*λ1*x1 + s2*λ2*x2)
    y3_O = dC1(λ1)+dC2(λ1)

    y1_axis.append( round(y1_w,2) )
    y2_axis.append( round(y2_R,2) )
    y3_axis.append( round(y3_O,2) )


print(f"\lambda = {x_axis}")
print(f"webster = {y1_axis}")
print(f"ROSCA  = {y2_axis}")
print(f"our  = {y3_axis}")
print(f"Vissim  = {y4_axis}")

import matplotlib.pyplot as plt
from sympy import *
import matplotlib.pyplot as plt 
from matplotlib import rcParams  
from matplotlib.patches import Rectangle










rcParams['font.family'] = 'Times New Roman'  
rcParams['text.usetex'] = True  
rcParams['axes.titlesize'] = 14      # 标题字号  
rcParams['axes.labelsize'] = 14      # 坐标轴标签字号  
rcParams['xtick.labelsize'] = 12     # x轴刻度字号  
rcParams['ytick.labelsize'] = 12     # y轴刻度字号  
rcParams['legend.fontsize'] = 12     # 图例字号  

# 创建画布和坐标轴
plt.figure(figsize=(10, 4))

# 绘制三条折线
plt.plot(x_axis, y1_axis, label='webster', marker='o', linestyle='-', color='b')
plt.plot(x_axis, y2_axis, label='ROSCA', marker='s', linestyle='--', color='g')
plt.plot(x_axis, y3_axis, label='Our model', marker='D', linestyle=':', color='r')
plt.plot(x_axis, y4_axis, label='Vissim Simulation', marker='D', linestyle=':', color='k')
# 添加图表元素
plt.xlabel('Green Time Ratio', fontsize=12)
plt.ylabel('Delay Time at Intersection (Seconds)', fontsize=12)
plt.xticks(x_axis)  # 确保显示所有x轴刻度
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(loc='upper center')

# 保存并显示图表
plt.tight_layout()
# plt.savefig('P2.png', dpi=600)
plt.show()






