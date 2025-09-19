import numpy as np
from scipy.optimize import minimize
from sympy import *



x0_1, x0_2 = 6, 8
f1 = [10.67, 10, 38]
f2 = [14, 13, 17]
u1, u2, T = 60, 60, 1.5

def dC1(s1,k):  
    if k == 0:  x1_k = x0_1
    else: x1_k = f1[k-1] * (1-s1[k-1]) * T
    return ( T*(f1[k]+u1)*(s1[k]-1)**2 )/ ( 2*u1 ) + (x1_k**2)/(2*T*u1*(u1 - f1[k]))

def dC2(s1,k):  
    if k == 0:  x2_k = x0_2
    else: x2_k = 0
    return (T*s1[k]*u2 + x2_k)**2/(2*T*u2*(u2 - f2[k]))

def SO_objective(s1):
    sum_total = 0
    for k in range(len(f1)):  # 可解的周期
        sum_k = dC1(s1,k)*f1[k] + dC2(s1,k)*f2[k]
        sum_total = sum_total + sum_k
    return sum_total 

def con0_r1(s):
    return s[0] - x0_1 / (T*(u1 - f1[0]))

def con0_r2(s):
    return ((T * (u2 - f2[0]) - x0_2) / (T * u2)) - s[0]

def con1_r1(s):
    return s[1] - (f1[0]*(1-s[0])*T) / (T*(u1 - f1[1]))

def con1_r2(s):
    return  ((T*(u2 - f2[1])) / (T*u2)) - s[1]

def con2_r1(s):
    return s[2] - (f1[1]*(1-s[1])*T) / (T*(u1 - f1[2]))

def con2_r2(s):
    return ((T*(u2 - f2[2])) / (T*u2)) - s[2]

constraints = [
    {'type': 'ineq', 'fun': con0_r1},
    {'type': 'ineq', 'fun': con0_r2},
    {'type': 'ineq', 'fun': con1_r1},
    {'type': 'ineq', 'fun': con1_r2},
    {'type': 'ineq', 'fun': con2_r1},
    {'type': 'ineq', 'fun': con2_r2}
]

def upper_level():
    Initial_guess_s1 = [0.3626,0.4673,0.716667]
    bo2 = [(0, 1), (0, 1), (0, 1)] 
    res_s1 = minimize(SO_objective, Initial_guess_s1, method='SLSQP', bounds=bo2,constraints=constraints, options={'maxiter': 5000, 'ftol': 1e-6, 'eps': 1e-6})
    if not res_s1.success:
        print(f"upper level Optimization failed: {res_s1.message}")
        return np.inf
    return res_s1

optimal_s = upper_level().x.tolist()
s = []
for s_i in optimal_s:
    s.append(round(s_i,9))
print( s )
print( upper_level().fun )

for k in range(len(s)+1):
    if k == 0:  x1_k = x0_1; x2_k = x0_2
    else: x1_k = f1[k-1] * (1-s[k-1]) * T; x2_k = 0
    print(f'x_1[{k}] = {x1_k},\t x_2[k] = {x2_k}')



# 符号模型
s0,s1,s2 = symbols('s_1[0] s_1[1] s_1[2]',real=True,nonzero=True,positive=True)

s = [s0,s1,s2]

f3 = f1
f4 = f2
total = 0
for k in range(len(f1)):
    integral_1 =  dC1(s,k)*f1[k]
    integral_2 =  dC2(s,k)*f2[k]
    sum_k = integral_1 + integral_2
    # print(sum_k)
    total = total + sum_k

print(f"\n目标函数: {simplify(total)}")

print(f'限制条件: {x0_1 / (T*(u1 - f1[0]))} <= s_1[0] <= {((T * (u2 - f2[0]) - x0_2) / (T * u2))}')
print(f'限制条件: {(f1[0]*(1-s[0])*T) / (T*(u1 - f1[1]))} <= s_1[1] <= {((T*(u2 - f2[1])) / (T*u2))}')
print(f'限制条件: {(f1[1]*(1-s[1])*T) / (T*(u1 - f1[2]))} <= s_1[2] <= {((T*(u2 - f2[2])) / (T*u2))}')





# print('----- 三维绘图 -----')
# ### 绘图 ###
# lo = []
# for _1 in range(11):
#     lo.append(_1/10)

# s0_set, s1_set, s2_set, value_set = [],[],[],[]
# for s0 in lo:
#     if x0_1 / (T*(u1 - f1[0])) <= s0 <= ((T * (u2 - f2[0]) - x0_2) / (T * u2)):
#         for s1 in lo:
#             if (f1[0]*(1-s1)*T) / (T*(u1 - f1[1])) <= s1 <= ((T*(u2 - f2[1])) / (T*u2)):
#                 for s2 in lo:
#                     if (f1[1]*(1-s1)*T) / (T*(u1 - f1[2])) <= s2 <= ((T*(u2 - f2[2])) / (T*u2)):
#                         s0_set.append(s0)
#                         s1_set.append(s1)
#                         s2_set.append(s2)
#                         value_set.append( SO_objective([s0,s1,s2]) )
            

# min_value = min(value_set)
# min_value_index = value_set.index(min_value)
# print(f'min_z^so = {min_value}')
# print(f's[0] = {s0_set[min_value_index]}')
# print(f's[1] = {s1_set[min_value_index]}')
# print(f's[2] = {s2_set[min_value_index]}')

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import rcParams  
# from matplotlib.patches import Rectangle


# # 设置全局字体和样式参数
# rcParams['font.family'] = 'Times New Roman'  
# rcParams['text.usetex'] = True  
# rcParams['axes.titlesize'] = 14      # 标题字号  
# rcParams['axes.labelsize'] = 14      # 坐标轴标签字号  
# rcParams['xtick.labelsize'] = 12     # x轴刻度字号  
# rcParams['ytick.labelsize'] = 12     # y轴刻度字号  
# rcParams['legend.fontsize'] = 14     # 图例字号  


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# sc = ax.scatter(s0_set, s1_set, s2_set, c=value_set, cmap='coolwarm')

# # 添加颜色条
# plt.colorbar(sc, label='$Z^{so}$',pad=0.05,ticks=[min(value_set), (min(value_set)+max(value_set))/2, max(value_set)],shrink=0.8)
# # pad影响与scatter绘图结果之间的距离, ticks颜色条刻度范围

# ax.set_xlabel('$s_1[0]$')
# ax.set_ylabel('$s_1[1]$')
# ax.set_zlabel('$s_1[2]$')
# # ax.view_init(elev=30, azim=30)
# # plt.savefig('P3.png', dpi=900)
# ax.view_init(elev=-25, azim=-125)       # 可修改参数以获取不同的视角
# #plt.savefig('P4.png', dpi=900)
# plt.show()

# print('----- 结束绘图 -----')



