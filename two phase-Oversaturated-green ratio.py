from sympy import *
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import minimize
from matplotlib.patches import Rectangle


x0_1, x0_2 = 10, 11
f1 = [35,35,25,20]
f2 = [30,24,20,16]
u1, u2, T = 60, 60, 1

def L1(s,k):
    sum_s = 0
    for j in range(k+1):
        sum_s += s[j]
    sum_f = 0
    for j in range(k):
        sum_f += f1[j]
    return x0_1 + f1[k]*s[k]*T - u1*T*sum_s + T*sum_f

def L2(s,k):
    sum_ = 0
    for j in range(k+1):
        sum_ += f2[j] + u2*s[j] - u2
    return x0_2 + T*sum_

def dC1(s1,k):  
    if k == 0:  x1_k = x0_1
    else: x1_k = f1[k-1] * (1-s1[k-1]) * T
    return ( T*(f1[k]+u1)*(s1[k]-1)**2 )/ ( 2*u1 ) + (x1_k**2)/(2*T*u1*(u1 - f1[k]))

def dC2(s1,k):  
    if k == 0:  x2_k = x0_2
    else: x2_k = 0
    return (T*s1[k]*u2 + x2_k)**2/(2*T*u2*(u2 - f2[k]))

def SO_objective(s,n):
    sum_x = 0
    for k in range(0,n+1):
        sum_x += L1(s,k)**2 + L2(s,k)**2 
    return sum_x + dC1(s,n+1) + dC2(s,n+1)

def upper_level(n):
    Initial_guess_s1 = [0.365508758772233, 0.613937319661386, 0.518018121791249, 0.494382022471910]
    bo2 = [(0, 1)] *(n+2)

    constraints = []
    for k in range(n+1):
        constraints.append({'type': 'ineq', 'fun': lambda s, k=k: L1(s, k)})
        constraints.append({'type': 'ineq', 'fun': lambda s, k=k: L2(s, k)})
    
    res_s1 = minimize(SO_objective, Initial_guess_s1,args=(n), method='SLSQP', bounds=bo2,constraints=constraints, options={ 'ftol': 1e-4, 'eps': 1e-2})
    if not res_s1.success:
        print(f"upper level Optimization failed: {res_s1.message}")
        return np.inf
    return res_s1


n = 2
s_set = upper_level(n).x
s_set = s_set.tolist()
print( f'\n[0,{n}]与[{n},{n+1}]周期内的最优绿信比: {s_set}' )




def x_1(s,k):
    x1_next = x0_1
    for j in range(k):
        x1_next = x1_next + f1[j]*T - u1*s[j]*T
    return x1_next

def x_2(s,k):
    x2_next = x0_2
    for j in range(k):
        x2_next = x2_next + f2[j]*T - u2*(1-s[j])*T
    return x2_next

# for k in range(len(s_set)):
#     print( f'x_1[{k}] = {x_1(s_set,k):.6g}, \t x_2[{k}] = {x_2(s_set,k):.6g}' )
#     print( f'L_1[{k}] = {L1(s_set,k):.6g}, \t L_2[{k}] = {L2(s_set,k):.6g}' )


# print( f'\n----- 符号表达式 -----' )
# s = symbols('s_{1}[0:11]')
# # n = 0
# print(f"\n目标函数: {latex(simplify( SO_objective(s,n) ))}")

# for k in range(len(s_set)):
#     if k <= 2:
#         print(f"\n限制条件: L_1[{k}] = {latex(simplify( L1(s,k) ))} > 0")
#         print(f"\n限制条件: L_2[{k}] = {latex(simplify( L2(s,k) ))} > 0")
#     else:
#         print(f"\n限制条件: 左端: {latex(simplify( x_1(s,k)/(T*(u1-f1[k])) ))} ")
#         print(f"\n限制条件: 右端: {latex(simplify( (T*(u2-f2[k])-x_2(s,k))/(T*u2)  ))}")


for s_k in s_set:
    print(round(s_k,4))










###### 不含可恢复周期 ######

# x0_1, x0_2 = 10, 11
# f1 = [41, 60, 65, 29, 14, 19]
# f2 = [28, 20, 15, 10, 10, 12]
# u1, u2, T = 60, 60, 1

# def L1(s, k):
#     sum_s = sum(s[:k+1])
#     sum_f = sum(f1[:k])
#     return x0_1 + f1[k] * s[k] * T - u1 * T * sum_s + T * sum_f

# def L2(s, k):
#     sum_ = sum(f2[j] + u2 * s[j] - u2 for j in range(k+1))
#     return x0_2 + T * sum_

# def dC1(s1,k):  
#     if k == 0:  x1_k = x0_1
#     else: x1_k = f1[k-1] * (1-s1[k-1]) * T
#     eq_1 = x1_k**2/(2*T*u1*(u1 - f1[k]))
#     eq_2 = T**2*(u1 + f1[k])*(s1[k]**3 - 3*s1[k]**2 + 3*s1[k] - 1) /(2*u1)
#     return eq_1 - eq_2 

# def dC2(s1,k):  
#     if k == 0:  x2_k = x0_2
#     else: x2_k = 0
#     return (T*s1[k]*u2 + x2_k)**2/(2*T*u2*(u2 - f2[k]))

# def SO_objective(s, n):
#     sum_x = sum(L1(s, k)**2 + L2(s, k)**2 for k in range(n+1))
#     return sum_x

# def upper_level(n):
#     Initial_guess_s1 = [0.5] * (n+1)
#     bo2 = [(0, 1)] * (n+1)

#     # Define constraints
#     constraints = []
#     for k in range(n+1):
#         constraints.append({'type': 'ineq', 'fun': lambda s, k=k: L1(s, k)})
#         constraints.append({'type': 'ineq', 'fun': lambda s, k=k: L2(s, k)})

#     res_s1 = minimize(SO_objective, Initial_guess_s1, args=(n), method='SLSQP', bounds=bo2, constraints=constraints, options={'ftol': 1e-4, 'eps': 1e-2})
    
#     if not res_s1.success:
#         print(f"Upper level optimization failed: {res_s1.message}")
#         return float('inf')
    
#     return res_s1


# n = 4
# s_set = upper_level(n).x
# s_set = s_set.tolist()
# print( f'\n[0,{n}]与[{n},{n+1}]周期内的最优绿信比: {s_set}' )

