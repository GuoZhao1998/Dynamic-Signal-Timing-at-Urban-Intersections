import numpy as np
from scipy.optimize import minimize

T = 2 # minute
Q0_1,Q0_5 = 19, 16
Q0_2,Q0_6 = 20, 18
Q0_3,Q0_7 = 3,  8
Q0_4,Q0_8 = 0,  0 
s1,s2,s3,s4,s5,s6,s7,s8 = [60]*8
f1 = [9,7,8]
f5 = [14,12,7]
f2 = [25,15,11]
f6 = [22,10,7]
f3 = [6,5,4]
f7 = [13,11,7]
f4 = [16,13,6]
f8 = [21,17,8]


def Tmin_15(k, Q1,Q5):
    return max([ Q1/(s1-f1[k]) , Q5/(s5-f5[k]) ])

def Tmin_26(k, Q2,Q6, Tmin_1_5):
    return max([ (Q2 + f2[k]*Tmin_1_5)/(s2-f2[k]) , (Q6 + f6[k]*Tmin_1_5)/(s6-f6[k]) ]) 

def Tmin_37(k, Q3,Q7, Tmin_1_5,Tmin_2_6):
    return max([ (Q3 + f3[k]*(Tmin_1_5+Tmin_2_6))/(s3-f3[k]), (Q7 + f7[k]*(Tmin_1_5+Tmin_2_6))/(s7-f7[k]) ]) 

def Tmin_48(k, Q4,Q8,Tmin_1_5,Tmin_2_6,Tmin_3_7):
    return max([ (Q4 + f4[k]*(Tmin_1_5+Tmin_2_6+Tmin_3_7))/(s4-f4[k]), (Q8 + f8[k]*(Tmin_1_5+Tmin_2_6+Tmin_3_7))/(s8-f8[k]) ]) 

def Theorem(k, Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8):
    Tmin_15_k = Tmin_15(k, Q1,Q5)
    Tmin_26_k = Tmin_26(k, Q2,Q6, Tmin_15_k)
    Tmin_37_k = Tmin_37(k, Q3,Q7, Tmin_15_k,Tmin_26_k)
    Tmin_48_k = Tmin_48(k, Q4,Q8, Tmin_15_k,Tmin_26_k,Tmin_37_k)
    if Tmin_15_k + Tmin_26_k + Tmin_37_k + Tmin_48_k > T:
        return True
    else:
        return False
    
def L15(lam, Q1, Q5, k):
    λ1,λ2,λ3,λ4 = lam
    L1 = Q1 + f1[k]*λ1*T - s1*λ1*T
    L5 = Q5 + f5[k]*λ1*T - s5*λ1*T
    return L1, L5

def L26(lam, Q2, Q6, k):
    λ1,λ2,λ3,λ4 = lam
    L2 = Q2 + f2[k]*(λ1+λ2)*T - s2*λ2*T
    L6 = Q6 + f6[k]*(λ1+λ2)*T - s6*λ2*T
    return L2, L6

def L37(lam, Q3, Q7, k):
    λ1,λ2,λ3,λ4 = lam
    L3 = Q3 + f3[k]*(λ1+λ2+λ3)*T - s3*λ3*T
    L7 = Q7 + f7[k]*(λ1+λ2+λ3)*T - s7*λ3*T
    return L3, L7

def L48(lam, Q4, Q8, k):
    λ1,λ2,λ3,λ4 = lam
    L4 = Q4 + f4[k]*T - s4*λ4*T
    L8 = Q8 + f8[k]*T - s8*λ4*T
    return L4, L8

def SO_objective(lam, Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8, k):
    L15_k = max(L15(lam,Q1,Q5,k))
    L26_k = max(L26(lam,Q2,Q6,k))
    L37_k = max(L37(lam,Q3,Q7,k))
    L48_k = max(L48(lam,Q4,Q8,k))
    return L15_k + L26_k + L37_k + L48_k

def cons_sum(lam):
    return lam.sum() - 1.0


# 初始化
Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8 = Q0_1,Q0_2,Q0_3,Q0_4,Q0_5,Q0_6,Q0_7,Q0_8

λset = []
for k in range(len(f1)):
    if Theorem(k, Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8):
        
        def obj(lam):  # 目标函数封装
            return SO_objective(lam, Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8, k)

        eps = 1e-5
        cons = [
            {'type':'ineq', 'fun': lambda x, Q1=Q1,Q5=Q5, k=k: max(L15(x,Q1,Q5,k)) - eps},
            {'type':'ineq', 'fun': lambda x, Q2=Q2,Q6=Q6, k=k: max(L26(x,Q2,Q6,k)) - eps},
            {'type':'ineq', 'fun': lambda x, Q3=Q3,Q7=Q7, k=k: max(L37(x,Q3,Q7,k)) - eps},
            {'type':'ineq', 'fun': lambda x, Q4=Q4,Q8=Q8, k=k: max(L48(x,Q4,Q8,k)) - eps},
            {'type':'eq',   'fun': lambda x: cons_sum(x)}
        ]
        # λ 的取值范围
        bounds = [(0.0,1.0)] * 4
        # 初始猜测
        x0 = np.array([0.25,0.25,0.25,0.25])

        res = minimize(obj, x0,method='SLSQP',bounds=bounds,constraints=cons,options={'ftol':1e-6})

        if not res.success:
            print(f"k={k} 时优化失败:", res.message)
            break

        lam_k = res.x.tolist()
        print(f"λ1[{k}] = {lam_k[0]:.4f},\t λ2[{k}] = {lam_k[1]:.4f},\t λ3[{k}] = {lam_k[2]:.4f},\t λ4[{k}] = {lam_k[3]:.4f}")
        print(f"g1[{k}] = {(lam_k[0]*2*60):.4f},\t g2[{k}] = {(lam_k[1]*2*60):.4f},\t g3[{k}] = {(lam_k[2]*2*60):.4f},\t g4[{k}] = {(lam_k[3]*2*60):.4f}")
        λset.append(lam_k)
        print(f"{(lam_k[0]*2*60+lam_k[1]*2*60+lam_k[2]*2*60+lam_k[3]*2*60):.4f}")
        print(f"Objective value = {res.fun}")
        # 计算实际的 Li_k 然后更新 Q1…Q8
        L1_k,L5_k = L15(lam_k, Q1,Q5, k)
        L2_k,L6_k = L26(lam_k, Q2,Q6, k)
        L3_k,L7_k = L37(lam_k, Q3,Q7, k)
        L4_k,L8_k = L48(lam_k, Q4,Q8, k)

        print(f'L1[{k}]={L1_k},\t L5[{k}]={L5_k},\t L2[{k}]={L2_k},\t L6[{k}]={L6_k}')
        print(f'L3[{k}]={L3_k},\t L7[{k}]={L7_k},\t L4[{k}]={L4_k},\t L8[{k}]={L8_k}\n')

        if L1_k > 0:
            Q1 = Q1 + f1[k]*T - s1*lam_k[0]*T
        else:
            Q1 = f1[k]*(lam_k[1]+lam_k[2]+lam_k[3])*T

        if L5_k > 0:
            Q5 = Q5 + f5[k]*T - s5*lam_k[0]*T
        else:
            Q5 = f5[k]*(lam_k[1]+lam_k[2]+lam_k[3])*T

        if L2_k > 0:
            Q2 = Q2 + f2[k]*T - s2*lam_k[1]*T
        else:
            Q2 = f2[k]*(lam_k[2]+lam_k[3])*T

        if L6_k > 0:
            Q6 = Q6 + f6[k]*T - s6*lam_k[1]*T
        else:
            Q6 = f6[k]*(lam_k[2]+lam_k[3])*T

        if L3_k > 0:
            Q3 = Q3 + f3[k]*T - s3*lam_k[2]*T
        else:
            Q3 = f3[k]*(lam_k[3])*T

        if L7_k > 0:
            Q7 = Q7 + f7[k]*T - s7*lam_k[2]*T
        else:
            Q7 = f7[k]*(lam_k[3])*T

        if L4_k > 0:
            Q4 = Q4 + f4[k]*T - s4*lam_k[3]*T
        else:
            Q4 = 0

        if L8_k > 0:
            Q8 = Q8 + f8[k]*T - s8*lam_k[3]*T
        else:
            Q8 = 0
        print(f'Q1[{k}]={Q1:.4f},\t Q5[{k}]={Q5:.4f},\t Q2[{k}]={Q2:.4f},\t Q6[{k}]={Q6:.4f}')
        print(f'Q3[{k}]={Q3:.4f},\t Q7[{k}]={Q7:.4f},\t Q4[{k}]={Q4:.4f},\t Q8[{k}]={Q8:.4f}\n')
    else:
        print(f"k = {k} This is a recoverable cycle")



for λi_k in λset:
    print(f"g1 = {(λi_k[0]*2*60):.4f},\t g2 = {(λi_k[1]*2*60):.4f},\t g3 = {(λi_k[2]*2*60):.4f},\t g4 = {(λi_k[3]*2*60):.4f}")





