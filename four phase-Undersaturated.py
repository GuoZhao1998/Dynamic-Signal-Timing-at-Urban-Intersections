
Q0_1,Q0_2,Q0_3,Q0_4,Q0_5,Q0_6,Q0_7,Q0_8 = [10]*8
s1,s2,s3,s4,s5,s6,s7,s8 = [60]*8
f1 = [8,12,20]
f2 = [24,19,34]
f3 = [4,4,8]
f4 = [12,14,18]
f5 = [6,10,16]
f6 = [20,17,22]
f7 = [7,11,15]
f8 = [10,15,24]
C = 2
def Tmin_15(k, Q1_k,Q5_k):
    return max([ Q1_k/(s1-f1[k]) , Q5_k/(s5-f5[k]) ])

def Tmin_26(k, Q2_k,Q6_k, Tmin_1_5):
    return max([(Q2_k + f2[k]*Tmin_1_5)/(s2-f2[k]), (Q6_k + f6[k]*Tmin_1_5)/(s6-f6[k])]) 

def Tmin_37(k, Q3_k,Q7_k, Tmin_1_5,Tmin_2_6):
    return max([(Q3_k + f3[k]*(Tmin_1_5+Tmin_2_6))/(s3-f3[k]), (Q7_k + f7[k]*(Tmin_1_5+Tmin_2_6))/(s7-f7[k])]) 

def Tmin_48(k, Q4_k,Q8_k,Tmin_1_5,Tmin_2_6,Tmin_3_7):
    return max([(Q4_k + f4[k]*(Tmin_1_5+Tmin_2_6+Tmin_3_7))/(s4-f4[k]), (Q8_k + f8[k]*(Tmin_1_5+Tmin_2_6+Tmin_3_7))/(s8-f8[k])]) 

def Qmin_i_k(k, T_15,T_26,T_37,T_48):
    if k==0: 
        Q1_k, Q2_k, Q3_k, Q4_k, Q5_k, Q6_k, Q7_k, Q8_k = Q0_1,Q0_2,Q0_3,Q0_4,Q0_5,Q0_6,Q0_7,Q0_8
    else: 
        Q1_k = f1[k-1]*(T_26 + T_37 + T_48)
        Q5_k = f5[k-1]*(T_26 + T_37 + T_48)
        Q2_k = f2[k-1]*(T_37 + T_48)
        Q6_k = f6[k-1]*(T_37 + T_48)
        Q3_k = f3[k-1]*(T_48)
        Q7_k = f7[k-1]*(T_48)
        Q4_k = 0
        Q8_k = 0
    return Q1_k, Q2_k, Q3_k, Q4_k, Q5_k, Q6_k, Q7_k, Q8_k

n = len(f1)
Q1_k, Q2_k, Q3_k, Q4_k, Q5_k, Q6_k, Q7_k, Q8_k = Q0_1,Q0_2,Q0_3,Q0_4,Q0_5,Q0_6,Q0_7,Q0_8
for k in range(n):
    Tmin_15k = Tmin_15(k, Q1_k,Q5_k)
    Tmin_26k = Tmin_26(k, Q2_k,Q6_k, Tmin_15k)
    Tmin_37k = Tmin_37(k, Q3_k,Q7_k, Tmin_15k,Tmin_26k)
    Tmin_48k = Tmin_48(k, Q4_k,Q8_k, Tmin_15k,Tmin_26k,Tmin_37k)
    Q1_k, Q2_k, Q3_k, Q4_k, Q5_k, Q6_k, Q7_k, Q8_k = Qmin_i_k(k+1, Tmin_15k,Tmin_26k,Tmin_37k,Tmin_48k)

    if (Tmin_15k + Tmin_26k + Tmin_37k + Tmin_48k) <= C:
        print(f'\nThe {k}-th cycle is a recoverable cycle.')
    else: 
        print(f'\nThe {k}th cycle is a unrecoverable cycle.')
