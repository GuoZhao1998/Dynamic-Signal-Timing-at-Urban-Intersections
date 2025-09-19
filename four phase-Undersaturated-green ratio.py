import numpy as np
from scipy.optimize import minimize

T,N=2,2; dim=4*N
f=np.array([[8,12,20],[24,19,34],[4,4,8],[12,14,18],
            [6,10,16],[20,17,22],[7,11,15],[10,15,24]])
s,Q0=np.full(8,60),np.full(8,10)
g=[0,1,2,3,0,1,2,3]  # "road number" corresponds to "Pattern mapping"

def split(x): return x.reshape(4,N)

def getQ(k,lam):
    # Calculate from Q_1[k] to Q_8[k]
    if k==0: return Q0
    l0,l1,l2,l3=lam[:,k-1]
    return np.array([f[0,k-1]*(l1+l2+l3)*T, f[1,k-1]*(l2+l3)*T, f[2,k-1]*l3*T, 0,
      f[4,k-1]*(l1+l2+l3)*T, f[5,k-1]*(l2+l3)*T, f[6,k-1]*l3*T, 0])

def Tc_all(k,Q,lam):
    # Tc_{i,i+4} = max (Tc_{i},Tc_{i+4})
    l0,l1,l2,l3=lam[:,k]
    return [
      max(Q[[0,4]]/(s[[0,4]]-f[[0,4],k])),
      max((Q[[1,5]]+f[[1,5],k]*l0*T)/(s[[1,5]]-f[[1,5],k])),
      max((Q[[2,6]]+f[[2,6],k]*(l0+l1)*T)/(s[[2,6]]-f[[2,6],k])),
      max((Q[[3,7]]+f[[3,7],k]*(l0+l1+l2)*T)/(s[[3,7]]-f[[3,7],k]))
    ]

def delay(i,k,lam,Tcs):
    l0,l1,l2,l3=lam[:,k]
    Qk,fi,si = getQ(k,lam)[i],f[i,k],s[i]
    Tc = Tcs[g[i]]
    # i=1,4 corresponds to mode 0, that is, the delay function d_{1,4}[k]
    if g[i]==0:
        return (Tc**2*(fi-si)+2*Tc*Qk+T**2*(fi+si)*(l0-1)**2)/(2*T*si)
    if g[i]==1:
        return (T*(fi+si)*(l0**2+(l2+l3)**2)/(2*si)+l0*(Tc*fi+Qk)/si+Tc**2*(fi-si)/(2*T*si)+Tc*Qk/(T*si))
    if g[i]==2:
        return (T*(fi+si)*((l3)**2+(l0+l1)**2)/(2*si)+(l0+l1)*(Tc*fi+Qk)/si
               +(Tc**2*(fi-si)+2*Tc*Qk)/(2*T*si))
    d=l0+l1+l2
    return ((fi+si)*d*T+Tc*(fi-si)+2*Qk)*(d*T+Tc)/(2*T*si)

def obj(x):
    lam=split(x); tot=0
    for k in range(N):
        Qk=getQ(k,lam)
        Tcs=Tc_all(k,Qk,lam)
        for i in range(8):
            tot += delay(i,k,lam,Tcs)*f[i,k]
    return tot

cons=[]
for k in range(N):
    for m in range(4):
        cons.append({'type':'ineq','fun': lambda x,k=k,m=m: split(x)[m,k]*T - Tc_all(k,getQ(k,split(x)),split(x))[m] })
    cons.append({'type':'eq','fun': lambda x,k=k: split(x)[:,k].sum() - 1})

bnds=[(0,1)]*dim
x0=np.full(dim,0.25)
res=minimize(obj,x0,method='SLSQP',bounds=bnds,constraints=cons,options={'ftol':1e-2,'eps':1e-2})

print("Optimization success:", res.success)
lam_opt=split(res.x)
for i in range(4):
    print(f"lambda {i+1} =", np.round(lam_opt[i],4))
print("Minimum objective value:", res.fun)
