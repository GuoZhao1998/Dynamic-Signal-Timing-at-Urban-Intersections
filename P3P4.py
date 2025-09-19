from sympy import *
import matplotlib.pyplot as plt 
from matplotlib import rcParams  
from matplotlib.patches import Rectangle


x0_1, x0_2 = 10, 11
f1 = [35,35,25,20]
f2 = [30,24,20,16]
u1, u2, T = 60, 60, 1
# s = [0.3847740265608951, 0.5807886454954859, 0.6177706554135172, 0.49002605348154915] # 我们的模型P3
s = [0.33515488, 0.36807507, 0.53660505, 0.64104805] # Chang的模型P4
print('----- 时间序列绘图(总量) -----')
def intersection(t0, t1, h0, h1, g0, g1):
    if t0 == t1:
        return False
    m1 = (h1 - h0) / (t1 - t0)    # 直线1的斜率
    m2 = (g1 - g0) / (t1 - t0)    # 直线2的斜率

    if m1 == m2:
        if h0 == g0:
            return False    # 重合
        return False  # 平行
    else:
        t_intersect = (m1*t0 - m2*t0 + g0 - h0) / (m1 - m2)
        if t0 <= t_intersect <= t1:
            N_intersect = m1 * (t_intersect - t0) + h0
            return t_intersect, N_intersect
        else: return False

def plot_1(x0_1,s,T):
    time,N_in,N_out = 0, x0_1, 0
    time_set1, in_set1, out_set1 = [time], [N_in], [N_out]

    for k in range(len(s)):
        time_0 = time
        Nin_0 = N_in
        Nout_0 = N_out

        time = time + s[k]*T
        N_in = N_in + f1[k]*s[k]*T
        N_out = N_out + u2*s[k]*T
        # print(f'k={k},\t L1[{k}]={N_in-N_out}')
        time_1 = time
        Nin_1 = N_in
        Nout_1 = N_out

        # print(f'time = {time:.6g}, \t N_in = {N_in:.6g}, \t N_out = {N_out:.6g}')
        point_01 = intersection(time_0, time_1, Nin_0, Nin_1, Nout_0, Nout_1)
        if point_01 == False:
            None
        else:
            time_set1.append(point_01[0])
            in_set1.append(point_01[1])
            out_set1.append(point_01[1])
            Nout_1 = Nin_1
            
        time_set1.append(time_1)
        in_set1.append(Nin_1)
        out_set1.append(Nout_1)

        time = time + (1-s[k])*T
        N_in = N_in + f1[k]*(1-s[k])*T
        N_out = Nout_1
        # print(f'time = {time:.6g}, \t N_in = {N_in:.6g}, \t N_out = {N_out:.6g}')

        time_set1.append(time)
        in_set1.append(N_in)
        out_set1.append(N_out)
    return time_set1, in_set1, out_set1



def plot_2(x0_2,s,T):
    time,N_in,N_out = 0, x0_2, 0
    time_set2, in_set2, out_set2 = [time], [N_in], [N_out]

    for k in range(len(s)):
        time = time + s[k]*T
        N_in = N_in + f2[k]*s[k]*T

        time_1 = time
        Nin_1 = N_in
        Nout_1 = N_out

        time_set2.append(time_1)
        in_set2.append(Nin_1)
        out_set2.append(Nout_1)
        
        time = time + (1-s[k])*T
        N_in = N_in + f2[k]*(1-s[k])*T
        N_out = N_out + u2*(1-s[k])*T
        # print(f'k={k},\t L2[{k}]={N_in-N_out}')
        time_2 = time
        Nin_2 = N_in
        Nout_2 = N_out

        point_12 = intersection(time_1, time_2, Nin_1, Nin_2, Nout_1, Nout_2)
        if point_12 != False:
            time_set2.append(point_12[0])
            in_set2.append(point_12[1])
            out_set2.append(point_12[1])
            N_out = Nin_2

        time_set2.append(time_2)
        in_set2.append(Nin_2)
        out_set2.append(N_out)

    return time_set2, in_set2, out_set2

time_set1, in_set1, out_set1 = plot_1(x0_1,s,T)
time_set2, in_set2, out_set2 = plot_2(x0_2,s,T)
print(out_set1)
print(out_set2)

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# 假设 time_set1, in_set1, out_set1, time_set2, in_set2, out_set2, s, T 已提前定义

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = True
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 14

fig, ax = plt.subplots(figsize=(12, 4))

# Road 1 进出曲线
ax.plot(time_set1, in_set1, 'b-', label='Road 1 Entries')
ax.plot(time_set1, out_set1, 'b--', label='Road 1 Exits')

# Road 1 灰色矩形（不同透明度方便区分）
for k in range(len(s)):
    rect1 = Rectangle((k*T+s[k]*T, 0), (1-s[k])*T, in_set1[-1], fill=True, color='0.7', alpha=0.4)
    ax.add_patch(rect1)



# Road 2 进出曲线
ax.plot(time_set2, in_set2, 'r-', label='Road 2 Entries')
ax.plot(time_set2, out_set2, 'r--', label='Road 2 Exits')


# 合并x刻度
custom_x_ticks = [0]
custom_x_labels = ["0"]
max_K = len(s)
for k in range(len(s)):
    tick1 = k*T + s[k]*T
    tick2 = (k+1)*T
    custom_x_ticks.extend([tick1, tick2])
    custom_x_labels.extend([f"{tick1:.4g}", f"{tick2:.4g}"])

ax.set_xticks(custom_x_ticks)
ax.set_xticklabels(custom_x_labels)
ax.set_xlabel('Time (minutes)', fontsize=16)
ax.set_ylabel('Number of vehicles', fontsize=16)

# y轴最大值取两条曲线的最大
ymax = max(in_set1[-1], in_set2[-1])
ax.set_xlim([0, max(time_set1[-1], time_set2[-1])])
ax.set_ylim([0, ymax])

# 图例去重复
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), loc='upper left', fontsize=12)

plt.tight_layout()
# plt.savefig('P4.png', dpi=600)
plt.show()

print('----- 结束绘图 -----')





