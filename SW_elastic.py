# -*- coding:utf-8 -*-
# datetime:2019/6/3 0003 15:09
# software: PyCharm

import cmath
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time


def generate(mx, n):       # mx:为列表形式的材料参数, n:为矩阵纬度
    rst = np.zeros([n, n], dtype=complex)
    for i in range(n):
        for j in range(n):
            g = G[i, :]-G[j, :]
            g = g[0]
            if abs(g) < 1e-5:
                rst[i, j] = mx[0]*pf + mx[1]*(1-pf)
            else:
                b1 = ii*(mx[0]+mx[1])
                b2 = np.exp(-ii*g*pf*a)-1
                rst[i, j] = np.mat(b1*b2/(g*a))
    return rst


def fun1(b, kk):     # 组装
    AA = -C44
    BB = np.zeros([NG, NG], dtype=complex)
    CC = np.zeros([NG, NG], dtype=complex)

    for i in range(NG):
        for j in range(NG):
            kGi = kk + G[i, :]
            kGj = kk + G[j, :]
            CC[i, j] = -C44[i, j] * kGi * kGj + LOU[i, j] * b ** 2
    UVA1 = np.hstack((BB, np.eye(NG, NG)))
    UVA2 = np.hstack((-CC, -BB))
    UVA = np.vstack((UVA1, UVA2))
    UVB1 = np.hstack((np.eye(NG, NG), BB))
    UVB2 = np.hstack((BB, AA))
    UVB = np.vstack((UVB1, UVB2))
    UVBF = np.mat(UVB).I
    res = np.dot(UVBF, UVA)
    k3, eek = np.linalg.eig(res)
    k3 = list(k3)
    eek1 = eek[0:NG, :]     # AG1
    # print(k3, '\n', eek1)
    # print(np.size(eek1))
    HH1 = np.zeros([NG, 2 * NG], dtype=complex)
    HH2 = np.zeros([NG, 2 * NG], dtype=complex)

    for i in range(0, 2 * NG):
        HH1[:, i:i + 1] = np.dot(C44, eek1[:, i]) * (ii * k3[i]) * np.exp(ii * k3[i] * h)
        HH2[:, i:i + 1] = np.dot(C44, eek1[:, i]) * (ii * k3[i]) * np.exp(-ii * k3[i] * h)
    HH = np.vstack((HH1, HH2))
    rst = np.linalg.det(HH)
    return rst


def switch2(b1,bstep,b2,e,n):
    bp1 = b1
    count = 1
    ANS = []
    while bp1 <= b2 and count <= n:
        bp4 = bp1 + bstep
        bp2 = bp4 - 0.618 * bstep
        bp3 = bp1 + 0.618 * bstep
        AA1 = abs(fun1(bp1,k))
        AA2 = abs(fun1(bp2,k))
        AA3 = abs(fun1(bp3,k))
        AA4 = abs(fun1(bp4,k))
        if (AA2 < AA1) & (AA2 < AA3):
            while abs(bp2-bp3) > e:
                if AA2 < AA3:
                    bp4 = bp3
                    bp3 = bp2
                    bp2 = bp1 + (1 - 0.618) * (bp4 - bp1)
                    AA2 = fun1(bp2,k)
                    AA3 = fun1(bp3,k)
                    #print(bp2,bp3)
                else:
                    bp1 = bp2
                    bp2 = bp3
                    bp3 = bp4 - (1 - 0.618) * (bp4 - bp1)
                    AA2 = fun1(bp2,k)
                    AA3 = fun1(bp3,k)
                    #print(bp2,bp3)
            ans1 = (bp2+bp3)/2
            print(ans1)
            ANS.append(ans1)
            count += 1
            bp1 = ans1 + bstep

        elif (AA3 < AA2) & (AA3 < AA4):
            while abs(bp2-bp3) > e:
                if AA2 < AA3:
                    bp4 = bp3
                    bp3 = bp2
                    bp2 = bp1 + (1 - 0.618) * (bp4 - bp1)
                    AA2 = fun1(bp2,k)
                    AA3 = fun1(bp3,k)
                    #print(bp2,bp3)
                else:
                    bp1 = bp2
                    bp2 = bp3
                    bp3 = bp4 - (1 - 0.618) * (bp4 - bp1)
                    AA2 = fun1(bp2,k)
                    AA3 = fun1(bp3,k)
                    #print(bp2,bp3)
            ans = (bp2+bp3)/2
            print(ans)
            ANS.append(ans)
            count += 1
            bp1 = ans + bstep

        elif (AA1<=AA2)&(AA2<=AA3)&(AA3<=AA4):
            bp1 = bp4
        elif (AA2 >= AA3) & (AA3 >= AA4):
            bp1 = bp3
        elif (AA1<=AA2)&(AA2<=AA3)&(AA3>AA4):
            bp1 = bp3
        else:
            bp1 = bp4
    ft = []
    lh = len(ANS)
    for each in ANS:
        ft.append(abs(fun1(each,k)))
    tj = sum(ft)/lh
    for each in ft:
        if each >= tj*8:
            num = ft.index(each)
            ANS[num] = 0
    return ANS


l1 = 0.01
l2 = 0.01
a = l1 + l2       # 晶格
pf = l1/a         # 填充率
c44 = [1.05e+10, 1.05e+10]
lou = [1980, 1980]
nrsquare = 5
ra1 = 2*np.pi/a
ii = cmath.sqrt(-1)
h = 5e-3
b = np.pi/a
kk = np.arange(b*1e-9,b+0.01,b/40)

###################
NG = 2*nrsquare+1
G = np.zeros([NG, 1])
count1 = 0
for l in np.arange(-nrsquare, nrsquare+1):
    G[count1, :] = l*ra1
    count1 += 1
# print(G)
C44 = generate(c44, NG)
LOU = generate(lou, NG)


####################

result1 = []
result2 = []
result3 = []
result4 = []
result5 = []
result6 = []
result7 = []
result8 = []
result9 = []
result0 = []
for k in kk:
    print('k='+str(k))
    a = switch2(1e+3, 1e+4, 15e+7, 1e+2, 6)
    try:
        result0.append(a[0])
    except Exception as e:
        print(e)
        result0.append(0)
    try:
        result1.append(a[1])
    except Exception as e:
        print(e)
        result1.append(0)
    try:
        result2.append(a[2])
    except Exception as e:
        print(e)
        result2.append(0)
    try:
        result3.append(a[3])
    except Exception as e:
        print(e)
        result3.append(0)
    try:
        result4.append(a[4])
    except Exception as e:
        print(e)
        result4.append(0)
    try:
        result5.append(a[5])
    except Exception as e:
        print(e)
        result5.append(0)
    try:
        result6.append(a[6])
    except Exception as e:
        print(e)
        result6.append(0)
    try:
        result7.append(a[7])
    except Exception as e:
        print(e)
        result7.append(0)
    try:
        result8.append(a[8])
    except Exception as e:
        print(e)
        result8.append(0)
    try:
        result9.append(a[9])
    except Exception as e:
        print(e)
        result9.append(0)



result123 = []
result123.append(result0)
result123.append(result1)
result123.append(result2)
result123.append(result3)
result123.append(result4)
result123.append(result5)
result123.append(result6)
result123.append(result7)
result123.append(result8)
result123.append(result9)
result = np.array(result123)
result = pd.DataFrame(result)
localtime = time.strftime("%Y-%m-%d %H-%M-%S", time.localtime())
ltime = str(localtime)
result.to_excel('./data/PWE(elastic)%s [H=%s, n=%s, l1=%s, l2=%s].xlsx' % (ltime, h, nrsquare, l1, l2), index=False)

x = np.arange(0, 1.001, 1/40)
plt.scatter(x, result0,s=15)
plt.scatter(x, result1,s=15)
plt.scatter(x, result2,s=15)
plt.scatter(x, result3,s=15)
plt.scatter(x, result4,s=15)
plt.scatter(x, result5,s=15)
plt.scatter(x, result6,s=15)
plt.scatter(x, result7,s=15)
plt.scatter(x, result8,s=15)
plt.scatter(x, result9,s=15)
# plt.xlim((1, 7))
# my_y_ticks = np.arange(0, 9.1, 1)

# plt.yticks(my_y_ticks)

# plt.xlabel(r'$\mathcal{kh}$', fontsize=15)
# plt.ylabel('Normalized frequency' + r'($\Omega$)', fontsize=10,)
plt.title('PWE(elastic) [H=%s, n=%s]' % (h, nrsquare))
plt.grid()
plt.savefig('./fig/PWE(elastic) [H=%s, n=%s, l1=%s, l2=%s].png' % (h, nrsquare, l1, l2))
plt.show()

