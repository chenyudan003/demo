# -*- coding:utf-8 -*-
# datetime:2019/6/3 0003 15:09
# software: PyCharm

import cmath
import numpy as np
import matplotlib.pyplot as plt


# import sys
# sys.setcheckinterval(10)


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
    AA1 = np.hstack((-C44, F12 - D44))
    AA2 = np.hstack((F12 - D44, -B4477))
    AA = np.vstack((AA1, AA2))
    BB = np.zeros([2 * NG, 2 * NG], dtype=complex)
    CC11 = np.zeros([NG, NG], dtype=complex)
    CC12 = np.zeros([NG, NG], dtype=complex)
    CC21 = np.zeros([NG, NG], dtype=complex)
    CC22 = np.zeros([NG, NG], dtype=complex)
    for i in range(NG):
        for j in range(NG):
            kGi = kk + G[i, :]
            kGj = kk + G[j, :]
            CC11[i, j] = -C44[i, j] * kGi * kGj + LOU[i, j] * b ** 2
            CC12[i, j] = -D44[i, j] * kGi * kGj + F12[i, j] * kGi * kGj
            CC21[i, j] = -D44[i, j] * kGi * kGj + F12[i, j] * kGj ** 2
            CC22[i, j] = -B4477[i, j] * kGi * kGj - A1[i, j]
    CC1 = np.hstack((CC11, CC12))
    CC2 = np.hstack((CC21, CC22))
    CC = np.vstack((CC1, CC2))
    UVA1 = np.hstack((BB, np.eye(2 * NG, 2 * NG)))
    UVA2 = np.hstack((-CC, -BB))
    UVA = np.vstack((UVA1, UVA2))
    UVB1 = np.hstack((np.eye(2 * NG, 2 * NG), BB))
    UVB2 = np.hstack((BB, AA))
    UVB = np.vstack((UVB1, UVB2))
    UVBF = np.mat(UVB).I
    res = np.dot(UVBF, UVA)
    k3, eek = np.linalg.eig(res)
    k3 = list(k3)
    eek1 = eek[0:NG, :]
    eek2 = eek[NG:2 * NG, :]
    # print(k3, '\n', eek1)
    # print(np.size(eek1))
    HH1 = np.zeros([NG, 4 * NG], dtype=complex)
    HH2 = np.zeros([NG, 4 * NG], dtype=complex)
    HH3 = np.zeros([NG, 4 * NG], dtype=complex)
    HH4 = np.zeros([NG, 4 * NG], dtype=complex)
    for i in range(0, 4 * NG):
        HH1[:, i:i + 1] = np.dot(C44, eek1[:, i]) * (ii * k3[i]) * np.exp(ii * k3[i] * h) + np.dot(D44, eek2[:, i]) * (
                    ii * k3[i]) * np.exp(
            ii * k3[i] * h) - np.dot(F44, eek2[:, i]) * (ii * k3[i]) * np.exp(ii * k3[i] * h)
        HH2[:, i:i + 1] = np.dot(C44, eek1[:, i]) * (ii * k3[i]) * np.exp(-ii * k3[i] * h) + np.dot(D44, eek2[:, i]) * (
                    ii * k3[i]) * np.exp(
            -ii * k3[i] * h) - np.dot(F44, eek2[:, i]) * (ii * k3[i]) * np.exp(-ii * k3[i] * h)
        HH3[:, i:i + 1] = np.dot(D44, eek1[:, i]) * (ii * k3[i]) * np.exp(ii * k3[i] * h) + np.dot(B4477,
                                                                                                   eek2[:, i]) * (
                                  ii * k3[i]) * np.exp(ii * k3[i] * h)
        HH4[:, i:i + 1] = np.dot(D44, eek1[:, i]) * (ii * k3[i]) * np.exp(-ii * k3[i] * h) + np.dot(B4477,
                                                                                                    eek2[:, i]) * (
                                  ii * k3[i]) * np.exp(-ii * k3[i] * h)
    HH = np.vstack((HH1, HH2, HH3, HH4))
    # HH = HH / (HH.max() - HH.min())
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
    return ANS


l1 = 0.012
l2 = 0.02
a = l1 + l2       # 晶格
pf = l1/a         # 填充率
a1 = [8.77e+10, 1.74e+10]
c44 = [3.25e+10, 1.26e+10]
d44 = [3.56, 2.22]
b4477 = [2*0.526e-9, 2*0.344e-9]
lou = [2.16e+3, 3.16e+3]
f44 = [-3.56, -2.56]
f12 = [-0.356, -0.356]

ra1 = 2*np.pi/a
ii = cmath.sqrt(-1)


nrsquare = 15
###################
NG = 2*nrsquare+1
G = np.zeros([NG, 1])
count1 = 0
for l in np.arange(-nrsquare, nrsquare+1):
    G[count1, :] = l*ra1
    count1 += 1
# print(G)
C44 = generate(c44, NG)
D44 = generate(d44, NG)
F12 = generate(f12, NG)
LOU = generate(lou, NG)
B4477 = generate(b4477, NG)
A1 = generate(a1, NG)
F44 = generate(f44, NG)


####################
h = 1e-9
b = np.pi/a
kk = np.arange(b*1e-7,b+0.1,b/40)
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
    a = switch2(0, 5e+4, 12e+7, 1e+2, 10)
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


plt.scatter(kk, result0,s=15)
plt.scatter(kk, result1,s=15)
plt.scatter(kk, result2,s=15)
plt.scatter(kk, result3,s=15)
plt.scatter(kk, result4,s=15)
plt.scatter(kk, result5,s=15)
plt.scatter(kk, result6,s=15)
plt.scatter(kk, result7,s=15)
plt.scatter(kk, result8,s=15)
plt.scatter(kk, result9,s=15)
plt.title('PWE(flex)2 dimension H=%s' % h)
plt.grid()
plt.savefig('PC(flex)2 [H=%s, n=%s].png' % (h, nrsquare))
plt.show()
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
np.savetxt('PC(flex)2 [H=%s, n=%s]' % (h, nrsquare), np.array(result123), delimiter=",")
