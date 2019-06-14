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




ee = 1e-9
l1 = 0.015
l2 = 0.025
a = l1 + l2       # 晶格
pf = l1/a         # 填充率
c44 = [3.25e+10, 1.26e+10]
lou = [2.16e+2, 3.16e+3]
nrsquare = 6
ra1 = 2*np.pi/a
ii = cmath.sqrt(-1)

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
h = 1e-9
result = []
result1 = []
result2 = []
q = 0.85
kk = ra1*0.5*q
x = np.arange(1e+4, 70e+4,2e+3)
for b in x:
    a = fun1(b,kk)
    # result.append(a.real)
    # result1.append(a.imag)
    result2.append(abs(a))

# plt.plot(x,result)
plt.plot(x,result2)
plt.title('det(flex) [H=%s, n=%s, k=%s]' % (h, nrsquare, q))
plt.grid()
plt.savefig('./fig/det(flex) [H=%s, n=%s, l1=%s, l2=%s, k=%s].png' % (h, nrsquare, l1, l2, q))
plt.show()

