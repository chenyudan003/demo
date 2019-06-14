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




ee = 1e-9
l1 = 0.005
l2 = 0.008
a = l1 + l2       # 晶格
pf = l1/a         # 填充率
a1 = [8.77e+10, 1.74e+10]
c44 = [3.25e+10, 1.05e+10]
d44 = [3.56, 8.26]
b4477 = [2*0.526e-9, 2*0.6e-9]
lou = [2160, 1980]
f44 = [-5.88, -10.56]
f12 = [-3.56, -8.26]
nrsquare = 5
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
D44 = generate(d44, NG)
F12 = generate(f12, NG)
LOU = generate(lou, NG)
B4477 = generate(b4477, NG)
A1 = generate(a1, NG)
F44 = generate(f44, NG)



####################
h = 1e-9
result = []
result1 = []
result2 = []
q = 0.6
kk = ra1*0.5*q
x = np.arange(5e+4, 220e+4,5e+3)
for b in x:
    a = fun1(b,kk)
    # result.append(a.real)
    # result1.append(a.imag)
    result2.append(abs(a))

# plt.plot(x,result)
# plt.plot(x,result1)
plt.plot(x,result2)
plt.title('det(flex) [H=%s, n=%s, k=%s]' % (h, nrsquare, q))
plt.grid()
plt.savefig('./fig/det(flex) [H=%s, n=%s, l1=%s, l2=%s, k=%s].png' % (h, nrsquare, l1, l2, q))
plt.show()

