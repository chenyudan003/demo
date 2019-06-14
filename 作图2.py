#-*- coding:utf-8 -*-
# datetime:2019/6/13 0013 13:43
# software: PyCharm
import numpy as np
import matplotlib.pyplot as plt
import time


h = 1e-09
z = np.arange(0,7.1,0.5)
x = np.arange(0,7.01,0.1)
a = np.loadtxt('PC(flex)1 [H=1e-09, n=6]',delimiter=',')
b = np.loadtxt('Flexoelectric H=1e-09 f12=0',delimiter=',')
d = np.loadtxt('Flexoelectric H=1e-09 f12=-0.356',delimiter=',')
c = np.loadtxt('Pure elastic 0.01',delimiter=',')
localtime = time.strftime("%Y-%m-%d %H-%M-%S", time.localtime())
ltime = str(localtime)

for i in range(6):
    y = a[i, :]
    x = np.arange(0,1,1/a.shape[1])
    plt.scatter(x, y, label='with flexoelectric',s=15)
# for i in range(a.shape[0]):
#     y = a[i, :]
#     if i == 0:
#         plt.plot(x, y, label='with flexoelectric', color='grey')
#     else:
#         plt.plot(x, y, color='grey')
#
# for i in range(b.shape[0]):
#     y = b[i, :]
#     if i == 0:
#         plt.plot(x, y,'-.', label='f12=f44=0',color='lightgreen')
#     else:
#         pass
#         plt.plot(x, y,'-.',color='lightgreen')
# for i in range(d.shape[0]):
#     y = d[i, :]
#     if i == 0:
#         plt.plot(x, y, '--', label='f12=0.1f44',color='black')
#     else:
#         pass
#         plt.plot(x, y, '--',color='black')

# for i in range(c.shape[0]):
#     y = c[i, :]
#     if i == 0:
#         plt.plot(x, y, '--',label='elastic',color='coral')
#     else:
#         plt.plot(x, y, '--', color='coral')

# for i in range(a.shape[0]):
#     each0 = []
#     y1 = c[i, :]
#     y2 = a[i, :]
#     for j in range(len(y1)):
#         each = (y1[j] - y2[j])
#         each0.append(each)
#     #plt.plot(x, each0)




plt.xlabel(r'$\mathcal{pi /a}$', fontsize=15)
plt.ylabel('Normalized frequency' + r'($\Omega$)', fontsize=10,)
# plt.xlim((0, 7))

#plt.ylim((0, 10))
# my_x_ticks = np.arange(1, 7.1, 1)
# my_y_ticks = np.arange(0, 9.1, 1)
# plt.yticks(my_y_ticks)
# plt.xticks(my_x_ticks)




#plt.title(' H=1e-9')
plt.grid(linestyle="--")
plt.title('PWE(flex)2 dimension H=%s' % h)

# plt.savefig('PC(flex)1 [H=%s, n=%s].png' % (h, nrsquare))
plt.show()