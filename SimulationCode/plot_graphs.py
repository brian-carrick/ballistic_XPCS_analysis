import io
from matplotlib import pyplot as plt

temps = []
with io.open("stress", mode="r") as f:
    next(f)
    for line in f:
        temps.append(line.split())


lam=[i[3] for i in temps]
lam=[float(i) for i in lam]


FE=[i[4] for i in temps]
FE=[float(i) for i in FE] #free energy stored in chain
deltaFE=[i[5] for i in temps]
deltaFE=[float(i) for i in deltaFE]

st0=[i[6] for i in temps]
st0=[float(i) for i in st0]
st1=[i[7] for i in temps]
st1=[float(i) for i in st1]
st2=[i[8] for i in temps]
st2=[float(i) for i in st2]
st3=[i[9] for i in temps]
st3=[float(i) for i in st3]
st4=[i[10] for i in temps]
st4=[float(i) for i in st4]
st5=[i[11] for i in temps]
st5=[float(i) for i in st5]


plt.plot(lam,st0,label='pxx')
plt.plot(lam,st1,label='pyy')
plt.plot(lam,st2,label='pzz')
plt.plot(lam,st3,label='pxy')
plt.plot(lam,st4,label='pyz')
plt.plot(lam,st5,label='pzx')

plt.xlabel('lambda')
plt.ylabel('pressure/stress components')

plt.figure()

plt.plot(lam,deltaFE,label='deltaFE')
plt.xlabel('lambda')
plt.ylabel('delta Free Energy')
plt.legend()
plt.show()
