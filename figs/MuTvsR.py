import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
#from matplotlib.ticker import ScalarFormatter
#sformatter=ScalarFormatter(useOffset=True,useMathText=True)
#sformatter.set_scientific(True)
#sformatter.set_powerlimits((-2,3))


tau=15
tau=int(input('Enter tau: '))
btype=0
btype=int(input('Enter btype: '))

filename='data/muTvsR_pi_tau'+str(tau)+'.txt'
print('filename=',filename)
pidata = np.loadtxt(filename,skiprows=0,unpack=True)
filename='data/muTvsR_K_tau'+str(tau)+'.txt'
print('filename=',filename)
Kdata = np.loadtxt(filename,skiprows=0,unpack=True)
filename='data/muTvsR_B'+str(btype)+'_tau'+str(tau)+'.txt'
print('filename=',filename)
Bdata = np.loadtxt(filename,skiprows=0,unpack=True)

plt.figure(figsize=(5,11))
fig=plt.figure(1)

rpi=pidata[0]
Npi=pidata[1]
Tpi=pidata[2]
Upi=pidata[3]
mupi=pidata[4]
rhopi=tau*pidata[5]
Upi_alt=pidata[7]

rK=Kdata[0]
NK=Kdata[1]
TK=Kdata[2]
UK=Kdata[3]
muK=Kdata[4]
rhoK=tau*Kdata[5]
UK_alt=Kdata[7]

rB=Bdata[0]
NB=Bdata[1]
TB=Bdata[2]
UB=Bdata[3]
muB=Bdata[4]
rhoB=tau*Bdata[5]
UB_alt=Bdata[7]

#######################################

ax = fig.add_axes([0.19,0.06,0.8,0.23])

plt.plot(rpi,1000*Tpi,linestyle='-',color='r',markersize=6,marker='o',markerfacecolor='r')
plt.plot(rK,1000*TK,linestyle='-',color='g',markersize=6,marker='o',markerfacecolor='g')
plt.plot(rB,1000*TB,linestyle='-',color='b',markersize=6,marker='o',markerfacecolor='b')

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,31,10), minor=False)
ax.set_xticklabels(np.arange(0,31,10), minor=False, family='serif')
ax.set_xticks(np.arange(0,31,5), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0f'))
plt.xlim(0,25)

ax.set_yticks(np.arange(0,200,40), minor=False)
ax.set_yticklabels(np.arange(0,200,40), minor=False, family='serif')
ax.set_yticks(np.arange(0,200,20), minor=True)
plt.ylim(20,170)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
#ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$r$ [fm]', fontsize=18, weight='normal')
plt.ylabel('$T$ [MeV]',fontsize=18)

#######################################

ax = fig.add_axes([0.19,0.29,0.8,0.23])

plt.plot(rpi,rhopi,linestyle='-',color='r',markersize=6,marker='o',markerfacecolor='r')
plt.plot(rK,rhoK,linestyle='-',color='g',markersize=6,marker='o',markerfacecolor='g')
plt.plot(rB,rhoB,linestyle='-',color='b',markersize=6,marker='o',markerfacecolor='b')

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,31,10), minor=False)
ax.set_xticklabels([])
ax.set_xticks(np.arange(0,31,5), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0f'))
plt.xlim(0,25)

ax.set_yticks(np.arange(0,4,0.5), minor=False)
ax.set_yticklabels(np.arange(0,4,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(0,4,0.1), minor=True)
plt.ylim(0,2.9999)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
#ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$r$ [fm]', fontsize=18, weight='normal')
plt.ylabel('$\\rho$ [fm$^{-2}$]',fontsize=18)

#######################################

ax = fig.add_axes([0.19,0.52,0.8,0.23])

x=np.array([],dtype=float)
y=np.array([],dtype=float)
z=np.array([],dtype=float)
s=np.prod(Npi.shape)
for i in range(0,s):
	if Npi[i]>4 :
		x=np.append(x,rpi[i])
		y=np.append(y,Upi[i])
		z=np.append(z,Upi_alt[i])
plt.plot(x,y,linestyle='-',color='r',markersize=6,marker='o',markerfacecolor='r')
#plt.plot(x,z,linestyle='-',color='r',markersize=6,marker='s',markerfacecolor='r')

x=np.array([],dtype=float)
y=np.array([],dtype=float)
z=np.array([],dtype=float)
s=np.prod(Npi.shape)
for i in range(0,s):
	if NK[i]>4 :
		x=np.append(x,rK[i])
		y=np.append(y,UK[i])
		z=np.append(z,UK_alt[i])
plt.plot(x,y,linestyle='-',color='g',markersize=6,marker='o',markerfacecolor='g')
#plt.plot(x,z,linestyle='-',color='g',markersize=6,marker='s',markerfacecolor='g')

x=np.array([],dtype=float)
y=np.array([],dtype=float)
z=np.array([],dtype=float)
s=np.prod(NB.shape)
for i in range(0,s):
	if NB[i]>4 :
		x=np.append(x,rB[i])
		y=np.append(y,UB[i])
		z=np.append(z,UB_alt[i])
plt.plot(x,y,linestyle='-',color='b',markersize=6,marker='o',markerfacecolor='b')
#plt.plot(x,z,linestyle='-',color='b',markersize=6,marker='s',markerfacecolor='b')

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,31,10), minor=False)
#ax.set_xticklabels(np.arange(0,31,10), minor=False, family='serif')
ax.set_xticklabels([])
ax.set_xticks(np.arange(0,31,5), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0f'))
plt.xlim(0,25)

ax.set_yticks(np.arange(0,6,1.0), minor=False)
ax.set_yticklabels(np.arange(0,6,1.0), minor=False, family='serif')
ax.set_yticks(np.arange(0,6,0.2), minor=True)
plt.ylim(0,2.999)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
#ax.yaxis.set_major_formatter(sformatter)

plt.xlabel(None)
plt.ylabel('$u_r$ ',fontsize=18)

#######################################

ax = fig.add_axes([0.19,0.75,0.80,0.23])

x=np.array([],dtype=float)
ypi=np.array([],dtype=float)
s=np.prod(Npi.shape)
for i in range(0,s):
	if Npi[i]>4 :
		x=np.append(x,rpi[i])
		ypi=np.append(ypi,mupi[i])
plt.plot(x,ypi,linestyle='-',color='r',markersize=6,marker='o',markerfacecolor='r')
#plt.plot(rpi,mupi*Tpi,linestyle='-',color='r',markersize=6,marker='o',markerfacecolor='r')

x=np.array([],dtype=float)
yK=np.array([],dtype=float)
s=np.prod(NK.shape)
for i in range(0,s):
	if NK[i]>4:
		x=np.append(x,rK[i])
		yK=np.append(yK,muK[i])
plt.plot(x,yK,linestyle='-',color='g',markersize=6,marker='o',markerfacecolor='g')

x=np.array([],dtype=float)
yB=np.array([],dtype=float)
s=np.prod(NB.shape)
for i in range(0,s):
	if NB[i]>4 :
		x=np.append(x,rB[i])
		yB=np.append(yB,muB[i])
plt.plot(x,yB,linestyle='-',color='b',markersize=6,marker='o',markerfacecolor='b')

x=np.array([],dtype=float)
DelY=np.array([],dtype=float)
for i in range(0,s):
	if NB[i]>4 :
		x=np.append(x,rB[i])
		DelY=np.append(DelY,2.0*yB[i]-5.0*ypi[i])
plt.plot(x,DelY,linestyle='-',color='k',markersize=6,marker='o',markerfacecolor='k')

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,31,10), minor=False)
#ax.set_xticklabels(np.arange(0,31,10), minor=False, family='serif')
ax.set_xticklabels([])
ax.set_xticks(np.arange(0,31,5), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0f'))
plt.xlim(0,25)

ax.set_yticks(np.arange(-5.0,20.0,5.0), minor=False)
ax.set_yticklabels(np.arange(-5.0,20.0,5.0), minor=False, family='serif')
ax.set_yticks(np.arange(-5.0,20.0,1.0), minor=True)
#plt.ylim(-25,1250)
#plt.ylim(-2,10.0)
plt.ylim(-1.5,10)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
#ax.yaxis.set_major_formatter(sformatter)

plt.xlabel(None)
plt.ylabel('$\mu/T$',fontsize=18)

#######################################


plt.savefig('MuTvsR.pdf')
os.system('open -a Preview MuTvsR.pdf')
#plt.show()
quit()
