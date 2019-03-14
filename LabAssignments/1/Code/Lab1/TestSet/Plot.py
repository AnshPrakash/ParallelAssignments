import sys
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np


data_set_size=[500,5000,50000,500000]
seq=[0.03298,0.258929,2.575969,25.511962]

omp2=[0.020707,0.159518,1.489727,15.88824]
omp4=[0.021177,0.152447,1.426512,14.767432]
omp8=[0.035059,0.178107,1.631854,14.435874]
omp16=[0.054059,0.185688,1.558936,14.593658]

OMP=[omp2,omp4,omp8,omp16]
pth2=[0.058816,0.384435,3.469053,34.836378]
pth4=[0.087848,0.602451,6.085788,63.299655]
pth8=[0.120449,0.594728,6.118216,63.606386]
pth16=[0.190338,0.67378,6.302639,64.929339]

PTHR=[pth2,pth4,pth8,pth16]

omp_fixed=[1.489727,1.426512,1.631854,1.558936]
pthr_fixed=[3.469053,6.085788,6.118216,6.302639]


sppomp=[[0]*4 for _ in range(4)]
spp_pth=[[0]*4 for _ in range(4)]

for i in range(4):
	for j in range(4):
		sppomp[i][j]=seq[j]/(OMP[i][j])
		spp_pth[i][j]=seq[j]/PTHR[i][j]
	# spp_pth[i]=seq[i]/(pth2[i])
# plt.plot(data_set_size,seq,'-p')
# plt.plot(data_set_size, pth2,'-^')
# plt.plot(data_set_size, pth4,'-*')
# plt.plot(data_set_size, pth8,'-o')
# plt.plot(data_set_size, pth16,'-h')

# plt.plot(data_set_size,seq,'-p')
# plt.plot(data_set_size, omp2,'-^')
# plt.plot(data_set_size, omp4,'-*')
# plt.plot(data_set_size, omp8,'-o')
# plt.plot(data_set_size, omp16,'-h')


# plt.plot(data_set_size,seq,'-p')
# plt.plot(data_set_size, omp2,'-^')
# plt.plot(data_set_size, pth2,'-*')

# plt.plot(data_set_size,seq,'-p')
# plt.plot(data_set_size, omp4,'-^')
# plt.plot(data_set_size, pth4,'-*')


# plt.plot(data_set_size,seq,'-p')
# plt.plot(data_set_size, omp8,'-^')
# plt.plot(data_set_size, pth8,'-*')


p=[1,2,4,8,16]
# # for fixed dataset
# omp=[1.489727,1.426512,1.631854,1.558936]
# pthread=[3.469053,6.085788,6.118216,6.302639]



# # plt.plot(data_set_size,seq,'-p')
# plt.plot(p, omp16,'-^')
# plt.plot(p, pth16,'-*')
spp_pth=list(map(list, zip(*spp_pth)))
sppomp=list(map(list, zip(*sppomp)))

plt.plot(p,[1]+spp_pth[0],'-^')
plt.plot(p,[1]+spp_pth[1],'-^')
plt.plot(p,[1]+spp_pth[2],'-^')
plt.plot(p,[1]+spp_pth[3],'-^')

# plt.plot(p,spp_pth,'-*')
plt.title('Speed Up of Pthread')
plt.ylabel('Speed Up')
plt.xlabel('Number of Processing Element')

# plt.legend(['seq', 'pth2', 'pth4', 'pth8','pth16'], loc='upper left')
plt.legend(['Data_500','Data_5000','Data_50000','Data_500000'], loc='upper left')


plt.show()