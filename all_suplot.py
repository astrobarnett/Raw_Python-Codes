#for velocity dispersion we used p[2] and replace -1 and do not minus the 1180
#for relative velocity we used p[1] and put -1 and minus the 1180

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from astropy import units as u
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors


##################H_alpha#############

a=6562 #observed wavelength here it is H_alpha
#following is for the H_alpha flux map
f=np.loadtxt('H_alpha_outnew.txt',usecols=(5),unpack=True) # we have extract the wavelength(flux) 
nsf=np.loadtxt('H_beta_outnew.txt',usecols=(4),unpack=True) # we have extract the wavelength 

#following is for the H_alpha relatve velocity map
v1=np.loadtxt('H_alpha_outnew.txt',usecols=(1),unpack=True) # we have extract the wavelength 
v1=((v1/a-1))*299793.4 #v=(l1/l0-1)*c
v1[:]=v1[:]-1180

#following is for the H_alpha velocity dispersion map
d1=np.loadtxt('H_alpha_outnew.txt',usecols=(2),unpack=True) # we have extract the wavelength 
d1=((d1/a))*299793.4 



##################H_beta#############

b=4861.31 #observed wavelength here it is H_beta
bf=np.loadtxt('H_beta_outnew.txt',usecols=(5),unpack=True) # we have extract the wavelength 
ns1=np.loadtxt('H_beta_outnew.txt',usecols=(4),unpack=True) # we have extract the wavelength 

#following is for the H_beta relatve velocity map
k=np.loadtxt('H_beta_out.txt',usecols=(1),unpack=True) # we have extract the wavelength 
k=((k/b-1))*299793.4 #v=(l1/l0-1)*c 
k[:]=k[:]-1180

#following is for the H_beta velocity dispersion map
j=np.loadtxt('H_beta_out.txt',usecols=(2),unpack=True) # we have extract the wavelength 
j=((j/b))*299793.4 

sn=bf/ns1

#condition to ignore the negative flux 
for i in range(len(bf)):
	if sn[i]<0.1:
		bf[i]=np.nan
		k[i]=np.nan
		j[i]=np.nan		


################## N_II #############
c=6583.41 #observed wavelength of N-II
nf=np.loadtxt('N-II_out.txt',usecols=(0),unpack=True) # we have extract the wavelength 
nf=abs(nf)
#following is for the N_II relatve velocity map
n=np.loadtxt('N-II_out.txt',usecols=(1),unpack=True) # we have extract the wavelength 
n=((n/c-1))*299793.4 #v=(l1/l0-1)*c 
n[:]=n[:]-1180

#following is for the N_II velocity dispersion map
q=np.loadtxt('N-II_out.txt',usecols=(2),unpack=True) # we have extract the wavelength 
q=((q/c))*299793.4 

################## O_II #############
d=5006.83 #observed wavelength of O-II
of=np.loadtxt('O-II_out.txt',usecols=(0),unpack=True) # we have extract the wavelength 
of=abs(of)
#following is for the O_II relatve velocity map
o=np.loadtxt('O-II_out.txt',usecols=(1),unpack=True) # we have extract the wavelength 
o=((o/d-1))*299793.4 #v=(l1/l0-1)*c 
o[:]=o[:]-1180

#following is for the O_II velocity dispersion map
r=np.loadtxt('O-II_out.txt',usecols=(2),unpack=True) # we have extract the wavelength 
r=((r/d))*299793.4 


################## [OIII]N1 #############
e=4958.92 #observed wavelength of O-II
gf=np.loadtxt('O-I_out.txt',usecols=(0),unpack=True) # we have extract the wavelength 
gf=abs(gf)
#following is for the O_I relatve velocity map
t=np.loadtxt('O-I_out.txt',usecols=(1),unpack=True) # we have extract the wavelength 
t=((t/e-1))*299793.4 #v=(l1/l0-1)*c 
t[:]=t[:]-1180

#following is for the O_I velocity dispersion map
w=np.loadtxt('O-I_out.txt',usecols=(2),unpack=True) # we have extract the wavelength 
w=((w/e))*299793.4 


################ Make the Grid ######################
#H_alpha
f1=np.zeros([580,580],float) #flux grid
v2=np.zeros([580,580],float) #velocity grid
d2=np.zeros([580,580],float) #dispersion grid
#H_beta
bf1=np.zeros([580,580],float)
k1=np.zeros([580,580],float)
j1=np.zeros([580,580],float)
#Nitrogen
c1=np.zeros([580,580],float)
n1=np.zeros([580,580],float) 
q1=np.zeros([580,580],float)
#oxygen_II
of1=np.zeros([580,580],float) 
o1=np.zeros([580,580],float) 
r1=np.zeros([580,580],float)
#oxygen_I
gf1=np.zeros([580,580],float)
t1=np.zeros([580,580],float)
w1=np.zeros([580,580],float)




x,y,bin_num=np.genfromtxt('out_step1_sq.txt10',skip_header=1,usecols=(0,1,2),unpack=True)


for i in range(len(x)):
	f1[int(x[i]-1),int(y[i]-1)]=f[int(bin_num[i]-1)]
	v2[int(x[i]-1),int(y[i]-1)]=v1[int(bin_num[i]-1)]
	d2[int(x[i]-1),int(y[i]-1)]=d1[int(bin_num[i]-1)]
	bf1[int(x[i]-1),int(y[i]-1)]=bf[int(bin_num[i]-1)]
	k1[int(x[i]-1),int(y[i]-1)]=k[int(bin_num[i]-1)]
	j1[int(x[i]-1),int(y[i]-1)]=j[int(bin_num[i]-1)]	
	c1[int(x[i]-1),int(y[i]-1)]=nf[int(bin_num[i]-1)]
	n1[int(x[i]-1),int(y[i]-1)]=n[int(bin_num[i]-1)]
	q1[int(x[i]-1),int(y[i]-1)]=q[int(bin_num[i]-1)]	
	of1[int(x[i]-1),int(y[i]-1)]=of[int(bin_num[i]-1)]
	o1[int(x[i]-1),int(y[i]-1)]=o[int(bin_num[i]-1)]
	r1[int(x[i]-1),int(y[i]-1)]=r[int(bin_num[i]-1)]	
	gf1[int(x[i]-1),int(y[i]-1)]=gf[int(bin_num[i]-1)]
	t1[int(x[i]-1),int(y[i]-1)]=t[int(bin_num[i]-1)]
	w1[int(x[i]-1),int(y[i]-1)]=w[int(bin_num[i]-1)]




#pdf = PdfPages("all_2D-maps.pdf")

fig,axs=plt.subplots(5,3)

#Halpha
im=axs[0,0].imshow(f1,norm=colors.LogNorm())
axs[0,0].set_ylabel('Size(")')
axs[0,0].set_yticklabels([])
axs[0,0].set_xticklabels([])
axs[0,0].xaxis.set_minor_locator(MultipleLocator(50))
axs[0,0].text(30,60,'(0.2")',color='r')
axs[0,0].set_title(r'$H\alpha \lambda6562$')
plt.colorbar(im,ax=axs[0,0])#,label='Flux (10$^{-16}$erg/s/cm$^{-2}$/A')


im=axs[0,1].imshow(v2,vmin=-300,vmax=300,cmap='RdBu')
axs[0,1].set_yticklabels([])
axs[0,1].set_xticklabels([])
axs[0,1].xaxis.set_minor_locator(MultipleLocator(50))
axs[0,1].text(30,60,'(0.2")',color='r')
plt.colorbar(im,ax=axs[0,1])#,label='Relative velocity ($Kms^{-1}$)')


im=axs[0,2].imshow(d2,vmin=0,vmax=200,cmap='hot_r')
axs[0,2].set_yticklabels([])
axs[0,2].set_xticklabels([])
axs[0,2].xaxis.set_minor_locator(MultipleLocator(50))
axs[0,2].text(30,60,'(0.2")',color='r')
plt.colorbar(im,ax=axs[0,2])#,label='Velocity dispersion ($Kms^{-1}$)')

#Hbeta
im=axs[1,0].imshow(bf1,norm=colors.LogNorm(),vmin=0.1,vmax=10000)
axs[1,0].set_ylabel('Size(")')
axs[1,0].set_yticklabels([])
axs[1,0].set_xticklabels([])
axs[1,0].xaxis.set_minor_locator(MultipleLocator(50))
axs[1,0].text(30,60,'(0.2")',color='r')
axs[1,0].set_title(r'$H\beta \lambda4861$')
plt.colorbar(im,ax=axs[1,0])#,label='Flux (10$^{-16}$erg/s/cm$^{-2}$/A')


im=axs[1,1].imshow(k1,vmin=-300,vmax=300,cmap='RdBu')
axs[1,1].set_yticklabels([])
axs[1,1].set_xticklabels([])
axs[1,1].xaxis.set_minor_locator(MultipleLocator(50))
axs[1,1].text(30,60,'(0.2")',color='r')
plt.colorbar(im,ax=axs[1,1])#,label='Relative velocity ($Kms^{-1}$)')


im=axs[1,2].imshow(j1,vmin=0,vmax=200,cmap='hot_r')
axs[1,2].set_yticklabels([])
axs[1,2].set_xticklabels([])
axs[1,2].xaxis.set_minor_locator(MultipleLocator(50))
axs[1,2].text(30,60,'(0.2")',color='r')
plt.colorbar(im,ax=axs[1,2])#,label='Velocity dispersion ($Kms^{-1}$)')

#Nitrogen
im=axs[2,0].imshow(c1,norm=colors.LogNorm())
axs[2,0].set_ylabel('Size(")')
axs[2,0].set_yticklabels([])
axs[2,0].set_xticklabels([])
axs[2,0].xaxis.set_minor_locator(MultipleLocator(50))
axs[2,0].text(30,60,'(0.2")',color='r')
axs[2,0].set_title(r'[NII] $\lambda6583$')
plt.colorbar(im,ax=axs[2,0])#,label='Flux (10$^{-16}$erg/s/cm$^{-2}$/A')


im=axs[2,1].imshow(n1,vmin=-300,vmax=300,cmap='RdBu')
axs[2,1].set_yticklabels([])
axs[2,1].set_xticklabels([])
axs[2,1].xaxis.set_minor_locator(MultipleLocator(50))
axs[2,1].text(30,60,'(0.2")',color='r')
plt.colorbar(im,ax=axs[2,1])#,label='Relative velocity ($Kms^{-1}$)')


im=axs[2,2].imshow(q1,vmin=0,vmax=200,cmap='hot_r')
axs[2,2].set_yticklabels([])
axs[2,2].set_xticklabels([])
axs[2,2].xaxis.set_minor_locator(MultipleLocator(50))
axs[2,2].text(30,60,'(0.2")',color='r')
plt.colorbar(im,ax=axs[2,2])#,label='Velocity dispersion ($Kms^{-1}$)')

#oxygen3-2 [O III]N2_5006.84 
im=axs[3,0].imshow(of1,norm=colors.LogNorm())
axs[3,0].set_ylabel('Size(")')
axs[3,0].set_yticklabels([])
axs[3,0].set_xticklabels([])
axs[3,0].xaxis.set_minor_locator(MultipleLocator(50))
axs[3,0].text(30,60,'(0.2")',color='r')
axs[3,0].set_title(r'[OIII] $\lambda5007$')
plt.colorbar(im,ax=axs[3,0])#,label='Flux (10$^{-16}$erg/s/cm$^{-2}$/A')


im=axs[3,1].imshow(o1,vmin=-300,vmax=300,cmap='RdBu')
axs[3,1].set_yticklabels([])
axs[3,1].set_xticklabels([])
axs[3,1].xaxis.set_minor_locator(MultipleLocator(50))
axs[3,1].text(30,60,'(0.2")',color='r')
plt.colorbar(im,ax=axs[3,1])#,label='Relative velocity ($Kms^{-1}$)')


im=axs[3,2].imshow(r1,vmin=0,vmax=200,cmap='hot_r')
axs[3,2].set_yticklabels([])
axs[3,2].set_xticklabels([])
axs[3,2].xaxis.set_minor_locator(MultipleLocator(50))
axs[3,2].text(30,60,'(0.2")',color='r')
plt.colorbar(im,ax=axs[3,2])#,label='Velocity dispersion ($Kms^{-1}$)')

#oxygen1 [OIII]N1 4958.92 
im=axs[4,0].imshow(gf1,norm=colors.LogNorm())
axs[4,0].set_ylabel('Size(")')
axs[4,0].set_yticklabels([])
axs[4,0].set_xticklabels([])
axs[4,0].xaxis.set_minor_locator(MultipleLocator(50))
axs[4,0].text(30,60,'(0.2")',color='r')
axs[4,0].set_title(r'[OIII] $\lambda4959$')
plt.colorbar(im,ax=axs[4,0])#,label='Flux (10$^{-16}$erg/s/cm$^{-2}$/A')


im=axs[4,1].imshow(t1,vmin=-300,vmax=300,cmap='RdBu')
axs[4,1].set_yticklabels([])
axs[4,1].set_xticklabels([])
axs[4,1].xaxis.set_minor_locator(MultipleLocator(50))
axs[4,1].text(30,60,'(0.2")',color='r')
plt.colorbar(im,ax=axs[4,1])#,label='Relative velocity ($Kms^{-1}$)')


im=axs[4,2].imshow(w1,vmin=0,vmax=200,cmap='hot_r')
axs[4,2].set_yticklabels([])
axs[4,2].set_xticklabels([])
axs[4,2].xaxis.set_minor_locator(MultipleLocator(50))
axs[4,2].text(30,60,'(0.2")',color='r')
plt.colorbar(im,ax=axs[4,2])#,label='Velocity dispersion ($Kms^{-1}$)')


#pdf.savefig(fig)
#plt.close()
#pdf.close()
fig.tight_layout()
plt.subplots_adjust(wspace=0.1, hspace=0.5)
#plt.tight_layout()
plt.show()
















