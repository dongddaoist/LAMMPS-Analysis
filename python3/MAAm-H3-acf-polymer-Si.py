from numpy import *
import math
import os.path
import numpy as np
f_inp=open("acf.inp",'r')
r_cut=float(f_inp.readline()) # cut-off distance for a pair to be defined
max_pair=5  #  number of atom type per
max_paird=max_pair*2
deltt=1.0  # time step in fs
n_dstep=2000  # dump frequency
dtime=deltt*n_dstep
max_frame=40000
max_cn=10  # 10 is the maximum coordination per central atom
nat_polymer=1016
n_paird=2
acf_index=np.zeros((2,max_paird)) # atom indices for RDF calculation
type_at_mat=np.zeros(2)
str1=f_inp.readline()
str1s=str1.split()
#print(str1s)
ni1=int(str1s[0])
ni2=int(str1s[1])
#print([ni1,ni2])
type_at_mat[0]=ni1
type_at_mat[1]=ni2
for jj in range(ni1):
   str1=f_inp.readline()
   str1s=str1.split()
   j1=jj*2+0
   j2=j1+1
   acf_index[0][j1]=float(str1s[j1])  
   acf_index[0][j2]=float(str1s[j2])
for jj in range(ni2):
   str1=f_inp.readline()
   str1s=str1.split()
   j1=jj*2+0
   j2=j1+1
   acf_index[1][j1]=float(str1s[j1])     
   acf_index[1][j2]=float(str1s[j2])
f_dat=open("lammps.dat",'r')
f_dat.readline()
str1=f_dat.readline()
str1s=str1.split()
nat=int(str1s[0])
ii=0
id_mat=np.zeros(())
while ii==0:
   str1=f_dat.readline()
   str1s=str1.split()
   if len(str1s)>0:
     if str1s[0]=="Atoms":
        ii=1
f_dat.readline()
iat_select=np.zeros((nat,n_paird))
kount_iat_select=np.zeros((n_paird))
fw=open("test",'w')
for ii in range(nat):
   str1=f_dat.readline()
   str1s=str1.split()
   typei=int(str1s[2])
   chgi=float(str1s[3])
   
   j1=0
   j2=1
   #for kk in range(ni1): #range(int(type_at_mat[0])):
   if ii<nat_polymer:
      tmp=int(kount_iat_select[j1])
      iat_select[tmp][j1]=ii+1 #jj+1
      kount_iat_select[j1]=kount_iat_select[j1]+1
   for kk in range(ni2): #range(int(type_at_mat[1])):
      k1=kk*2
      k2=k1+1
      diff_chg=abs(acf_index[j2][k2]-chgi)
      if int(acf_index[j2][k1])==typei and  diff_chg< 0.0001:
         tmp=int(kount_iat_select[j2])
         iat_select[tmp][j2]=ii+1
         kount_iat_select[j2]=kount_iat_select[j2]+1 

for ii in range(nat):
   for jj in range(n_paird): 
      fw.write('%8d'%iat_select[ii][jj])
   fw.write('\n')
kount_iat_select.astype(int)
#print(kount_iat_select)
f_files=open("files",'r')
str1=f_files.readline()
nfile=int(str1)
ntrj=0
ipos=0.0
nfreq=1
n_eval=3000 #int(min(kount_iat_select)*max_cn)  

pair_indx=np.zeros([n_eval,2,max_frame])
kount_pair=np.zeros(max_frame)
acf=np.zeros(max_frame-1)
acf_knt=np.zeros(max_frame-1)
for zz in range(nfile):
   filei=f_files.readline()
   filei.split()
   #filei.astype.string()
   #print(filei)
   li=len(filei)
   fff=open(filei[0:li-1],'r')
   while True:
     str2=fff.readline()
     if str2=='':
       break
     ipos=ipos+1.0   
     x_mat=np.zeros((nat,3))
     stepi=fff.readline()
     fff.readline()
     fff.readline()
     fff.readline()
     box0=np.zeros((3,2))
     box_size=np.zeros(3)
     hbox_size=np.zeros(3)
     vi=1.0
     for ii in range(3):
         str2=fff.readline()
         str2s=str2.split()
         for kk in range(2):
            box0[ii][kk]=float(str2s[kk])
         box_size[ii]=box0[ii][1]-box0[ii][0]
         vi=vi*box_size[ii]
         hbox_size[ii]=box_size[ii]*0.5
     fff.readline()
     for ii in range(nat):
         str2=fff.readline()
         str2s=str2.split()
         iat=int(str2s[0])
         imol=int(str2s[1])
         for kk in range(3):
            x_mat[iat-1][kk]=float(str2s[kk+2])-box0[kk][0] 
     if int(np.mod(ipos,nfreq))==0:
       #print(['ntrj',ntrj])
#       if int(ntrj)>100:
#           break
       for ii in range(int(kount_iat_select[0])):
           for jj in range(int(kount_iat_select[1])):
               dr=np.zeros(3)
               indi=int(iat_select[ii][0])
               indj=int(iat_select[jj][1])
              # print([indi,indj])
               for kk in range(3):
                   dr[kk]=abs(x_mat[indi][kk]-x_mat[indj][kk])
                   if dr[kk]>hbox_size[kk]:
                       dr[kk]=box_size[kk]-dr[kk]
               if dr[0]< r_cut and dr[1]< r_cut and dr[2]< r_cut: 
                   disij=sqrt(dr[0]**2+dr[1]**2+dr[2]**2)
                   #print(disij)
                   if disij<r_cut:
                      tmp=int(kount_pair[int(ntrj)])
                      if ii<jj:
                         pair_indx[tmp][0][ntrj]=ii
                         pair_indx[tmp][1][ntrj]=jj
                         kount_pair[int(ntrj)]=tmp+1
                      elif ii>jj:
                         pair_indx[tmp][0][ntrj]=jj
                         pair_indx[tmp][1][ntrj]=ii
                         kount_pair[int(ntrj)]=tmp+1
       ntrj=ntrj+1               
   fff.close()
f3=open("acf-tot.dat","w")
imax=int(ntrj-1)
acf_knt=np.zeros(imax)
for ii in range(imax):
    nj=int(ntrj-(ii+1))
#    print(ii)
    for jj in range(nj):
       j1=jj
       j2=jj+(ii+1)
       n0=kount_pair[int(j1)]
       if n0>=1:
          ni=0
          n1=int(kount_pair[j1])
          n2=int(kount_pair[j2])
          if n1>0 and n2>0:
            for k1 in range(n1):
              for k2 in range(n2):
                 if pair_indx[k1][0][j1]==pair_indx[k2][0][j2] and pair_indx[k1][1][j1]==pair_indx[k2][1][j2]:
                    ni=ni+1
                    break
          ni=ni/kount_pair[int(j1)]
          acf[ii]=acf[ii]+ni         
          acf_knt[ii]= acf_knt[ii]+1.0
    if acf_knt[ii] !=0:
       acf[ii]=acf[ii]/acf_knt[ii]
       timex=dtime*nfreq*(ii+1)/1000.0 # unit in ps
       f3.write('%10.4f'%timex)
       f3.write('%10.4f'%acf[ii])
       f3.write('\n')
f3.close()   
#return()
