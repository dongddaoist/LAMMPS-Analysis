from numpy import *
import math
import os.path
import numpy as np
from numpy import linalg as la
f_files=open("files",'r')
str1=f_files.readline()
nfile=int(str1)
# system specifit
g1= 9 #np.array([6]) # starting of vector
g2=11  #np.array([8]) # ending of vector
g3=15 #np.array([12]) # ending of vector
 
v_ref=np.array([0.0,0.0,1.0])
r_ref=math.sqrt(v_ref[0]**2+v_ref[1]**2+v_ref[2]**2)
#f_trj=open("1.gro",'r')
itrj=0
max_mol_type=4
max_N_ptype=20000
kount_mol=np.zeros(max_mol_type)
mol_type=[]
indx_mol=np.zeros((max_N_ptype,2*max_mol_type)) # itart and iend for teach mol
n_mol_type=0
cmol=[]
r_min=0.0
r_bin=0.05
r_max=20.0

a_min=-0.5
a_bin=0.05
a_max=1.0

iselect=1 # which mol type
#nfreqi=input('how many trajectories for one analysis: ')
nfreq=1 #int(nfreqi)
n_bin=math.ceil((r_max-r_min)/r_bin)
n_a_bin=math.ceil((a_max-a_min)/a_bin)
r_array=np.zeros((n_bin,1))
a_array=np.zeros((n_a_bin,1))
for ii in range(n_bin):
   rf=r_min+float(ii)*r_bin
   r_array[ii]=rf+0.5*r_bin
for ii in range(n_a_bin):
   af=a_min+float(ii)*a_bin
   a_array[ii]=af+0.5*a_bin

angle_array=np.zeros(n_bin)
array2=np.zeros(n_bin)
array3=np.zeros(n_bin)
array4=np.zeros(n_bin)
sumi=np.zeros(n_bin) 
sumi2=np.zeros(n_bin)
ftest=open('test2','w')
for zz in range(nfile):
  filei=f_files.readline()
  filei.split()
  li=len(filei)
  f_trj=open(filei[0:li-1],'r')  
  while True:
   str2=f_trj.readline()
   if str2=='':
       break
   itrj=itrj+1
   print(['id of this trj is: ',itrj])
   str2=f_trj.readline()
   nat=int(str2)
   if itrj==1:
      for ii in range(nat):
          str1=f_trj.readline()
          res_num=int(str1[0:5])
          res_name=str1[5:10]
          at_name=str1[10:15]
          at_num=int(str1[15:20])
          if n_mol_type==0:
              mol_type.append(res_name)
              #kount_mol[n_mol_type]=kount_mol[n_mol_type]+1
              istart=int(2*n_mol_type)
              #iend=istart+1
              tmp=int(kount_mol[n_mol_type])
              indx_mol[tmp][istart]=ii
              indx_mol[tmp][istart+1]=ii
              cmol=res_num
              n_mol_type=n_mol_type+1
          else:
              old_type=0  #indicate whether it's new mol type
              for jj in range(n_mol_type):
                  if mol_type[jj]==res_name and res_num==cmol:
                      #kount_mol[jj]=kount_mol[jj]+1
                      old_type=1
                      iend=2*jj+1
                      tmp=int(kount_mol[jj])
                      indx_mol[tmp][iend]=indx_mol[tmp][iend]+1
                  elif mol_type[jj]==res_name and res_num!=cmol:
                      old_type=1
                      istart=int(2*jj)
                      iend=istart+1
                      kount_mol[jj]=kount_mol[jj]+1
                      tmp=int(kount_mol[jj])
                      indx_mol[tmp][istart]=ii
                      indx_mol[tmp][iend]=ii
                      cmol=res_num
              if old_type==0:
                  mol_type.append(res_name)
                  istart=int(2*n_mol_type)
                  tmp=int(kount_mol[n_mol_type])
                  indx_mol[tmp][istart]=ii
                  indx_mol[tmp][istart+1]=ii
                  cmol=res_num
                  n_mol_type=n_mol_type+1
      for ii in range(max_N_ptype):
          for jj in range(n_mol_type):
              ftest.write('%6d'%indx_mol[ii][2*jj])
              ftest.write('%6d'%indx_mol[ii][2*jj+1])
          ftest.write('\n')   
      for ii in range(n_mol_type):
          kount_mol[ii]=kount_mol[ii]+1
      kount_mol.astype(int)    
   else:
       x_mat=np.zeros((nat,3))
       for ii in range(nat):
          str1=f_trj.readline()
          for jj in range(3):
              jstart=20+jj*8
              jend=20+(jj+1)*8
              x_mat[ii][jj]=float(str1[jstart:jend])
       #   x_mat[ii][2]=x_mat[ii][2]+10.0
       #   if x_mat[ii][2]> 20.0:
       #       x_mat[ii][2]=x_mat[ii][2]-20.0
       box0=np.zeros(3)
       hbox=np.zeros(3)
       for kk in range(3):
           box0[kk]=float(str2s[kk])
           hbox[kk]=box0[kk]*0.5
       if np.mod(itrj,nfreq)==0:    
         for ii in range(int(kount_mol[iselect-1])):
             ave_v1=np.zeros(3)
             jstart=indx_mol[ii][2*iselect-2]
             dr=np.zeros(3)
             jj1=int(jstart+g1-1)
             for kk in range(3):
                 ave_v1[kk]=x_mat[jj1][kk]
             print('here')
             print(ave_v1)
             jj2=int(jstart+g2-1)
             ave_v2=np.zeros(3)
             for kk in range(3):
                  dk=x_mat[jj2][kk]-x_mat[jj1][kk]
                  ave_v2[kk]=x_mat[jj2][kk]
                  if dk>hbox[kk]:
                      ave_v2[kk]=ave_v2[kk]-box0[kk]
                  elif  dk<-hbox[kk]:
                      ave_v2[kk]=ave_v2[kk]+box0[kk]     
             print(ave_v2)
             jj3=int(jstart+g3-1)
             ave_v3=np.zeros(3)
             for kk in range(3):
                  dk=x_mat[jj3][kk]-x_mat[jj1][kk]
                  ave_v3[kk]=x_mat[jj3][kk]
                  if dk>hbox[kk]:
                      ave_v3[kk]=ave_v3[kk]-box0[kk]
                  elif  dk<-hbox[kk]:
                      ave_v3[kk]=ave_v3[kk]+box0[kk]
             print(ave_v3)
             v1=np.zeros(3)
             v2=np.zeros(3)
             r1=0
             r2=0
             for kk in range(3):
                 v1[kk]=ave_v2[kk]-ave_v1[kk]
                 v2[kk]=ave_v3[kk]-ave_v1[kk]
             vc=np.zeros(3)
             for kk in range(3):
                 vc[kk]=(ave_v2[kk]+ave_v3[kk])/2.0
             vf1=np.zeros(3)
             vf2=np.zeros(3)
             for kk in range(3):
                 vf1[kk]=v1[kk]
                 vf2[kk]=v2[kk]
             vf1[2]=vc[2]
             vf2[2]=vc[2]

             vvv=np.cross(v1,v2)
             vvv2=np.zeros(3)
             if vvv[2]<0:
                 for kk in range(3):
                     vvv2[kk]=-vvv[kk]
             else:
                 for kk in range(3):
                     vvv2[kk]=vvv[kk]
             rvv2=math.sqrt(vvv2[0]**2+vvv2[1]**2+vvv2[2]**2) 

           
             vvv21=np.cross(vf1,vf2)
             vvv3=np.zeros(3)
             if vvv21[2]<0:
                 for kk in range(3):
                     vvv3[kk]=-vvv21[kk]
             else:
                 for kk in range(3):
                     vvv3[kk]=vvv21[kk]
             rvv3=math.sqrt(vvv3[0]**2+vvv3[1]**2+vvv3[2]**2)
             if rvv2==0 or rvv3==0:
                 continue
             cosj= np.dot(vvv2,vvv3)/(rvv2*rvv3) #/(la.norm(vc1)*la.norm(vc2))
             if cosj<0:
                 cosj=-cosj
             if cosj>1:
                 print('error')
                 print(cosj)
                 cosj=1.0
             if cosj<0:
                 print('error2')
                 print(cosj)
                 cosj=0.0
             z_ave=(ave_v1[2]+ave_v2[2]+ave_v3[2])/3
             loc1=int(round((z_ave-r_min)/r_bin))
             angletrad=np.arccos(cosj)
             if angletrad>np.pi/2.0:
                 angletrad=np.pi-angletrad

             anglet=angletrad/np.pi*180.0
             #anglej=90-anglej
             cost=math.cos(angletrad)
             sint=math.sin(angletrad)  
             cos2t=(2*cost*cost-1)
             sin2t=(2*sint*cost)

             if abs(cosj)>10**-6:
                angle_array[loc1-1]=angle_array[loc1-1]+anglet#/abs(sint) #np.arccos(cosj)/np.pi*180.0
                sumi[loc1-1]=sumi[loc1-1]+1 #sin2t #1/abs(cosj)
                array2[loc1-1]=array2[loc1-1]+cost#*sin2t  # cosj 
                array3[loc1-1]=array3[loc1-1]+cos2t#*sin2t#/cosj #cos(2*j)
                sumi2[loc1-1]=sumi2[loc1-1]+1
                array4[loc1-1]=array4[loc1-1]+0.5*(3*cos2t**2-1)#*sin2t
            # loc2=int(round((aj-a_min)/a_bin))
            # if aj0!=0:
   if np.mod(itrj,10*nfreq)==0: 
      angle_arrayi=np.zeros(n_bin)
      array2i=np.zeros(n_bin)
      array3i=np.zeros(n_bin)
      array4i=np.zeros(n_bin)
      for ii in range(n_bin):
          angle_arrayi[ii]=angle_array[ii]
          array2i[ii]=array2[ii]
          array3i[ii]=array3[ii]
          array4i[ii]=array4[ii]
      f1=open("alignment-i.dat","w")
      for ii in range(n_bin):
          tmp=0
          tmp2=0
          tmp3=0
          tmp4=0
          tmp5=0
          tmp6=0
          if sumi[ii]!=0:
                  tki=sumi[ii]
                  tmp2=angle_arrayi[ii]/tki
                  tmp3=sumi[ii]/(itrj/nfreq)
                  tmp4=array2i[ii]/tki
                  tmp5=array3i[ii]/tki
                  tmp6=array4i[ii]/tki
          f1.write('%12.4f'%r_array[ii])    
          f1.write('%12.4f'%tmp2)
          f1.write('%12.4f'%tmp3)
          f1.write('%12.4f'%tmp4)
          f1.write('%12.4f'%tmp5)
          f1.write('%12.4f'%tmp6)  #p2
          f1.write('\n')
      f1.write('\n')   
      f1.close()
   str2=f_trj.readline()
   str2s=str2.split()

# output
f1=open("density-profile.dat","w")
angle_arrayi=angle_array
f1=open("alignment.dat","w")
for ii in range(n_bin):
    tmp=0
    tmp2=0
    tmp3=0
    meani=sum(sumi[150:250])/len(sumi[150:250])
    if sumi[ii]!=0:
            tmp2=90-angle_arrayi[ii] /sumi[ii]
            tmp3=sumi[ii]/(itrj/nfreq)#/sumi[ii]
        
    f1.write('%10.4f'%r_array[ii])
    f1.write('%12.4f'%tmp2)
    f1.write('%12.4f'%tmp3)
    f1.write('\n')
f1.write('\n')
f1.close()

