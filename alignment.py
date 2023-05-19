from numpy import *
import math
import os.path
import numpy as np


f_files=open("files",'r')
str1=f_files.readline()
nfile=int(str1)
# system specifit
g1=np.array([9]) # starting of vector
n_g1=len(g1)
g2=np.array([11,15]) # ending of vector
n_g2=len(g2)

v_ref=np.array([0.0,0.0,1.0])
r_ref=sqrt(v_ref[0]**2+v_ref[1]**2+v_ref[2]**2)
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

density_array=np.zeros(n_bin) # total rdfs
density_array1=np.zeros(n_bin) # total rdfs
density_array2=np.zeros(n_bin)
density_array3=np.zeros(n_bin)
angle_array=np.zeros(n_bin)
sumi=np.zeros(n_bin) 
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
       if np.mod(itrj,nfreq)==0 and itrj>6000:    
         for ii in range(int(kount_mol[iselect-1])):
             ave_v1=np.zeros(3)
             jstart=indx_mol[ii][2*iselect-2]
             dr=np.zeros(3)
             avez=0.0
             for j1 in range(n_g1):
                 if j1==0:
                     jj1=int(jstart+g1[j1]-1)
                     for kk in range(3):
                        ave_v1[kk]=ave_v1[kk]+x_mat[jj1][kk]
                     avez=avez+x_mat[jj1][2]   
                 else:
                     jj2=int(jstart+g1[j1]-1)
                     for kk in range(3):
                         dk=x_mat[jj2][kk]-x_mat[jj1][kk]
                         ave_v1[kk]=ave_v1[kk]+x_mat[jj2][kk]
                         if dk>hbox[kk]:
                            ave_v1[kk]=ave_v1[kk]-box0[kk]
                         elif dk<-hbox[kk]:
                            ave_v1[kk]=ave_v1[kk]+box0[kk]
                     avez=avez+x_mat[jj1][2]        
             for kk in range(3):
                 ave_v1[kk]=ave_v1[kk]/float(n_g1)
             ave_v2=np.zeros(3)
             for j2 in range(n_g2):
                 jj2=int(jstart+g2[j2]-1)
                 for kk in range(3):
                     dk=x_mat[jj2][kk]-x_mat[jj1][kk]
                     ave_v2[kk]=ave_v2[kk]+x_mat[jj2][kk]
                     if dk>hbox[kk]:
                         ave_v2[kk]=ave_v2[kk]-box0[kk]
                     elif  dk<-hbox[kk]:
                         ave_v2[kk]=ave_v2[kk]+box0[kk]
             r=0.0               
             dr=np.zeros(3)
             for kk in range(3):
                 ave_v2[kk]=ave_v2[kk]/float(n_g2)                   
                 dr[kk]= ave_v2[kk]- ave_v1[kk]
                 r=r+dr[kk]**2
             r=sqrt(r)
             cosj= dr[0]*v_ref[0]+dr[1]*v_ref[1]+dr[2]*v_ref[2]
             cosj= cosj/(r*r_ref)
             p2=0.5*(3*cosj**2-1)
#             aj0=np.arccos(cosj)
#             aj=aj0/np.pi*180.0
             zj=0    
             jend=int(indx_mol[ii][2*iselect-1])
             nk=jend-jstart+1
             for kk in range(int(nk)):
                 tmp=int(jstart+kk)
                 zj=zj+x_mat[tmp][2]
             zj=zj/float(nk)
             loc1=int(round((avez-r_min)/r_bin))
             density_array[loc1-1]=density_array[loc1-1]+p2
             density_array1[loc1-1]=density_array1[loc1-1]+cosj
             density_array2[loc1-1]=density_array2[loc1-1]+cosj*cosj
             density_array3[loc1-1]=density_array3[loc1-1]+np.power(cosj,3)
             angle_array[loc1-1]=angle_array[loc1-1]+np.arccos(cosj)/np.pi*180.0
             sumi[loc1-1]=sumi[loc1-1]+1
            # loc2=int(round((aj-a_min)/a_bin))
            # if aj0!=0:
            #   density_array[loc1-1][loc2-1]=density_array[loc1-1][loc2-1]+ 1 #/abs(sin(aj0))
   if np.mod(itrj,10*nfreq)==0:
      density_arrayi=np.zeros(n_bin) 
      density_array1i=np.zeros(n_bin)
      density_array2i=np.zeros(n_bin)
      density_array3i=np.zeros(n_bin)

      density_arrayi=density_array  
      density_array1i=density_array1
      density_array2i=density_array2
      density_array3i=density_array3

      angle_arrayi=angle_array
      f1=open("alignment-i.dat","w")
      for ii in range(n_bin):
#          meani=sum(sumi[150:250])/len(sumi[150:250])

          tmp=0
          tmp1=0
          tmp2=0
          tmp3=0
          tmp4=0
          tki=sumi[ii]
          if tki !=0: 
              tmp=density_arrayi[ii]/tki
              tmp1=density_array1i[ii]/tki
              tmp2=density_array2i[ii]/tki
              tmp3=density_array3i[ii]/tki
          #tmp4=density_array2i[ii]/tki
                    
          f1.write('%18.6f'%r_array[ii])   
          f1.write('%18.6f'%tmp)
          f1.write('%18.6f'%tmp1)
          f1.write('%18.6f'%tmp2)
          f1.write('%18.6f'%tmp3)
          tmp4=tki/(itrj/nfreq)
          f1.write('%18.6f'%tmp4)
          f1.write('\n')
      f1.write('\n')   
      f1.close()
   str2=f_trj.readline()
   str2s=str2.split()

# output
f1=open("density-profile.dat","w")
density_arrayi=density_array
angle_arrayi=angle_array
density_array2i=density_array2
f1=open("alignment.dat","w")
for ii in range(n_bin):
    tmp=0
    tmp2=0
    tmp3=0
    tmp4=0
    meani=sum(sumi[150:250])/len(sumi[150:250])
    
    if sumi[ii]!=0:
        tki=sumi[ii]#/meani
        tmp=density_arrayi[ii] /sumi[ii]
        tmp2=angle_arrayi[ii] /sumi[ii]
        tmp3=sumi[ii]#/sumi[ii]
        tmp4=density_array2i[ii]/tki
    f1.write('%10.4f'%r_array[ii])
    f1.write('%12.6f'%tmp)
    f1.write('%12.6f'%tmp2)
    f1.write('%12.6f'%tmp3)
    f1.write('%12.6f'%tmp4)
    f1.write('\n')
f1.write('\n')
f1.close()

