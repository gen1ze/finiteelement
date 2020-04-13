import numpy as np
import matplotlib.pyplot as plt
no_elem=int(input('give no of elements'))
decision=int(input('is elements of equal length,if true give 1 otherwise 0'))
if decision==1:
    length=float(input('give whole beam length'))
    temp=float(length/no_elem)
    length_of_elems=[]
    i=0
    while i<no_elem:
        length_of_elems.append(temp)
        i+=1
    print('length of elements starting from element 1 is:')
    print(length_of_elems)
    print(len(length_of_elems))
else:
    length_of_elems=list(map(float,input('give lengths of each element seperated by spaces').split()))
    length=sum(length_of_elems)
    print('length of elements starting from element 1 is:')
    print(length_of_elems)
    print(len(length_of_elems))
#getting element level stiffness
l=length_of_elems[:]
k=np.ndarray((no_elem,4,4))
k1=np.ndarray((no_elem,4,4))
ei=float(input('give value of E*I'))

i=0
while i<no_elem:
    k1[i,:,:]=np.array([[12,6*l[i],-12,6*l[i]],
                       [6*l[i],4*(l[i]**2),-6*l[i],2*(l[i]**2)],
                       [-12,-6*l[i],12,-6*l[i]],
                       [6*l[i],2*(l[i]**2),-6*l[i],4*(l[i]**2)]])
    k[i,:,:]=(ei*k1[i,:,:])/((l[i]**2)*l[i])
    i+=1
print('element level stiffness are\n',k)
print(k.shape)

#getting loads(considering only uniformly distributed load)
qo=float(input('give load/length'))
print('length of beam must be 5x or more than the 2 dimention for validity of Euler-Bernoulli theory')
z=float(input('give distance of outermost fiber fron neutral axis'))
q=np.ndarray((no_elem,4,1))
i=0
while i<no_elem:
    q[i,:,:]=np.array([qo*l[i]/2,qo*(l[i]**2)/12,qo*l[i]/2,-qo*(l[i]**2)/12]).reshape(4,1)
    i+=1
print('force martices are',q)
print(q.shape)

#getting global stiffness matrix
no_nodes=no_elem+1
K=np.zeros((no_nodes*2,no_nodes*2))
j=0
while j<=((no_nodes*2)-4):
    K_temp=np.zeros((no_nodes*2,no_nodes*2))
    i=int((j/2))
    K_temp[j:j+4,j:j+4]=k[i,:,:]
    K+=K_temp
    j+=2

print('global stiffness matrix is\n',K)
print(K.shape)

#getting global force matrix
Q=np.zeros((no_nodes*2,1))
i=0
while i<=2*(no_elem-1):
    Q_temp=np.zeros((no_nodes*2,1))
    j=int(i/2)
    Q_temp[i:i+4,0:]=q[j,:,:]
    Q+=Q_temp
    i+=2

print('global force matrix is\n',Q)
print(Q.shape)

#solving KU=Q
#by default deflection and slope at node 1 is assumed zero(True for a cantilever beam)
U_temp=np.linalg.solve(K[2:,2:],Q[2:,:])
U1=np.insert(U_temp,0,0,axis=0)
U=np.insert(U1,0,0,axis=0)
print('deflections are\n',U)
print(U.shape)
print('maximum deflection(through fea) is\n',max(U))
print('minimum deflection(through fea) is\n',min(U))

#plotting graphs
a=[]
i=0
while i<(no_nodes*2):
    a.append(U[i,0])
    i+=2
#print(a)
node_num=np.linspace(start=1,stop=no_nodes,num=no_nodes,endpoint=True)
b=[]
i=1
while i<(no_nodes*2):
    b.append(U[i,0])
    i+=2
#print(b)
plt.subplot(2,1,1)
plt.scatter(node_num,b,label='node_slope')
plt.scatter(node_num,a,label='node_displacement')
plt.legend()
plt.xlabel='node location'
plt.grid()


A=np.ndarray((no_elem,4,1))
M=np.array([[2,-2,1,1],
            [-3,3,-2,-1],
            [0,0,1,0],
            [1,0,0,0]])
B=np.ndarray((no_elem,4,1))
i=0
while i<no_elem:
    B[i,:,:]=np.array([a[i],a[i+1],b[i],b[i+1]]).reshape(4,1)
    A[i,:,:]=np.dot(M,B[i,:,:])
    i+=1
#print(A)
#start
A1=A[:,:2,:]
#end
u=np.linspace(start=0,stop=1,endpoint=True,num=50)
r=np.zeros((50,no_elem))
i=0
while i<no_elem:
    j=0
    while j<50:
        r[j,i]=u[j]*(u[j]**2)*A[i,0,:]+(u[j]**2)*A[i,1,:]+u[j]*A[i,2,:]+A[i,3,:]
        j+=1
    i+=1
#print(r)
#print(r.shape)
#start
r2=np.zeros((50,no_elem))
i=0
while i<no_elem:
    j=0
    while j<50:
        r2[j,i]=6*u[j]*A1[i,0,:]+2*A1[i,1,:]
        j+=1
    i+=1
#end
plt.subplot(2,1,2)
plt.title('element wise deformed beam(shape functions)')
plt.ylabel('displacement')
plt.plot(u,r)


R=np.zeros((50*no_elem,1))
i=0
while i<50*no_elem:
    R1 = np.zeros((50 * no_elem, 1))
    k=int(i/50)
    R1[i:i+50,0]=r[:,k]
    R+=R1
    i+=50
#print(R)
#print(R.shape)
#start
R2=np.zeros((50*no_elem,1))
i=0
while i<50*no_elem:
    R3 = np.zeros((50 * no_elem, 1))
    k=int(i/50)
    R3[i:i+50,0]=r2[:,k]
    R2+=R3
    i+=50
#stop

plt.figure(2)
plt.ylabel('displacement')
x=np.linspace(start=0,stop=no_elem,endpoint=True,num=50*no_elem)
x1=np.linspace(start=0,stop=length,endpoint=True,num=50*no_elem)
exact_solution=(((qo*(length**2)*(x1**2))/4)-(((qo*length)*((x1**2)*x1))/6)+((qo*(x1**2)*(x1**2))/24))/ei
plt.title('deformed beam')
plt.grid()
#plt.annotate(('max deflection',max(R)),xy=(no_elem,max(R)),xytext=(1,1),arrowprops=dict())
plt.plot((x1*no_elem)/length,exact_solution,label='exact solution')
plt.plot(x,R,label='fea solution')
plt.legend()
#start
#print(z*R2)
#print(R2.shape)
plt.figure(3)
plt.subplot(2,1,1)
x3=np.linspace(start=0,stop=length,endpoint=True,num=50*no_elem)
exact_strain=(2*z*qo*(length-x3)**2)/ei
#print(exact_strain)
#print(exact_strain.shape)
plt.plot((x3*no_elem)/length,exact_strain,label='exact strain')
plt.plot(x,z*R2,label='fea strain')
plt.title('comparing axial strains')
plt.legend()
plt.subplot(2,1,2)
plt.plot((x3*no_elem)/length,exact_strain)
plt.title('exact strain')
plt.grid()

#start
plt.figure(4)
plt.plot((x1*no_elem)/length,exact_solution,label='exact solution')
plt.scatter(node_num-1,a,label='node_displacement(through fea)')
plt.plot(x,R,label='fea solution')
plt.legend()
plt.grid()

print('maximum deflection (by bending equation)\n',max(exact_solution))
print('minimum deflection (by bending equation)\n',min(exact_solution))
print('maximum axial strain (by fea)\n',max(z*R2))
print('minimum axial strain (by fea)\n',min(z*R2))
print('maximum axial strain (by bending equation)\n',max(exact_strain))
print('minimum axial strain (by bending equation)\n',min(exact_strain))
#stop
plt.show()