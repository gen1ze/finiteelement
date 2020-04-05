import numpy  as np
n=int(input('give no of elements'))
l=list(map(float,input('give k1 k2 etc upto kn seperated by spaces').split()))
k=np.array([[1,-1],[-1,1]])
K=np.zeros((n+1,n+1))
for i in range(n):
    k1=k*l[i]
    K_temp=np.zeros((n+1,n+1))
    K_temp[i:i+2,i:i+2]=k1
    K+=K_temp
print('global stiffness matrix is\n',K)
print(K.shape)


#getting global force matrix( considering point forces only)
print('no of force to be given is equal to n+1')
print('give force at first node as -R')
f=list(map(float,input('give forces at node seperated by spaces').split()))
F=np.array(f).reshape(n+1,1)
print('global force matrix is\n',F)
print(F.shape)


u=np.linalg.solve(K[1:,1:],F[1:,:])
U=np.insert(u,0,0,axis=0)
print('deformation is\n',U)
print(U.shape)
