import matplotlib.pyplot as plt
import numpy as np
import os
import re

# coordinates of element
x = np.empty([3,3])
# x[:,0] = np.array([1.0,5.0,-3.0])
# x[:,1] = np.array([10.0,4.0,3.0])
# x[:,2] = np.array([-3.0,4.0,2.0])

x[:,0] = np.array([0.0,2.0,2.0])
x[:,1] = np.array([0.0,7.0,2.0])
x[:,2] = np.array([0.0,7.0,7.0])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# shift
xshift = np.empty([3,3])
for i in [0,1,2]:
    xshift[:,i] = x[:,i] - x[:,0]

print(x[:,0])
print(x[:,1])
print(x[:,2])

print(xshift[:,0])
print(xshift[:,1])
print(xshift[:,2])

# transformation matrix
e1loc_g = x[:,1] - x[:,0]
e1loc_g = e1loc_g/np.linalg.norm(e1loc_g)
e3loc_g = np.cross(xshift[:,1], xshift[:,2])
e3loc_g = e3loc_g/np.linalg.norm(e3loc_g)
e2loc_g = np.cross(e3loc_g, e1loc_g)

Ql2g = np.empty([3,3])
Ql2g[:,0] = e1loc_g
Ql2g[:,1] = e2loc_g
Ql2g[:,2] = e3loc_g
Qg2l = np.transpose(Ql2g)

triangle_x = np.empty([3,4])
triangle_x[:,0] = xshift[:,0]
triangle_x[:,1] = xshift[:,1]
triangle_x[:,2] = xshift[:,2]
triangle_x[:,3] = triangle_x[:,0]
ax.plot(triangle_x[0,:], triangle_x[1,:], triangle_x[2,:])
triangle_x = np.empty([3,4])
triangle_x[:,0] = x[:,0]
triangle_x[:,1] = x[:,1]
triangle_x[:,2] = x[:,2]
triangle_x[:,3] = triangle_x[:,0]
ax.plot(triangle_x[0,:], triangle_x[1,:], triangle_x[2,:])
ax.text(triangle_x[0,0], triangle_x[1,0], triangle_x[2,0], '1', zorder = 1)
ax.text(triangle_x[0,1], triangle_x[1,1], triangle_x[2,1], '2', zorder = 1)
ax.text(triangle_x[0,2], triangle_x[1,2], triangle_x[2,2], '3', zorder = 1)

ax.quiver(x[0,0], x[1,0], x[2,0], e1loc_g[0], e1loc_g[1], e1loc_g[2], length=5.0, normalize=True,  color = 'red')
ax.quiver(x[0,0], x[1,0], x[2,0], e2loc_g[0], e2loc_g[1], e2loc_g[2], length=5.0, normalize=True,  color = 'blue')
ax.quiver(x[0,0], x[1,0], x[2,0], e3loc_g[0], e3loc_g[1], e3loc_g[2], length=5.0, normalize=True,  color = 'green')
ax.axis('equal')

# local
xshift_loc = np.empty([3,3])
xshift_loc[:,0] = np.matmul(Qg2l, xshift[:,0])
xshift_loc[:,1] = np.matmul(Qg2l, xshift[:,1])
xshift_loc[:,2] = np.matmul(Qg2l, xshift[:,2])

# coordinates in reference space
xhat = np.empty([2,3])
xhat[:,0] = np.array([0.0,0.0])
xhat[:,1] = np.array([2.0,0.0])
xhat[:,2] = np.array([2.0,2.0])

print(xshift_loc[0:2,1:3])
A = np.matmul(xhat[:,1:3], np.linalg.inv(xshift_loc[0:2,1:3]))
Ainv = np.linalg.inv(A)

# test the transformation matrix A
xhat_test = np.empty([2,3])
xhat_test[:,0] = np.matmul(A, xshift_loc[0:2,0])
xhat_test[:,1] = np.matmul(A, xshift_loc[0:2,1])
xhat_test[:,2] = np.matmul(A, xshift_loc[0:2,2])

# ax.scatter(x4[0,:], x4[1,:], x4[2,:], marker='o', color = 'r')

# A = np.matmul(xhat234, np.linalg.inv(x234shift))

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

triangle_x = np.empty([2,4])
triangle_x[:,0] = xshift_loc[0:2,0]
triangle_x[:,1] = xshift_loc[0:2,1]
triangle_x[:,2] = xshift_loc[0:2,2]
triangle_x[:,3] = triangle_x[:,0]
ax1.plot(triangle_x[0,:], triangle_x[1,:])
ax1.text(triangle_x[0,0], triangle_x[1,0], '1', zorder = 1)
ax1.text(triangle_x[0,1], triangle_x[1,1], '2', zorder = 1)
ax1.text(triangle_x[0,2], triangle_x[1,2], '3', zorder = 1)
triangle_x[:,0] = xhat_test[0:2,0]
triangle_x[:,1] = xhat_test[0:2,1]
triangle_x[:,2] = xhat_test[0:2,2]
triangle_x[:,3] = triangle_x[:,0]
ax1.plot(triangle_x[0,:], triangle_x[1,:])
ax1.text(triangle_x[0,0], triangle_x[1,0], '1', zorder = 1)
ax1.text(triangle_x[0,1], triangle_x[1,1], '2', zorder = 1)
ax1.text(triangle_x[0,2], triangle_x[1,2], '3', zorder = 1)
ax1.axis('equal')

dx_dxhat = np.linalg.inv(A)
jacob = np.linalg.det(dx_dxhat)
print(jacob)

# 1D Gaussian quadrature for transformed space
a = np.array([-np.sqrt(0.6),0.0,np.sqrt(0.6)])
w = np.array([5.0/9.0,8.0/9.0,5.0/9.0])

xintpt_bar = np.empty([2,9])
wtot = np.empty([9])
ic = -1
for i in range(0,3):
    for j in range(0,3):
        ic = ic + 1
        xintpt_bar[0,ic] = a[i]
        xintpt_bar[1,ic] = a[j]
        wtot[ic] = jacob*0.5*(a[i] + 1.0)*w[i]*w[j]

ax1.scatter(xintpt_bar[0,:], xintpt_bar[1,:])

# integration points to standard triangle
xintpt_hat = np.empty([2,9])
xintpt_hat[0,:] = xintpt_bar[0,:] + 1.0
xintpt_hat[1,:] = 0.5*(xintpt_bar[0,:] + 1.0)*(xintpt_bar[1,:] + 1.0)

ax1.scatter(xintpt_hat[0,:], xintpt_hat[1,:])

# integration points 2D triangle
xintpt_2D = np.empty([2,9])
for i in range(0,9):
    xintpt_2D[:,i] = np.matmul(Ainv, xintpt_hat[:,i])

ax1.scatter(xintpt_2D[0,:], xintpt_2D[1,:])

# integration points in 3D
xintpt = np.empty([3,9])
xintpt[0:2,:] = xintpt_2D[0:2,:]
xintpt[2,:] = np.zeros([1,9])
for i in range(0,9):
    xintpt[:,i] = np.matmul(Ql2g, xintpt[:,i]) + x[:,0]

ax.scatter(xintpt[0,:], xintpt[1,:], xintpt[2,:])

# evaluate integral 1/r
intg = 0.0
for i in range(0,9):
    intg = intg + 1.0/np.linalg.norm(xintpt[:,i] - x[:,0])*wtot[i]

print(wtot)
print(intg)


plt.show()