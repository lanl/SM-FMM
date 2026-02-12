import matplotlib.pyplot as plt
import numpy as np
import os
import re

# coordinates of element
x = np.empty([3,3])
# x[:,0] = np.array([1.0,5.0,-3.0])
# x[:,1] = np.array([10.0,4.0,3.0])
# x[:,2] = np.array([-3.0,4.0,2.0])

x[:,0] = np.array([0.0,0.0,0.0])
x[:,1] = np.array([5.0,0.0,0.0])
x[:,2] = np.array([5.0,5.0,0.0])

x0 = np.array([4.0,1.0,0.5])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# shift
xshift = np.empty([3,3])
for i in [0,1,2]:
    xshift[:,i] = x[:,i] - x[:,0]
x0shift = x0 - x[:,0]

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
ax.scatter(x0[0], x0[1], x0[2])
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
x0shift_loc = np.matmul(Qg2l, x0shift)

# coordinates in reference space
xhat = np.empty([2,3])
xhat[:,0] = np.array([0.0,0.0])
xhat[:,1] = np.array([2.0,0.0])
xhat[:,2] = np.array([2.0,2.0])

print(xshift_loc[0:2,1:3])
A = np.matmul(xhat[:,1:3], np.linalg.inv(xshift_loc[0:2,1:3]))
Ainv = np.linalg.inv(A)

x0hat = np.matmul(A, x0shift_loc[0:2])
b = x0shift_loc[2] # distance of source point from plane of triangle

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

# sinh transformation of triangle - stays triangle
eta = np.empty(2)
eta = - np.arcsinh(-x0hat/b)
mu = np.empty(2)
mu = 0.5*(np.arcsinh((2.0-x0hat)/b) + eta)
x_grid = np.empty([20*20,2])
x_grid_trans = np.empty([20*20,2])
ic = -1
for i in range(0,20):
    for j in range(0,20):
        xtmp = np.array([float(i)*0.1,float(j)*0.1])
        if (xtmp[1] <= xtmp[0]):
            ic = ic + 1
            x_grid[ic,:] = xtmp
            x_grid_trans[ic,0] = x0hat[0] + b*np.sinh(xtmp[0]*mu[0] - eta[0])
            x_grid_trans[ic,1] = x0hat[1] + b*np.sinh(xtmp[1]*mu[1] - eta[1])
            
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.scatter(x_grid[0:ic+1,0],x_grid[0:ic+1,1])
ax2.scatter(x_grid_trans[0:ic+1,0],x_grid_trans[0:ic+1,1])
ax2.axis('equal')

# 2D Gaussian quadrature for transformed space (triangle)
xbar = np.empty([2,3])
for i in range(0,3):
    xbar[:,i] = x0hat + b*np.sinh(xhat[:,i]*mu - eta)

xintpt_bar = np.empty([2,4])
w = np.array([-27.0/48.0,25.0/48.0,25.0/48.0,25.0/48.0])
wtot = np.empty(4)
traingular_coord = np.empty([3,4])
traingular_coord[:,0] = np.array([1.0/3.0,1.0/3.0,1.0/3.0])
traingular_coord[:,1] = np.array([0.6,0.2,0.2])
traingular_coord[:,2] = np.array([0.2,0.6,0.2])
traingular_coord[:,3] = np.array([0.2,0.2,0.6])

area_hat = 2.0
for i in range(0,4):
    xintpt_bar[:,i] = traingular_coord[0,i]*xbar[:,0] + traingular_coord[1,i]*xbar[:,1] + traingular_coord[2,i]*xbar[:,2]
    dx1hat_dx1bar = 1.0/mu[0]/b/np.sqrt(1.0 + np.power((xintpt_bar[0,i] - x0hat[0])/b,2))
    dx2hat_dx2bar = 1.0/mu[1]/b/np.sqrt(1.0 + np.power((xintpt_bar[1,i] - x0hat[1])/b,2))
    wtot[i] = w[i]*abs(jacob)*abs(dx1hat_dx1bar*dx2hat_dx2bar)*area_hat

print(wtot)

triangle_x[:,0] = xbar[0:2,0]
triangle_x[:,1] = xbar[0:2,1]
triangle_x[:,2] = xbar[0:2,2]
triangle_x[:,3] = triangle_x[:,0]
ax1.plot(triangle_x[0,:], triangle_x[1,:])
ax1.text(triangle_x[0,0], triangle_x[1,0], '1', zorder = 1)
ax1.text(triangle_x[0,1], triangle_x[1,1], '2', zorder = 1)
ax1.text(triangle_x[0,2], triangle_x[1,2], '3', zorder = 1)
ax1.scatter(xintpt_bar[0,:],xintpt_bar[1,:])

# integration points to standard triangle
xintpt_hat = np.empty([2,4])
for i in range(0,4):
    xintpt_hat[:,i] = (np.arcsinh((xintpt_bar[:,i] - x0hat)/b) + eta)/mu

print(xintpt_bar)
print(x0hat)
print(b)

ax1.scatter(xintpt_hat[0,:], xintpt_hat[1,:])

# integration points 2D triangle
xintpt_2D = np.empty([2,4])
for i in range(0,4):
    xintpt_2D[:,i] = np.matmul(Ainv, xintpt_hat[:,i])

ax1.scatter(xintpt_2D[0,:], xintpt_2D[1,:])

# integration points in 3D
xintpt = np.empty([3,4])
xintpt[0:2,:] = xintpt_2D[0:2,:]
xintpt[2,:] = np.zeros([1,4])
for i in range(0,4):
    xintpt[:,i] = np.matmul(Ql2g, xintpt[:,i]) + x[:,0]

ax.scatter(xintpt[0,:], xintpt[1,:], xintpt[2,:])


print(xintpt_hat)
print(xintpt)

# evaluate integral 1/r
intg = 0.0
for i in range(0,4):
    intg = intg + 1.0/np.linalg.norm(xintpt[:,i] - x0)*wtot[i]

print(intg)
plt.show()