import matplotlib.pyplot as plt
import numpy as np
import os
import re
import scipy as sp
# from scipy.optimize import fsolve

# coordinates of element
x = np.empty([3,4])
x[:,0] = np.array([1.0,5.0,-3.0])
x[:,1] = np.array([10.0,4.0,3.0])
x[:,2] = np.array([-3.0,4.0,2.0])
x[:,3] = x[:,0] + (x[:,2] - x[:,1] - 0.3*(x[:,1]-x[:,0]))
x0 = np.array([4.0,4.0,0.05])

# x[:,0] = np.array([0.0,0.0,0.0])
# x[:,1] = np.array([5.0,0.0,0.0])
# x[:,2] = np.array([5.0,5.0,0.0])
# x[:,3] = np.array([0.0,5.0,0.0])
# 
# x0 = np.array([4.0,4.0,0.05])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# shift
xcenter = 0.25*(x[:,0] + x[:,1] + x[:,2] + x[:,3])
xshift = np.empty([3,4])
for i in [0,1,2,3]:
    xshift[:,i] = x[:,i] - xcenter
x0shift = x0 - xcenter

# transformation matrix
e1loc_g = x[:,1] - x[:,0]
e1loc_g = e1loc_g/np.linalg.norm(e1loc_g)
e3loc_g = np.cross(x[:,1] - x[:,0], x[:,3] - x[:,0])
e3loc_g = e3loc_g/np.linalg.norm(e3loc_g)
e2loc_g = np.cross(e3loc_g, e1loc_g)

Ql2g = np.empty([3,3])
Ql2g[:,0] = e1loc_g
Ql2g[:,1] = e2loc_g
Ql2g[:,2] = e3loc_g
Qg2l = np.transpose(Ql2g)

quad_x = np.empty([3,5])
quad_x[:,0] = xshift[:,0]
quad_x[:,1] = xshift[:,1]
quad_x[:,2] = xshift[:,2]
quad_x[:,3] = xshift[:,3]
quad_x[:,4] = quad_x[:,0]
ax.plot(quad_x[0,:], quad_x[1,:], quad_x[2,:])
ax.scatter(x0[0], x0[1], x0[2])
quad_x = np.empty([3,5])
quad_x[:,0] = x[:,0]
quad_x[:,1] = x[:,1]
quad_x[:,2] = x[:,2]
quad_x[:,3] = x[:,3]
quad_x[:,4] = quad_x[:,0]
ax.plot(quad_x[0,:], quad_x[1,:], quad_x[2,:])
ax.text(quad_x[0,0], quad_x[1,0], quad_x[2,0], '1', zorder = 1)
ax.text(quad_x[0,1], quad_x[1,1], quad_x[2,1], '2', zorder = 1)
ax.text(quad_x[0,2], quad_x[1,2], quad_x[2,2], '3', zorder = 1)
ax.text(quad_x[0,3], quad_x[1,3], quad_x[2,3], '4', zorder = 1)

ax.quiver(xcenter[0], xcenter[1], xcenter[2], e1loc_g[0], e1loc_g[1], e1loc_g[2], length=5.0, normalize=True,  color = 'red')
ax.quiver(xcenter[0], xcenter[1], xcenter[2], e2loc_g[0], e2loc_g[1], e2loc_g[2], length=5.0, normalize=True,  color = 'blue')
ax.quiver(xcenter[0], xcenter[1], xcenter[2], e3loc_g[0], e3loc_g[1], e3loc_g[2], length=5.0, normalize=True,  color = 'green')
ax.axis('equal')

# local
xshift_loc = np.empty([3,4])
xshift_loc[:,0] = np.matmul(Qg2l, xshift[:,0])
xshift_loc[:,1] = np.matmul(Qg2l, xshift[:,1])
xshift_loc[:,2] = np.matmul(Qg2l, xshift[:,2])
xshift_loc[:,3] = np.matmul(Qg2l, xshift[:,3])
x0shift_loc = np.matmul(Qg2l, x0shift)

# callculate auxiliary vectors
xa_loc = 0.25*(xshift_loc[:,0] + xshift_loc[:,1] + xshift_loc[:,2] + xshift_loc[:,3])
xb_loc = 0.25*(-xshift_loc[:,0] + xshift_loc[:,1] + xshift_loc[:,2] - xshift_loc[:,3])
xc_loc = 0.25*(-xshift_loc[:,0] - xshift_loc[:,1] + xshift_loc[:,2] + xshift_loc[:,3])
xd_loc = 0.25*(xshift_loc[:,0] - xshift_loc[:,1] + xshift_loc[:,2] - xshift_loc[:,3])
xa = 0.25*(x[:,0] + x[:,1] + x[:,2] + x[:,3])
xb = 0.25*(-x[:,0] + x[:,1] + x[:,2] - x[:,3])
xc = 0.25*(-x[:,0] - x[:,1] + x[:,2] + x[:,3])
xd = 0.25*(x[:,0] - x[:,1] + x[:,2] - x[:,3])

# transformations
def hat2true2(xhat, xa, xb, xc, xd):
    x =  [xa[0] + xb[0]*xhat[0] + xc[0]*xhat[1] + xd[0]*xhat[0]*xhat[1],
          xa[1] + xb[1]*xhat[0] + xc[1]*xhat[1] + xd[1]*xhat[0]*xhat[1]]
    return x
def hat2true3(xhat, xa, xb, xc, xd):
    x =  [xa[0] + xb[0]*xhat[0] + xc[0]*xhat[1] + xd[0]*xhat[0]*xhat[1],
          xa[1] + xb[1]*xhat[0] + xc[1]*xhat[1] + xd[1]*xhat[0]*xhat[1],
          xa[2] + xb[2]*xhat[0] + xc[2]*xhat[1] + xd[2]*xhat[0]*xhat[1]]
    return x
def bar2hat(xbar, x0hat, mu, eta, b):
    xhat = [x0hat[0] + b*np.sinh(mu[0]*xbar[0] - eta[0]),
            x0hat[1] + b*np.sinh(mu[1]*xbar[1] - eta[1])]
    # xhat = xbar # force same
    return xhat
def dtrue_dhat(xhat, xa, xb, xc, xd):
    dx_dxhat = np.empty([2,2])
    dx_dxhat[:,0] = xb[0:2] + xd[0:2]*xhat[1]
    dx_dxhat[:,1] = xc[0:2] + xd[0:2]*xhat[0]
    return dx_dxhat
def dhat_dbar(xbar, x0hat, mu, eta, b):
    dxhat_dxbar = np.zeros([2,2])
    dxhat_dxbar[0,0] = b*np.cosh(mu[0]*xbar[0] - eta[0])*mu[0]
    dxhat_dxbar[1,1] = b*np.cosh(mu[1]*xbar[1] - eta[1])*mu[1]
    # dxhat_dxbar[0,0] = 1.0 # force same
    # dxhat_dxbar[1,1] = 1.0
    return dxhat_dxbar

# find x0hat and set b, and mu and eta
def func(xhat, xa_loc, xb_loc, xc_loc, xd_loc, x0shift_loc):
    R = hat2true2(xhat, xa_loc, xb_loc, xc_loc, xd_loc) - x0shift_loc[0:2]
    return R

x0hat = sp.optimize.fsolve(func, [0.0, 0.0], args=(xa_loc, xb_loc, xc_loc, xd_loc, x0shift_loc))
b = x0shift_loc[2]
mu = 0.5*(np.arcsinh((1.0 + x0hat)/b) + np.arcsinh((1.0 - x0hat)/b))
eta = 0.5*(np.arcsinh((1.0 + x0hat)/b) - np.arcsinh((1.0 - x0hat)/b))

print('x0shift_loc',x0shift_loc)
print('x0hat',x0hat)
print('b',b)
print('mu',mu)
print('eta',eta)
# exit()

print(hat2true2(x0hat, xa_loc, xb_loc, xc_loc, xd_loc))
print(x0shift_loc)

# 1D Gaussian quadrature for transformed space
a = np.array([-np.sqrt(0.6),0.0,np.sqrt(0.6)])
w = np.array([5.0/9.0,8.0/9.0,5.0/9.0])

# integration points in bar space
xintpt_bar = np.empty([2,9])
ic = -1
for i in range(0,3):
    for j in range(0,3):
        ic = ic + 1
        xintpt_bar[0,ic] = a[i]
        xintpt_bar[1,ic] = a[j]

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.scatter(xintpt_bar[0,:], xintpt_bar[1,:])
ax1.scatter(x0hat[0], x0hat[1])
quad_x = np.empty([2,5])
quad_x[:,0] = [-1.0,-1.0]
quad_x[:,1] = [1.0,-1.0]
quad_x[:,2] = [1.0,1.0]
quad_x[:,3] = [-1.0,1.0]
quad_x[:,4] = quad_x[:,0]
ax1.plot(quad_x[0,:], quad_x[1,:])
ax1.axis('equal')

# integration points in hat space
ic = -1
xintpt_hat = np.empty([2,9])
for i in range(0,3):
    for j in range(0,3):
        ic = ic + 1
        xintpt_hat[:,ic] = bar2hat(xintpt_bar[:,ic], x0hat, mu, eta, b)

ax1.scatter(xintpt_hat[0,:], xintpt_hat[1,:])

# integration points in true space (global)
ic = -1
xintpt = np.empty([3,9])
for i in range(0,3):
    for j in range(0,3):
        ic = ic + 1
        xintpt[:,ic] = hat2true3(xintpt_hat[:,ic], xa, xb, xc, xd)
        print('xintpt[:,ic]',xintpt[:,ic])

ax.scatter(xintpt[0,:], xintpt[1,:], xintpt[2,:])

# total weights
wtot = np.empty([9])
ic = -1
for i in range(0,3):
    for j in range(0,3):
        ic = ic + 1
        dx_dxhat = dtrue_dhat(xintpt_hat[:,ic], xa_loc, xb_loc, xc_loc, xd_loc)
        dxhat_dxbar = dhat_dbar(xintpt_bar[:,ic], x0hat, mu, eta, b)
        # print(abs(np.linalg.det(dx_dxhat)))
        wtot[ic] = abs(np.linalg.det(dx_dxhat))*abs(np.linalg.det(dxhat_dxbar))*w[i]*w[j]
        print('jacob_ref2quad',abs(np.linalg.det(dx_dxhat)))
        print('jacob_sinh',abs(np.linalg.det(dxhat_dxbar)))
        print('wtot[ic]',wtot[ic])


#A = 0.5*[(x1y2  + x2y3 + x3y4 + x4y1) – (x2y1 + x3y2 + x4y3 + x1y4)]
A = 0.5*((xshift_loc[0,0]*xshift_loc[1,1] + xshift_loc[0,1]*xshift_loc[1,2] + xshift_loc[0,2]*xshift_loc[1,3] + xshift_loc[0,3]*xshift_loc[1,0]) - 
   (xshift_loc[0,1]*xshift_loc[1,0] + xshift_loc[0,2]*xshift_loc[1,1] + xshift_loc[0,3]*xshift_loc[1,2] + xshift_loc[0,0]*xshift_loc[1,3]))
print(np.sum(wtot)) # sum of weights should be close to area
print(A)

# evaluate integral 1/r
intg = 0.0
ic = -1
for i in range(0,3):
    for j in range(0,3):
        ic = ic + 1
        intg = intg + 1.0/np.linalg.norm(x0 - xintpt[:,ic])*wtot[ic]

print(intg)


plt.show()


exit()
