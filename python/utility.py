import numpy as np
from mpmath import mp

def Rx(phi):
    data = np.zeros([3,3])
    data[0][0] = 1
    data[1][1] = mp.cos(phi)
    data[1][2] = mp.sin(phi)*(-1)
    data[2][1] = mp.sin(phi)
    data[2][2] = mp.cos(phi)
    return np.mat(data)

def Ry(phi):
    data = np.zeros([3,3])
    data[0][0] = mp.cos(phi)
    data[0][2] = mp.sin(phi)
    data[1][1] = 1
    data[2][0] = (-1)*mp.sin(phi)
    data[2][2] = mp.cos(phi)
    return np.mat(data)

def Rz(phi):
    data = np.zeros([3,3])
    data[0][0] = mp.cos(phi)
    data[0][1] = mp.sin(phi)*(-1)
    data[1][0] = mp.sin(phi)
    data[1][1] = mp.cos(phi)
    data[2][2] = 1
    return np.mat(data)

Tx = np.mat([[-1,0,0],[0,1,0],[0,0,1]])
Ty = np.mat([[1,0,0],[0,-1,0],[0,0,1]])
Tz = np.mat([[1,0,0],[0,1,0],[0,0,-1]])
Txz = np.multiply(Tx,Tz)
Tyz = np.multiply(Ty,Tz)
Txy = np.multiply(Tx,Ty)

def weierstrassFunction(tao)->complex:
    return complex(1,0)/mp.sqrt(tao*tao*tao*tao*tao*tao*tao*tao-complex(14,0)*tao*tao*tao*tao+1)
