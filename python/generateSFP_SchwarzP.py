import numpy as np
from mpmath import mp
import argparse
import utility

if __name__ == "__main__":
    file2 = open('../data/uv/uv_SchwarzP.off', 'r')
    file2.readline()
    file2.readline()
    text = file2.readline()
    uvList = []
    while text:
        x,y,z = text.split(" ")
        text = file2.readline()
        temp = []
        temp.append(float(x))
        temp.append(float(y))
        temp.append(float(z))
        uvList.append(temp)
    file2.close()

    temp = []
    for item in uvList:
        temp.append([item[0],(-1)*item[1]])
    for item in temp:
        uvList.append(item)
    
    K = mp.ellipk(0.25)
    K_d = mp.ellipk(0.75)
    ka = 2.0/K_d
    ttt = []
    for item in uvList:
        complexuv = mp.mpc(item[0], item[1])
        x_i_j = mp.asin((2.0*mp.sqrt(2)*complexuv)/mp.sqrt(1.0+4.0*complexuv*complexuv+complexuv*complexuv*complexuv*complexuv))
        y_i_j = mp.asin(((-1.0)*2.0*mp.sqrt(2)*complexuv)/mp.sqrt(1.0+4.0*complexuv*complexuv+complexuv*complexuv*complexuv*complexuv))
        z_i_j = mp.asin((4.0*complexuv*complexuv)/(1.0+complexuv*complexuv*complexuv*complexuv))
        ellipfx = mp.ellipf(x_i_j, 0.25)
        ellipfy = mp.ellipf(y_i_j, 0.75)
        ellipfz = mp.ellipf(z_i_j, 97.2465-56*mp.sqrt(3))
        fxuv = (-1)*ka*(ellipfx.imag)/(2.0*mp.sqrt(2))
        fyuv = ka*(ellipfy.real)/(2.0*mp.sqrt(2))
        fzuv = (-1)*ka*(ellipfz.imag)*0.25
        if fxuv>(-1)*1e-5 and fxuv<1e-5:
            fxuv = 0
        if fyuv>(-1)*1e-5 and fyuv<1e-5:
            fyuv = 0
        if fzuv>(-1)*1e-5 and fzuv<1e-5:
            fzuv = 0
        temp_fp = []
        temp_fp.append(fxuv)
        temp_fp.append(fyuv)
        temp_fp.append(fzuv)
        ttt.append(temp_fp)

# calculate offset
    new_ttt = []
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--offset" , help="offset")
    args = parser.parse_args()
    d = float(args.offset)*(-1)
    cursor = 0
    for item in uvList:
        tao = mp.mpc(item[0], item[1])
        r = utility.weierstrassFunction(tao)
        phi = []
        phi.append(complex(1,0)*(complex(1,0)-tao*tao)*r)
        phi.append(complex(0,1)*(complex(1,0)+tao*tao)*r)
        phi.append(2*tao*r)
        nx = phi[1]*mp.conj(phi[2])
        nx = nx.imag*2
        ny = phi[2]*mp.conj(phi[0])
        ny = ny.imag*2
        nz = phi[0]*mp.conj(phi[1])
        nz = nz.imag*2
        length = mp.sqrt(nx*nx+ny*ny+nz*nz)
        new_item = []
        if length!=0:
            new_item.append(ttt[cursor][0]+d*nx/length)
            new_item.append(ttt[cursor][1]+d*ny/length)
            new_item.append(ttt[cursor][2]+d*nz/length)
        else:
            new_item.append(0)
            new_item.append(0)
            new_item.append(0)
        new_ttt.append(new_item)
        cursor = cursor+1
# mapping
    arr = np.array(new_ttt)
    matrix = np.mat(arr)
    fp_1 = matrix

    n = int(mp.sqrt(len(uvList)))-1
    cursor = 0
    data = []
    for i in range(n+1):
        temp = []
        for j in range(n+1):
            temp.append([])
        data.append(temp)
    for i in range(n+1):
        for j in range(i,n+1):
            data[i][j].append(fp_1[cursor+int(len(uvList)/2)][0,0])
            data[i][j].append(fp_1[cursor+int(len(uvList)/2)][0,1])
            data[i][j].append(fp_1[cursor+int(len(uvList)/2)][0,2])
            if i != j:
                data[j][i].append(fp_1[cursor][0,0])
                data[j][i].append(fp_1[cursor][0,1])
                data[j][i].append(fp_1[cursor][0,2])
            cursor = cursor + 1

    fileName = "../data/SFP/SchwarzP/SFP_SchwarzP_offset_"+ str(int(100*d)) +".off"
    with open(fileName, 'w') as fp:
        fp.write('OFF\n'+str((n+1)*((n+1)))+' 0 0\n')
        for i in range(n+1):
            for j in range(n+1):
                fp.write(str(data[i][j][0]))
                fp.write(" ")
                fp.write(str(data[i][j][1]))
                fp.write(" ")
                fp.write(str(data[i][j][2]))
                fp.write('\n')
        fp.close()
