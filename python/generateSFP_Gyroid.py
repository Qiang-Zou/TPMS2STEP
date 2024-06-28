import numpy as np
from mpmath import mp
import argparse
import utility

def spaceLength(a, b):
    length = mp.sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))
    return length

if __name__ == "__main__":
    file = open('../data/uv/uv_Gyroid.off', 'r')
    file.readline()
    file.readline()
    text = file.readline()
    uvList = []
    while text:
        x,y,z = text.split(" ")
        text = file.readline()
        temp = []
        temp.append(float(x))
        temp.append(float(y))
        temp.append(float(z))
        uvList.append(temp)
    uvArray = np.array(uvList)
    uvParam = np.mat(uvArray)
    file.close()

    ttt = []
    K = mp.ellipk(0.25)
    K_d = mp.ellipk(0.75)
    K_d_d = K*K+K_d*K_d
    ka = mp.sqrt(K_d_d)/(K*K_d)
    for item in uvList:
        complexuv = mp.mpc(item[0], item[1])
        x_i_j = mp.asin((2.0*mp.sqrt(2)*complexuv)/mp.sqrt(1.0+4.0*complexuv*complexuv+complexuv*complexuv*complexuv*complexuv))
        y_i_j = mp.asin(((-1.0)*2.0*mp.sqrt(2)*complexuv)/mp.sqrt(1.0+4.0*complexuv*complexuv+complexuv*complexuv*complexuv*complexuv))
        z_i_j = mp.asin((4.0*complexuv*complexuv)/(1.0+complexuv*complexuv*complexuv*complexuv))
        theta = mp.pi/2.0-mp.atan(K_d/K)
        complex_theta = mp.mpc(mp.cos(theta), mp.sin(theta))
        ellipfx = mp.ellipf(x_i_j, 0.25)*complex_theta
        ellipfy = mp.ellipf(y_i_j, 0.75)*complex_theta
        ellipfz = mp.ellipf(z_i_j, 97.2465-56*mp.sqrt(3))*complex_theta
        # ellipfz = mp.ellipf(z_i_j, 97-56*mp.sqrt(3))*complex_theta
        fxuv = ka*(ellipfx.real)/(2.0*mp.sqrt(2))
        fyuv = ka*(ellipfy.imag)/(2.0*mp.sqrt(2))
        fzuv = ka*(ellipfz.real)*0.25
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

    n = int(mp.sqrt(len(uvList)))-1
    # save the file for uniform sampling
    n2 = n+1
    with open("../data/uv/sampleParams.txt", 'w') as fp:
        for i in range(n2-1):
            fp.write(str(spaceLength(ttt[n2*n2-1-(i+1)*(i+1)-2*(i+1)],ttt[n2*n2-1-i*i-2*i])))
            fp.write("\n")
        for i in range(n2-1):
            fp.write(str(spaceLength(ttt[n2*n2-1-(i+1)*(i+1)-(i+1)],ttt[n2*n2-1-i*i-i])))
            fp.write("\n")
        fp.close()

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--offset" , help="offset")
    args = parser.parse_args()
    d = float(args.offset)*(-1)
    offset_1 = []
    offset_2 = []
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
            new_item_1 = []
            new_item_2 = []
            if length!=0:
                new_item_1.append(ttt[cursor][0]+d*nx/length)
                new_item_1.append(ttt[cursor][1]+d*ny/length)
                new_item_1.append(ttt[cursor][2]+d*nz/length)
                new_item_2.append(ttt[cursor][0]-d*nx/length)
                new_item_2.append(ttt[cursor][1]-d*ny/length)
                new_item_2.append(ttt[cursor][2]-d*nz/length)
            else:
                new_item_1.append(0)
                new_item_1.append(0)
                new_item_1.append(0)
                new_item_2.append(0)
                new_item_2.append(0)
                new_item_2.append(0)
            offset_1.append(new_item_1)
            offset_2.append(new_item_2)
            cursor = cursor+1
    arr_1 = np.array(offset_1)
    arr_2 = np.array(offset_2)
    matrix_1 = np.mat(arr_1)
    matrix_2 = np.mat(arr_2)

    fp_1 = matrix_1*utility.Rz(mp.pi*(-0.25))*utility.Ty + np.array([1,0.5,0.25])
    fp_2 = matrix_2*utility.Rz(mp.pi*0.25)*utility.Tyz + np.array([1,0.5,0.25])
    fp_3 = matrix_1*utility.Rz(mp.pi*(-0.25))*utility.Tx + np.array([1,0.5,0.25])
    fp_4 = matrix_2*utility.Rz(mp.pi*0.25)*utility.Txz + np.array([1,0.5,0.25])

    #print(n)
    cursor = 0
    data = []
    for i in range(2*n+1):
        temp = []
        temp2 = []
        for j in range(2*n+1):
            temp.append([])
            temp2.append([])
        data.append(temp)
    for item in range(n+1):
        i = n-item
        for j in range(2*i+1):
            points = [fp_1[cursor], fp_2[cursor], fp_3[cursor], fp_4[cursor]]
            cursor = cursor + 1
            if j != 2*i:
                x0 = n-(i-j)
                y0 = n-i
                data[x0][y0].append(points[0][0,0])
                data[x0][y0].append(points[0][0,1])
                data[x0][y0].append(points[0][0,2])
                x1 = n+i
                y1 = n-(i-j)
                data[x1][y1].append(points[1][0,0])
                data[x1][y1].append(points[1][0,1])
                data[x1][y1].append(points[1][0,2])
                x2 = n-(j-i)
                y2 = n+i
                data[x2][y2].append(points[2][0,0])
                data[x2][y2].append(points[2][0,1])
                data[x2][y2].append(points[2][0,2])
                x3 = n-i
                y3 = n-(j-i)
                data[x3][y3].append(points[3][0,0])
                data[x3][y3].append(points[3][0,1])
                data[x3][y3].append(points[3][0,2])
    data[n][n] = np.array([1,0.5,0.25-d])

    dat = []
    for i in range(len(data)):
        for j in range(len(data[i])):
            temp = data[i][j]
            temp[0] = temp[0]-1
            temp[1] = temp[1]-0.5
            dat.append(temp)
    arr = np.array(dat)
    matrix = np.mat(arr)

    fileName = "../data/SFP/Gyroid/SFP_Gyroid_offset_"+ str(int(100*d)) +".off"
    fileHeader = 'OFF\n'+str((2*n+1)*((2*n+1)))+' 0 0'
    with open(fileName, "w") as fp:
        np.savetxt(fp, matrix, header=fileHeader, comments='')
        fp.close()
