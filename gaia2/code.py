import math
import numpy as np


def calibration(x,y,ksi,eta,n,Nmin,rmax):
    Q = int((n+1)*(n+2)/2)
    C = np.ones((np.size(x), Q))
    q = 0
    for i in range(n + 1):
        for j in range(n + 1):
            if (i + j <= n):
                C[:, q] = (x ** i) * (y ** j)
                q += 1
    We = np.diag(np.ones(np.size(x)))
    flag = 0
    it = 0
    while (flag == 0):
        Zx = np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(C), We), C)),\
                    np.dot(np.dot(np.transpose(C), We), ksi))
        Zy = np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(C), We), C)),\
                    np.dot(np.dot(np.transpose(C), We), eta))
        rx = np.dot(We, ksi - np.dot(C, Zx))
        ry = np.dot(We, eta - np.dot(C, Zy))
        r = np.sqrt(rx ** 2 + ry ** 2)
        kmax = np.argmax(r)
        flag = 1
        if (np.size(ksi) - it <= Nmin):
            break
        if (r[kmax] > rmax):
            We[kmax, kmax] = 0
            flag = 0
            it += 1
    uwex = np.dot(np.dot(np.transpose(rx), We), rx) / (np.size(x) - it - Q)
    uwey = np.dot(np.dot(np.transpose(ry), We), ry) / (np.size(y) - it - Q)
    return Zx,Zy,math.sqrt(uwex),math.sqrt(uwey),np.size(x)-it


def transform(x,y,Z):
    Q = np.size(Z)
    n = int(np.max(np.roots(np.array([1,3,-2*(Q-1)]))))
    C = np.ones((np.size(x), Q))
    q = 0
    for i in range(n+1):
        for j in range(n+1):
            if (i + j <= n):
                C[:, q] = (x ** i) * (y ** j)
                q += 1
    return np.dot(C,Z)