import numpy as np
import time
class GMRES:

    def __init__(self):
        self.h = []
        self.v = []
        self.beta = [] 
        self.j = []

    def GMREShook(self, m, del_value):
        # Called by GMRESm.  Return hookstep vector y s.t. approx |y|=del.
        #  c.f. Viswanath (2008) arXiv:0809.1498
        a = self.h[:self.j+1, :self.j]

        u, s, v = np.linalg.svd(a)
        #u, s, v = svd(a)
        #s = np.diag(s)
        p = self.beta * u[0, :self.j].T
        #print("size(s) = " + str(s.shape) + ", size(u) = " + str(u.shape) + ", size(v) = " + str(v.shape) + ", size(p) = " + str(p.shape))
   #     print("s = " + str(s))
        mu = max(s[self.j-1] * s[self.j-1] * 1e-6, 1e-99)
  #      print("mu = " + str(mu))
        qn = 1e99

        while qn > del_value:
            mu = mu * 1.1
            q = p * s / (mu + s * s)
            qn = np.sqrt(np.sum(q * q))
           # print(q)
           # print(qn)

        y = np.matmul(v.T, q)
    #    print("size(v) = " + str(v.shape) + ", v = " + str(v))

        p = -np.dot(self.h[:self.j+1, :self.j], y)
        p[0] = p[0] + self.beta
        del_value = np.sqrt(np.sum(p * p))
   #     print("y = " + str(y) + ", p = " + str(p) + ", del_value = " + str(del_value))
   #     time.sleep(6)

        return y, del_value


    def GMRESm(self, m, n, x, b, matvec, psolve, dotprd, res, del_value, its, info):

        if info == 2:
            y, del_value = self.GMREShook(m, del_value)
            z = np.dot(self.v[:, :self.j], y[:self.j])
            x = psolve(n, z)
            info = 0
            return x, res, del_value, its, info

        tol = res
        imx = its
        its = 0
        self.v = np.zeros((n, m + 1))

        while True:  # Restart
            res_ = 1e99
            stgn = 1.0 - 1e-14

            self.beta = np.sqrt(dotprd(n, x, x))
            if self.beta == 0.0:
                w = np.zeros(n)
            else:
                w = matvec(n, x)
            w = b - w
            self.beta = np.sqrt(dotprd(n, w, w))
            self.v[:, 0] = w / self.beta

            self.h = np.zeros((m + 1, m))
            for j in range(m):
                
                self.j = j+1
                its += 1
                z = self.v[:, j]
                #print("z = " +str(z))
                z = psolve(n, z)
                #print("z = " +str(z))
                w = matvec(n, z)
                #print("w = " +str(w))

                for i in range(j+1):
                    self.h[i, j] = dotprd(n, w, self.v[:, i])
                 #   print("h[i,j] = " +str(self.h[i,j]))
                    w = w - self.h[i, j] * self.v[:, i]
                 #   print("w = " +str(w))

                self.h[j + 1, j] = np.sqrt(dotprd(n, w, w))
                #print("h[j+1,j] = " +str(self.h[j+1,j]))
                self.v[:, j + 1] = w / self.h[j + 1, j]
                #print("v[:,j+1] = " +str(self.v[:,j+1]))

                p = np.zeros(self.j + 1)
                #print("p" + str(p))
                p[0] = self.beta
                #print("p" + str(p))
                h_ = self.h[0:self.j + 1, 0:self.j]
                #print("h_" + str(h_))
                y = np.dot(np.linalg.pinv(h_), p)
                #print("y" + str(y))

                p = -np.dot(self.h[0:self.j + 1, 0:self.j], y)
                #print("p" + str(p))
                p[0] = p[0] + self.beta
                #print("p" + str(p))
                res = np.sqrt(np.sum(p * p))
                #print("res" + str(res))
                if info == 1:
                    print(f'gmresm: it={its}  res={res}')
                #time.sleep(1)
                #time.sleep(10)
      #          print("Res: " + str(res) + ", tol: " + str(tol) + ", its: " + str(its) + ", imx: " + str(imx) + ", res_: " + str(res_) + ", del_value: " + str(del_value) + ", j = " + str(j))
                #time.sleep(2)               
                done = (res <= tol) or (its == imx) or (res > res_)
      #          print(done)
                if done or self.j == m:
                    #print("done: " + done + ", j = " + self.j + ", m = " + m + ", del_value = " + del_value + "its")
                    if del_value > 0.0:

#                        print("info: " + str(info) + ", m: " + str(m) + ", j: " + str(j))

                        y, del_value = self.GMREShook(m, del_value)
                    z = np.dot(self.v[:, :j], y[:j])
                    z = psolve(n, z)
                    x = x + z
                    if its == imx:
                        info = 2
                    if res > res_:
                        info = 1
                    if res <= tol:
                        info = 0
                    if done:
                        return x, res, del_value, its, info
                    if del_value > 0.0:
                        print('gmres: WARNING: m too small. restart affects hookstep.')

                res_ = res * stgn
        
            