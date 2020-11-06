import numpy as np
class solver():
    def  __init__(self,ap,aw,ae,ad,au,ab,af,ap0,b,c):
        self.ap=ap
        self.aw=aw
        self.ae=ae
        self.ad=ad
        self.au=au
        self.ab=ab
        self.af=af
        self.ap0=ap0
        self.b=b
        self.c=c


    def fun(self,x,y,z):
        return x*self.ap.shape[1]*self.ap.shape[2]+y*self.ap.shape[2]+z

    def solve1(self):
        print("Using Gauss seedlel method")
        runtime=1
        dt=0.1

        self.phi=np.zeros((int(runtime/dt)+1,self.ap.shape[0],self.ap.shape[1],self.ap.shape[2]))

        self.phi[:,0,:,:]=10
#         print(self.phi[1,0,:,:])

        for t in range(1,self.phi.shape[0]):
            xi=self.phi[t,:,:,:].flatten()
            xf=np.copy(xi)
            itr=1

#             print("calculating for t=",t*dt," sec")
            while True:
                for x in range(1,self.phi.shape[1]):
                    for y in range(1,self.phi.shape[2]):
                        for z in range(1,self.phi.shape[3]):
                            i=self.fun(x,y,z)
                            xf[i]=self.b[x,y,z]+self.ap0[x,y,z]*self.phi[t-1,x,y,z]
                            if(x>0):
                                xf[i]+=self.aw[x,y,z]*xf[self.fun(x-1,y,z)]
                            if(y>0):
                                xf[i]+=self.ad[x,y,z]*xf[self.fun(x,y-1,z)]
                            if(z>0):
                                xf[i]+=self.ab[x,y,z]*xf[self.fun(x,y,z-1)]
                            if x < self.phi.shape[1]-1:
                                 xf[i]+=self.ae[x,y,z]*xi[self.fun(x+1,y,z)]
                            if y < self.phi.shape[2]-1:
                                 xf[i]+=self.au[x,y,z]*xi[self.fun(x,y+1,z)]
                            if z < self.phi.shape[3]-1:
                                 xf[i]+=self.af[x,y,z]*xi[self.fun(x,y,z+1)]
                            xf[i]/=self.ap[x,y,z]
                if np.max(np.absolute(xf-xi))>1e-6:
                    xi=np.copy(xf)
                    itr+=1
                    if itr>100:
                        break
                else:
                    print("for time t:",t," Itr:",itr)
                    break
            self.phi[t,:,:,:]=np.resize(xf,(self.phi.shape[1],self.phi.shape[2],self.phi.shape[3]))
        return self.phi




    def solve2(self):
        print("TDMA Gauss seedle")
        runtime=1
        dt=0.1

        self.phi=np.zeros((int(runtime/dt)+1,self.ap.shape[0]+1,self.ap.shape[1]+1,self.ap.shape[2]+1))
        P,Q=np.array([np.zeros(self.ap.shape[2],) for _ in range(2)])
        error=np.zeros((self.ap.shape[0]-1,self.ap.shape[1]-1,self.ap.shape[2]-1),)
        self.phi[:,0,:,:]=10
#         print(self.phi[1,0,:,:])

        for t in range(1,self.phi.shape[0]):
            # print("calculating for t=",t*dt," sec")
            itr=1
            error_init=0
            # if(t>1):
            #     break
            while True:
                for x in range(1,self.phi.shape[1]-1):
                    for y in range(1,self.phi.shape[2]-1):
                        P[0]=0
                        Q[0]=self.phi[t,x,y,0]
                        for z in range(1,self.phi.shape[3]-1):
                            tp=self.ap[x,y,z]-self.ab[x,y,z]*P[z-1]
                            P[z]=self.af[x,y,z]/tp
                            Q[z]=self.b[x,y,x]+self.ap0[x,y,z]*self.phi[t-1,x,y,z]+self.ab[x,y,z]*Q[z-1]
                            Q[z]+=self.aw[x,y,x]*self.phi[t,x-1,y,z]+self.ae[x,y,z]*self.phi[t,x+1,y,z]
                            Q[z]+=self.ad[x,y,x]*self.phi[t,x,y-1,z]+self.au[x,y,z]*self.phi[t,x,y+1,z]
                            Q[z]/=tp
                        for z in range(self.phi.shape[3]-2,0,-1):
                            self.phi[t,x,y,z]=P[z]*self.phi[t,x,y,z+1]+Q[z]
                error=self.b[1:,1:,1:]+np.multiply(self.ap0[1:,1:,1],self.phi[t-1,1:-1,1:-1,1:-1])+np.multiply(self.aw[1:,1:,1:],self.phi[t,:-2,1:-1,1:-1])+np.multiply(self.ae[1,:1:,1:],self.phi[t,2:,1:-1,1:-1])
                error+=np.multiply(self.ab[1:,1:,1:],self.phi[t,1:-1,1:-1,:-2])+np.multiply(self.af[1:,1:,1:],self.phi[t,1:-1,1:-1,2:])
                error+=np.multiply(self.ad[1:,1:,1:],self.phi[t,1:-1,:-2,1:-1])+np.multiply(self.au[1:,1:,1:],self.phi[t,1:-1,2:,1:-1])
                error=self.phi[t,1:-1,1:-1,1:-1]-np.divide(error,self.ap[1:,1:,1:])
                # print(self.phi[t,1,1:-1,1:-1])
                err=np.max(np.absolute(error))
                if err>1e-6:
                    print("At itr=",itr," Error:",err)
                    itr+=1
                    if(np.absolute(err-error_init)<1e-12):
                        break
                    error_init=err
                    if itr>100:
                        break
                else:
                    print("for time t:",t," Itr:",itr)
                    break
        return self.phi
