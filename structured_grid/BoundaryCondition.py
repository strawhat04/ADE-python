import numpy as np
'''
d>dirichlet
f>neumann flux
a>advective

''' 

def Bcond(phi, Sx, Sy, Sz,t, nx, ny, nz, ap, aw, ae, ad, au, ab, af, ap0, b, Q,R):
    # R=np.zeros((nz,ny,nx+2),)
    # Q=np.zeros((nz,ny,nx+2),)
    
    #these Bc need to implement for TDMA specifically
    Q[:,:,0]=R[:,:,0]=0

    bc=[]

    bc.append({'type':'d', 'val':20})   #x0 face
    bc.append({'type':'d', 'val':15})   #y0 face
    bc.append({'type':'d', 'val':10})   #z0 face
    bc.append({'type':'d', 'val':5})    #xn face
    bc.append({'type':'d', 'val':20})   #yn face
    bc.append({'type':'d', 'val':0})    #zn face
    #Computing flux coefficient for dirichlet condition

    #for x0 face boundary
    if bc[0]['type']=='d':
        phi[t,:,:,1]=bc[0]['val']
        Q[:,:,1]=phi[t,1:-1,1:-1,1]
        R[:,:,1]=0
    if bc[0]['type']=='a':
        pass
    if bc[0]['type']=='n':
        aw[:,:,0]=0
        b[:,:,0]=bc[0]['val']*Sx


    #for y0 face boundary
    if bc[1]['type']=='d':
        phi[t,:,1,:]=bc[1]['val']
    if bc[1]['type']=='a':
        pass
    if bc[1]['type']=='f':
        ad[:,1,:]=0
        b[:,0,:]=bc[1]['val']*Sy

    #for z0 face bcoundary
    if bc[2]['type']=='d':
        phi[t,1,:,:]=bc[2]['val']
    if bc[2]['type']=='a':
        pass
    if bc[2]['type']=='f':
        ab[2,:,:]=0
        b[0,:,:]=bc[2]['val']*Sz

    #for xn face boundary
    if bc[3]['type']=='d':
        phi[t,:,:,nx]=bc[3]['val']
        Q[:,:,nx]=phi[t,1:-1,1:-1,nx]
        R[:,:,nx]=0

    if bc[3]['type']=='a':
        pass
    if bc[3]['type']=='f':
        ae[:,:,nx-1]=0
        b[:,:,nx-1]=bc[3]['val']*Sx
        # Q[:,:,nx]=b[:,:,nx-1]/ap[:,:,nx-1]

    #for yn face boundary
    if bc[4]['type']=='d':
        phi[t,:,ny,:]=bc[4]['val']
    if bc[4]['type']=='a':
        pass
    if bc[4]['type']=='f':
        au[:,ny-1,:]=0
        b[:,ny-1,:]=bc[4]['val']*Sy

    #for zn face boundary
    if bc[5]['type']=='d':
        phi[t,nz,:,:]=bc[5]['val']
    if bc[5]['type']=='a':
        pass
    if bc[5]['type']=='f':
        af[nz-1,:,:]=0
        b[nz-1,:,:]=bc[5]['val']*Sz

    return phi, b, Q, R