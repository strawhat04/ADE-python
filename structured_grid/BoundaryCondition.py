import numpy as np

'''
d>dirichlet
n>newman
a>advective

''' 

def Bcond(phi, b, Sx, Sy, Sz,t, nx, ny, nz):
    bc=[]

    bc.append({'type':'d', 'val':7})	#x0 face
    bc.append({'type':'a', 'val':7})	#y0 face
    bc.append({'type':'d', 'val':7})	#z0 face
    bc.append({'type':'a', 'val':7})	#xn face
    bc.append({'type':'a', 'val':7})	#yn face
    bc.append({'type':'a', 'val':7})	#zn face
    #Computing flux coefficient for dirichlet condition

    #for x0 face boundary
    if bc[0]['type']=='d':
    	phi[t,:,:,1]=bc[0]['val']
    if bc[0]['type']=='a':
    	pass
    if bc[0]['type']=='n':
    	aw[:,:,0]=0
    	b[:,:,0]+=bc[0]['val']*Sx

    #for y0 face boundary
    if bc[1]['type']=='d':
    	phi[t,:,1,:]=bc[1]['val']
    if bc[1]['type']=='a':
    	pass
    if bc[1]['type']=='n':
    	ad[:,1,:]=0
    	b[:,1,:]+=bc[1]['val']*Sy

    #for z0 face bcoundary
    if bc[2]['type']=='d':
    	phi[t,1,:,:]=bc[2]['val']
    if bc[2]['type']=='a':
    	pass
    if bc[2]['type']=='n':
    	abc[2,:,:]=0
    	b[2,:,:]+=bc[2]['val']*Sz

    #for xn face boundary
    if bc[3]['type']=='d':
    	phi[t,:,:,nx]=bc[3]['val']
    if bc[3]['type']=='a':
    	pass
    if bc[3]['type']=='n':
    	aw[:,:,nx]=0
    	b[:,:,nx]+=bc[3]['val']*Sx

    #for yn face boundary
    if bc[4]['type']=='d':
    	phi[t,:,ny,:]=bc[4]['val']
    if bc[4]['type']=='a':
    	pass
    if bc[4]['type']=='n':
    	aw[:,ny,:]=0
    	b[:,ny,:]+=bc[4]['val']*Sy

    #for zn face boundary
    if bc[5]['type']=='d':
    	phi[t,nz,:,:]=bc[5]['val']
    if bc[5]['type']=='a':
    	pass
    if bc[5]['type']=='n':
    	aw[nz,:,:]=0
    	b[nz,:,:]+=bc[5]['val']*Sz

    return phi, b	