import numpy as np
from scipy import sparse




def compute_dos(cheb_moments, A,B,sigma, num_points=5000):
    def jackson_kernel(M):
        """ Computes the Jackson kernel coefficients for given number of moments M. """
        k = np.arange(M)
        ans= ((M - k) * np.cos(np.pi * k / M) + np.sin(np.pi * k / M) / np.tan(np.pi / M)) / M
        ans[0]*=0.5
        return ans
    """ Computes the density of states using Chebyshev moments and the Jackson kernel. """
    M = len(cheb_moments)
    kernel = jackson_kernel(M)
    x = np.linspace(-0.99, 0.99, num_points)  # energy values where DOS is computed
    T_k = np.cos(np.arccos(x)[:, None] * np.arange(M))  # Chebyshev polynomials T_k(x)
    dos = T_k @ (cheb_moments * kernel)  # DOS computation using matrix multiplication
    dos*=2/np.pi/np.sqrt(1-x**2)
    
    dos/=B
    x*=B
    x+=A
    dx=(x[1]-x[0])
    npoints=int(round(5*sigma/dx))
    ker_x=np.arange(-npoints,npoints)*dx 
    ker_y=1/(2*np.pi*sigma**2)**0.5*np.exp(-ker_x**2/2/sigma**2)
    dos=np.convolve(dos,ker_y,"same")*dx
    return x, dos

def compute_dos_exact(hiidx,hjidx,hvalue,num_points=5000,sigma=0.05):
    H=sparse.coo_matrix((hvalue,(hiidx,hjidx))).todense()
    e,_=np.linalg.eigh(H)
    e=np.sort(e)
    emin=e[0]-5*sigma
    emax=e[-1]+5*sigma 
    Elist=np.linspace(emin,emax,num_points)
    DE=np.zeros_like(Elist)
    dE=Elist[1]-Elist[0]
    for i in range(e.shape[0]):
        idx_aux=int(round((e[i]-emin)/dE))
        DE[idx_aux]+=1/dE 
    npoints=int(round(5*sigma/dE))
    ker_x=np.arange(-npoints,npoints)*dE
    ker_y=1/(2*np.pi*sigma**2)**0.5*np.exp(-ker_x**2/2/sigma**2)
    DE=np.convolve(DE,ker_y,"same")*dE


    return Elist,DE

    
def compute_eh_list(mumn, Nx, Ny, hw, A, B, Fermi, beta, gph):
    #All units should be in Hartree right? 

    xlist = np.linspace(-0.99, 0.99, Nx)
    ylist = np.linspace(-0.99, 0.99, Ny)
    NTx, NTy = mumn.shape
    NTxlist = np.arange(NTx)
    NTylist = np.arange(NTy)
    Jx = 1.0 / NTx * ((NTx - NTxlist) * np.cos(np.pi * NTxlist / NTx) + np.sin(np.pi * NTxlist / NTx) * np.cos(np.pi / NTx) / np.sin(np.pi / NTx))
    Jy = 1.0 / NTy * ((NTy - NTylist) * np.cos(np.pi * NTylist / NTy) + np.sin(np.pi * NTylist / NTy) * np.cos(np.pi / NTy) / np.sin(np.pi / NTy))
    Jx[0] *= 0.5
    Jy[0] *= 0.5
    Jxx, Jyy = np.meshgrid(Jx, Jy,indexing="ij")
    mutmn = mumn * Jxx * Jyy 
    NTyy, yy = np.meshgrid(NTylist, ylist,indexing="ij")
    Dny = np.cos(NTyy * np.arccos(yy)) / np.sqrt(1.0 - yy**2)
    xx, NTxx = np.meshgrid(xlist, NTxlist,indexing="ij")
    Dxm = np.cos(NTxx * np.arccos(xx)) / np.sqrt(1.0 - xx**2)
    phi = Dxm @ mutmn @ Dny 
    
    #factor = np.pi**2 / B**2 / 4
    phi *= 4/B**2/np.pi**2

    Ei = ((xlist * B) + A)
    Ej = ((ylist * B) + A)
    factor=np.zeros((Nx,Ny))
    gij=gph*2
    Eii,Ejj=np.meshgrid(Ei,Ej,indexing="ij")
    factor=2*np.pi*2.0/np.sqrt(2*np.pi*gij**2)*np.exp(-(hw-(Ejj-Eii) )**2/2.0/gij**2)*1.0/(1.0+np.exp( beta*(Eii-Fermi) ) )*(1.0-(1.0/( 1.0+np.exp(beta*(Ejj-Fermi )))))
    hole_factor=2*np.pi*2.0/np.sqrt(2.0*np.pi*gij**2)*np.exp(-(hw-(Eii-Ejj))**2 / 2.0 / gij**2) * 1.0 / (1.0 + np.exp(beta * (Ejj - Fermi))) * (1.0 - (1.0 / (1.0 + np.exp(beta * (Eii - Fermi)))))

    factor=factor*phi
    hole_factor=hole_factor*phi
    dx = Ei[1] - Ei[0]
    Nelist=np.sum(factor,axis=0)*dx
    Nhlist=np.sum(hole_factor,axis=0)*dx

    return Ej, Nelist, Nhlist

def compute_exact_eh_dis(hiidx,hjidx,havalue, phiidx,phijdix,phivalue, Nx, Ny, A, B, gph, Fermi, beta, hw):
    """
    This function computes the electron hole distribution based on a given sparse matrix H and Phi,
    it converts the matrix to real, and then performs the computation via diagonalization.
    """
    pi=np.pi
    HD=sparse.coo_matrix((havalue,(hiidx,hjidx))).todense()
    Phi=sparse.coo_matrix((phivalue,(phiidx,phijdix))).todense()
    E, C = np.linalg.eigh(HD)
    np.save("E_py",E)
    np.save("C_py",C)
    #E = np.real(E)
    Cdagger=np.conjugate(C).T
    np.save("Cdagger_py",Cdagger)
    Intermediate=Phi.dot(C)
    np.save("Intermediate_py",Intermediate)
    Melement = Cdagger.dot(Intermediate)
    np.save("Melement_nosquare_py",Melement)
    Intermediate=Melement.copy()
    Melement = Intermediate.real**2+Intermediate.imag**2
    np.save("Melement_py",Melement)
    print("Melement done!")
    
    xlist = np.linspace(-0.99, 0.99, Nx)
    ylist = np.linspace(-0.99, 0.99, Ny)
    Ei = xlist * B + A
    Ef = ylist * B + A
    
    phi = np.zeros((Nx, Ny))
    DEi = Ei[1] - Ei[0]
    DEf = Ef[1] - Ef[0]
    
    for ii in range(len(E)):
        Eii_aux = E[ii]
        Ei_idx = int(np.floor((Eii_aux - Ei[0]) / DEi))
        for ff in range(len(E)):
            Eff_aux = E[ff]
            Ef_idx = int(np.floor((Eff_aux - Ef[0]) / DEf))
            phi[Ei_idx, Ef_idx] += Melement[ii, ff]
    np.save("phi_py",phi)
    
    factor = np.zeros((Nx, Ny))
    
    for i in range(Nx):
        for j in range(Ny):
            gif = gph * 2
            factor[i, j] = (2 * pi * 2.0 / np.sqrt(2.0 * pi * gif**2) *
                            np.exp(-(hw - (Ef[j] - Ei[i]))**2 / (2.0 * gif**2)) *
                            1.0 / (1.0 + np.exp(beta * (Ei[i] - Fermi))) *
                            (1.0 - (1.0 / (1.0 + np.exp(beta * (Ef[j] - Fermi))))))
    
    hole_factor = np.zeros((Nx, Ny))
    
    for i in range(Nx):
        for j in range(Ny):
            gif = gph * 2
            hole_factor[i, j] = (2 * pi * 2.0 / np.sqrt(2.0 * pi * gif**2) *
                                 np.exp(-(hw - (Ei[i] - Ef[j]))**2 / (2.0 * gif**2)) *
                                 1.0 / (1.0 + np.exp(beta * (Ef[j] - Fermi))) *
                                 (1.0 - (1.0 / (1.0 + np.exp(beta * (Ei[i] - Fermi))))))
    
    Nelist = np.zeros(Ny)
    Nhlist = np.zeros(Ny)
    
    dx = Ei[1] - Ei[0]
    dy = Ef[1] - Ef[0]
    
    for i in range(Ny):
        for j in range(Nx):
            Nelist[i] += phi[j, i] * factor[j, i]
            Nhlist[i] += phi[j, i] * hole_factor[j, i]
    
    Nelist /= dy
    Nhlist /= dy
    
    return Ef, Nelist, Nhlist

def compute_eh_list_test_version(mumn, Nx, Ny, hw, A, B, Fermi, beta, gph):
    #All units should be in Hartree right? 
    """
    This function computes the electron hole distribution but leaves all the npy and npz file for comparison. 
    """

    xlist = np.linspace(-0.99, 0.99, Nx)
    ylist = np.linspace(-0.99, 0.99, Ny)
    NTx, NTy = mumn.shape
    NTxlist = np.arange(NTx)
    NTylist = np.arange(NTy)
    Jx = 1.0 / NTx * ((NTx - NTxlist) * np.cos(np.pi * NTxlist / NTx) + np.sin(np.pi * NTxlist / NTx) * np.cos(np.pi / NTx) / np.sin(np.pi / NTx))
    Jy = 1.0 / NTy * ((NTy - NTylist) * np.cos(np.pi * NTylist / NTy) + np.sin(np.pi * NTylist / NTy) * np.cos(np.pi / NTy) / np.sin(np.pi / NTy))
    Jx[0] *= 0.5
    Jy[0] *= 0.5
    Jxx, Jyy = np.meshgrid(Jx, Jy,indexing="ij")
    mutmn = mumn * Jxx * Jyy 
    NTyy, yy = np.meshgrid(NTylist, ylist,indexing="ij")
    Dny = np.cos(NTyy * np.arccos(yy)) / np.sqrt(1.0 - yy**2)
    xx, NTxx = np.meshgrid(xlist, NTxlist,indexing="ij")
    Dxm = np.cos(NTxx * np.arccos(xx)) / np.sqrt(1.0 - xx**2)
    phi = Dxm @ mutmn @ Dny 
    
    #factor = np.pi**2 / B**2 / 4
    phi *= np.pi**2 / B**2 / 4
    np.save("phi_python.npy",phi)

    Ei = ((xlist * B) + A)
    Ej = ((ylist * B) + A)
    factor=np.zeros((Nx,Ny))
    gij=gph*2
    Eii,Ejj=np.meshgrid(Ei,Ej,indexing="ij")
    factor=2*np.pi*2.0/np.sqrt(2*np.pi*gij**2)*np.exp(-(hw-(Ejj-Eii) )**2/2.0/gij**2)*1.0/(1.0+np.exp( beta*(Eii-Fermi) ) )*(1.0-(1.0/( 1.0+np.exp(beta*(Ejj-Fermi )))))
    np.save("factor_python.npy",factor)
    hole_factor=2*np.pi*2.0/np.sqrt(2.0*np.pi*gij**2)*np.exp(-(hw-(Eii-Ejj))**2 / 2.0 / gij**2) * 1.0 / (1.0 + np.exp(beta * (Ejj - Fermi))) * (1.0 - (1.0 / (1.0 + np.exp(beta * (Eii - Fermi)))))
    np.save("hole_factor_python.npy",hole_factor)


    

    
    factor=factor*phi
    hole_factor=hole_factor*phi
    dx = Ei[1] - Ei[0]
    Nelist=np.sum(factor,axis=0)*dx
    Nhlist=np.sum(hole_factor,axis=0)*dx
    
    np.save("Nelist.npy",Nelist)
    np.save("Nhlist.npy",Nhlist)
    return Ej, Nelist, Nhlist



if __name__=="__main__":
    import numpy as np
    from shapelib import Sphere
    """
    np.random.seed(114514)
    epsr=-3.2096-1.86j
    epsm=1
    Ball=Sphere(10.5*4,epsw=epsr,epsm=epsm)
    from ase.build import bulk
    Ag=bulk("Ag","fcc",4.13)
    Ag_SK=np.load(r"Ag_SK15.npz")
    Ag_SK=Ag_SK["SK"]
    onsite_args=Ag_SK[20:].tolist()
    hopping_args=Ag_SK[:20].tolist()
    hopping_args.append(4.13)
    
    ase_nanoparticle=construct_ase_nanoparticle(Ag,Ball,extend=(120,120,120))
    iidx,jidx,value=construct_Hamiltonian(ase_nanoparticle,{"Ag":4.13/2**0.5/2+0.01},spd_hopping,onsite_model,hopping_args,onsite_args,A=13.38281572788353,B=11.886513341034474+6)
    Hdim=np.max(iidx)+1
    v=np.random.rand(Hdim)-0.5 
    v=v/np.sqrt(np.linalg.norm(v)**2/Hdim)*(1+0j)
    
    
    #A=13.38281572788353
    #B=11.886513341034474(+ 1)
    phiiidx,phijidx,phivalue=Ball.compute_phi(ase_nanoparticle.get_positions())
    np.savez("phi.npz",iidx=phiiidx,jidx=phijidx,value=phivalue)
    np.savez("H.npz",iidx=iidx,jidx=jidx,value=value)
    np.save("rvector.npy",v)
    from KPM import compute_mumn
    NNR=500
    from time import time
    t0=time()
    #print("Now let's compute mumn")
    mumn=compute_mumn(iidx,jidx,value,phiiidx,phijidx,phivalue,NNR,25,seed_number=114514)
    print(mumn)
    #mumn=np.reshape(mumn,(NNR,NNR))
    np.save("mumn.npy",mumn)"""
    
    

