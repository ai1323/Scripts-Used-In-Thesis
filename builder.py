from ase import Atoms
from kpm.shapelib import Shape
import numpy as np 
from scipy.spatial import KDTree
from typing import Callable
from kpm.evaltool import eval_ase
from kpm.materialslib import Materials
from ase import io


def construct_ase_nanoparticle(shape:Shape,materials:Materials,offset=[[0.0,0.0,0.0]]):
    #Need an undate of unit cell method. 
    #Create the bulk cell for slicing 
    for ii in range(shape.nregion):
        offset[ii]=np.array(offset[ii])


    natoms_converged=False
    extend=[10,10,10]

    for i in range(shape.nregion):
        
        while not natoms_converged:
            super_cell=materials.cells[i].repeat(extend)
            #print(len(super_cell))
            offset2=np.array([extend[k]//2 for k in range(3)]).dot(materials.cells[i].get_cell())
            super_cell.set_cell(None)
            super_cell.translate(-offset2+offset[i])
            inshape_atoms=shape.inshape(super_cell.get_positions(),i)
            inshape_lists=[i for i in range(len(super_cell)) if inshape_atoms[i]]
            result=super_cell[inshape_lists]
            natoms=sum(inshape_atoms)

            #extend_test=np.array([extend[i]+4 for i in range(3)]).dot(shape.materials[i].get_cell())
            extend_test=[extend[k]+8 for k in range(3)]
            super_cell_test=materials.cells[i].repeat(extend_test)

            offset2=np.array([extend_test[k]//2 for k in range(3)]).dot(materials.cells[i].get_cell())
            super_cell_test.set_cell(None)
            super_cell_test.translate(-offset2+offset[i])
            inshape_atoms_test=shape.inshape(super_cell_test.get_positions(),i)
            natoms_test=sum(inshape_atoms_test)

            natoms_converged=(natoms==natoms_test)
            if not natoms_converged:
                extend=[extend[i]*2 for i in range(3)]


        
        if i==0:
            ase_nanoparticle=result
        else:
            ase_nanoparticle.extend(result)
        natoms_converged=False
        extend=[10,10,10]



    return ase_nanoparticle

class nanoparticle_nl:
    def __init__(self, positions: np.ndarray, cutoffs: list, skin=0.01):
        self.max_cutoff = max(cutoffs) * 2+skin
        self.positions = positions
        self.Tree = KDTree(positions)
        self.cutoffs = cutoffs
        self.skin = skin

    def get_neighbors(self, idx):
        Ri = self.positions[idx, :]
        neighbor_list = []
        neighbors_list_aux = self.Tree.query_ball_point(Ri, self.max_cutoff)
        
        for j in neighbors_list_aux:
            if np.linalg.norm(self.positions[j, :] - Ri) < self.cutoffs[idx] + self.cutoffs[j] + self.skin:
                if j != idx:
                    neighbor_list.append(j)
        return neighbor_list



def spd_hopping(R:np.ndarray,
             ss_sig1:float,pp_sig1:float,pp_pi1:float,dd_sig1:float,dd_pi1:float,dd_del1:float,sp_sig1:float,sd_sig1:float,pd_sig1:float,pd_pi1:float,
             ss_sig2:float,pp_sig2:float,pp_pi2:float,dd_sig2:float,dd_pi2:float,dd_del2:float,sp_sig2:float,sd_sig2:float,pd_sig2:float,pd_pi2:float,
             Es:float,Ep:float,Ed:float,
             R0:float,nn1_skin=0.01):
    """
    Comment: klist[Nk,3], where Nk is the number of k-points, the k-vector provided should be in unit of 2pi/Angstrom
    R: 3 dimensional array.
    ss_sig to a_pd_pi2: fitting parameters
    R0:The unstrained distance.
    R: the hopping position
    """
    def getHoppingParameters(R):
        r=np.linalg.norm(R)
        if r==0:
            raise Exception("Error, r=0")
        if r<1/2**0.5*R0+nn1_skin:
            ss_sig=ss_sig1
            pp_sig=pp_sig1
            pp_pi=pp_pi1
            dd_sig=dd_sig1
            dd_pi=dd_pi1
            dd_del=dd_del1
            sp_sig=sp_sig1
            sd_sig=sd_sig1
            pd_sig=pd_sig1
            pd_pi=pd_pi1
        elif r>3*R0:
            raise Exception("Error, the radius is too big, please check the input. ")
        else:
            ss_sig=ss_sig2
            pp_sig=pp_sig2
            pp_pi=pp_pi2
            dd_sig=dd_sig2
            dd_pi=dd_pi2
            dd_del=dd_del2
            sp_sig=sp_sig2
            sd_sig=sd_sig2
            pd_sig=pd_sig2
            pd_pi=pd_pi2
            pd_pi=pd_pi2
        return sp_sig,ss_sig,pp_sig,pp_pi,sd_sig,pd_sig,pd_pi,dd_sig,dd_pi,dd_del
    VR=np.zeros([9,9],dtype=np.complex128)



    r=np.linalg.norm(R)
    l,m,n=R/r
    sp_sig,ss_sig,pp_sig,pp_pi,sd_sig,pd_sig,pd_pi,dd_sig,dd_pi,dd_del=getHoppingParameters(R)
    #s and d interaction
    #Same orbital interaction
    VR[0,0]+=ss_sig
    VR[1,1]+=( 3*l**2*m**2*dd_sig + ( l**2 + m**2 - 4*l**2*m**2)*dd_pi +(n**2+l**2*m**2)*dd_del)
    VR[2,2]+=( 3*l**2*n**2*dd_sig + ( l**2 + n**2 - 4*l**2*n**2)*dd_pi +(m**2+l**2*n**2)*dd_del)
    VR[3,3]+=( 3*m**2*n**2*dd_sig + ( m**2 + n**2 - 4*m**2*n**2)*dd_pi +(l**2+m**2*n**2)*dd_del)
    VR[4,4]+=(3/4*(l**2-m**2)**2*dd_sig+(l**2+m**2-(l**2-m**2)**2)*dd_pi+(n**2+(l**2-m**2)**2/4)*dd_del)
    VR[5,5]+=((n**2-.5*(l**2+m**2))**2*dd_sig+3*n**2*(l**2+m**2)*dd_pi +3/4*(l**2+m**2)**2*dd_del)
    VR[6,6]+=(l**2*pp_sig+(1-l**2)*pp_pi)
    VR[7,7]+=(m**2*pp_sig+(1-m**2)*pp_pi)
    VR[8,8]+=(n**2*pp_sig+(1-n**2)*pp_pi)
    #p-p interaction
    VR[6,7]+=(l*m*pp_sig-l*m*pp_pi)
    VR[6,8]+=(l*n*pp_sig-l*n*pp_pi)
    VR[7,8]+=(m*n*pp_sig-m*n*pp_pi)
    #dd interaction
    VR[1,2]+=(3*l**2*m*n*dd_sig+m*n*(1-4*l**2)*dd_pi+m*n*(l**2-1)*dd_del) #xy-xz
    VR[1,3]+=(3*l*m**2*n*dd_sig+l*n*(1-4*m**2)*dd_pi+l*n*(m**2-1)*dd_del) #xy-yz
    VR[1,4]+=(1.5*l*m*(l**2-m**2)*dd_sig+2*l*m*(m**2-l**2)*dd_pi+.5*l*m*(l**2-m**2)*dd_del) #xy - x2-y2
    VR[1,5]+=(3**.5*l*m*(n**2-.5*(l**2+m**2))*dd_sig-2*3**.5*l*m*n**2*dd_pi+.5*3**.5*l*m*(1+n**2)*dd_del) #xy - 3z**2-r**2
    VR[2,3]+=(3*n**2*m*l*dd_sig+m*l*(1-4*n**2)*dd_pi+l*m*(n**2-1)*dd_del) #xz-yz
    VR[2,4]+=(1.5*n*l*(l**2-m**2)*dd_sig+n*l*(1-2*(l**2-m**2))*dd_pi-n*l*(1-.5*(l**2-m**2))*dd_del) #xz ->x**2-y**2
    VR[2,5]+=(3**.5*l*n*(n**2-.5*(l**2+m**2))*dd_sig+3**.5*l*n*(l**2+m**2-n**2)*dd_pi-.5*3**.5*l*n*(l**2+m**2)*dd_del) #xz ->z**2
    VR[3,4]+=(1.5*m*n*(l**2-m**2)*dd_sig-m*n*(1+2*(l**2-m**2))*dd_pi+m*n*(1+(l**2-m**2)/2)*dd_del) #yz ->x**2-y**2
    VR[3,5]+=(3**.5*m*n*(n**2-.5*(l**2+m**2))*dd_sig+3**.5*m*n*(l**2+m**2-n**2)*dd_pi-.5*3**.5*m*n*(l**2+m**2)*dd_del)#yz->z**2
    VR[4,5]+=(.5*3**.5*(l**2-m**2)*(n**2-.5*(l**2+m**2))*dd_sig+3**.5*n**2*(m**2-l**2)*dd_pi+3**.5*(1+n**2)*(l**2-m**2)/4*dd_del)
    #s-d interaction
    VR[0,1]+=3**.5 * l*m*sd_sig
    VR[0,2]+=3**.5*l*n*sd_sig
    VR[0,3]+=3**.5*n*m*sd_sig
    VR[0,4]+=3**.5/2*(l**2-m**2)*sd_sig
    VR[0,5]+=(n**2-.5*(l**2+m**2))*sd_sig
    #sp interaction
    VR[0,6]+=(l*sp_sig)
    VR[0,7]+=(m*sp_sig)
    VR[0,8]+=(n*sp_sig)
    #pd interaction
    VR[6,1]+=(3**.5*l**2*m*pd_sig+m*(1-2*l**2)*pd_pi)#x->xy
    VR[7,1]+=(3**0.5*m**2*l*pd_sig+l*(1-2*m**2)*pd_pi)#y->xy
    VR[8,1]+=(3**0.5*l*m*n*pd_sig-2*l*m*n*pd_pi)#z->xy
    VR[6,2]+=(3**.5*l**2*n*pd_sig+n*(1-2*l**2)*pd_pi) #x->xz
    VR[7,2]+=(3**0.5*l*m*n*pd_sig-2*l*m*n*pd_pi) #xz->y
    VR[8,2]+=(3**0.5*n**2*l*pd_sig+l*(1-2*n**2)*pd_pi) #xz->z
    VR[6,3]+=(3**0.5*l*m*n*pd_sig-2*l*m*n*pd_pi) #yz->x
    VR[7,3]+=(3**0.5*m**2*n*pd_sig+n*(1-2*m**2)*pd_pi) #yz->y
    VR[8,3]+=(3**0.5*n**2*m*pd_sig+m*(1-2*n**2)*pd_pi) #yz->z
    VR[6,4]+=(3**0.5/2*l*(l**2-m**2)*pd_sig+l*(1-l**2+m**2)*pd_pi) #x**2-y**2->x
    VR[7,4]+=(3**0.5/2*m*(l**2-m**2)*pd_sig-m*(1+l**2-m**2)*pd_pi) #x**2-y**2->y
    VR[8,4]+=(3**0.5/2*n*(l**2-m**2)*pd_sig-n*(l**2-m**2)*pd_pi) #x**2-y**2->z
    VR[6,5]+=(l*(n**2-(l**2+m**2)/2)*pd_sig-3**0.5*l*n**2*pd_pi) #3z**2-r**2->x
    VR[7,5]+=(m*(n**2-(l**2+m**2)/2)*pd_sig-3**0.5*m*n**2*pd_pi) #3z**2-r**2->y
    VR[8,5]+=(n*(n**2-(l**2+m**2)/2)*pd_sig+3**0.5*n*(l**2+m**2)*pd_pi)#3z**2-r**2->z
    #For those hopping not listed in SK table
    #Compute the opposite hopping and use Hermiticity of Hamiltonian to
    #derive these elements.
    #Now flip the direction:
    l,m,n=-R/r
    #dd interaction
    VR[2,1]+=(3*l**2*m*n*dd_sig+m*n*(1-4*l**2)*dd_pi+m*n*(l**2-1)*dd_del) #xy-xz
    VR[3,1]+=(3*l*m**2*n*dd_sig+l*n*(1-4*m**2)*dd_pi+l*n*(m**2-1)*dd_del) #xy-yz
    VR[4,1]+=(1.5*l*m*(l**2-m**2)*dd_sig+2*l*m*(m**2-l**2)*dd_pi+.5*l*m*(l**2-m**2)*dd_del) #xy - x2-y2
    VR[5,1]+=(3**.5*l*m*(n**2-.5*(l**2+m**2))*dd_sig-2*3**.5*l*m*n**2*dd_pi+.5*3**.5*l*m*(1+n**2)*dd_del) #xy - 3z**2-r**2
    VR[3,2]+=(3*n**2*m*l*dd_sig+m*l*(1-4*n**2)*dd_pi+l*m*(n**2-1)*dd_del) #xz-yz
    VR[4,2]+=(1.5*n*l*(l**2-m**2)*dd_sig+n*l*(1-2*(l**2-m**2))*dd_pi-n*l*(1-.5*(l**2-m**2))*dd_del) #xz ->x**2-y**2
    VR[5,2]+=(3**.5*l*n*(n**2-.5*(l**2+m**2))*dd_sig+3**.5*l*n*(l**2+m**2-n**2)*dd_pi-.5*3**.5*l*n*(l**2+m**2)*dd_del) #xz ->z**2
    VR[4,3]+=(1.5*m*n*(l**2-m**2)*dd_sig-m*n*(1+2*(l**2-m**2))*dd_pi+m*n*(1+(l**2-m**2)/2)*dd_del) #yz ->x**2-y**2
    VR[5,3]+=(3**.5*m*n*(n**2-.5*(l**2+m**2))*dd_sig+3**.5*m*n*(l**2+m**2-n**2)*dd_pi-.5*3**.5*m*n*(l**2+m**2)*dd_del)#yz->z**2
    VR[5,4]+=(.5*3**.5*(l**2-m**2)*(n**2-.5*(l**2+m**2))*dd_sig+3**.5*n**2*(m**2-l**2)*dd_pi+3**.5*(1+n**2)*(l**2-m**2)/4*dd_del)
    #p-p interaction
    VR[7,6]+=(l*m*pp_sig-l*m*pp_pi)
    VR[8,6]+=(l*n*pp_sig-l*n*pp_pi)
    VR[8,7]+=(m*n*pp_sig-m*n*pp_pi)
    #ds interaction
    VR[1,0]+=3**.5 * l*m*sd_sig
    VR[2,0]+=3**.5*l*n*sd_sig
    VR[3,0]+=3**.5*n*m*sd_sig
    VR[4,0]+=3**.5/2*(l**2-m**2)*sd_sig
    VR[5,0]+=(n**2-.5*(l**2+m**2))*sd_sig
    #ps interaction
    VR[6,0]+=(l*sp_sig)
    VR[7,0]+=(m*sp_sig)
    VR[8,0]+=(n*sp_sig)
    #dp interaction
    VR[1,6]+=(3**.5*l**2*m*pd_sig+m*(1-2*l**2)*pd_pi)#x->xy
    VR[1,7]+=(3**0.5*m**2*l*pd_sig+l*(1-2*m**2)*pd_pi)#y->xy
    VR[1,8]+=(3**0.5*l*m*n*pd_sig-2*l*m*n*pd_pi)#z->xy
    VR[2,6]+=(3**.5*l**2*n*pd_sig+n*(1-2*l**2)*pd_pi) #x->xz
    VR[2,7]+=(3**0.5*l*m*n*pd_sig-2*l*m*n*pd_pi) #xz->y
    VR[2,8]+=(3**0.5*n**2*l*pd_sig+l*(1-2*n**2)*pd_pi) #xz->z
    VR[3,6]+=(3**0.5*l*m*n*pd_sig-2*l*m*n*pd_pi) #yz->x
    VR[3,7]+=(3**0.5*m**2*n*pd_sig+n*(1-2*m**2)*pd_pi) #yz->y
    VR[3,8]+=(3**0.5*n**2*m*pd_sig+m*(1-2*n**2)*pd_pi) #yz->z
    VR[4,6]+=(3**0.5/2*l*(l**2-m**2)*pd_sig+l*(1-l**2+m**2)*pd_pi) #x**2-y**2->x
    VR[4,7]+=(3**0.5/2*m*(l**2-m**2)*pd_sig-m*(1+l**2-m**2)*pd_pi) #x**2-y**2->y
    VR[4,8]+=(3**0.5/2*n*(l**2-m**2)*pd_sig-n*(l**2-m**2)*pd_pi) #x**2-y**2->z
    VR[5,6]+=(l*(n**2-(l**2+m**2)/2)*pd_sig-3**0.5*l*n**2*pd_pi) #3z**2-r**2->x
    VR[5,7]+=(m*(n**2-(l**2+m**2)/2)*pd_sig-3**0.5*m*n**2*pd_pi) #3z**2-r**2->y
    VR[5,8]+=(n*(n**2-(l**2+m**2)/2)*pd_sig+3**0.5*n*(l**2+m**2)*pd_pi)#3z**2-r**2->z
    return VR

def spd_onsite(ss_sig1:float,pp_sig1:float,pp_pi1:float,dd_sig1:float,dd_pi1:float,dd_del1:float,sp_sig1:float,sd_sig1:float,pd_sig1:float,pd_pi1:float,
        ss_sig2:float,pp_sig2:float,pp_pi2:float,dd_sig2:float,dd_pi2:float,dd_del2:float,sp_sig2:float,sd_sig2:float,pd_sig2:float,pd_pi2:float,
        Es:float,Ep:float,Ed:float,
        R0:float):
    return np.diag([Es,Ed,Ed,Ed,Ed,Ed,Ep,Ep,Ep])


#I think this function requires python 3.8 or above. 
def construct_Hamiltonian(ase_nanoparticle:Atoms,materials:Materials,hopping_model=spd_hopping,onsite_model=spd_onsite,
                          onsite_kwargs={},hopping_kwargs={},nl_skin=0.01,A=0,B=1,testing=False):
    """
    SK_dict[str]: containing all args for sdp_hopping. 

    """
    def addHR(iidx,jidx,value,HR,iptr,jptr):
        Ni,Nj=HR.shape
        for i in range(Ni):
            for j in range(Nj):
                Hij=HR[i,j]
                if np.abs(Hij)>1e-10:
                    iidx.append(iptr+i)
                    jidx.append(jptr+j)
                    value.append(Hij)
    cutoffs=[]
    if testing:
        f=open("neighbor_test.txt","w")
    for i in range(len(ase_nanoparticle)):
        try:
            cutoff=materials.SK_dict[ase_nanoparticle[i].symbol][23]/2
        except:
            materials.SK_dict[ase_nanoparticle[i].symbol]=np.loadtxt(ase_nanoparticle[i].symbol+"_SK.txt")
            cutoff=materials.SK_dict[ase_nanoparticle[i].symbol][23]/2
        cutoffs.append(cutoff)


    #Construct the neighborlist 
    nl=nanoparticle_nl(ase_nanoparticle.get_positions(),cutoffs,skin=nl_skin)
    
    
    #ase_nanoparticle.set_positions(ase_nanoparticle.get_positions()+np.array([big_num,big_num,big_num]))
    #Need to set cell for nl to work
    #ase_nanoparticle.set_cell(np.array([[big_num,0,0],[0,big_num,0],[0,0,big_num]]))
    iidx=[]
    jidx=[]
    value=[]

    #Before running the actual test, let's test the validity of nl.get_neighbors():

    
    for i in range(len(ase_nanoparticle)):
        
        try:
            Atom_i_SK=np.array(materials.SK_dict[ase_nanoparticle[i].symbol])
        except:
            materials.SK_dict[ase_nanoparticle[i].symbol]=np.load(ase_nanoparticle[i].symbol+"_SK.npz")
        
        Ri=ase_nanoparticle.get_positions()[i,:]
        js=nl.get_neighbors(i)
        #onsite_args=Atom_i_SK[20:23]
        Onsite=onsite_model(*Atom_i_SK,**onsite_kwargs)
        Onsite-=np.eye(9)*A
        addHR(iidx,jidx,value,Onsite,i*9,i*9)
        if testing:
            f.write(f"Atom {i} has {len(js)} neighbors with following distances:\n")
        for j in (js):
            if j==i: 
                raise Exception("Error: j==i")
            try:
                Atom_j_SK=np.array(materials.SK_dict[ase_nanoparticle[j].symbol])
            except:
                materials.SK_dict[ase_nanoparticle[j].symbol]=np.load(ase_nanoparticle[j].symbol+"_SK.npz")
                Atom_j_SK=materials.SK_dict[ase_nanoparticle[j].symbol]
            
            Rj=ase_nanoparticle.get_positions()[j,:]
            
            #For now let's assume that i and j both have 9 orbitals
            Rij=Ri-Rj
            #hopping_args=(Atom_i_SK[:20]+Atom_j_SK[:20]).tolist().append((Atom_i_SK[23]+Atom_j_SK[23])/2)
            distance=np.linalg.norm(Rij)
            if testing:
                f.write(f"\tNeighbor: {j}, Distance: {distance:.3f}\n")
            HR=hopping_model(Rij,*((Atom_i_SK+Atom_j_SK)/2).tolist(),**hopping_kwargs)#-np.eye(9)*A
            addHR(iidx,jidx,value,HR,i*9,j*9)
        if testing:
            f.write("\n")
    iidx=np.array(iidx,dtype=np.int64)
    jidx=np.array(jidx,dtype=np.int64)
    value=(np.array(value,dtype=np.complex128)/B)
    if testing:
        f.close()
    return iidx,jidx,value


