import numpy as np 
from ase import Atoms 
from scipy import interpolate
from kpm.evaltool import eval_ase
from kpm.materialslib import Materials
import os 


def read_comsol_output(filename):
    """
    Reads a COMSOL output file with complex numbers.

    Parameters:
    - filename: str, the path to the COMSOL output file.

    Returns:
    - x, y, z: NumPy arrays containing the x, y, and z coordinates.
    - V: A NumPy array containing the complex potential values.
    """
    from kpm.units import Angstrom

    # Define a function to parse the complex numbers correctly
    def parse_complex_number(complex_str):
        return complex(complex_str.replace('i', 'j'))

    # Initialize lists to store your data
    x, y, z, V = [], [], [], []

    with open(filename, 'r') as file:
        for line in file:
            # Skip lines that do not start with a number
            if not line[0].isdigit() and not line[0] == '-':
                continue
            
            # Split the line into components
            parts = line.split()
            
            # Append the real numbers directly, parse the complex number
            x.append(float(parts[0]))
            y.append(float(parts[1]))
            z.append(float(parts[2]))
            V.append(parse_complex_number(parts[3]))

    # Convert lists to NumPy arrays
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    V = np.array(V)*Angstrom
    V-=np.mean(V)

    return x, y, z, V


class Shape:
    def inshape(self,region):
        return True 
    def __init__(self,nregion=1,epsm=1,omega=2.4,compute_phi_kwargs={}):
        
        #self.materials = [eval_ase(materials[i]) for i in range(nregion)]
        self.compute_phi_kwargs=compute_phi_kwargs
        #self._materials=materials
        self.nregion=nregion
        self.epsm=epsm
        self.compute_phi_kwargs=compute_phi_kwargs
        self.omega=omega
        
        
        pass 
    def inshape(self,R, region):
        return False

class Sphere(Shape):
    def __init__(self,radius,nregion=1,*args,**kwargs):
        super().__init__(nregion=1,**kwargs)
        self.radius=radius

    def inshape(self,R,region):
        import numpy as np 
        #Requires numpy >1.8.0 to use axis
        return np.linalg.norm(R,axis=1)<self.radius

    def compute_phi(self,positions,materials):
        if materials.epsw==None:
            raise Exception("Computation of phir requires epsw value of the sphereical nanoparticle")
        z=positions[:,2]
        phir= -3*self.epsm/(materials.epsw[0]+2*self.epsm)*z
        phiiidx=[]
        phijidx=[]
        phivalue=[]
        for i in range(z.shape[0]):
            for j in range(9):
                phiiidx.append(i*9+j)
                phijidx.append(i*9+j)
                phivalue.append(phir[i])
        return np.array(phiiidx,dtype=np.int64),np.array(phijidx,dtype=np.int64),np.array(phivalue,dtype=np.complex128)

class TouchedSphere(Shape):
    def __init__(self, radius=2, d=0.5, nregion=2, *args, **kwargs):
        super().__init__(nregion=nregion,**kwargs)  # Corrected super() call
        self.radius = float(radius)
        self.d = float(d)
        #We need to obtain the epsw data if it is absent 
    
    def inshape(self,R:np.ndarray,region:int):
        match region:
            case 0:
                Rcopy=R.copy()
                Rcopy[:,2]-=self.d 
                cond1=np.linalg.norm(Rcopy,axis=1)<self.radius 
                cond2=R[:,2]>=0
                return np.logical_and(cond1,cond2)
            case 1:
                Rcopy=R.copy()
                Rcopy[:,2]+=self.d 
                cond1=np.linalg.norm(Rcopy,axis=1)<self.radius
                cond2=R[:,2]<0
                return np.logical_and(cond1,cond2)

    def compute_phi(self,positions:np.ndarray,materials:Materials,Rmax=None,angle=0.0,mesh_size=1,pg="mul",V_th=0.05):
        if materials.epsw==None:
            raise Exception("Computation of phir requires epsw value of a touching Sphere")
        if Rmax==None:
            Rmax=20*self.radius

        if os.path.isfile("comsol_output.txt"):
            _,_,_,V=read_comsol_output("comsol_output.txt")
        
            phiiidx=[]
            phijidx=[]
            phivalue=[]
            for i in range(len(V)):
                phiiidx.extend([i*9+j for j in range(9)])
                phijidx.extend([i*9+j for j in range(9)])
                phivalue.extend([V[i]]*9)
            return np.array(phiiidx,dtype=np.int64),np.array(phijidx,dtype=np.int64),np.array(phivalue,dtype=np.complex128)

        else:
            import mph 
            client = mph.start()
            pymodel = client.create('Model')

            model = pymodel.java
            model.component().create("comp1", True);
            model.component("comp1").geom().create("geom1", 3);
            model.component('comp1').geom('geom1').create('sph1', 'Sphere');
            #create the upper sphere 
            model.component('comp1').geom('geom1').feature("sph1").set("pos", [0.0,0.0,self.d]);
            model.component('comp1').geom('geom1').feature("sph1").set("r", float(self.radius));

            model.component('comp1').geom('geom1').create('blk1', 'Block');
            model.component('comp1').geom('geom1').feature('blk1').set('base', 'center');
            model.component("comp1").geom("geom1").feature("blk1").set("pos",  [0.,0.,-self.radius]);
            model.component('comp1').geom('geom1').feature('blk1').set('size', [float(2*self.radius),float(2*self.radius),float(2*self.radius)]);
            model.component('comp1').geom('geom1').run('sph1');
            model.component('comp1').geom('geom1').run('blk1');
            model.component("comp1").geom("geom1").create("dif1", "Difference");
            model.component("comp1").geom("geom1").feature("dif1").selection("input").set("sph1");
            model.component("comp1").geom("geom1").feature("dif1").selection("input2").set("blk1");
            model.component("comp1").geom("geom1").run("dif1");

            #Now create the lower half
            model.component('comp1').geom('geom1').create('sph2', 'Sphere');
            #create the upper sphere 
            model.component('comp1').geom('geom1').feature("sph2").set("pos", [0,0,-self.d]);
            model.component('comp1').geom('geom1').feature("sph2").set("r", float(self.radius));

            model.component('comp1').geom('geom1').create('blk2', 'Block');
            model.component('comp1').geom('geom1').feature('blk2').set('base', 'center');
            model.component("comp1").geom("geom1").feature("blk2").set("pos",  [0.,0.,float(self.radius)]);
            model.component('comp1').geom('geom1').feature('blk2').set('size', [float(2*self.radius),float(2*self.radius),float(2*self.radius)]);
            model.component('comp1').geom('geom1').run('sph2');
            model.component('comp1').geom('geom1').run('blk2');
            model.component("comp1").geom("geom1").create("dif2", "Difference");
            model.component("comp1").geom("geom1").feature("dif2").selection("input").set("sph2");
            model.component("comp1").geom("geom1").feature("dif2").selection("input2").set("blk2");
            model.component("comp1").geom("geom1").run("dif2");

            model.component('comp1').geom('geom1').create('blk3', 'Block');
            model.component('comp1').geom('geom1').feature('blk3').set('base', 'center');
            model.component('comp1').geom('geom1').feature('blk3').set('size', [Rmax,Rmax,Rmax]);
            model.component('comp1').geom('geom1').run('blk3');

            model.component("comp1").material().create("mat1", "Common", 'comp1');
            model.component("comp1").material("mat1").propertyGroup("def").set("relpermittivity", str(self.epsm));

            model.component("comp1").material().create("mat2", "Common", 'comp1');
            model.component("comp1").material("mat2").propertyGroup("def").set("relpermittivity", str(materials.epsw[0]));

            model.component("comp1").material().create("mat3", "Common", 'comp1');
            model.component("comp1").material("mat3").propertyGroup("def").set("relpermittivity", str(materials.epsw[1]));

            model.component("comp1").geom("geom1").run();

            

            model.component("comp1").material("mat1").selection().set(1);
            model.component("comp1").material("mat2").selection().set(3);
            model.component("comp1").material("mat3").selection().set(2);


            #Create the selection
            model.component("comp1").geom("geom1").create("ballsel1", "BallSelection");
            model.component("comp1").geom("geom1").feature("ballsel1").set("entitydim", "2");
            model.component("comp1").geom("geom1").feature("ballsel1").set("posx", Rmax/2);
            model.component("comp1").geom("geom1").feature("ballsel1").set("r", 0.005);
            model.component("comp1").geom("geom1").run("ballsel1");
            model.component("comp1").geom("geom1").create("ballsel2", "BallSelection");
            model.component("comp1").geom("geom1").run("ballsel2");
            model.component("comp1").geom("geom1").feature("ballsel2").set("entitydim", "2");
            model.component("comp1").geom("geom1").feature("ballsel2").set("posx", -Rmax/2);
            model.component("comp1").geom("geom1").feature("ballsel2").set("r", 0.005);
            model.component("comp1").geom("geom1").run("ballsel2");

            model.component("comp1").geom("geom1").create("ballsel3", "BallSelection");
            model.component("comp1").geom("geom1").run("ballsel3");
            model.component("comp1").geom("geom1").feature("ballsel3").set("entitydim", "2");
            model.component("comp1").geom("geom1").feature("ballsel3").set("posy", Rmax/2);
            model.component("comp1").geom("geom1").feature("ballsel3").set("r", 0.005);
            model.component("comp1").geom("geom1").run("ballsel3");
            model.component("comp1").geom("geom1").create("ballsel4", "BallSelection");
            model.component("comp1").geom("geom1").run("ballsel4");
            model.component("comp1").geom("geom1").feature("ballsel4").set("entitydim", "2");
            model.component("comp1").geom("geom1").feature("ballsel4").set("posy", -Rmax/2);
            model.component("comp1").geom("geom1").feature("ballsel4").set("r", 0.005);
            model.component("comp1").geom("geom1").run("ballsel4");

            model.component("comp1").geom("geom1").create("ballsel5", "BallSelection");
            model.component("comp1").geom("geom1").run("ballsel5");
            model.component("comp1").geom("geom1").feature("ballsel5").set("entitydim", "2");
            model.component("comp1").geom("geom1").feature("ballsel5").set("posz", Rmax/2);
            model.component("comp1").geom("geom1").feature("ballsel5").set("r", 0.005);
            model.component("comp1").geom("geom1").run("ballsel5");
            model.component("comp1").geom("geom1").create("ballsel6", "BallSelection");
            model.component("comp1").geom("geom1").run("ballsel6");
            model.component("comp1").geom("geom1").feature("ballsel6").set("entitydim", "2");
            model.component("comp1").geom("geom1").feature("ballsel6").set("posz", -Rmax/2);
            model.component("comp1").geom("geom1").feature("ballsel6").set("r", 0.005);
            model.component("comp1").geom("geom1").run("ballsel6");

            model.component("comp1").geom("geom1").create("unisel1", "UnionSelection");
            model.component("comp1").geom("geom1").feature("unisel1").set("entitydim", "2");
            model.component("comp1").geom("geom1").feature("unisel1").set("input", ["ballsel1", "ballsel2", "ballsel3", "ballsel4", "ballsel5", "ballsel6"]);
            model.component("comp1").geom("geom1").run();

            model.component("comp1").physics().create("es", "Electrostatics", "geom2");
            model.component("comp1").physics("es").create("pot1", "ElectricPotential", 2);

            model.component("comp1").physics("es").feature("pot1").selection().named("geom1_unisel1");
            #model.component("comp1").physics("es").feature("pot1").selection().set([1, 2, 3, 4, 5,23]);
            model.component("comp1").physics("es").feature("pot1").set("V0", "(z*cos({}*pi/180)+x*sin({}*pi/180))[V/m]".format(angle,angle));

            model.component("comp1").mesh().create("mesh1");
            model.component("comp1").mesh("mesh1").autoMeshSize(mesh_size);
            model.component("comp1").mesh("mesh1").run();
            
            model.study().create("std1");

            model.study("std1").create("stat", "Stationary");
            model.study("std1").feature("stat").activate("es", True);


            model.sol().create("sol1");
            model.sol("sol1").study("std1");

            model.study("std1").feature("stat").set("notlistsolnum", "1");
            model.study("std1").feature("stat").set("notsolnum", "1");
            model.study("std1").feature("stat").set("listsolnum", "1");
            model.study("std1").feature("stat").set("solnum", "1");

            model.sol("sol1").create("st1", "StudyStep");
            model.sol("sol1").feature("st1").set("study", "std1");
            model.sol("sol1").feature("st1").set("studystep", "stat");
            model.sol("sol1").create("v1", "Variables");
            model.sol("sol1").feature("v1").set("control", "stat");
            model.sol("sol1").create("s1", "Stationary");
            model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
            model.sol("sol1").feature("s1").create("i1", "Iterative");
            model.sol("sol1").feature("s1").feature("i1").set("linsolver", "cg");
            model.sol("sol1").feature("s1").feature("i1").create("mg1", "Multigrid");
            model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("prefun", "amg");
            model.sol("sol1").feature("s1").feature("fc1").set("linsolver", "i1");
            model.sol("sol1").feature("s1").feature().remove("fcDef");
            model.sol("sol1").feature("s1").feature("i1").set("linsolver", "gmres");
            model.sol("sol1").feature("s1").feature("i1").set("itrestart", "200");
            model.sol("sol1").attach("std1");

            if pg=="iso":
                model.result().create("pg1", "PlotGroup3D");
                model.result("pg1").create("iso1", "Isosurface");
                model.result().param().set("V_th", "{}".format(V_th));
                
                model.result("pg1").feature("iso1").set("expr", "V-(z*cos({}*pi/180)+x*sin({}*pi/180))[V/m]".format(angle,angle));
                model.result("pg1").feature("iso1").set("levelmethod", "levels");
                model.result("pg1").feature("iso1").set("levels", "V_th");
                model.result("pg1").feature().duplicate("iso2", "iso1");
                model.result("pg1").feature("iso2").set("levels", "-V_th");
            else:
                model.result().create("pg1", "PlotGroup3D");
                model.result("pg1").label("Electric Potential (es)");
                model.result("pg1").set("frametype", "spatial");
                model.result("pg1").set("data", "dset1");
                model.result("pg1").feature().create("mslc1", "Multislice");
                model.result("pg1").feature("mslc1").set("showsolutionparams", "on");
                model.result("pg1").feature("mslc1").set("colortable", "RainbowLight");
                model.result("pg1").feature("mslc1").set("showsolutionparams", "on");
                model.result("pg1").feature("mslc1").set("data", "parent");
            

            

            model.sol("sol1").runAll();
            
            np.savetxt("comsol_input.txt",positions)
            model.result().export().create("data1", "Data");
            model.result().export("data1").setIndex("expr", "V", 0);
            model.result().export("data1").set("location", "file");
            model.result().export("data1").set("filename", "comsol_output.txt");
            model.result().export("data1").set("coordfilename", "comsol_input.txt");
            model.result().export("data1").run();

            model.save("comsol_output.mph")
            _,_,_,V=read_comsol_output("comsol_output.txt")
            
            phiiidx=[]
            phijidx=[]
            phivalue=[]
            for i in range(len(V)):
                phiiidx.extend([i*9+j for j in range(9)])
                phijidx.extend([i*9+j for j in range(9)])
                phivalue.extend([V[i]]*9)
            return np.array(phiiidx,dtype=np.int64),np.array(phijidx,dtype=np.int64),np.array(phivalue,dtype=np.complex128)
    def __eq__(self,other):
        ans=np.allclose(self.d,other.d)
        ans=ans&np.allclose(self.radius,other.radius)
        ans=ans&(self.nregion==other.nregion)
        ans=ans&np.allclose(self.epsm,other.epsm)
        ans=ans&np.allclose(self.omega,other.omega)
        return ans


        

 

        

if __name__=="__main__":
    pass