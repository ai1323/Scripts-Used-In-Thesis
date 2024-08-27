def eval_ase(expression):
    from ase import Atom 
    from ase import Atoms 
    from ase import build 
    from ase import io
    from ase import collections
    from ase import db
    import numpy as np
    fccmaterial=["Au","Ag","Cu"]
    if (expression[:2] in fccmaterial):
        a0=float(expression[3:])
        return Atoms(expression[:2],cell=np.array([[1/2,1/2,0],[-1/2,1/2,0],[0,1/2,1/2]])*a0)
    else:
        return eval(expression)
