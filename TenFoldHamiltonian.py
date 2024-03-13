import numpy as np
def OccupiedStatesCount(Eigen):#count the number of negative eigen-values.
    number=0
    for iindex in range(len(Eigen)):
        if(Eigen[iindex]<0):
            number=number+1
        else:
            pass
    return number
class TenFoldHamiltonian:
    dimension=0#real space dimension of the system
    lattice_size=[]#numbers of lattice sites in each dimension
    spin=0#internal degree of freedom,total spin j
    symmetry_vec=[]#symmetry vectors for [T^2,C^2,S]
    ###the input argument of a random model

    internal_dim = 0  # internal degree of freedom in each lattice site
    internal_pair_num = 0 #the number of the spin pairs in internal space
    total_lattice_sites_number=0#the total number of space lattice sites
    total_hamilt_dim=0#the total dimension of the hamiltonian metrics with internal degree
    ###some information that are usually used.

    unitary_time_reversal=[]#the overall time reversal operator
    unitary_particle_hole=[]#the overall particle hole operator
    unitary_chiral=[]#the overall chiral operator
    ###the symmetry operators

    tenfold_hamilt=[]#the random hamiltonian with the given symmetry
    topology_invariant=0#the corresponding topological invariant
    ###the random tenfold hamiltonian and the corresponding topological invariant

    def __init__(self,dim,lat_size,sp,sym):#initialization function for the simple parameters
        #the basic parameters
        self.dimension=dim
        self.lattice_size=lat_size
        self.spin=sp
        self.symmetry_vec=sym

        #the usually used parameters
        self.internal_dim=int(self.spin*2+1)
        self.internal_pair_num=int(self.internal_dim/2)

        total_sites_number = 0#calculate the total number of lattice sites
        if (self.dimension == 0):
            total_sites_number = 1
        else:
            total_sites_number = 1
            for iindex in range(len(self.lattice_size)):
                total_sites_number = total_sites_number * self.lattice_size[iindex]
        self.total_lattice_sites_number = total_sites_number

        #the total dimension of the hamiltonian including internal degree
        self.total_hamilt_dim=int(self.total_lattice_sites_number*self.internal_dim)
    def ConstructSymmetryOperators(self):#construct the symmetry operators
        internal_dim=self.internal_dim
        internal_pair_num=self.internal_pair_num
        one_site_T=[]
        one_site_P=[]
        one_site_S=[]

        ###construct the one site time reversal operator
        if (self.symmetry_vec[0]==1):
            one_site_T=np.kron(np.array([[0,1],[1,0]]),np.eye(internal_pair_num))
        elif (self.symmetry_vec[0]==0):
            one_site_T=np.kron(np.array([[1, 0], [0, 1]]), np.eye(internal_pair_num))
        else:
            one_site_T=np.kron(np.array([[0, 1], [-1, 0]]), np.eye(internal_pair_num))

        ###construct the one site particle hole operator
        if (self.symmetry_vec[1]==1):
            one_site_P=np.kron(np.array([[1,0],[0,-1]]),np.eye(internal_pair_num))
        elif (self.symmetry_vec[1]==0):
            one_site_P = np.kron(np.array([[1, 0], [0, 1]]), np.eye(internal_pair_num))
        else:
            one_site_P = np.kron(np.array([[0, 1], [-1, 0]]), np.eye(internal_pair_num))

        ###construct the one site chiral operator
        if (self.symmetry_vec[0]*self.symmetry_vec[1]==0):#
            if (self.symmetry_vec[2]==0):
                one_site_S=np.kron(np.array([[1, 0], [0, 1]]), np.eye(internal_pair_num))##witout chiral symmetry
            else:
                one_site_S=np.kron(np.array([[1, 0], [0, -1]]), np.eye(internal_pair_num))##only with chiral symmetry
        elif (self.symmetry_vec[0]==1 and self.symmetry_vec[1]==1):
            one_site_S =1j*np.dot(one_site_P,one_site_T)##BDI class should add a phase factor to chiral operator so that S^2==1
        else:
            one_site_S=np.dot(one_site_P,one_site_T)

        ##constuct the overall symmetry operators
        self.unitary_time_reversal=np.kron(np.eye(self.total_lattice_sites_number),one_site_T)
        self.unitary_particle_hole=np.kron(np.eye(self.total_lattice_sites_number),one_site_P)
        self.unitary_chiral=np.kron(np.eye(self.total_lattice_sites_number),one_site_S)

    ####functions to connstruct the random ten fold hamiltonian
    def ConstructTenfoldHamiltonian(self):
        #read out the global parameters
        internal_dim = self.internal_dim
        internal_pair_num = self.internal_pair_num
        total_hamilt_dim = self.total_hamilt_dim
        hamilt = []#the hamiltonian
        sym_vec = self.symmetry_vec#symmetry vector
        if(self.dimension==0):#the zero dimensional case
            ###create a random hamiltonian without any symmetry
            rd = np.random.RandomState()
            hamilt = rd.uniform(-1, 1, (total_hamilt_dim,total_hamilt_dim)) + 1j * rd.uniform(-1, 1,(total_hamilt_dim,total_hamilt_dim))
            hamilt =0.5*(hamilt + hamilt.T.conjugate())
            if (sym_vec[0] == 0 and sym_vec[1] == 0 and sym_vec[2] == 0):# class A
                self.tenfold_hamilt = hamilt
            elif (sym_vec[0] == 0 and sym_vec[1] == 0 and sym_vec[2] == 1):#class AIII
                hamilt = hamilt - self.unitary_chiral @ hamilt @ self.unitary_chiral
                self.tenfold_hamilt = hamilt
            elif (sym_vec[0] == 1 and sym_vec[1] == 0):# class AI
                pass
            elif (sym_vec[0] == 1 and sym_vec[1] == 1):# class BDI
                pass
            elif (sym_vec[0] == 0 and sym_vec[1] == 1):# class D
                pass
            elif (sym_vec[0] == -1 and sym_vec[1] == 1):# class DIII
                pass
            elif (sym_vec[0] == -1 and sym_vec[1] == 0):# class AII
                pass
            elif (sym_vec[0] == -1 and sym_vec[1] == -1):# class CII
                pass
            elif (sym_vec[0] == 0 and sym_vec[1] == -1):# class C
                pass
            else:# class CI
                pass
        elif(self.dimension==1):#the one dimensional case
            if (sym_vec[0] == 0 and sym_vec[1] == 0 and sym_vec[2] == 0):# class A
                pass
            elif (sym_vec[0] == 0 and sym_vec[1] == 0 and sym_vec[2] == 1):# class AIII
                pass
            elif (sym_vec[0] == 1 and sym_vec[1] == 0):# class AI
                pass
            elif (sym_vec[0] == 1 and sym_vec[1] == 1): # class BDI
                pass
            elif (sym_vec[0] == 0 and sym_vec[1] == 1): # class D
                pass
            elif (sym_vec[0] == -1 and sym_vec[1] == 1):# class DIII
                pass
            elif (sym_vec[0] == -1 and sym_vec[1] == 0):# class AII
                pass
            elif (sym_vec[0] == -1 and sym_vec[1] == -1):# class CII
                pass
            elif (sym_vec[0] == 0 and sym_vec[1] == -1): # class C
                pass
            else:
                pass  # class CI
        else:#other case remained further consideration
            pass

    ####functions to calculate the topological invariant for the previous hamiltonian
    def CalculateTopologicalInvariant(self):
        hamilt=self.tenfold_hamilt
        sym_vec = self.symmetry_vec  # symmetry vector
        eigen_val=[]
        eigen_vec=[]
        if (self.dimension == 0):  # the zero dimensional case
            if (sym_vec[0] == 0 and sym_vec[1] == 0 and sym_vec[2] == 0):  # class A
                [eigen_val,eigen_vec] = np.linalg.eig(hamilt)
                self.topology_invariant = OccupiedStatesCount(eigen_val)
            elif (sym_vec[0] == 0 and sym_vec[1] == 0 and sym_vec[2] == 1):  # class AIII
                [eigen_val,eigen_vec] = np.linalg.eig(hamilt)
                self.topology_invariant = OccupiedStatesCount(eigen_val) - int(self.total_hamilt_dim / 2)
            elif (sym_vec[0] == 1 and sym_vec[1] == 0):  # class AI
                pass
            elif (sym_vec[0] == 1 and sym_vec[1] == 1):  # class BDI
                pass
            elif (sym_vec[0] == 0 and sym_vec[1] == 1):  # class D
                pass
            elif (sym_vec[0] == -1 and sym_vec[1] == 1):  # class DIII
                pass
            elif (sym_vec[0] == -1 and sym_vec[1] == 0):  # class AII
                pass
            elif (sym_vec[0] == -1 and sym_vec[1] == -1):  # class CII
                pass
            elif (sym_vec[0] == 0 and sym_vec[1] == -1):  # class C
                pass
            else:  # class CI
                pass
        elif (self.dimension == 1):  # the one dimensional case
            if (sym_vec[0] == 0 and sym_vec[1] == 0 and sym_vec[2] == 0):  # class A
                pass
            elif (sym_vec[0] == 0 and sym_vec[1] == 0 and sym_vec[2] == 1):  # class AIII
                pass
            elif (sym_vec[0] == 1 and sym_vec[1] == 0):  # class AI
                pass
            elif (sym_vec[0] == 1 and sym_vec[1] == 1):  # class BDI
                pass
            elif (sym_vec[0] == 0 and sym_vec[1] == 1):  # class D
                pass
            elif (sym_vec[0] == -1 and sym_vec[1] == 1):  # class DIII
                pass
            elif (sym_vec[0] == -1 and sym_vec[1] == 0):  # class AII
                pass
            elif (sym_vec[0] == -1 and sym_vec[1] == -1):  # class CII
                pass
            elif (sym_vec[0] == 0 and sym_vec[1] == -1):  # class C
                pass
            else:
                pass  # class CI
        else:  # other case remained further consideration
            pass