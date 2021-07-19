from quspin.operators import hamiltonian # Hamiltonians and operators
from quspin.basis import spinless_fermion_basis_general # Hilbert space for fermions
from quspin.basis.user import user_basis # Hilbert space user basis
from quspin.basis.user import next_state_sig_32,op_sig_32,map_sig_32,count_particles_sig_32 # user basis data types signatures
from numba import carray,cfunc,jit # numba helper functions
from numba import uint32,int32 # numba data types
import numpy as np
from scipy.special import comb
from gmpy2 import popcount



@jit
def comb_states(N_tot,N_imp,s_imp,up_states,dn_states):
	states = np.zeros(up_states.size*dn_states.size,dtype=np.uint32)
	i = 0
	
	for s_dn in dn_states:
		for s_up in up_states:
			states[i] = (s_imp << (2*N_tot)) + (s_dn << N_tot) + s_up
			i += 1

	return states

def generate_Sz_basis_Kondo(N_imp,N_tot,Nf,Sz):
	Ns_imp = 2**N_imp

	sub_sectors = []

	for s_imp in range(Ns_imp):
		Sup = popcount(s_imp)
		Sz_imp = 2*Sup - N_imp


		for Nup in range(Nf+1):
			Ndn = Nf - Nup

			if (Sz_imp + (Nup - Ndn)) == Sz:
				basis_up = spinless_fermion_basis_general(N_tot,Nf=Nup)
				basis_dn = spinless_fermion_basis_general(N_tot,Nf=Ndn)


				ferm_states = comb_states(N_tot,N_imp,s_imp,basis_up.states,basis_dn.states)

				sub_sectors.append(ferm_states)


	if sub_sectors:
		return np.hstack(sub_sectors)
	else:
		raise ValueError("combination of Sz and Nf not feasible")


######  function to read user-imported basis into QuSpin 

@cfunc(next_state_sig_32)
def next_state(s,counter,N,args):
    # return pre-calculated basis state.
    # add one to counter because the first state is already checked.
    return args[counter+1] # = basis

class function_wrapper(object):
    """
    This class provides a wrapper for the user-imported basis,
    as well as the functions required for the `user_basis` functionality.
    #
    This is needed to easily pass parameters (defined as class attributes) to the
    functions `get_so_pcon()` and `get_Ns_pcon`.
    """
    def __init__(self,basis):
        self.basis = basis
    #
    # python function to calculate the starting state to generate the particle conserving basis
    def get_s0_pcon(self,N,Np):
        """ calculates the starting state to generate the particle conserving basis. """
        # ignore input arguments as basis is already calculated.
        return self.basis[0]
    # 
    # python function to calculate the size of the particle-conserved basis, 
    # i.e. BEFORE applying pre_check_state and symmetry maps
    def get_Ns_pcon(self,N,Np):
        """ calculates the size of the particle conservation basis (ignoring symmetries at this stage). """
        # ignore input arguments as basis is already calculated.
        return self.basis.size


@jit(uint32(uint32,uint32),locals=dict(f_count=uint32,),nopython=True,nogil=True)
def _count_particles_32(state,site_ind):
	# auxiliary function to count number of fermions, i.e. 1's in bit configuration of the state, up to site site_ind
	# CAUTION: 32-bit integers code only!
	f_count = state & ((0x7FFFFFFF) >> (31-site_ind));
	f_count = f_count - ((f_count >> 1) & 0x55555555);
	f_count = (f_count & 0x33333333) + ((f_count >> 2) & 0x33333333);
	return (((f_count + (f_count >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24

@cfunc(op_sig_32,
	locals=dict(s=int32,sign=int32,n=int32,b=uint32,f_count=uint32), )
def op(op_struct_ptr,op_str,site_ind,N,args):
	# using struct pointer to pass op_struct_ptr back to C++ see numba Records
	op_struct = carray(op_struct_ptr,1)[0]
	err = 0

	if site_ind >= args[0]:
		site_ind = N - site_ind - 1
		f_count = _count_particles_32(op_struct.state,site_ind)
		sign = -1 if f_count&1 else 1
	else:
		site_ind = N - site_ind - 1
		sign = 1


	# convention, site_ind = N - bit_ind - 1

	n = (op_struct.state>>site_ind)&1 # either 0 or 1
	b = (1<<site_ind)

	if op_str==43: # "+" is integer value 43 = ord("+")
		op_struct.matrix_ele *= (0.0 if n else sign)
		op_struct.state ^= b 

	elif op_str==45: # "-" is integer value 45 = ord("-")
		op_struct.matrix_ele *= (sign if n else 0.0)
		op_struct.state ^= b 
		
	elif op_str==110: # "n" is integer value 110 = ord("n")
		op_struct.matrix_ele *= n

	elif op_str==122: # "z" is integer value 122 = ord("n")
		op_struct.matrix_ele *= (n-0.5)

	elif op_str==73: # "I" is integer value 73 = ord("I")
		pass

	else:
		op_struct.matrix_ele = 0
		err = -1

	return err

class site_mapping(object):
	def __init__(self,N_tot,N_imp):
		self._N_tot = N_tot
		self._N_imp = N_imp


	def site(self,site,part):
		if part == "u":
			if site >= 0 and site < self._N_tot:
				return self._N_imp+self._N_tot+site
			else:
				raise ValueError("site falls outside of system")
		elif part == "d":
			if site >= 0 and site < self._N_tot:
				return self._N_imp+site
			else:
				raise ValueError("site falls outside of system")
		elif part == "i":
			if site >= 0 and site < self._N_imp:
				return site
			else:
				raise ValueError("site falls outside of system")
		else:
			raise ValueError

def get_Kondo_basis(N_imp,N_tot,Nf,Sz):
	N = 2*N_tot+N_imp

	if N > 32:
		raise ValueError("System size too large")

	op_args=np.array([N_imp],dtype=np.uint32)
	kondo_states = generate_Sz_basis_Kondo(N_imp,N_tot,Nf,Sz)

	FW = function_wrapper(kondo_states)
	pcon_dict = dict(Np=(),next_state=next_state,next_state_args=kondo_states,
	                 get_Ns_pcon=FW.get_Ns_pcon,get_s0_pcon=FW.get_s0_pcon)
	op_dict = dict(op=op,op_args=op_args)
	noncommuting_bits = [(np.arange(2*N_tot),-1)]
	basis = user_basis(np.uint32,N,op_dict,allowed_ops=set("+-znI"),
		sps=2,pcon_dict=pcon_dict,noncommuting_bits=noncommuting_bits)

	return basis,site_mapping(N_tot,N_imp)

# np.set_printoptions(linewidth=1000000)

# L = 3
# J = 0.1

# basis,m = get_Kondo_basis(1,L,L,0)

# print(basis)


# hop_pm = [[ 1.0,m.site(i,p),m.site(i+1,p)] for p in ["d","u"] for i in range(L-1)]
# hop_mp = [[-1.0,m.site(i,p),m.site(i+1,p)] for p in ["d","u"] for i in range(L-1)]

# static = [["+-",hop_pm],["-+",hop_mp]]

# kwargs = dict(basis=basis,dtype=np.float64,check_symm=False,check_pcon=False,check_herm=False)


# H = hamiltonian(static,[],**kwargs)

# print(H-H.T)



