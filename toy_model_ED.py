from quspin.basis import spin_basis_general,spinful_fermion_basis_general,spinless_fermion_basis_general,tensor_basis
from quspin.operators import hamiltonian,quantum_operator
from quspin.tools.evolution import expm_multiply_parallel
import numpy as np
from scipy.linalg import eigh_tridiagonal,block_diag
from scipy.sparse.linalg import eigsh
from scipy.sparse import eye
from Kondo_basis import get_Kondo_basis
import sys


import matplotlib
matplotlib.rcParams['lines.linewidth'] = 1
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['legend.frameon'] = False


import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

fig_ratio = 1.3


def get_graphene(radgr,mapping_file):
	d = np.loadtxt(mapping_file)
	if radgr == 0:
		return np.array([],dtype=np.float64),np.array([],dtype=np.float64)
	else:
		return d[:radgr,0],d[:radgr-1,1]

def solve_single_particle(a,b):
	E,V = eigh_tridiagonal(a,b)
	C = np.hstack(([0],E.cumsum()))
	i = C.argmin()
	return int(i),C[i]


def get_H(a,b,Nf_tot=None,Sz=None):
	L = a.size

	if Nf_tot is None:
		nf,E_sp = solve_single_particle(a,b)
		Nf_tot = 2*nf

	if Sz is None:
		Sz = 0


	basis,m = get_Kondo_basis(2,L,Nf_tot,Sz)


	Jzz_list = [[ 0.5,m.site(0,"i"),m.site( 0 ,"u")],[-0.5,m.site(0,"i"),m.site( 0 ,"d")],
				[ 0.5,m.site(1,"i"),m.site(L-1,"u")],[-0.5,m.site(1,"i"),m.site(L-1,"d")]
				]

	Jpm_list = [
				[+0.5,m.site(0,"i"),m.site( 0 ,"u"),m.site( 0 ,"d")],
				[+0.5,m.site(1,"i"),m.site(L-1,"u"),m.site(L-1,"d")]
				]

	Jmp_list = [
				[-0.5,m.site(0,"i"),m.site( 0 ,"u"),m.site( 0 ,"d")],
				[-0.5,m.site(1,"i"),m.site(L-1,"u"),m.site(L-1,"d")]
				]

	hop_pm = ([[ b[i],m.site(i,"u"),m.site(i+1,"u")] for i in range(L-1)]+
			  [[ b[i],m.site(i,"d"),m.site(i+1,"d")] for i in range(L-1)])

	hop_mp = ([[-b[i],m.site(i,"u"),m.site(i+1,"u")] for i in range(L-1)]+
			  [[-b[i],m.site(i,"d"),m.site(i+1,"d")] for i in range(L-1)])

	mu_list = ([[a[i],m.site(i,"u")] for i in range(L)]+
			   [[a[i],m.site(i,"d")] for i in range(L)])

	n_list  = ([[1,m.site(i,"u")] for i in range(L)]+
			   [[1,m.site(i,"d")] for i in range(L)])
	static_SP = [
				["+-",hop_pm],
				["-+",hop_mp],
				]

	static_V  = [
				["zn",Jzz_list],
				["-+-",Jpm_list],
				["+-+",Jmp_list],
				]

	static_mu = [["n",mu_list],]
	kwargs = dict(check_symm=False,check_pcon=False,
		check_herm=False,basis=basis,dtype=np.float64)

	H = quantum_operator(dict(H0=static_SP,J=static_V,mu=static_mu),**kwargs)

	Czz_list = [[1.0,m.site(0,"i"),m.site(1,"i")]]
	Cxy_list = [[0.5,m.site(0,"i"),m.site(1,"i")]]
	Sz_list = [[1.0,m.site(0,"i")]]
	static_Czz = [["zz",Czz_list]]
	static_Cxy = [["+-",Cxy_list],["-+",Cxy_list]]
	static_Sz = [["z",Sz_list]]

	Czz = hamiltonian(static_Czz,[],**kwargs)
	Cxy = hamiltonian(static_Cxy,[],**kwargs)
	Sz = hamiltonian(static_Sz,[],**kwargs)

	rho_op = {"uu":{},"dd":{},"SS":{}}
	# print(Nf_tot,Sz)
	# print(basis)

	Jzz = [[1.0,m.site(0,"i"),m.site(1,"i")]]
	Jxy = [[0.5,m.site(0,"i"),m.site(1,"i")]]
	static_SS = [["zz",Jzz],["+-",Jxy],["-+",Jxy]]

	SS = hamiltonian(static_SS,[],**kwargs)
	PS = (0.25*eye(basis.Ns,basis.Ns)-SS)
	PT = (0.75*eye(basis.Ns,basis.Ns)+SS)

	for i in range(L):
		for j in range(L):
			for s1,s2 in [("u","u"),("d","d")]:

				static = [["+-",[[1.0,m.site(i,s1),m.site(j,s2)]]]]
				op =  hamiltonian(static,[],**kwargs)

				rho_op[s1+s2][(i,j)] = op


			Jzz = [[ 0.5,m.site(0,"i"),m.site(i,"u"),m.site(i,"u"),m.site(j,"u")],
				   [-0.5,m.site(0,"i"),m.site(i,"d"),m.site(i,"u"),m.site(j,"u")],
				   ]

			static = [
					  ["zn+-",Jzz],
					  ]
			op =  hamiltonian(static,[],**kwargs)
			rho_op["SS"][(i,j)] = op

	K_dict = {"K{}-{}-{}".format(spin,x,y):(spin,x,y)for x in range(L) for y in range(L) for spin in ["u","d"]}
	ops_dict = {}

	for key,(spin,x,y) in K_dict.items():
			K_list = [[1.0,m.site(x,spin),m.site(y,spin)]]
			ops_dict[key] = [["+-",K_list]]


	K = quantum_operator(ops_dict,**kwargs)


	return dict(H=H,m=m,Sz=Sz,Czz=Czz,Cxy=Cxy,rho_op=rho_op,PS=PS,PT=PT,K=K,K_dict=K_dict)

def mu_plots_NPT(Lc,radgr,savefig,ext="pdf"):
	"""
	This function produces Figures. XX and XX 
	"""

	mu=0.0
	tc=1.0
	mugr=0.0
	tprime=1.0
	J=0.5
	gap = 0

	L = 2*radgr+Lc

	sim = dict(Lc=Lc,radgr=radgr,tc=tc,tprime=tprime,J=J,gap=gap)

	mapping_file = "mappings/single_site_graphene.txt"


	np.set_printoptions(linewidth=1000000)
	a_gr,b_gr = get_graphene(radgr,mapping_file)
	a_ch = np.zeros(Lc,dtype=np.float64)
	b_ch = np.zeros(Lc-1,dtype=np.float64)
	a_ch[:] = 1.0
	b_ch[:] = tc + (gap/4)*(-1)**np.arange(Lc-1)
	A = np.hstack((a_gr,a_ch,a_gr[::-1]))
	B = -np.hstack((b_gr,[tprime],b_ch,[tprime],b_gr[::-1]))


	systems = {}

	for nf in range(1,2*L,1):
		systems[nf] = get_H(A,B,Nf_tot=nf,Sz=nf%2)


	# mu_list = np.around(np.linspace(-5,5,501),10)
	# slices = {}
	# mu_slice = [-0.5,0,0.1]

	mu_list = np.around(np.linspace(-4,4,401),10)
	slices = {}
	mu_slice = [1.0,-0.5]

	for mu in mu_slice:
		if mu not in mu_list:
			i = np.searchsorted(mu_list,mu)
			mu_list = np.insert(mu_list,i,mu)


	labels= {mu_slice[0]:{"a":"$(a)$","b":"$(c)$","c":"$(e)$"},
			 mu_slice[1]:{"a":"$(b)$","b":"$(d)$","c":"$(f)$"}
	}

	C_list = []
	Nf_list = []
	E_hyb_list = []
	E_nhyb_list = []

	for mu in mu_list:
		a_gr,b_gr = get_graphene(radgr,mapping_file)
		a_ch = np.zeros(Lc,dtype=np.float64)
		b_ch = np.zeros(Lc-1,dtype=np.float64)
		a_ch[:] = -mu
		b_ch[:] = tc + (gap/4)*(-1)**np.arange(Lc-1)
		A = np.hstack((a_gr,a_ch,a_gr[::-1]))
		B1 = -np.hstack((b_gr,[tprime],b_ch,[tprime],b_gr[::-1]))
		B2 = -np.hstack((b_gr,[0],b_ch,[0],b_gr[::-1]))

		E_hyb, V_hyb  = eigh_tridiagonal(A,B1)
		E_nhyb, V_nhyb = eigh_tridiagonal(A,B2)
		E_hyb_list.append(E_hyb)
		E_nhyb_list.append(E_nhyb)

		h_SP = (V_hyb*E_hyb).dot(V_hyb.T)


		H_SP = hamiltonian([h_SP],[],dtype=np.float64)


		C = np.hstack(([0],E_hyb.cumsum()))
		i = int(C.argmin())

		L = A.size

		E = np.inf
		Nf = 2*i-1
		psi0 = None

		pars = dict(H0=1,J=J,mu=-mu)

		for new_Nf in range(2*i-4,2*i+5,1):
			if new_Nf in systems:

				H = systems[new_Nf]["H"].aslinearoperator(pars=pars)
				[new_E],new_psi0 = eigsh(H,k=1,which="SA")

				if new_E < E:
					E = new_E
					psi0 = new_psi0.ravel()
					Nf = new_Nf


		Nf_list.append(Nf)
		rho = np.zeros((L,L))

		for (i,j),O in systems[Nf]["rho_op"]["dd"].items():
			rho[i,j] = O.expt_value(psi0)

		N_no_dn,V_no_dn = np.linalg.eigh(rho)

		V_no_dn = V_no_dn[:,::-1].copy()
		N_no_dn = N_no_dn[::-1].copy()

		for (i,j),O in systems[Nf]["rho_op"]["uu"].items():
			rho[i,j] = O.expt_value(psi0)

		N_no_up,V_no_up = np.linalg.eigh(rho)

		N_no = N_no_up[::-1].copy()
		V_no = V_no_up[:,::-1].copy()

		N_no_up = N_no_up[::-1].copy()
		V_no_up = V_no_up[:,::-1].copy()


		Czz = systems[Nf]["Czz"].expt_value(psi0)
		Cxy = systems[Nf]["Cxy"].expt_value(psi0)

		C_list.append(Czz)
		print(mu,Czz)

		if mu in mu_slice:
			K_dict = systems[Nf]["K_dict"]
			K = systems[Nf]["K"]
			m = systems[Nf]["m"]

			e_up,u_up = np.linalg.eig(V_no_up.T)
			e_dn,u_dn = np.linalg.eig(V_no_dn.T)

			k = {}
			k["u"] = ((np.log(e_up)*u_up).dot(u_up.T.conj()))
			k["d"] = ((np.log(e_dn)*u_dn).dot(u_dn.T.conj()))

			pars = {key:k[spin][x,y] for key,(spin,x,y) in K_dict.items()}
			psi1 = expm_multiply_parallel(K.tocsr(pars=pars)).dot(psi0)
			
			sents = np.zeros((L,L))
			sents_imp = np.zeros((L,))

			for i in range(L):
				# sub_sys_A = [m.site(0,"i"),m.site(i,"u"),m.site(i,"d")]
				# sents_imp[i] = K.basis.ent_entropy(psi1,sub_sys_A=sub_sys_A,density=False)["Sent_A"]

				for j in range(L):
					if i==j:
						sub_sys_A = ([m.site(o,"i") for o in range(2)]+
									 [m.site(i,"u"),m.site(i,"d")])
						
					else:

						sub_sys_A = ([m.site(o,"i") for o in range(2)]+
									 [m.site(i,"u"),m.site(i,"d")]+
									 [m.site(j,"u"),m.site(j,"d")])

					sents[i,j] = K.basis.ent_entropy(psi1,sub_sys_A=sub_sys_A,density=False)["Sent_A"]/np.log(2)

			slices[mu] = (sents,N_no,N_no_dn,V_no,sents_imp)


	C_list = np.array(C_list)
	Nf_list = np.array(Nf_list)
	E_hyb_list = np.array(E_hyb_list)
	E_nhyb_list = np.array(E_nhyb_list)

	print()
	w = 3.5

	f,(ax2,ax1) = plt.subplots(2,1,figsize=(w,1.7*w/fig_ratio),sharex="col")


	for mu in mu_slice:
		(sents,N_no,N_no_dn,V_no,sents_imp) = slices[mu]
		ax1.axvline(-mu,ymin=-3,ymax=3,color="red",linestyle=":")
		# ax1.text(mu+0.1,6.1,label,fontsize=12,color="black")

	ln, = ax1.plot(-mu_list,np.zeros_like(mu_list),linestyle=":")
	ax1.plot(-mu_list,E_hyb_list,color=ln.get_color())
	
	# ax1.set_xlabel(r"$\mu$",fontsize=12)
	ax1.set_ylabel(r"$E_{\rm electron}$",fontsize=12,labelpad=0)
	ax1.text(-0.2,0.95,"$(b)$",fontsize=12,transform=ax1.transAxes)
	# ax1.grid(linestyle=":")


	for label,mu in zip([r"$\mathrm{FM}$",r"$\mathrm{AFM}$"],mu_slice):
		(sents,N_no,N_no_dn,V_no,sents_imp) = slices[mu]
		ax2.axvline(-mu,ymin=-3,ymax=3,color="red",linestyle=":")
		x = (mu_list[-1]-mu)/(mu_list[-1]-mu_list[0])
		ax2.text(x,1.01,label,fontsize=12,color="black",ha="center",va="bottom",transform=ax2.transAxes)

	ax2.plot(-mu_list,C_list)
	ax1.set_xlabel(r"$V_g$",fontsize=12)
	ax2.set_ylabel(r"$\langle 0|S^z_1S^z_2|0\rangle$",fontsize=12,labelpad=0)

	# ax2.set_ylim((-0.8,0.5))
	ax2.text(-0.2,0.95,"$(a)$",fontsize=12,transform=ax2.transAxes)
	# ax2.grid(linestyle=":")

	ax1.xaxis.set_major_locator(ticker.MultipleLocator(1))
	ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
	ax1.yaxis.set_major_locator(ticker.MultipleLocator(2))
	ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
	ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
	ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.02))
	f.tight_layout()
	f.subplots_adjust(hspace=0.03)
	if savefig:
		f.savefig("../corr_Lc_{}_J_{}.{}".format(Lc,J,ext),bbox_inches="tight",pad_inches=0.025)

	w = 4.5
	h = 6.0
	f,Axs = plt.subplots(3,len(mu_slice),figsize=(w,h),sharey="row")
	Axs = np.asarray(Axs)

	Axs[0,0].set_title(r"$\mathrm{FM}$")
	Axs[0,1].set_title(r"$\mathrm{AFM}$")
	
	Axs[[1,2],:] = Axs[[2,1],:]

	for mu,axs in zip(mu_slice,Axs.T[:]):
		(sents,N_no,N_no_dn,V_no,sents_imp) = slices[mu]
		ma = np.logical_not(np.isinf(sents))

		sent_min = sents[ma].min()
		sent_max = sents[ma].max()

		ind = np.argwhere(sents==sent_min)

		i,j = ind[0]

		if N_no[i] < N_no[j]:
			axs[0].plot(V_no[:,i],marker="x",label=r"$|\phi_{}\rangle$".format(i))
			if j!=i:
				axs[0].plot(V_no[:,j],marker="+",label=r"$|\phi_{}\rangle$".format(j))
		else:
			axs[0].plot(V_no[:,j],marker="x",label=r"$|\phi_{}\rangle$".format(j))
			if j!=i:
				axs[0].plot(V_no[:,i],marker="+",label=r"$|\phi_{}\rangle$".format(i))			

		axs[0].text(0.125,0.075,labels[mu]["a"],fontsize=10,transform=axs[0].transAxes)
		axs[2].text(0.125,0.075,labels[mu]["b"],fontsize=10,transform=axs[2].transAxes)
		axs[1].text(0.125,0.075,labels[mu]["c"],fontsize=10,transform=axs[1].transAxes)

		alpha = [r"$\phi_{}$".format(i) for i in range(L)]
		im = axs[1].matshow(sents)

		axs[0].set_xlabel("$x$",labelpad=0)
		axs[0].set_ylabel(r"$\langle x|\phi_i\rangle$")
		axs[0].set_xticks(list(range(L)))
		axs[0].legend(ncol=2,fontsize=8,loc="upper center")
		axs[0].set_ylim((-0.65,0.75))

		axs[1].set_xticklabels(alpha,fontsize=10)
		axs[1].set_xticks(list(range(L)))
		axs[1].set_yticklabels(['']+alpha,fontsize=10)


		divider = make_axes_locatable(axs[1])
		# cax = divider.append_axes('right', size='5%', pad=0.05)
		# f.colorbar(im, cax=cax, orientation='vertical')
		cax = divider.new_vertical(size='5%', pad=0.05, pack_start=True)
		f.add_axes(cax)
		f.colorbar(im, cax=cax, orientation='horizontal',ticks = [sents.min(),sents.max()])
		axs[1].set_xticks(list(range(L)))
		axs[1].xaxis.tick_top()
		axs[1].xaxis.set_ticks_position('both')



		axs[2].plot(N_no,marker=".",label=r"$\sigma=\uparrow$")
		axs[2].plot(N_no_dn,marker=".",label=r"$\sigma=\downarrow$")
		axs[2].legend(fontsize=8,loc="upper right")
		axs[2].set_ylabel(r"$n_{i\sigma}$")
		axs[2].set_xticks(list(range(L)))
		axs[2].set_xticklabels(alpha,fontsize=10)


	f.tight_layout()
	f.subplots_adjust(hspace=0.25)
	if savefig:
		f.savefig("../slices_Lc_{}_J_{}..{}".format(Lc,J,ext),bbox_inches="tight",pad_inches = 0.03)
	else:
		plt.show()



mu_plots_NPT(7,1,True,"pdf")
