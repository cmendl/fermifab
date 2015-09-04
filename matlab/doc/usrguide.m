%% FermiFab User's Guide (Matlab version)
%
% <html><img style="float: right;" border: 0px; src="fermifab.png" alt="FermiFab" /></html>
%
% Christian B. Mendl, 2008-2015


%% 1.1  What is FermiFab?
%
% FermiFab stands for "fermion laboratory" and is a quantum physics toolbox
% for small numbers of fermions. It is mainly concerned with the
% representation of (symbolic) fermionic wavefunctions and the calculation
% of their reduced density matrices (RDMs). The toolbox transparently
% handles the inherent antisymmetrization of wavefunctions and incorporates
% the creation and annihilation formalism.
%
% FermiFab is available at <http://sourceforge.net/projects/fermifab>.


%% 1.2 Fermi states
% A fundamental building block of multi-fermion quantum mechanics are Slater
% determinants, which you can think of as a collection of "orbitals" (or slots),
% some of which are occupied by a fermionic particle (e.g., an electron):
%
% <<slater.png>>
%
% In this illustration, each circle represents an orbital, and the
% filled circles denote occupied orbitals.
% In mathematical terms, the available number of orbitals 'orbs' is the
% dimension of the underlying single-particle Hilbert space $\mathcal{H}$
% and the number of occupied orbitals is the particle number $N$.
% Thus there are altogether
%
% $${\mathrm{orbs} \choose N}$$
%
% Slater determinants. Each quantum mechanical Fermi state is a linear
% combination of these Slater determinants with complex coefficients.
% Formally, the quantum Fermi state is an element of the $N$-particle
% _Fock space_ $\wedge^N \mathcal{H}$.
%
% In FermiFab, you can define, say, a $N = 4$ particle state (or wavefunction)
% $\psi$ and 6 available orbitals with the |fermistate| command:
orbs = 6; N = 4;
x = zeros(nchoosek(orbs,N),1); x(1)=1/sqrt(2); x(2)=1i/sqrt(2);
psi = fermistate(orbs,N,x)

%%
% The vector $x$ contains the Slater determinant coefficients of $\psi$ in
% lexicographical order. Let's assign names to the orbitals of $\psi$:
psi = set(psi,'orbnames',{'1s' '1s~' '2s' '2s~' '2p' '2p~'})

%%
% In this example, the names denote electronic subshells of atoms,
% with the tilde ~ representing spin-down, otherwise spin-up.

%%
% The rank-one projector in bra-ket notation
%
% $$\vert\psi\rangle \langle\psi\vert$$
%
% or "density matrix" of $\psi$ can be calculated intuitively by
psi*psi'

%%
% Note that the result is now an _operator_ acting on the 4-particle
% antisymmetrized Hilbert space $\wedge^4 \mathcal{H}$ !


%% 1.3 Reduced density matrices (RDMs)
% The core feature of the toolbox is the efficient calculation of reduced density matrices.
% For example, the two-body reduced density matrix $\Gamma_{\psi}$ defined by
%
% $$\langle ij\vert \Gamma_{\psi} \vert kl\rangle
% = \langle\psi\vert a^{\dagger}_k a^{\dagger}_l a_j a_i \vert\psi\rangle$$
%
% can be calculated by the |rdm| command as
rdm(psi,2)

%%
% In the above formula, $a^{\dagger}_i$ and $a_i$ denote creation and
% annihilation operators, respectively.

%%
% More generally, for wavefunctions $\psi_k \in \wedge^{N_k}\mathcal{H}$
% and orthonormal basis sets $(\varphi_{ki})_i$ of $\wedge^{p_k}\mathcal{H}$
% with $1 \le p_k \le N_k$ and $k = 1,2$,
% define the reduced density matrix
%
% $$\gamma_{\vert\psi_1\rangle\langle\psi_2\vert}:
% \wedge^{p_2}\mathcal{H} \to \wedge^{p_1}\mathcal{H}$$
%
% by
%
% $$\langle\varphi_{1j}\vert \gamma_{\vert\psi_1\rangle\langle\psi_2\vert} \vert\varphi_{2i}\rangle
% := \langle\psi_2\vert a^{\dagger}_{\varphi_{2i}} a_{\varphi_{1j}} \vert\psi_1\rangle
%  = \langle a_{\varphi_{2i}} \psi_2 \vert a_{\varphi_{1j}} \psi_1\rangle.$$
%
% The significance of this definition can be seen as follows.
% Any linear map $b: \wedge^{p_1}\mathcal{H} \to \wedge^{p_2}\mathcal{H}$
% with matrix representation $\left(b_{ij}\right)$ may be "lifted" to
% an operator $B: \wedge^{N_1}\mathcal{H} \to \wedge^{N_2}\mathcal{H}$ by
%
% $$B = \sum_{i,j} b_{ij}\,a^\dagger_{\varphi_{2i}} a_{\varphi_{1j}},$$
%
% implemented as |p2N| in FermiFab.
% A prominent example is the Coulomb operator with $p_1 = p_2 = 2$,
% which describes the _pairwise_ interaction between charged particles.
% Now, the $B$ expectation value with respect to $\vert\psi_1\rangle\langle\psi_2\vert$ equals
%
% $$\langle\psi_2\vert B \vert\psi_1\rangle
% \equiv \sum_{i,j} b_{ij} \langle\psi_2\vert a^{\dagger}_{\varphi_{2i}} a_{\varphi_{1j}} \vert\psi_1\rangle
% = \mathrm{tr}\big[ b\,\gamma_{\vert\psi_1\rangle\langle\psi_2\vert}\big].$$
%
% In other words, this equation switches from $\wedge^{N_k}\mathcal{H}$ to
% $\wedge^{p_k}\mathcal{H}$ for $k = 1,2$.
% For many applications, this is the only possibility to avoid the
% "curse of dimensionality" induced by the large $N_1$, $N_2$-particle systems.


%% 1.4 Symbolic computations
% FermiFab is compatible with the "Symbolic Math Toolbox".
% So you may insert symbolic variables as well in the
% above examples (assuming the Symbolic Math Toolbox is installed)
syms a b c
y = sym(zeros(1,nchoosek(orbs,N)));
y(1) = a; y(3) = 1i*b^2; y(4) = 1/c;
psi = set(psi,'data',y)
rdm(psi,2)


%% 1.5 Tensor products of operators
% For any linear operator (i.e., matrix) $A: \mathcal{H} \to \mathcal{H}$,
% it holds that
%
% $$\langle j_1,\dots,j_N \vert (A \otimes\cdots\otimes A) \vert i_1,\dots,i_N \rangle
% = \det \langle j_k \vert A \vert i_\ell\rangle_{k,\ell}.$$
%
% That is, we obtain a matrix representation of $A \otimes\cdots\otimes A$
% acting on $\wedge^N \mathcal{H}$.
% The |tensor_op| command implements precisely this operation.
% The following code lines are taken from the "natural orbitals" example
% in test/norbs.m:
orbs = 6; N = 4;
psi = fermistate(orbs,N,crand(nchoosek(orbs,N),1));
[U,D] = eig(rdm(psi,1));

%%
% |crand| generates pseudorandom complex numbers (similar to |rand|),
% and |eig| computes eigenvalues and -vectors. Thus, the eigenvectors of the
% one-body reduced density matrix of $\psi$ are stored in $U$.
% Performing a corresponding base change on $\wedge^N \mathcal{H}$ using these eigenvectors
% should result in a diagonal one-body RDM [Löwdin 1955]:
psi = tensor_op(U,N)'*psi;
G = get(rdm(psi,1),'data');
err = norm(G-diag(diag(G)))

%%
% In many physical applications, one can take advantage of unitary base changes
% on $\mathcal{H}$ such that subsequent computations are simplified.
% The above code shows how to implement the according base change on $\wedge^N \mathcal{H}$.


%% 1.6 State configurations
% For performance and memory efficiency reasons, the concept of "configurations"
% are built into FermiFab, i.e. you can partition the orbitals
% into several groups, each containing a fixed number of particles.
% Let's say your system involves 3 particles in 9 orbitals, with exactly 2 particles
% in the first 5 orbitals and 1 particle in the remaining 4 orbitals. Then a
% |fermistate| reflecting this configuration is specified by
orbs = [5,4]; N = [2,1];
psi = fermistate(orbs,N)

%%
% Note that $\vert126\rangle$ is the lexicographically first basis vector
% respecting the configuration constraints, and that $\psi$ requires only
%
% $${5 \choose 2} \cdot {4 \choose 1} = 40$$
%
% rather than ${9 \choose 3} = 84$ complex entries:
length(psi)

%%
% The command |rdm| works transparently for all configurations,
% so $\psi$ behaves like a standard 9-orbital 3-particle state.

%%
% What happens if you add two |fermistate| wavefunctions with different but compatible
% configurations (i.e., the total number of orbitals and particles agrees)?
orbs = [2,7]; N = [1,2];
phi = fermistate(orbs,N)
%%
length(phi)
%%
% Adding $\psi$ and $\varphi$ gives
chi = psi+phi
%%
% as expected - so how is this accomplished? FermiFab has detected that
% it needs to merge the two configurations, resulting in the full-fledged
% 9-orbital 3-particle state, which is reflected by
length(chi)


%% 1.7 References
%
% * [Mendl 2011] C. B. Mendl.
%   _The FermiFab toolbox for fermionic many-particle quantum systems_.
%   Comput. Phys. Commun. 182, 1327-1337 (2011)
%   <http://dx.doi.org/10.1016/j.cpc.2011.01.028 DOI>,
%   preprint <http://arxiv.org/abs/1103.0872>
% * [Mazziotti 2007] D. A. Mazziotti.
%   _Reduced-density-matrix mechanics: with application to many-electron atoms and molecules_.
%   Advances in Chemical Physics 134, Wiley (2007)
% * [Coleman & Yukalov 2000] A. J. Coleman and V. I. Yukalov.
%   _Reduced density matrices: Coulson's challenge_.
%   Lecture Notes in Chemistry 72, Springer (2000)
% * [Ando 1963] T. Ando.
%   _Properties of Fermion density matrices_.
%   Rev. Mod. Phys. 35, 690–702 (1963)
%   <http://dx.doi.org/10.1103/RevModPhys.35.690 DOI>
% * [Löwdin 1955] P.-O. Löwdin.
%   _Quantum theory of many-particle systems I_.
%   Phys. Rev. 97, 1474–1489 (1955)
%   <http://dx.doi.org/10.1103/PhysRev.97.1474 DOI>
