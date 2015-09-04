function err = rdm_single_test(orbs,N)
%RDM_SINGLE_TEST - Single-particle reduced density matrix test file
%
%   err = RDM_SINGLE_TEST(orbs,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(mod(orbs,2) == 0,'Number of orbitals must be even.');

% fictitious operator acting on spatial modes only
h = bosonop(orbs/2,1,1,crand(orbs/2));

% 1-body Hamiltonian acting on orbitals
% containing position and spin
H = fermiop(orbs,1,1,kron(get(h,'data'),eye(2)));

psi = fermistate(orbs,N,crand(nchoosek(orbs,N),1));

err = ...
	abs(hprod(psi,p2N(H,N)*psi) - trace(H*rdm(psi,1))) + ...
	abs(trace(H*rdm(psi,1)) - trace(h*spintrace_sing(rdm(psi,1))));
