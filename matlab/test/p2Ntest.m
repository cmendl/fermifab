function err = p2Ntest(orbs,p,N)
%P2NTEST - Verify that <psi|H|psi> = trace(h rdm(psi)) for a given H = sum_ij h_ij a^dagger_i a_j
%
%   P2NTEST(orbs,p,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(length(orbs)==1 && length(p)==1 && length(N)==1,'Only simple configurations supported.');

h = fermiop(orbs,p,p,crand(nchoosek(orbs,p)));

% random wavefunction 'psi'
psi = fermistate(orbs,N,crand(nchoosek(orbs,N),1));

err = abs(hprod(psi,p2N(h,N)*psi) - trace(h*rdm(psi,p)));
