function err = config_test(orbs1,orbs2,N1,N2)
%CONFIG_TEST - Particle configuration test file
%
%   err = CONFIG_TEST(orbs1,orbs2,N1,N2)
%
%   Example: CONFIG_TEST([4,5],[3,6],[3,2],[1,4])
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(sum(orbs1)==sum(orbs2) && sum(N1)==sum(N2),'Total number of orbitals and particles must be equal.');

psi1 = fermistate(orbs1,N1);
psi1 = set(psi1,'data',crand(length(psi1),1));
chi1 = fermistate(sum(orbs1),sum(N1),gen_comprise(orbs1,N1,true(length(orbs1)-1,1))*get(psi1,'data'));

psi2 = fermistate(orbs2,N2);
psi2 = set(psi2,'data',crand(length(psi2),1));
chi2 = fermistate(sum(orbs2),sum(N2),gen_comprise(orbs2,N2,true(length(orbs2)-1,1))*get(psi2,'data'));

% plus uses 'gen_comprise' if necessary
psi = psi1+psi2;
orbs = get(psi,'orbs'); N = get(psi,'N');
chi = fermistate(sum(orbs),sum(N),gen_comprise(orbs,N,true(length(orbs)-1,1))*get(psi,'data'));

err = norm(chi1+chi2-chi);
