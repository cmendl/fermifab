function err = calcGtest(orbs,N)
%
%   err = calcGtest(orbs,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

% generate random fermi state
psi = fermistate(orbs,N);
x = crand(size(psi)); x = x/norm(x);
psi = set(psi,'data',x);

% calculate G operator
G = calcG(rdm(psi,2),rdm(psi,1));

%%
% reference implementation
Gref = calcGref(rdm(psi,2),rdm(psi,1));
% store error
err = norm(G-Gref);

%%
% G operator based solely on two-body RDM
G2 = calcG(rdm(psi,2),get(psi,'N'));
% accumulate error
err = err + norm(G-G2);

%%
% G operator directly from N-particle state
GN = calcGN(psi);
% accumulate error
err = err + norm(G-GN);
