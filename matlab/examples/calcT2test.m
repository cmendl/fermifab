function err = calcT2test(orbs,N)
%
%   err = calcT2test(orbs,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

% generate random fermi state
psi = fermistate(orbs,N);
x = crand(size(psi)); x = x/norm(x);
psi = set(psi,'data',x);

% calculate T2 operator
T = calcT2(rdm(psi,2),rdm(psi,1));

%%
% T2 operator based solely on two-body RDM
Td = calcT2(rdm(psi,2),get(psi,'N'));
% store error
err = norm(T-Td);

%%
% T operator directly from N-particle state psi
TN = calcT2N(psi);
% accumulate error
err = err + norm(T-TN);
