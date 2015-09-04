function err = calcT1test(orbs,N)
%
%   err = calcT1test(orbs,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

% generate random fermi state
psi = fermistate(orbs,N);
x = crand(size(psi)); x = x/norm(x);
psi = set(psi,'data',x);

% calculate T1 operator
T = calcT1(rdm(psi,2),rdm(psi,1));

%%
% T1 operator based solely on two-body RDM
Td = calcT1(rdm(psi,2),get(psi,'N'));
% store error
err = norm(T-Td);

%%
% T operator directly from N-particle state psi
TN = calcT1N(psi);
% accumulate error
err = err + norm(T-TN);

%%
% generalized T1 operator for p = 2 should be gamma^2(psi) + Q
g1 = rdm(psi,1);	% one-body RDM
g2 = rdm(psi,2);	% two-body RDM
Q = calcQ(g2,g1);
% accumulate error
err = err + norm(calcT1N(psi,2) - (g2+Q));
