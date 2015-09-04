function err = calcQtest(orbs,N)
%
%   err = calcQtest(orbs,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

% generate random fermi state
psi = fermistate(orbs,N);
x = crand(size(psi)); x = x/norm(x);
psi = set(psi,'data',x);

% calculate Q operator
Q = calcQ(rdm(psi,2),rdm(psi,1));

%%
% Q operator based solely on two-body RDM
Q2 = calcQ(rdm(psi,2),get(psi,'N'));
% accumulate error
err = norm(Q-Q2);

%%
% Q operator directly from N-particle state
QN = calcQN(psi);
% accumulate error
err = err + norm(Q-QN);

%%
% generalized Q operator for p = 1 should be id - gamma^1(psi)
g1 = rdm(psi,1);	% one-body RDM
id = set(g1,'data',eye(size(g1)));
% accumulate error
err = err + norm(calcQN(psi,1) - (id-g1));
