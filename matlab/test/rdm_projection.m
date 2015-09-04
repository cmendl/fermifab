function psi = rdm_projection(orbs,p,N)
%RDM_PROJECTION - Find a wavefunction psi such that rdm(psi,p) is an orthogonal projection
%
%   psi = RDM_PROJECTION(orbs,p,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

fprintf('Trying to find a %i-body wavefunction psi such that rdm(psi,%i) is an orthogonal projection...\n',N,p);
options.MaxIter = 1e4; options.MaxFunEvals = 1e4;
[x,m] = fminsearch(@(x)f(x,orbs,p,N),crand(nchoosek(orbs,N),1),options);
psi = fermistate(orbs,N,x/norm(x));
fprintf('Deviation from projection: %g\n',m);

% eigenvalues of 1-body density matrix
disp('Eigenvalues of 1-body density matrix:');
disp(eig(get(rdm(psi,1),'data')));


%%
function y = f(x,orbs,p,N)

G = rdm(fermistate(orbs,N,x/norm(x)),p);
y = norm(G^2-G,'fro');
