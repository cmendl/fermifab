function [psi,m] = maxlambda(orbs,p,N)
%MAXLAMBDA - Maximize the spectral radius rho(rdm(psi)) over all wavefunctions psi
%
%   [psi,m] = MAXLAMBDA(orbs,p,N)
%
%   Maximize the spectral radius rho(rdm(psi)) over all psi,
%   should be 0.5*N*(orbs-N+2)/orbs if p = 2 and N, orbs even.
%   see C. Yang:  Concept of Off-Diagonal Long-Range Order and
%   the Quantum Phases of Liquid He and of Superconductors,
%   Rev. Mod. Phys. 34, 694-704 (1962)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

fprintf('Trying to maximize the spectral radius of the %i-body reduced density matrix for N == %i and orbs == %i...\n',p,N,orbs);
options.MaxIter = 1e4; options.MaxFunEvals = 1e4;
[x,m] = fminsearch(@(x)calc_lambda(x,orbs,p,N),crand(nchoosek(orbs,N),1),options);
psi = fermistate(orbs,N,x/norm(x));
m = -m;
fprintf('Maximal spectral radius estimate: %g\n',m);
if p==2 && ~mod(orbs,2) && ~mod(N,2)
	fprintf('Absolute difference to analytic value: %g\n',abs(m-0.5*N*(orbs-N+2)/orbs));
end


%% maximal eigenvalue of rdm(psi)
function y = calc_lambda(x,orbs,p,N)

options.disp = 0;
% use 'real' to remove numerical artefacts
y = -real(eigs(get(rdm(fermistate(orbs,N,x/norm(x)),p),'data'),1,'lm',options));
