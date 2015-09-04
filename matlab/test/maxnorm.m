function [psi,m] = maxnorm(orbs,p,N)
%MAXNORM - Try to maximize trace(rdm(psi)'*rdm(psi)) over all psi.
%          Seems that the maximizers are exactly the Slater determinants.
%
%   [psi,m] = MAXNORM(orbs,p,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

[x,m] = fminsearch(@(x)calc_norm(x,orbs,p,N),crand(nchoosek(orbs,N),1),optimset('MaxFunEvals',5000));
psi = fermistate(orbs,N,x/norm(x));
m = -m;


%%
function y = calc_norm(x,orbs,p,N)

y = -norm(rdm(fermistate(orbs,N,x/norm(x)),p),'fro')^2;
