function [err,A] = kraus_test(orbs,p,N,varargin)
%KRAUS_TEST - Kraus operator test
%
%   [err,A] = KRAUS_TEST(orbs,p,N,...)
%
%   Verify that the p-body RDM of psi equals
%   sum_j a_{chi_j} |psi><psi| a^dagger_{chi_j}
%   where (chi_j) is an orthonormal basis of wedge^{N-p} H
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if ~isempty(varargin) && strcmp(varargin{1},'pair')
	pairflag = 1;
else
	pairflag = 0;
end

if ~pairflag
	psi = fermistate(orbs,N,crand(nchoosek(orbs,N),1));

	G = fermiop(orbs,p,p);
	chi = fermistate(orbs,N-p);
	A = cell(length(chi),1);	% Kraus operators
	for j=1:length(chi)
		x = zeros(length(chi),1); x(j) = 1;
		chi = set(chi,'data',x);
		A{j} = annihil(chi,N);
		phi = A{j}*psi;
		G = G + phi*phi';
		A{j} = sparse(get(A{j},'data'));
	end

	err = norm(G - rdm(psi,p));
else
	assert(mod(orbs,2)==0 && mod(p,2)==0 && mod(N,2)==0);

	psi = fermistate(orbs,N,crand(nchoosek(orbs/2,N/2),1),'pair');

	G = fermiop(orbs,p,p,'pair');
	chi = fermistate(orbs,N-p,'pair');
	A = cell(length(chi),1);	% Kraus operators
	for j=1:length(chi)
		x = zeros(length(chi),1); x(j) = 1;
		chi = set(chi,'data',x);
		A{j} = annihil(chi,N);
		phi = A{j}*psi;
		G = G + phi*phi';
		A{j} = sparse(get(A{j},'data'));
	end

	err = norm(G - rdm(psi,p));
end
