function G = calcGN(psi,p)
%
%   G = calcGN(psi)
%   G = calcGN(psi,p)
%
%   Calculate the generalized G operator of the G-representability condition
%   directly from the N-particle state psi; i.e., the matrix
%   <psi| a^dagger_j a_i a^dagger_k a_ell |psi>
%
%   Reference: "The reduced density matrix method for electronic structure
%   calculations and the role of three-index representability conditions"
%   Journal of Chemical Physics 120, 2095 (2004); doi:10.1063/1.1636721
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

% default p value
if nargin<2, p = 1; end

assert(isa(psi,'fermistate'));
assert(p < get(psi,'N'));

% kernel matrices
K = gen_rdm(get(psi,'orbs'),p,get(psi,'N'));

% raw input data
x = get(psi,'data');

% raw output data
G = zeros(length(K)*[1,1,1,1]);

for i=1:length(K)
	for j=1:length(K)
		for k=1:length(K)
			for l=1:length(K)
				% K{l,k} equals a^dagger_k a_l
				G(i,j,k,l) = (K{j,i}*x)'*(K{l,k}*x);
			end
		end
	end
end

G = reshape(G,length(K)^2*[1,1]);

% embed G into a 'fermiop'
G = fermiop(get(psi,'orbs')*[1,1],[p,p],[p,p],G);
