function T = calcT1N(psi,p)
%
%   T = calcT1N(psi)
%   T = calcT1N(psi,p)
%
%   Calculate the generalized T1 operator of the T1-representability condition
%   directly from the N-particle state psi; i.e., the matrix
%   <psi| A'*A + A*A' |psi> with A_ijk = a_|ijk>
%
%   Reference: "The reduced density matrix method for electronic structure
%   calculations and the role of three-index representability conditions"
%   Journal of Chemical Physics 120, 2095 (2004); doi:10.1063/1.1636721
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

% default p value
if nargin<2, p = 3; end

assert(isa(psi,'fermistate'));
assert(p < get(psi,'N'));

orbs = get(psi,'orbs');

% raw input data
x = get(psi,'data');

%%
% kernel matrices: annihilation operators
Ka = gen_annihil(orbs,p,get(psi,'N'));
assert(length(x)==size(Ka,2));

A = zeros(size(Ka,1),length(Ka{1,1}));
for i=1:size(Ka,1)
	for j=1:size(Ka,2)
		A(i,:) = A(i,:) + x(j)*Ka{i,j}';
	end
end

%%
% kernel matrices: creation operators (conjugate transpose of annihilation operators)
Kc = gen_annihil(orbs,p,get(psi,'N')+p);
assert(length(x)==size(Kc,1));

C = zeros(size(Kc,2),length(Kc{1,1}));
for i=1:size(Kc,2)
	for j=1:size(Kc,1)
		% conjugate transpose of annihilation: switch i and j in Kc{i,j}
		C(i,:) = C(i,:) + x(j)*Kc{j,i}';
	end
end

%%
% positive semidefinite T1 operator
T = (A'*A).' + C'*C;

% embed T1 into a 'fermiop'
T = fermiop(orbs,p,p,T);
