function T = calcT2N(psi)
%
%   T = calcT2N(psi)
%
%   Calculate the T2 operator of the T2-representability condition
%   directly from the N-particle state psi; i.e., the matrix
%   <psi| A'*A + A*A' |psi> with A_ijk = a^\dagger_i a_|jk>
%
%   Reference: "The reduced density matrix method for electronic structure
%   calculations and the role of three-index representability conditions"
%   Journal of Chemical Physics 120, 2095 (2004); doi:10.1063/1.1636721
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(psi,'fermistate'));

orbs = get(psi,'orbs');

% raw input data
x = get(psi,'data');

%%
% annihilation operators N -> N-2
Ka = gen_annihil(orbs,2,get(psi,'N'));
assert(length(x)==size(Ka,2));

% creation operators N-2 -> N-1 (conjugate transpose of annihilation operators)
Kc = gen_annihil(orbs,1,get(psi,'N')-1);
assert(size(Ka,1)==size(Kc,1));

A = zeros(size(Kc,2),length(Ka{1,1})*length(Kc{1,1}));
for i=1:size(Kc,2)
	for j=1:size(Ka,1)
		for k=1:size(Ka,2)
			% conjugate transpose of annihilation: switch i and j in Kc{i,j}
			A(i,:) = A(i,:) + kron(Ka{j,k}*x(k),Kc{j,i}).';
		end
	end
end

%%
% annihilation operators N -> N-1
Ka = gen_annihil(orbs,1,get(psi,'N'));
assert(length(x)==size(Ka,2));

% creation operators N-1 -> N+1 (conjugate transpose of annihilation operators)
Kc = gen_annihil(orbs,2,get(psi,'N')+1);
assert(size(Ka,1)==size(Kc,1));

C = zeros(size(Kc,2),length(Ka{1,1})*length(Kc{1,1}));
for i=1:size(Kc,2)
	for j=1:size(Ka,1)
		for k=1:size(Ka,2)
			% conjugate transpose of annihilation: switch i and j in Kc{i,j}
			C(i,:) = C(i,:) + kron(Kc{j,i},Ka{j,k}*x(k)).';
		end
	end
end

%%
% positive semidefinite T2 operator
T = A'*A + (C'*C).';

% embed T2 into a 'fermiop'
T = fermiop([orbs,orbs],[1,2],[1,2],T);
