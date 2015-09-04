function Q = calcQN(psi,p)
%
%   Q = calcQN(psi)
%   Q = calcQN(psi,p)
%
%   Calculate the generalized Q operator of the Q-representability condition
%   directly from the N-particle state psi; i.e., the matrix
%   <psi| a_|ij> a^dagger_|kl> |psi>
%
%   Reference: "The reduced density matrix method for electronic structure
%   calculations and the role of three-index representability conditions"
%   Journal of Chemical Physics 120, 2095 (2004); doi:10.1063/1.1636721
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

% default p value
if nargin<2, p = 2; end

assert(isa(psi,'fermistate'));
assert(p < get(psi,'N'));

orbs = get(psi,'orbs');

% kernel matrices
K = gen_annihil(orbs,p,get(psi,'N')+p);

x = get(psi,'data');
assert(length(x)==size(K,1));

% creation operators (conjugate transpose of annihilation operators)
C = zeros(size(K,2),length(K{1,1}));
for i=1:size(K,2)
	for j=1:size(K,1)
		% conjugate transpose of annihilation: switch i and j in K{i,j}
		C(i,:) = C(i,:) + x(j)*K{j,i}';
	end
end

% Q operator
Q = C'*C;

% embed Q into a 'fermiop'
Q = fermiop(orbs,p,p,Q);
