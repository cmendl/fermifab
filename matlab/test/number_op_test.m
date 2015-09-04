function err = number_op_test(orbs,N)
%NUMBER_OP_TEST - Test pair structure commutator relations
%
%   err = NUMBER_OP_TEST(orbs,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if mod(orbs,2)~=0
	error('''orbs'' must be even.');
end

disp('Testing pair structure commutator relations ...');

% random Hermitian matrix
A = crand(orbs/2);
A = A+A';
% set coefficients of h corresponding to the |2j-1,2j> basis vectors
h = zeros(nchoosek(orbs,2));
for j=1:orbs/2
	% index of |2j-1,2j>
	jInd = coords2fermi([2*j-1,2*j],orbs);
	for k=1:orbs/2
		kInd = coords2fermi([2*k-1,2*k],orbs);
		h(jInd,kInd) = A(j,k);
	end
end
h = fermiop(orbs,2,2,h);
% N-body Hamiltonian
H = p2N(h,N);

m = cell(orbs/2,1);
for j=1:orbs/2
	% n_{2j-1}
	n2j1 = zeros(orbs); n2j1(2*j-1,2*j-1) = 1;
	n2j1 = p2N(fermiop(orbs,1,1,n2j1),N);
	% n_{2j}
	n2j = zeros(orbs); n2j(2*j,2*j) = 1;
	n2j = p2N(fermiop(orbs,1,1,n2j),N);
	% identity
	id = fermiop(orbs,N,N,eye(nchoosek(orbs,N)));
	% special number operator
	% m{j} = K{2*j-1,2*j-1}*(eye(size(K{1,1}))-K{2*j,2*j}) + (eye(size(K{1,1}))-K{2*j-1,2*j-1})*K{2*j,2*j};
	m{j} = n2j1*(id-n2j) + (id-n2j1)*n2j;
end

err = 0;
for j=1:orbs/2
	err = err + norm(com(m{j},H));	% commutator should be zero
	for k=j:orbs/2
		err = err + norm(com(m{j},m{k})) + norm(com(m{j}*m{k},H));
	end
end
fprintf('error: %g\n',err);
