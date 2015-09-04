function err = annihil_test(orbs,p,N)
%ANNIHIL_TEST - Annihilation operator test file
%
%   err = ANNIHIL_TEST(orbs,p,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

err = 0;

% kernel matrices for current configuration
K = gen_rdm(orbs,p,N);

for j=1:size(K,1)
	x = zeros(size(K,1),1); x(j) = 1;
	chi_j = fermistate(orbs,p,x);
	for k=1:size(K,2)
		x = zeros(size(K,2),1); x(k) = 1;
		chi_k = fermistate(orbs,p,x);
		err = err + norm(K{j,k} - get(create(chi_k,N)*annihil(chi_j,N),'data'));
	end
end
