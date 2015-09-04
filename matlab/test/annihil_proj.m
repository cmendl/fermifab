function err = annihil_proj(orbs,p,N)
%ANNIHIL_PROJ - sum_j a_{phi_j} |Psi><Psi| a^dagger_{phi_j} should be independent
%               of the particular choice of the p-body basis phi_j for any N-body |Psi>
%
%   err = ANNIHIL_PROJ(orbs,p,N)
%
%   Copyright (c) 2008, Christian B. Mendl
%   All rights reserved.

psi = fermistate(orbs,N,crand(nchoosek(orbs,N),1));

% standard orthonormal basis
B1 = eye(nchoosek(orbs,p));
P1 = fermiop(orbs,N-p,N-p);
for j=1:size(B1,2)
	phi = annihil(fermistate(orbs,p,B1(:,j)),N)*psi;
	P1 = P1 + phi*phi';
end

% random orthonormal basis
[B2,R] = qr(crand(nchoosek(orbs,p)));	% random unitary matrix
P2 = fermiop(orbs,N-p,N-p);
for j=1:size(B2,2)
	phi = annihil(fermistate(orbs,p,B2(:,j)),N)*psi;
	P2 = P2 + phi*phi';
end

err = norm(P1-P2);
