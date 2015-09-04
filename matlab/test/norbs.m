function [psi,err] = norbs(orbs,N)
%NORBS - "natural orbitals": eigenstates of the 1-body RDM
%
%   [psi,err] = NORBS(orbs,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

% random wavefunction
psi = fermistate(orbs,N,crand(nchoosek(orbs,N),1));

% "natural orbitals"
[U,~] = eig(rdm(psi,1));
psi = tensor_op(U,N)'*psi;
% the 1-body RDM should now be diagonal
G = get(rdm(psi,1),'data');
err = norm(G-diag(diag(G)));
