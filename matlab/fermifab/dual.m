function chi = dual(psi)
%DUAL - Antilinear dual operator: chi = a_psi |12...orbs>
%
%   chi = DUAL(psi)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

chi = annihil(psi,get(psi,'orbs'))*fermistate(get(psi,'orbs'),get(psi,'orbs'),1,get(psi,'orbnames'));
