function C = create(chi,N)
%CREATE - Fermion creation operator
%
%   C = CREATE(chi,N)
%
%   Calculate the fermion creation operator
%   a^dagger_chi: wedge^N H -> wedge^(N-p) H
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

C = annihil(chi,N)';
