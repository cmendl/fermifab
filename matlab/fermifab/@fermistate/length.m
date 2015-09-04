function n = length(psi)
%LENGTH - Fermi state internal data length
%
%   n = LENGTH(psi)
%
%   Same as
%       nchoosekprod(psi.orbs,psi.N)		if pairflag is 0
%       nchoosekprod(psi.orbs/2,psi.N/2)	if pairflag is 1
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

n = length(psi.data);
