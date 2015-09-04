function z = kron(x,y)
%KRON - Kronecker tensor product of two Fermi operators
%
%   z = KRON(x,y)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(x,'fermiop') && isa(y,'fermiop'));

z = fermiop([x.orbs,y.orbs],[x.pFrom,y.pFrom],[x.pTo,y.pTo],kron(x.data,y.data),[x.orbnames,y.orbnames]);
