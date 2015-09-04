function z = kron(x,y)
%KRON - Kronecker tensor product of two Fermi states
%
%   z = KRON(x,y)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(x,'fermistate') && isa(y,'fermistate'));

z = fermistate([x.orbs,y.orbs],[x.N,y.N],kron(x.data,y.data),[x.orbnames,y.orbnames]);
