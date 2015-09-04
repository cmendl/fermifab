function z = plus(x,y)
%PLUS - Add Fermi states
%
%   z = PLUS(x,y)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(x,'fermistate') && isa(y,'fermistate'));
assert(sum(x.orbs)==sum(y.orbs) && sum(x.N)==sum(y.N) && x.pairflag==y.pairflag);

% comprise Hilbert space dimensions
[m1,m2,orbs,N] = comprise_config(x.orbs,y.orbs,x.N,y.N);

z = fermistate(orbs,N,gen_comprise(x.orbs,x.N,m1)*x.data+gen_comprise(y.orbs,y.N,m2)*y.data,x.orbnames);
