function z = plus(x,y)
%PLUS - Add two Fermi operators
%
%   z = PLUS(x,y)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(x,'fermiop') && isa(y,'fermiop'));
assert(all(x.orbs == y.orbs & x.pFrom == y.pFrom & x.pTo == y.pTo & x.pairflag == y.pairflag));

z = x;	% copy all parameters
z.data = x.data+y.data;
