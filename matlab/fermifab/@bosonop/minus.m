function z = minus(x,y)
%MINUS - Subtract two Bose operators
%
%   z = MINUS(x,y)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(x,'bosonop') && isa(y,'bosonop'));
assert(x.modes == y.modes && x.pFrom == y.pFrom && x.pTo == y.pTo);

z = x;	% copy all parameters
z.data = x.data-y.data;
