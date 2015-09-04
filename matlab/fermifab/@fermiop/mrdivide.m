function z = mrdivide(x,y)
%MRDIVIDE - Matrix right division for Fermi operators
%
%   z = MRDIVIDE(x,y)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(x,'fermiop'));

if isnumeric(y)
	z = x;
	z.data = z.data/y;
elseif isa(y,'fermiop')
	assert(all(x.orbs == y.orbs & x.pFrom == y.pFrom & x.pTo == y.pTo));
	z = x;
	z.data = z.data/y.data;
else
	error('Invalid input arguments.');
end
