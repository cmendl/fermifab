function z = mtimes(x,y)
%MTIMES - Product of two Bose operators
%
%   z = MTIMES(x,y)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if isa(x,'bosonop')
	if isnumeric(y)
		z = x;
		z.data = y*z.data;
	elseif isa(y,'bosonstate')
		assert(all(x.modes==get(y,'modes') & x.pFrom==get(y,'N')));
		z = bosonstate(x.modes,x.pTo,x.data*get(y,'data'));
	elseif isa(y,'bosonop')
		assert(all(x.modes==y.modes & x.pFrom==y.pTo));
		z = bosonop(x.modes,y.pFrom,x.pTo,x.data*y.data,x.modenames);
	else
		error('Invalid input arguments.');
	end
elseif isnumeric(x)	% cannot use 'switch' here
	assert(isa(y,'bosonop'));
	z = y;
	z.data = x*y.data;
else
	error('Invalid input arguments.');
end
