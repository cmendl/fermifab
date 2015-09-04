function z = mtimes(x,y)
%MTIMES - Product of a Fermi operator and Fermi state, or two Fermi operators
%
%   z = MTIMES(x,y)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if isa(x,'fermiop')
	if isnumeric(y)
		assert(size(y,1)==size(y,2),'Dimensions must not change.');
		z = x;
		z.data = y*z.data;
	elseif isa(y,'fermistate')
		assert(all(x.orbs==get(y,'orbs') & x.pFrom==get(y,'N')) && x.pairflag==get(y,'pairflag'));
		if ~x.pairflag
			z = fermistate(x.orbs,x.pTo,x.data*get(y,'data'));
		else
			z = fermistate(x.orbs,x.pTo,x.data*get(y,'data'),'pair');
		end
	elseif isa(y,'fermiop')
		assert(all(x.orbs==y.orbs & x.pFrom==y.pTo) && x.pairflag==y.pairflag);
		if ~x.pairflag
			z = fermiop(x.orbs,y.pFrom,x.pTo,x.data*y.data,x.orbnames);
		else
			z = fermiop(x.orbs,y.pFrom,x.pTo,x.data*y.data,x.orbnames,'pair');
		end
	else
		error('Invalid input arguments.');
	end
elseif isnumeric(x)	% cannot use 'switch' here
	assert(size(x,1)==size(x,2),'Dimensions must not change.');
	z = y;
	z.data = x*y.data;
else
	error('Invalid input arguments.');
end
