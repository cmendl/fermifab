function z = mtimes(x,y)
%MTIMES - Multiply Fermi operators and states
%
%   z = MTIMES(x,y)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if isa(y,'fermistate')
	assert(length(x) == 1);
	z = y;
	z.data = x*y.data;
elseif isa(y,'fermiop')
	assert(x.pairflag == get(y,'pairflag') && size(y,1) == 1 && all(x.orbs==get(y,'orbs')));
	if ~x.pairflag
		z = fermiop(x.orbs,get(y,'pFrom'),x.N,x.data*get(y,'data'),x.orbnames);
	else
		z = fermiop(x.orbs,get(y,'pFrom'),x.N,x.data*get(y,'data'),x.orbnames,'pair');
	end
else
	error('Invalid input arguments.');
end
