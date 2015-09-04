function y = ctranspose(x)
%CTRANSPOSE - Complex conjugate transpose of a Fermi state, resulting in a Fermi operator
%
%   y = CTRANSPOSE(x)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if ~x.pairflag
	y = fermiop(x.orbs,x.N,zeros(size(x.orbs)),x.data',x.orbnames);
else
	y = fermiop(x.orbs,x.N,zeros(size(x.orbs)),x.data',x.orbnames,'pair');
end
