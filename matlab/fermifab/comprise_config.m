function [m1,m2,orbs,N] = comprise_config(orbs1,orbs2,N1,N2)
%COMPRISE_CONFIG - Comprise configurations
%
%   [m1,m2,orbs,N] = COMPRISE_CONFIG(orbs1,orbs2,N1,N2)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(sum(orbs1)==sum(orbs2) && sum(N1)==sum(N2),'Total number of orbitals and particles must be equal.');
assert(all(orbs1>0) && all(orbs2>0),'Number of orbitals must be positive.');

% cumulative number of orbitals
cumorbs1 = cumsum(orbs1);
cumorbs2 = cumsum(orbs2);
% cumulative particle numbers
n1 = cumsum(N1);
n2 = cumsum(N2);

% merge tables
m1 = true(size(cumorbs1));
m2 = true(size(cumorbs2));
orbs = orbs1(1); N = N1(1);
j1 = 1; j2 = 1;
while j1<length(cumorbs1)
	if cumorbs1(j1)<cumorbs2(j2)
		orbs(end) = orbs(end)+orbs1(j1+1);
		N(end) = N(end)+N1(j1+1);
		j1 = j1+1;
	elseif cumorbs2(j2)<cumorbs1(j1)
		j2 = j2+1;
	else
		if n1(j1)==n2(j2)
			% both mode and particle numbers match
			m1(j1) = false;
			m2(j2) = false;
			orbs(end+1) = orbs1(j1+1);
			N(end+1) = N1(j1+1);
		else
			orbs(end) = orbs(end)+orbs1(j1+1);
			N(end) = N(end)+N1(j1+1);
		end
		j1 = j1+1;
		j2 = j2+1;
	end
end
assert(sum(orbs)==cumorbs1(end) && sum(N)==n1(end));
m1 = m1(1:end-1);
m2 = m2(1:end-1);
