function G = calcGref(g2,g1)
%
%   G = calcGref(g2,g1)
%
%   Calculate the G operator of the G-representability condition
%   from the one- and two-body RDMs g1 and g2, respectively
%
%   Reference: "The reduced density matrix method for electronic structure
%   calculations and the role of three-index representability conditions"
%   Journal of Chemical Physics 120, 2095 (2004); doi:10.1063/1.1636721
%
%   Reference implementation
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

orbs = get(g2,'orbs');

g1 = get(g1,'data');
g2 = get(g2,'data');

G = zeros(orbs*[1,1,1,1]);
for i=1:orbs
	for j=1:orbs
		for k=1:orbs
			for l=1:orbs
				[x,f1] = fermiIX(i,l,orbs);
				[y,f2] = fermiIX(j,k,orbs);
				G(i,j,k,l) = f1*f2*g2(x,y) + (i==k)*g1(l,j);
			end
		end
	end
end
G = reshape(G,orbs^2*[1,1]);

% embed G into a 'fermiop'
G = fermiop([orbs,orbs],[1,1],[1,1],G);


%%
function [x,fac] = fermiIX(i,j,K)

if i==j
	x = 1;
	fac = 0;
	return;
end

if i<j
	fac = 1;
else % i>j
	fac = -1;
	% swap i and j
	t = j;
	j = i;
	i = t;
end

% convert to [12,13,14,...,23,24,...] enumeration index
x = coords2fermi([i,j],K);
