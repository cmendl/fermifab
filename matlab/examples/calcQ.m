function Q = calcQ(g2,g1)
%
%   Q = calcQ(g2,g1)
%   Q = calcQ(g2,N)
%
%   Calculate the Q operator of the Q-representability condition
%   from the one- and two-body RDMs g1 and g2, respectively
%
%   Reference: "The reduced density matrix method for electronic structure
%   calculations and the role of three-index representability conditions"
%   Journal of Chemical Physics 120, 2095 (2004); doi:10.1063/1.1636721
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if isnumeric(g1)
	% interpret argument 'g1' as the particle number N
	N = g1;
	% calculate the one-body RDM from the two-body RDM
	g1 = rdm(g2,1)/(N-1);
end

id = set(g2,'data',eye(size(g2)));

Q = g2 - p2N(g1,2) + id;
