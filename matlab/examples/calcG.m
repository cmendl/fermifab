function G = calcG(g2,g1)
%
%   G = calcG(g2,g1)
%   G = calcG(g2,N)
%
%   Calculate the G operator of the G-representability condition
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

orbs = get(g2,'orbs');

g1 = get(g1,'data');
g2 = get(g2,'data');

% assemble matrix G with entries (ij,kl) equal to Gamma(|i,j>,|k,l>)
G = zeros(orbs^2);

X = fermi2coords(orbs,2);
% entries ij with i < j
IXs = X(1,:) + (X(2,:)-1)*orbs;
% entries ij with i > j -> inverse sign
IXg = X(2,:) + (X(1,:)-1)*orbs;

G(IXs,IXs) =  g2;
G(IXg,IXg) =  g2;
G(IXs,IXg) = -g2;
G(IXg,IXs) = -g2;

% reorder indices
G = reshape(G,orbs*[1,1,1,1]);
G = permute(G,[1,3,4,2]);
G = reshape(G,orbs^2*[1,1]);

% add contribution of 1-body RDM
G = G + kron(g1.',eye(orbs));

% embed G into a 'fermiop'
G = fermiop([orbs,orbs],[1,1],[1,1],G);
