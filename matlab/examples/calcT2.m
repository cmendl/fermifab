function T = calcT2(g2,g1)
%
%   T = calcT2(g2,g1)
%   T = calcT2(g2,N)
%
%   Calculate the T2 operator of the T2-representability condition
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

%%

% assemble matrix T with entries (ij,kl) equal to Gamma(|ij>,|kl>)
T = zeros(orbs^2);

X = fermi2coords(orbs,2);
% entries ij with i < j
IXs = X(1,:) + (X(2,:)-1)*orbs;
% entries ij with i > j -> inverse sign
IXg = X(2,:) + (X(1,:)-1)*orbs;

T(IXs,IXs) =  g2;
T(IXg,IXg) =  g2;
T(IXs,IXg) = -g2;
T(IXg,IXs) = -g2;

% tensor identity matrix
T = kron(T,eye(orbs));

% reorder indices
T = reshape(T,orbs*[1,1,1,1,1,1]);
T = permute(T,[2,1,6,5,4,3]);
T = reshape(T,orbs*[1,orbs,1,orbs]);

% collect antisymmetrized entries
T = T(:,IXs,:,IXs) + T(:,IXg,:,IXg) - T(:,IXs,:,IXg) - T(:,IXg,:,IXs);

T = -reshape(T,orbs*length(g2)*[1,1]);

%%

T = T + kron(eye(size(g2)),g1) + kron(g2.',eye(size(g1)));

%%
% embed T2 into a 'fermiop'
T = fermiop([orbs,orbs],[1,2],[1,2],T);
