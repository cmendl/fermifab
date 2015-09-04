function [C,lambda] = com_rdm(varargin)
%COM_RDM - Test if some 1-body RDMs of the eigenstates of a 2-body RDM commute
%
%   [C,lambda] = COM_RDM(varargin)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

switch nargin
	case 0
		orbs = 6;
		N = 4;
		x = crand(nchoosek(orbs,N),1); x = x/norm(x);
		psi = fermistate(orbs,N,x);
	case 1
		psi = varargin{1};
		orbs = get(psi,'orbs');
	case 2
		orbs = varargin{1};
		N = varargin{2};
		x = crand(nchoosek(orbs,N),1); x = x/norm(x);
		psi = fermistate(orbs,N,x);
	otherwise
		error('Invalid number of input arguments.');
end

[U,T] = schur(get(rdm(psi,2),'data'));
lambda = real(diag(T));

g = cell(size(U,2),1);
for j=1:length(g)
	g{j} = rdm(fermistate(orbs,2,U(:,j)),1);
end

% A = zeros(length(g));
C = zeros(length(g));
for j=1:size(C,1)
	for k=1:size(C,2)
		% A(j,k) = norm(get(acom(g{j},g{k}),'data'));
		C(j,k) = norm(get(com(g{j},g{k}),'data'));
	end
end
