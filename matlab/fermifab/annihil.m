function A = annihil(chi,N)
%ANNIHIL - Fermion annihilation operator
%
%   A = ANNIHIL(chi,N)
%
%   Calculate the fermion annihilation operator
%   a_chi: wedge^N H -> wedge^(N-p) H
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(chi,'fermistate'));
assert(N > get(chi,'N'));

% annihilation matrix kernels
persistent Kann;

if ~isfield(Kann,'data')
	Kann.data = [];
end
d = ['d' num2str(get(chi,'orbs')) '_' num2str(get(chi,'N')) '_' num2str(N)];
if ~get(chi,'pairflag')
	if ~isfield(Kann.data,d)
		Kann.data.(d) = gen_annihil(get(chi,'orbs'),get(chi,'N'),N);
	end
else
	assert(mod(N,2)==0);	% N must be even
	d = [d '_pair'];
	if ~isfield(Kann.data,d)
		Kann.data.(d) = gen_annihil_pairs(get(chi,'orbs')/2,get(chi,'N')/2,N/2);
	end
end
% annihilation matrices for current configuration
K = Kann.data.(d);

v = get(chi,'data');
for j=1:size(K,1)
	for k=1:size(K,2)
		A(j,k) = v'*K{j,k};
	end
end
% useful for symbolic computations
if isa(A,'sym'), A = simplify(A); end

if ~get(chi,'pairflag')
	A = fermiop(get(chi,'orbs'),N,N-get(chi,'N'),A);
else
	A = fermiop(get(chi,'orbs'),N,N-get(chi,'N'),A,'pair');
end
