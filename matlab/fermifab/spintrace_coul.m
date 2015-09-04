function G = spintrace_coul(state,varargin)
%SPINTRACE_COUL - Trace out spin coordinates of a two-body Fermi state or operator
%                 for mapping it to Coulomb integrals (ab|cd)
%
%   G = SPINTRACE_COUL(state)
%   G = SPINTRACE_COUL(state,orbnames)
%   G = SPINTRACE_COUL(...,options)
%
%   Trace out spin coordinates, i.e., calculate the spatial-bosonic
%   reduced density matrix of a two-body Fermi state which maps to
%   Coulomb integrals (ab|cd).
%   Assuming alternating spin up and spin down orbitals.
%
%   Input:      state       two-body 'fermistate' or 'fermiop' fermionic state
%               orbnames    cell array of target orbital names
%               options     options in the form 'param',value
%                               'cache' <true|false> cache kernel matrices
%                               (possibly faster, but uses more memory)
%
%   Output:     G           spatial reduced density matrix
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

% check number of orbitals
assert(length(get(state,'orbs'))==1,'Only simple configuration supported.');
assert(mod(get(state,'orbs'),2)==0,'Number of orbitals must be even due to spin structure.');

% simple 2-body configuration allowed only
if isa(state,'fermistate')
	assert(get(state,'N')==2,'Only 2-body configuration allowed.');
else
	assert(get(state,'pFrom')==2 && get(state,'pTo')==2,'Only 2-body configuration allowed.');
end

% orbital names
orbnames = {};
if ~isempty(varargin) && iscell(varargin{1})
	orbnames = varargin{1};
	varargin = varargin(2:end);
end

% process options
options = processopts(struct('cache',true),varargin{:});

% store kernel matrices in a persistent
% variable to save computation time
persistent cacheK;
% unique kernel denominator
d = sprintf('orbs%u',get(state,'orbs')/2);
% check whether requested kernel is cached already
if ~isfield(cacheK,d)
	% kernel matrices for current configuration
	K = gen_spintrace_coul(get(state,'orbs')/2);
	% cache kernel matrices
	if options.cache, cacheK.(d) = K; end
else
	K = cacheK.(d);
end

x = get(state,'data');
G = zeros(size(K));
if ~isa(x,'sym')
	% try sparse matrix
	if issparse(x), G = sparse(G); end
else
	% symbolic variables - 'sparse' doesn't work
	G = sym(G);
end
if isa(state,'fermistate')
	for j=1:size(K,1)
		for k=1:size(K,2)
			G(j,k) = x'*K{j,k}*x;
		end
	end
else
	% 'state' is a 'fermiop'
	for j=1:size(K,1)
		for k=1:size(K,2)
			G(j,k) = trace_prod(K{j,k},x);
		end
	end
end

if isa(G,'sym'), G = simplify(G); end

% bosonic operator on spatial orbitals
if isempty(orbnames)
	G = bosonop(get(state,'orbs')/2,2,2,G);
else
	G = bosonop(get(state,'orbs')/2,2,2,G,orbnames);
end
