function H = p2N(h,N,varargin)
%P2N - Calculate N-body from p-body operator,
%      H = sum_ij h_ij a^dagger_i a_j
%
%   H = P2N(h,N,options)
%
%   Input:      h           p-body operator
%               N           target particle number
%               options     options in the form 'param',value
%                               'cache' <true|false> cache kernel matrices
%                               (possibly faster, but uses more memory)
%
%   Output:     H           N-body operator generated from h
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(h,'fermiop'));
assert(length(get(h,'pFrom'))==1 && length(get(h,'pTo'))==1,'Only simple configurations supported.');
assert(get(h,'pFrom')==get(h,'pTo'),'Input operator ''h'' must be symmetric.');

% process options
options = processopts(struct('cache',true),varargin{:});

% store kernel matrices in a persistent
% variable to save computation time
persistent cacheK;
% unique kernel denominator
d = [sprintf('orbs%u',get(h,'orbs')) '_' sprintf('p%u',get(h,'pFrom')) '_' sprintf('N%u',N)];
% check whether requested kernel is cached already
if ~isfield(cacheK,d)
	% kernel matrices for current configuration
	K = gen_rdm(get(h,'orbs'),get(h,'pFrom'),N);
	% cache kernel matrices
	if options.cache, cacheK.(d) = K; end
else
	K = cacheK.(d);
end

A = get(h,'data');
assert(~nnz(size(A)-size(K)));
H = zeros(size(K{1,1}));
if ~isa(A,'sym')
	% try sparse matrices
	if issparse(A), H = sparse(H); end
	for j=1:size(K,1)
		for k=1:size(K,2)
			H = H + A(j,k)*K{k,j};
		end
	end
else
	% some optimizations for symbolic variables
	H = sym(H);
	for j=1:size(K,1)
		for k=1:size(K,2)
			% K{k,j} is a sparse matrix
			H(K{k,j} ==  1) = H(K{k,j} ==  1) + A(j,k);
			H(K{k,j} == -1) = H(K{k,j} == -1) - A(j,k);
		end
	end
	H = simplify(H);
end

H = fermiop(get(h,'orbs'),N,N,H,get(h,'orbnames'));
