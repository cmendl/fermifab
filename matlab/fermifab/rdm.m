function G = rdm(state,p)
%RDM - Calculate reduced density matrices
%
%   G = RDM(state,p,options)
%
%   Calculate the p-body reduced density matrix
%   of a N-body pure or mixed state
%
%   Input:      state       'fermistate' or 'fermiop' fermionic state
%               p           target particle number
%
%   Output:     G           reduced density matrix
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

% check input
assert(length(p)==1,'Input argument ''p'' must be a scalar.');
if get(state,'pairflag')
	assert(mod(p,2)==0,'Expecting even ''p'' for pair structure.');
end

% number of orbitals
orbs = get(state,'orbs');

% particle numbers
if isa(state,'fermistate')
	N1 = get(state,'N'); N2 = N1;
else
	N1 = get(state,'pFrom');
	N2 = get(state,'pTo');
	assert(sum(N1)==sum(N2),'Input and output particle numbers of the density matrix must agree.');
end

% effective configuration index
IX = N1<orbs | N2<orbs;

% raw input data
x = get(state,'data');

% raw output data
if ~get(state,'pairflag')
	G = zeros(nchoosek(sum(orbs),p));
	% generalized 'configuration generator'
	gen_config_switch = @gen_config;
else
	G = zeros(nchoosek(sum(orbs)/2,p/2));
	% pair structure
	gen_config_switch = @(orbs,N,maxC) 2*gen_config(orbs,N/2,maxC/2);
end
if ~isa(x,'sym')
	% try sparse matrix
	if issparse(x), G = sparse(G); end
else
	% 'sparse sym' not supported by Matlab
	G = sym(G);	
end

% for each configuration of the reduced particle space...
for pc1=gen_config_switch(length(N2),p,N2)
	% N1 and N2 can store different configurations, although sum(N1)==sum(N2)
	pc = [pc1';N1-N2+pc1'];	% N2-p1==N1-p2
	% no negative numbers in configuration of reduced particle space allowed
	if any(pc(2,:)<0), continue; end
	if sum(pc1(IX))>0
		% kernel matrices for current effective configuration
		if ~get(state,'pairflag')
			K = gen_rdm(orbs(IX),pc1(IX),N1(IX),N2(IX));
		else
			K = gen_rdm_pairs(orbs(IX)/2,pc1(IX)/2,N1(IX)/2,N2(IX)/2);
		end

		% actually calculate RDM for current configuration
		if ~isa(x,'sym')
			if issparse(x), Gc = sparse(size(K,1),size(K,2));
			else Gc = zeros(size(K)); end
			% switch between 'fermistate' and 'fermiop'
			if isa(state,'fermistate')
				for j=1:numel(K)
					Gc(j) = x'*K{j}*x;
				end
			else	% 'state' is a 'fermiop'
				for j=1:numel(K)
					Gc(j) = trace_prod(K{j},x);
				end
			end
		else
			Gc = sym(zeros(size(K)));
			% optimization for symbolic variables - 'sparse' doesn't work
			if isa(state,'fermistate')
				for j=1:numel(K)
					[r,c] = find(K{j} == 1);
					if ~isempty(r), Gc(j) = x(r)'*x(c); end
					[r,c] = find(K{j} == -1);
					if ~isempty(r), Gc(j) = Gc(j)-x(r)'*x(c); end
				end
			else
				% 'state' is a 'fermiop'
				for j=1:numel(K)
					[r,c] = find(K{j} == 1);
					if ~isempty(r), Gc(j) = trace(x(c,r)); end
					[r,c] = find(K{j} == -1);
					if ~isempty(r), Gc(j) = Gc(j)-trace(x(c,r)); end
				end
			end
			Gc = simplify(Gc);
		end
	else
		if isa(state,'fermistate')
			if ~isa(x,'sym'), Gc = norm(x,2)^2;
			else Gc = sum(x.*conj(x)); end
		else
			Gc = trace(x);
		end
	end

	% take into account Hilbert spaces completely filled by the input state
	if ~get(state,'pairflag')
		dim = nchoosekprod(orbs(~IX),pc1(~IX));
	else
		dim = nchoosekprod(orbs(~IX)/2,pc1(~IX)/2);
	end
	if issparse(x)
		Gc = kron(speye(dim),Gc);
	else
		Gc = kron(eye(dim),Gc);
	end
	% rearrange to original order
	dim = zeros(2,length(orbs));	% dimensions of reduced particle spaces
	if ~get(state,'pairflag')
		for j=1:2
			for k=1:length(orbs), dim(j,k) = nchoosek(orbs(k),pc(j,k)); end
		end
	else
		for j=1:2
			for k=1:length(orbs), dim(j,k) = nchoosek(orbs(k)/2,pc(j,k)/2); end
		end
	end
	dim = [dim(:,~IX),dim(:,IX)];
	% completely filled spaces at the beginning
	order = 1:length(orbs);
	order = [order(~IX),order(IX)];
	Gc = ipermOp(Gc,order,dim(1,:),dim(2,:));

	% populate output state
	X = cell(2,1);
	if ~get(state,'pairflag')
		for j=1:2, X{j} = gen_comprise(orbs,pc(j,:),true(length(orbs)-1,1)); end
	else
		for j=1:2, X{j} = gen_comprise(orbs/2,pc(j,:)/2,true(length(orbs)-1,1)); end
	end
	G = G + X{1}*Gc*X{2}';
end

% put raw data inside a 'fermiop'
if ~get(state,'pairflag')
	G = fermiop(sum(orbs),p,p,G,get(state,'orbnames'));
else
	G = fermiop(sum(orbs),p,p,G,get(state,'orbnames'),'pair');
end
