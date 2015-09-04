function psi = fermistate(varargin)
%FERMISTATE - Construct a Fermi state
%
%   psi = FERMISTATE()
%   psi = FERMISTATE(phi)
%   psi = FERMISTATE(orbs,N)
%   psi = FERMISTATE(orbs,N,data)
%   psi = FERMISTATE(orbs,N,data,orbnames)
%   psi = FERMISTATE(...,'pair')
%       assumes pair structure; effectively (orbs/2, N/2)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if ~isempty(varargin) && ischar(varargin{end})
	assert(strcmp(varargin{end},'pair'));
	varargin = varargin(1:end-1);
	psi.pairflag = 1;
else
	% set 'pair flag' to false (default)
	psi.pairflag = 0;
end

switch length(varargin)
	case 0
		% standard values
		psi.orbs = 6;
		psi.N = 4;
		% arbitrary state
		if ~psi.pairflag
			psi.data = zeros(nchoosekprod(psi.orbs,psi.N),1);
		else
			psi.data = zeros(nchoosekprod(psi.orbs/2,psi.N/2),1);
		end
		psi.data(1) = 1;
		% custom orbital names, change with 'set' command
		psi.orbnames = cell(1,psi.orbs);
		for j=1:psi.orbs, psi.orbnames{j} = num2str(j); end;
		% register class
		psi = class(psi,'fermistate');
		% check whether the symbolic math toolbox is installed
		if ~isempty(ver('symbolic'))
			superiorto('sym');
		end
	case 1
		% copy constructor
		assert(isa(varargin{1},'fermistate') && psi.pairflag==varargin{1}.pairflag);
		psi = varargin{1};
		% check whether the symbolic math toolbox is installed
		if ~isempty(ver('symbolic'))
			superiorto('sym');
		end
	case 2
		assert(isa(varargin{1},'numeric'));
		assert(isa(varargin{2},'numeric'));
		psi.orbs = round(varargin{1});
		psi.N = round(varargin{2});
		assert(~nnz(psi.N < 0 | psi.N > psi.orbs));
		% arbitrary state
		if ~psi.pairflag
			psi.data = zeros(nchoosekprod(psi.orbs,psi.N),1);
		else
			assert(all(mod(psi.orbs,2)==0 & mod(psi.N,2)==0));
			psi.data = zeros(nchoosekprod(psi.orbs/2,psi.N/2),1);
		end
		psi.data(1) = 1;
		% custom orbital names
		psi.orbnames = cell(1,sum(psi.orbs));
		for j=1:length(psi.orbnames), psi.orbnames{j} = num2str(j); end;
		% register class
		psi = class(psi,'fermistate');
		% check whether the symbolic math toolbox is installed
		if ~isempty(ver('symbolic'))
			superiorto('sym');
		end
	case 3
		% set 'orbs' and 'N'
		assert(isa(varargin{1},'numeric'));
		assert(isa(varargin{2},'numeric'));
		psi.orbs = round(varargin{1});
		psi.N = round(varargin{2});
		assert(all(psi.N >= 0 & psi.N <= psi.orbs));
		% copy Slater coefficients
		psi.data = reshape(varargin{3},[],1);
		assert(size(psi.data,2)==1);
		if isa(psi.data,'sym'), psi.data = simplify(psi.data); end
		if ~psi.pairflag
			assert(size(psi.data,1)==nchoosekprod(psi.orbs,psi.N),'Input data length does not equal prod_j nchoosek(orbs(j),N(j)).');
		else
			assert(mod(psi.orbs,2)==0 && mod(psi.N,2)==0 && size(psi.data,1)==nchoosekprod(psi.orbs/2,psi.N/2));
		end
		% custom orbital names
		psi.orbnames = cell(1,sum(psi.orbs));
		for j=1:length(psi.orbnames), psi.orbnames{j} = num2str(j); end;
		% register class
		psi = class(psi,'fermistate');
		% check whether the symbolic math toolbox is installed
		if ~isempty(ver('symbolic'))
			superiorto('sym');
		end
	case 4
		% set 'orbs' and 'N'
		assert(isa(varargin{1},'numeric'));
		assert(isa(varargin{2},'numeric'));
		psi.orbs = round(varargin{1});
		psi.N = round(varargin{2});
		assert(~nnz(psi.N < 0 | psi.N > psi.orbs));
		% copy Slater coefficients
		psi.data = reshape(varargin{3},[],1);
		assert(size(psi.data,2)==1);
		if isa(psi.data,'sym'), psi.data = simplify(psi.data); end
		if ~psi.pairflag
			assert(size(psi.data,1)==nchoosekprod(psi.orbs,psi.N),'Input data length does not equal prod_j nchoosek(orbs(j),N(j)).');
		else
			assert(mod(psi.orbs,2)==0 && mod(psi.N,2)==0 && size(psi.data,1)==nchoosekprod(psi.orbs/2,psi.N/2));
		end
		% orbital names
		assert(iscell(varargin{4}));
		psi.orbnames = reshape(varargin{4},1,[]);
		assert(length(psi.orbnames)==sum(psi.orbs));
		% register class
		psi = class(psi,'fermistate');
		% check whether the symbolic math toolbox is installed
		if ~isempty(ver('symbolic'))
			superiorto('sym');
		end
	otherwise
		error('Invalid number of input arguments.');
end
