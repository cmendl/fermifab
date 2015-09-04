function op = fermiop(varargin)
%FERMIOP - Construct a Fermi operator
%
%	op = FERMIOP();
%	op = FERMIOP(orbs,pFrom,pTo);
%	op = FERMIOP(orbs,pFrom,pTo,data);
%	op = FERMIOP(orbs,pFrom,pTo,data,orbnames);
%	op = FERMIOP(...,'pair');
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if ~isempty(varargin) && ischar(varargin{end})
	assert(strcmp(varargin{end},'pair'));
	varargin = varargin(1:end-1);
	op.pairflag = 1;
else
	% set 'pair flag' to false (default)
	op.pairflag = 0;
end

switch length(varargin)
	case 0
		% standard values
		op.orbs = 6;
		op.pFrom = 2;
		op.pTo = 2;
		if ~op.pairflag
			op.data = zeros(nchoosekprod(op.orbs,op.pTo),nchoosekprod(op.orbs,op.pFrom));
		else
			op.data = zeros(nchoosekprod(op.orbs/2,op.pTo/2),nchoosekprod(op.orbs/2,op.pFrom/2));
		end
		% custom orbital names
		op.orbnames = cell(1,sum(op.orbs));
		for j=1:length(op.orbnames), op.orbnames{j} = num2str(j); end;
		% register class
		op = class(op,'fermiop');
	case 1
		% copy constructor
		assert(isa(varargin{1},'fermiop') && op.pairflag==varargin{1}.pairflag);
		op = varargin{1};
	case 3
		for j=1:3
			assert(isa(varargin{j},'numeric'));
		end
		op.orbs  = reshape(round(varargin{1}),1,[]);
		op.pFrom = reshape(round(varargin{2}),1,[]);
		op.pTo   = reshape(round(varargin{3}),1,[]);
		assert(~nnz(op.pFrom < 0 | op.pFrom > op.orbs | op.orbs < 1 | ...
			op.pTo < 0 | op.pTo > op.orbs));
		% fill 'data' with zeros
		if ~op.pairflag
			op.data = sparse(nchoosekprod(op.orbs,op.pTo),nchoosekprod(op.orbs,op.pFrom));
		else
			assert(~nnz(mod(op.orbs,2) | mod(op.pFrom,2) | mod(op.pTo,2)));
			op.data = sparse(nchoosekprod(op.orbs/2,op.pTo/2),nchoosekprod(op.orbs/2,op.pFrom/2));
		end
		% custom orbital names
		op.orbnames = cell(1,sum(op.orbs));
		for j=1:length(op.orbnames), op.orbnames{j} = num2str(j); end;
		% register class
		op = class(op,'fermiop');
	case 4
		for j=1:3
			assert(isa(varargin{j},'numeric'));
		end
		op.orbs  = reshape(round(varargin{1}),1,[]);
		op.pFrom = reshape(round(varargin{2}),1,[]);
		op.pTo   = reshape(round(varargin{3}),1,[]);
		assert(~nnz(op.pFrom < 0 | op.pFrom > op.orbs | op.orbs < 1 | ...
			op.pTo < 0 | op.pTo > op.orbs),'Invalid orbital or particle number argument.');
		% copy data
		op.data = varargin{4};
		if ~op.pairflag
			assert(size(op.data,1)==nchoosekprod(op.orbs,op.pTo) && size(op.data,2)==nchoosekprod(op.orbs,op.pFrom),'Input data dimension mismatch.');
		else
			assert(~nnz(mod(op.orbs,2)~=0 | mod(op.pFrom,2)~=0 | mod(op.pTo,2)~=0));
			assert(size(op.data,1)==nchoosekprod(op.orbs/2,op.pTo/2) && size(op.data,2)==nchoosekprod(op.orbs/2,op.pFrom/2));
		end
		if isa(op.data,'sym'), op.data = simplify(op.data); end
		% custom orbital names
		op.orbnames = cell(1,sum(op.orbs));
		for j=1:length(op.orbnames), op.orbnames{j} = num2str(j); end;
		% register class
		op = class(op,'fermiop');
	case 5
		for j=1:3
			assert(isa(varargin{j},'numeric'));
		end
		op.orbs  = reshape(round(varargin{1}),1,[]);
		op.pFrom = reshape(round(varargin{2}),1,[]);
		op.pTo   = reshape(round(varargin{3}),1,[]);
		assert(~nnz(op.pFrom < 0 | op.pFrom > op.orbs | op.orbs < 1 | ...
			op.pTo < 0 | op.pTo > op.orbs),'Invalid orbital or particle number argument.');
		% copy data
		op.data = varargin{4};
		if ~op.pairflag
			assert(size(op.data,1)==nchoosekprod(op.orbs,op.pTo) && size(op.data,2)==nchoosekprod(op.orbs,op.pFrom),'Input data dimension mismatch.');
		else
			assert(~nnz(mod(op.orbs,2)~=0 | mod(op.pFrom,2)~=0 | mod(op.pTo,2)~=0));
			assert(size(op.data,1)==nchoosekprod(op.orbs/2,op.pTo/2) && size(op.data,2)==nchoosekprod(op.orbs/2,op.pFrom/2));
		end
		if isa(op.data,'sym'), op.data = simplify(op.data); end
		% orbital names
		assert(iscell(varargin{5}));
		op.orbnames = reshape(varargin{5},1,[]);
		assert(length(op.orbnames)==sum(op.orbs),'Wrong number of orbital names.');
		% register class
		op = class(op,'fermiop');
	otherwise
		error('Invalid number of input arguments.');
end
