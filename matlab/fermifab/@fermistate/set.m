function psi = set(psi,varargin)
%SET - Set Fermi state information and data
%
%   psi = SET(psi,varargin)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

while ~isempty(varargin)
	switch varargin{1}
		case 'orbs'
			assert(isa(varargin{2},'numeric'));
			psi.orbs = round(varargin{2});
			varargin = varargin(3:end);
		case 'N'
			assert(isa(varargin{2},'numeric'));
			psi.N = round(varargin{2});
			varargin = varargin(3:end);
		case 'data'
			if size(varargin{2},1) < size(varargin{2},2)
				varargin{2} = varargin{2}.';
			end
			psi.data = varargin{2};
			if isa(psi.data,'sym')
				psi.data = simplify(psi.data);
			end
			varargin = varargin(3:end);
		case 'base'
			assert(isa(varargin{2},'numeric') && length(varargin{2})==1);
			varargin{2} = round(varargin{2});
			assert(1 <= varargin{2} && varargin{2} <= length(psi.data));
			psi.data = zeros(size(psi.data));
			psi.data(varargin{2}) = 1;
			varargin = varargin(3:end);
		case 'orbname'
			assert(ischar(varargin{3}))
			psi.orbnames{varargin{2}} = varargin{3};
			varargin = varargin(4:end);
		case 'orbnames'
			assert(iscell(varargin{2}));
			psi.orbnames = reshape(varargin{2},1,[]);
			varargin = varargin(3:end);
		case 'pairflag'
			if ~varargin{2}
				psi.pairflag = 0;
			else
				psi.pairflag = 1;
			end
			varargin = varargin(3:end);
		otherwise
			error('Invalid property ''%s''.',varargin{1});
	end
end

assert(~nnz(psi.N < 0 | psi.N > psi.orbs));
if ~psi.pairflag
	assert(size(psi.data,1)==nchoosekprod(psi.orbs,psi.N) && size(psi.data,2)==1,'Invalid dimension of the data vector.');
else
	assert(~nnz(mod(psi.orbs,2) | mod(psi.N,2)));
	assert(size(psi.data,1)==nchoosekprod(psi.orbs/2,psi.N/2) && size(psi.data,2)==1);
end
% make sure we have the right number of orbital names
if length(psi.orbnames) < sum(psi.orbs)
	psi.orbnames{sum(psi.orbs)} = num2str(sum(psi.orbs));
elseif length(psi.orbnames) > sum(psi.orbs)
	psi.orbnames = psi.orbnames(1:sum(psi.orbs));
end
for j=1:sum(psi.orbs)
	if isempty(psi.orbnames{j}), psi.orbnames{j} = num2str(j); end
end
