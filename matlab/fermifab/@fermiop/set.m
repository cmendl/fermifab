function op = set(op,varargin)
%SET - Set Fermi operator information and data
%
%   op = SET(op,varargin)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

while ~isempty(varargin)
	switch varargin{1}
		case 'orbs'
			assert(isa(varargin{2},'numeric') && length(varargin{2})==1);
			op.orbs = round(varargin{2});
		case 'pFrom'
			assert(isa(varargin{2},'numeric') && length(varargin{2})==1);
			op.pFrom = round(varargin{2});
		case 'pTo'
			assert(isa(varargin{2},'numeric') && length(varargin{2})==1);
			op.pTo = round(varargin{2});
		case 'data'
			op.data = varargin{2};
		case 'orbnames'
			op.orbnames = varargin{2};
		case 'pairflag'
			if ~varargin{2}
				psi.pairflag = 0;
			else
				psi.pairflag = 1;
			end
		otherwise
			error('Invalid property.');
	end
	varargin = varargin(3:end);
end

if isa(op.data,'sym'), op.data = simplify(op.data); end

assert(~nnz(op.pFrom < 0 | op.pFrom > op.orbs | op.pTo < 0 | op.pTo > op.orbs | op.orbs < 1));
if ~op.pairflag
	assert(size(op.data,1)==nchoosekprod(op.orbs,op.pTo) && size(op.data,2)==nchoosekprod(op.orbs,op.pFrom));
else
	assert(size(op.data,1)==nchoosekprod(op.orbs/2,op.pTo/2) && size(op.data,2)==nchoosekprod(op.orbs/2,op.pFrom/2));
end
