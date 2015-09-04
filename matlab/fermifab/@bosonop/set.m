function op = set(op,varargin)
%SET - Set Bose operator information and data
%
%   op = SET(op,varargin)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

while ~isempty(varargin)
	switch varargin{1}
		case 'modes'
			assert(isa(varargin{2},'numeric') && length(varargin{2})==1);
			op.modes = round(varargin{2});
		case 'pFrom'
			assert(isa(varargin{2},'numeric') && length(varargin{2})==1);
			op.pFrom = round(varargin{2});
		case 'pTo'
			assert(isa(varargin{2},'numeric') && length(varargin{2})==1);
			op.pTo = round(varargin{2});
		case 'data'
			op.data = varargin{2};
		case 'modenames'
			op.modenames = varargin{2};
		otherwise
			error('Invalid property.');
	end
	varargin = varargin(3:end);
end

if isa(op.data,'sym'), op.data = simplify(op.data); end

assert(op.modes > 0 && op.pFrom >= 0 && op.pTo >= 0 && op.modes == length(op.modenames));
assert(size(op.data,1)==nchoosek(op.modes+op.pTo-1,op.pTo) && size(op.data,2)==nchoosek(op.modes+op.pFrom-1,op.pFrom));
