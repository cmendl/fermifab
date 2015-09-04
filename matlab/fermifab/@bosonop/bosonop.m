function op = bosonop(varargin)
%BOSONOP - Construct a Bose operator
%
%   op = BOSONOP();
%   op = BOSONOP(modes,pFrom,pTo);
%   op = BOSONOP(modes,pFrom,pTo,data);
%   op = BOSONOP(modes,pFrom,pTo,data,modenames);
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

switch length(varargin)
	case 0
		% standard values
		op.modes = 6;
		op.pFrom = 2;
		op.pTo = 2;
		op.data = zeros(nchoosek(op.modes+op.pTo-1,op.pTo),nchoosek(op.modes+op.pFrom-1,op.pFrom));
		% default mode names, change with 'set' command
		op.modenames = cell(1,op.modes);
		for j=1:op.modes, op.modenames{j} = num2str(j); end;
		% register class
		op = class(op,'bosonop');
	case 1
		% copy constructor
		assert(isa(varargin{1},'bosonop'));
		op = varargin{1};
	case 3
		for j=1:3
			assert(isa(varargin{j},'numeric') && length(varargin{j})==1);
		end
		op.modes = round(varargin{1});
		op.pFrom = round(varargin{2});
		op.pTo = round(varargin{3});
		assert(op.modes > 0 && op.pFrom >= 0 && op.pTo >= 0);
		% fill 'data' with zeros
		op.data = zeros(nchoosek(op.modes+op.pTo-1,op.pTo),nchoosek(op.modes+op.pFrom-1,op.pFrom));
		% default mode names
		op.modenames = cell(1,op.modes);
		for j=1:op.modes, op.modenames{j} = num2str(j); end;
		% register class
		op = class(op,'bosonop');
	case 4
		for j=1:3
			assert(isa(varargin{j},'numeric') && length(varargin{j})==1);
		end
		op.modes = round(varargin{1});
		op.pFrom = round(varargin{2});
		op.pTo = round(varargin{3});
		assert(op.modes > 0 && op.pFrom >= 0 && op.pTo >= 0);
		% copy data
		op.data = varargin{4};
		assert(size(op.data,1)==nchoosek(op.modes+op.pTo-1,op.pTo) && ...
			size(op.data,2)==nchoosek(op.modes+op.pFrom-1,op.pFrom),'Input data dimension mismatch.');
		if isa(op.data,'sym'), op.data = simplify(op.data); end
		% default mode names
		op.modenames = cell(1,op.modes);
		for j=1:op.modes, op.modenames{j} = num2str(j); end;
		% register class
		op = class(op,'bosonop');
	case 5
		for j=1:3
			assert(isa(varargin{j},'numeric') && length(varargin{j})==1);
		end
		op.modes = round(varargin{1});
		op.pFrom = round(varargin{2});
		op.pTo = round(varargin{3});
		assert(op.modes > 0 && op.pFrom >= 0 && op.pTo >= 0);
		% copy data
		op.data = varargin{4};
		assert(size(op.data,1)==nchoosek(op.modes+op.pTo-1,op.pTo) && ...
			size(op.data,2)==nchoosek(op.modes+op.pFrom-1,op.pFrom),'Input data dimension mismatch.');
		if isa(op.data,'sym'), op.data = simplify(op.data); end
		% mode names
		assert(iscell(varargin{5}));
		op.modenames = reshape(varargin{5},1,[]);
		assert(length(op.modenames)==op.modes,'Wrong number of mode names.');
		% register class
		op = class(op,'bosonop');
	otherwise
		error('Invalid number of input arguments.');
end
