function op = sparse(op)
%SPARSE - Convert internal data to 'sparse' format
%
%   op = SPARSE(op)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if ~isa(op.data,'sym')
	op.data = sparse(op.data);
else
	% sparse 'sym' not supported;
	% convert non-zero entries only to speed up computations
	IX = op.data~=0;
	d = double(op.data(IX));
	op.data = sparse(size(op.data,1),size(op.data,2));
	op.data(IX) = d;
end
