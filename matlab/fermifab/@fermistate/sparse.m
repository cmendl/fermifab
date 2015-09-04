function psi = sparse(psi)
%SPARSE - Convert internal data to 'sparse' format
%
%   psi = SPARSE(psi)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if ~isa(psi.data,'sym')
	psi.data = sparse(psi.data);
else
	% sparse 'sym' not supported
	psi.data = sparse(double(psi.data));
end
