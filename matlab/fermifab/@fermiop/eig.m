function varargout = eig(op,varargin)
%EIG - Eigenvalues and eigenvectors of Fermi operators
%
%   d = EIG(op,varargin)
%   [V,D] = EIG(op,varargin)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(op,'fermiop'));

varargout = cell(1,max(nargout,1));
[varargout{:}] = eig(op.data,varargin{:});

if nargout==2
	op.data = varargout{1};
	varargout{1} = op;
end
