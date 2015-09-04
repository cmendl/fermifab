function n = norm(op,varargin)
%NORM - Norm of a Fermi operator
%
%   n = NORM(op,varargin)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

n = norm(op.data,varargin{:});
