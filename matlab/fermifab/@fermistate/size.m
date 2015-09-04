function d = size(psi,varargin)
%SIZE - Size of internal data
%
%   d = SIZE(psi,varargin)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

d = size(psi.data,varargin{:});
