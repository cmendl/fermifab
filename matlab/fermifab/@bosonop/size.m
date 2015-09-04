function d = size(op,varargin)
%SIZE - Size of internal data
%
%   d = SIZE(op,varargin)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

d = size(op.data,varargin{:});
