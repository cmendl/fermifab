function n = numel(x,varargin)
%NUMEL - Number of array elements of a Bose operator
%
%   n = NUMEL(x,varargin)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

n = numel(x.data,varargin{:});
