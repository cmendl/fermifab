function n = numel(op,varargin)
%NUMEL - Number of array elements of a Fermi operator
%
%   n = NUMEL(op,varargin)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

n = numel(op.data,varargin{:});
