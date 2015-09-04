%GEN_CONFIG - Generate configurations
%
%   C = GEN_CONFIG(nc,N)
%   C = GEN_CONFIG(nc,N,maxC)
%
%   Enumerate all configurations of 'N' particles in 'nc' slots,
%   such that not more than 'maxC(i)' particles are in slot 'i'.
%
%   Examples:
%
%       GEN_CONFIG(3,2) ->
%           2     1     0     1     0     0
%           0     1     2     0     1     0
%           0     0     0     1     1     2
%
%       GEN_CONFIG(3,2,[2,2,1]) ->
%           2     1     0     1     0
%           0     1     2     0     1
%           0     0     0     1     1
%
%   C mex function
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.
