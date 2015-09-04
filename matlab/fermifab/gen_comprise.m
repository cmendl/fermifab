%GEN_COMPRISE - Comprise input configuration scope
%
%   X = GEN_COMPRISE(orbs,N,merge_next)
%
%   Comprise input configuration scope, i.e., create index mapping from input to target configuration
%   defined by merging subspace neighbors in the input configuration. For example, from
%   orbs = [3,4,6], N = [1,2,4] via merge_next = [true,false] to orbs = [7,6], N = [3,4]
%
%   C mex function
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.
