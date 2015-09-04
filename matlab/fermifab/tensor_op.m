function y = tensor_op(x,N)
%TENSOR_OP - Operator tensor product
%
%   y = TENSOR_OP(x,N)
%
%   Calculate the tensor product operator
%   (x \otimes x ... \otimes x): wedge^N H -> wedge^N H
%   for a 1-body operator x: H -> H
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(x,'fermiop'));
assert(sum(get(x,'pFrom'))==1 && sum(get(x,'pTo'))==1,'Input must be a single-particle operator.');

y = fermiop(get(x,'orbs'),N,N,calc_tensor_op(full(get(x,'data')),N),get(x,'orbnames'));
