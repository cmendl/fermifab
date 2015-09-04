function y = trace_prod(A,B)
%TRACE_PROD - Calculate trace(A*B) efficiently
%
%   y = TRACE_PROD(A,B)
%
%   Input:      A, B    fermiop
%                           internal data: m x n and n x m matrices, respectively
%
%   Output:     y       trace(A*B)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(A,'fermiop') && isa(B,'fermiop'),'Both input arguments must be ''fermiop'' class');
assert(all(A.orbs==B.orbs & A.pFrom==B.pTo & A.pTo==B.pFrom),'Dimension mismatch.');
assert(A.pairflag==B.pairflag,'Pair structure flags must agree.');

y = trace_prod(A.data,B.data);
