function y = trace(op)
%TRACE - Matrix trace
%
%   y = TRACE(op)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(all(op.pFrom==op.pTo));
y = trace(op.data);
