function y = mpower(x,n)
%MPOWER - Matrix power a Fermi operator
%
%   y = MPOWER(x,n)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(x.pFrom == x.pTo);

y = x;
y.data = x.data^n;
