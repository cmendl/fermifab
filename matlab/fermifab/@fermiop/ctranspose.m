function y = ctranspose(x)
%CTRANSPOSE - Complex conjugate transpose of a Fermi operator
%
%   y = CTRANSPOSE(x)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

y = x;

% switch pFrom and pTo
y.pFrom = x.pTo;
y.pTo = x.pFrom;

y.data = x.data';
