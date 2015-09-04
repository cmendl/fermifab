function z = mrdivide(x,y)
%MRDIVIDE - Matrix right division for Fermi states and operators
%
%   z = MRDIVIDE(x,y)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(x,'fermistate'));

z = x;
z.data = x.data/y;
