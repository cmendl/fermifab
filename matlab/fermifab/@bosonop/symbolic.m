function op = symbolic(op)
%SYMBOLIC - Convert internal data to symbolic format
%
%   op = SYMBOLIC(op)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

op.data = sym(op.data);
