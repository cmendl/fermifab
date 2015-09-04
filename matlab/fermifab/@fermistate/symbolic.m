function psi = symbolic(psi)
%SYMBOLIC - Convert internal data to symbolic format
%
%   psi = SYMBOLIC(psi)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

psi.data = sym(psi.data);
