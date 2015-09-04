function op = simplify(op)
%SIMPLIFY - Algebraic simplification of internal data enries
%
%   op = SIMPLIFY(op)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

op.data = simplify(op.data);
