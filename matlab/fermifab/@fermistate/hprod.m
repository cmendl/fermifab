function y = hprod(chi,psi)
%HPROD - Hilbert space inner product
%
%   y = HPROD(chi,psi)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(isa(chi,'fermistate') && isa(psi,'fermistate'));
assert(chi.orbs == psi.orbs && chi.N == psi.N && chi.pairflag == psi.pairflag);

y = chi.data'*psi.data;
if isa(y,'sym'), y = simplify(y); end
