function c = nchoosekprod(v,k)
%NCHOOSEKPROD - Compute the product of binomial coefficients
%
%   c = NCHOOSEKPROD(v,k)
%
%   Compute prod_j nchoosek(v(j),k(j))
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(length(v) == length(k));

c = 1;
for j=1:length(v)
	c = c*nchoosek(v(j),k(j));
end
