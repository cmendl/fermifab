function T = permOp(A,order,dim1,dim2)
%PERMOP - Reorder operator tensor product
%
%   T = PERMOP(A,order,dim)
%   T = PERMOP(A,order,dim1,dim2)
%
%   Rearrange operator A according to 'order';
%   A must be a prod(dim) x prod(dim) or prod(dim1) x prod(dim2) matrix
%
%   Example:
%       PERMOP(kron(kron(A1,A2),A3),[3,1,2],dim) == kron(kron(A3,A1),A2)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if nargin==3
	dim2 = dim1;
end

% check input arguments
assert(length(dim1)==length(dim2),'dim1 and dim2 must be of the same length.');
assert(length(dim1)==length(order),'dim and order must be of the same length.');

order = order(end:-1:1);
p = [length(dim1)+1-order,2*length(dim1)+1-order];

% workaround for sparse matrices - reshape only works on full matrices
if issparse(A)
	A = full(A);
	sparsemode = true;
else
	sparsemode = false;
end

T = reshape(permute(reshape(A,[dim1(end:-1:1),dim2(end:-1:1)]),p),[prod(dim1),prod(dim2)]);

if sparsemode
	T = sparse(T);
end
