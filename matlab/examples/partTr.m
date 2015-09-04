function A = partTr(A,sys,dim)
%
%   A = partTr(A,sys,dim)
%
%   Partial trace, trace out subsystem 'sys';
%   'A' must be a prod(dim) x prod(dim) matrix
%
%   Example:
%       A = kron(kron(A1,A2),A3) -> partTr(A,2,dim) == trace(A2)*kron(A1,A3)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

A = reshape(A,[dim(end:-1:1),dim(end:-1:1)]);

% move system to be traced out to the end
k = length(dim);
A = permute(A,[1:k-sys,k-sys+2:2*k-sys,2*k-sys+2:2*k,k-sys+1,2*k-sys+1]);

d = dim(sys);
A = reshape(A,(prod(dim)/d)^2,d^2);

A = sum(A(:,diag(reshape(1:d^2,d,d))),2);

A = reshape(A,prod(dim)/d,prod(dim)/d);
