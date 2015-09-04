function T = partTp(A,sys,dim)
%
%   T = partTp(A,sys,dim)
%
%   Partial transpose with respect to subsystem 'sys';
%   'A' must be a prod(dim) x prod(dim) matrix
%
%   Example:
%       A = kron(kron(A1,A2),A3) -> partTp(A,2,dim) == kron(kron(A1,A2.'),A3)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

p = 1:2*length(dim);
p(length(dim)-sys+1) = 2*length(dim)-sys+1;
p(2*length(dim)-sys+1) = length(dim)-sys+1;

T = reshape(permute(reshape(A,[dim(end:-1:1),dim(end:-1:1)]),p),[prod(dim),prod(dim)]);
