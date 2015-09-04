function A = elim_zeros(A,varargin)
%ELIM_ZEROS - Set small entries to zero
%
%	A = ELIM_ZEROS(A)
%	A = ELIM_ZEROS(A,tol)
%
%   Set all entries in A with absolute value smaller
%   than a tolerance 'tol' to zero. If no tolerance is
%   provided, a small relative absolute value is used.
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if ~isnumeric(A), return; end

if nargin > 1
	tol = varargin{1};
else
	tol = max(size(A))*eps(normest(A));
end

X = real(A);
Y = imag(A);

X = X.*(abs(X)>tol);
Y = Y.*(abs(Y)>tol);

A = X + 1i*Y;
