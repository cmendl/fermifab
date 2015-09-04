function y = trace_prod(A,B)
%TRACE_PROD - Calculate trace(A*B) efficiently
%
%   y = TRACE_PROD(A,B)
%
%   Input:      A, B    m x n and n x m matrices, respectively
%
%   Output:     y       trace(A*B)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if issparse(A)&&issparse(B)
	% rely on sparse structure
	y = trace(A*B);
else
	y = 0;
	for j=1:size(A,1)
		y = y + A(j,:)*B(:,j);
	end
end

if isa(y,'sym')
	y = simplify(y);
end
