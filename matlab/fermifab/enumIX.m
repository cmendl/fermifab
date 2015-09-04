function IX = enumIX(vardims)
%ENUMIX - Enumerate index combinations
%
%   IX = ENUMIX(vardims)
%
%   Example:
%
%   IX = ENUMIX([2,3]) ->
%
%       1     1
%       1     2
%       1     3
%       2     1
%       2     2
%       2     3
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

IX = zeros(prod(vardims),length(vardims));
for j=1:size(IX,2)
	IX(:,j) = kron(ones(1,prod(vardims(1:j-1))),kron(1:vardims(j),ones(1,prod(vardims(j+1:end)))))';
end
