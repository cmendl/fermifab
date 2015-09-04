function err = config_test2(nc,N,maxC)
%CONFIG_TEST2 - Particle configuration test file
%
%   err = CONFIG_TEST2(nc,N,maxC)
%
%   Example: CONFIG_TEST2(3,5,[3,2,4])
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

maxC = reshape(maxC,1,[]);
assert(length(maxC)==nc);

% 'implements' configurations
X = boson2coords(nc,N);

% unrestricted (bosonic) case
C = gen_config(nc,N);
err = 0;
for j=1:size(C,2)
	c = zeros(nc,1);
	for k=1:nc
		c(k) = sum(X(:,j)==k);
	end
	err = err + norm(c-C(:,j));
end

% including restrictions
C = gen_config(nc,N,maxC);
IX = true(1,size(X,2));
for j=1:size(X,2)
	for k=1:nc
		if sum(X(:,j)==k)>maxC(k)
			IX(j) = false;
			break;
		end
	end
end
X = X(:,IX);
if size(C,2)~=size(X,2)
	err = 1;
	return;
end
for j=1:size(C,2)
	c = zeros(nc,1);
	for k=1:nc
		c(k) = sum(X(:,j)==k);
	end
	err = err + norm(c-C(:,j));
end
