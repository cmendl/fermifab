function options = processopts(defopts,varargin)
%PROCESSOPTS - Process options in the form ('param',value)
%
%   options = PROCESSOPTS(defopts,varargin)
%
%   Input:      defopts         default options
%               varargin        'param', value, ...
%
%   Output:     options         options structure
%
%   Example:
%       PROCESSOPTS(struct('tol',1e-3,'maxiter',5,'method','fancy'),'tol',1e-5,'maxiter',42)
%       -> struct with fields
%               tol: 1e-5
%           maxiter: 42
%            method: 'fancy'
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

% default options
options = defopts;

% iterate through arguments
while ~isempty(varargin)
	if isfield(options,varargin{1})
		options.(varargin{1}) = varargin{2};
	else
		warning('Unrecognized option ''%s''',varargin{1});
	end
	varargin = varargin(3:end);
end
