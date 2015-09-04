function Y = crand(varargin)
%CRAND - Complex random numbers
%
%   Y = CRAND(varargin)
%
%   Uniformly distributed pseudorandom complex
%   numbers from the interval [-1/2,1/2] x i[-1/2,1/2],
%   similar to 'rand'
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

Y = rand(varargin{:})-ones(varargin{:})/2+1i*(rand(varargin{:})-ones(varargin{:})/2);
