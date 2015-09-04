function op = elim_zeros(op,varargin)
%ELIM_ZEROS - Set small entries of internal Fermi operator data to zero
%
%   op = ELIM_ZEROS(op,...)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

op.data = elim_zeros(op.data,varargin{:});
