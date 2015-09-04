function psi = elim_zeros(psi,varargin)
%ELIM_ZEROS - Set small entries of internal Fermi state data to zero
%
%   psi = ELIM_ZEROS(psi,...)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

psi.data = elim_zeros(psi.data,varargin{:});
