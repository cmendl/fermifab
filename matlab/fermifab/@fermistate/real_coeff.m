function psi = real_coeff(psi)
%REAL_COEFF - Apply a complex phase rotation to obtain real coefficients
%
%   psi = REAL_COEFF(psi)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

% 'max' compares the absolute values of its arguments
psi.data = elim_zeros(exp(-1i*angle(max(psi.data)))*psi.data);
