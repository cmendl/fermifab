function err = tensor_op_test(orbs,p,N)
%TENSOR_OP_TEST - Operator tensor product test
%
%   err = TENSOR_OP_TEST(orbs,p,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

err = 0;

% random "wave function"
psi = fermistate(orbs,N,crand(nchoosek(orbs,N),1));

% random unitary matrix
[U,R] = qr(crand(orbs));
U = fermiop(orbs,1,1,U);

Up = tensor_op(U,p);
UN = tensor_op(U,N);

% check if 'UN' is actually unitary
err = err + norm(UN'*UN - fermiop(orbs,N,N,eye(size(UN))));

err = err + norm(Up'*rdm(UN*psi,p)*Up - rdm(psi,p));
