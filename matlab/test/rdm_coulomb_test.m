function err = rdm_coulomb_test(orbs,N)
%RDM_COULOMB_TEST - Two-particle reduced density matrix for the Coulomb operator test
%
%   err = RDM_COULOMB_TEST(orbs,N)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

assert(mod(orbs,2)==0,'Number of orbitals must be even.');
assert(N>1,'Particle number must at least be 2.');

% fictitious integral values
ints = crand(orbs/2,orbs/2,orbs/2,orbs/2);

% 2-body Coulomb operator on orbitals
% containing position and spin
X = fermi2coords(orbs,2);
V = zeros(size(X,2));
for j=1:length(V)
	x = round(X(:,j)/2);	% position
	s = mod(X(:,j),2);		% spin
	for k=1:length(V)
		y = round(X(:,k)/2);
		t = mod(X(:,k),2);
		p = symmetry([x(1),y(1),x(2),y(2)]);
		q = symmetry([x(1),y(2),x(2),y(1)]);
		V(j,k) = ...
			~xor(s(1),t(1))*~xor(s(2),t(2))*ints(p(1),p(2),p(3),p(4)) - ...
			~xor(s(1),t(2))*~xor(s(2),t(1))*ints(q(1),q(2),q(3),q(4));
	end
end
V = fermiop(orbs,2,2,V);

% 2-body Coulomb operator, stores spatial Coulomb integrals (ab|cd),
% so effectively divides number of orbitals by two
X = boson2coords(orbs/2,2);
W = zeros(size(X,2));
for j=1:size(W,1)
	for k=1:j	% need to fill lower triangular part only
		p = symmetry([X(:,j)',X(:,k)']);
		W(j,k) = ints(p(1),p(2),p(3),p(4));
	end
end
W = bosonop(orbs/2,2,2,W);

% pure states
psi = fermistate(orbs,N,crand(nchoosek(orbs,N),1));

err = ...
	abs(hprod(psi,p2N(V,N)*psi) - trace(V*rdm(psi,2))) + ...
	abs(trace(V*rdm(psi,2)) - trace(W*spintrace_coul(rdm(psi,2))));

% mixed states
rho = fermiop(orbs,N,N,crand(nchoosek(orbs,N)));

err = err + ...
	abs(trace(p2N(V,N)*rho) - trace(V*rdm(rho,2))) + ...
	abs(trace(V*rdm(rho,2)) - trace(W*spintrace_coul(rdm(rho,2))));


%%
function p = symmetry(p)

% assume that all orbitals are real
if p(1)>p(2), p([1,2]) = p([2,1]); end
if p(3)>p(4), p([3,4]) = p([4,3]); end

% use reverse-lexicographical Bose map ordering (last entry more significant than first)
if p(2)>p(4) || (p(2)==p(4) && p(1)>p(3))
	p = p([3,4,1,2]);
end
