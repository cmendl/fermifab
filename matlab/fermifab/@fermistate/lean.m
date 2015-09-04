function psi = lean(psi,varargin)
%LEAN - Delete references to completely filled orbitals
%
%   psi = LEAN(psi)
%   psi = LEAN(psi,IX)  considers orbitals indexed by 'IX' only
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if nargin==1
	IX = 1:length(psi.orbs);
else
	IX = varargin{1};
end
if isempty(IX), return; end;
IX = unique(IX);
assert(1<=IX(1) && IX(end)<=length(psi.orbs),'Index out of range.');

% calculate complement of 'IX', i.e., non-affected orbitals
IY = true(1,length(psi.orbs));
IY(IX) = false;

% retained orbitals
IK = find(psi.N < psi.orbs | IY);

% orbital names
cumH = [0,cumsum(psi.orbs)];
IY = [];
for j=IK, IY = [IY,cumH(j)+1:cumH(j+1)]; end
psi.orbnames = psi.orbnames(IY);

% orbitals and particles
psi.orbs = psi.orbs(IK);
psi.N = psi.N(IK);

% vacuum state
if isempty(psi.orbs)
	psi.orbs = 1;
	psi.N = 0;
	psi.orbnames = {1};
end

assert(length(psi.data)==nchoosekprod(psi.orbs,psi.N));
