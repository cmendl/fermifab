function display(psi)
%DISPLAY - Display a Fermi state
%
%   DISPLAY(psi)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

fprintf('\n%s = \n',inputname(1));
if ~psi.pairflag
	fprintf('\tFermi State (orbs == %d, N == %d)\n',sum(psi.orbs),sum(psi.N));
else
	fprintf('\tFermi State (orbs == %d, N == %d, pair structure)\n',sum(psi.orbs),sum(psi.N));
end
fprintf('\t%s\n\n',char(psi));
