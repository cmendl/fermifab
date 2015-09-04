function val = get(psi,propName)
%GET - Access Fermi state information and data
%
%   val = GET(psi,propName)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

switch propName
	case 'orbs'
		val = psi.orbs;
	case 'N'
		val = psi.N;
	case 'data'
		val = psi.data;
	case 'orbnames'
		val = psi.orbnames;
	case 'pairflag'
		val = psi.pairflag;
	otherwise
		error(['Unknown property ''' propName '''.']);
end
