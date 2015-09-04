function val = get(op,propName)
%GET - Access Fermi operator information and data
%
%   val = GET(op,propName)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

switch propName
	case 'orbs'
		val = op.orbs;
	case 'pFrom'
		val = op.pFrom;
	case 'pTo'
		val = op.pTo;
	case 'data'
		val = op.data;
	case 'orbnames'
		val = op.orbnames;
	case 'pairflag'
		val = op.pairflag;
	otherwise
		error(['Unknown property ''' propName '''.']);
end
