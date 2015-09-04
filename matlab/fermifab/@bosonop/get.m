function val = get(op,propName)
%GET - Access Bose operator information and data
%
%   val = GET(op,propName)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

switch propName
	case 'modes'
		val = op.modes;
	case 'pFrom'
		val = op.pFrom;
	case 'pTo'
		val = op.pTo;
	case 'data'
		val = op.data;
	case 'modenames'
		val = op.modenames;
	otherwise
		error(['Unknown property ''' propName '''.']);
end
