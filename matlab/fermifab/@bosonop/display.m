function display(op)
%DISPLAY - Display a Bose operator
%
%   DISPLAY(op)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

fprintf('\n%s = \n\tBoson operator %u-body space -> %u-body space (modes == %d)\n\n',...
	inputname(1),op.pFrom,op.pTo,op.modes);
% output matrix
if nnz(op.data~=0)
	fprintf('Matrix representation w.r.t. ordered bosonic basis\n%s -> %s:\n\n',...
		print_base(op.modes,op.pFrom,op.modenames),print_base(op.modes,op.pTo,op.modenames));
	disp(op.data);
else
	fprintf('\t<zero matrix>\n\n');
end


function s = print_base(modes,N,modenames)

assert(modes==length(modenames));

if sum(N)==0
	s = '/vac>';
else
	X = boson2coords(modes,N);

	% determine maximum mode name length
	lg_max = 0;
	for j=1:length(modenames)
		lg = length(modenames{j});
		if lg > lg_max, lg_max = lg; end;
	end
	% choose delimiter accordingly
	if lg_max>1, delim = ' ';
	else delim = ''; end

	if size(X,2) < 5
		s = '(';
		for j=1:size(X,2)
			s = [s,'/'];
			for k=1:size(X,1)
				s = [s,modenames{X(k,j)},delim];
			end
			% remove trailing delimiter
			s = [s(1:end-length(delim)),'>, '];
		end
		% remove trailing ', '
		s = [s(1:end-2),')'];
	else
		s = '(/';
		for k=1:size(X,1)
			s = [s,modenames{X(k,1)},delim];
		end
		s = [s(1:end-length(delim)),'>, /'];
		for k=1:size(X,1)
			s = [s,modenames{X(k,2)},delim];
		end
		s = [s(1:end-length(delim)),'>, ... /'];
		for k=1:size(X,1)
			s = [s,modenames{X(k,end)},delim];
		end
		s = [s(1:end-length(delim)),'>)'];
	end
end
