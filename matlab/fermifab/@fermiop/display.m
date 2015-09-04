function display(op)
%DISPLAY - Display a Fermi operator
%
%   DISPLAY(op)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

fprintf('\n%s = \n',inputname(1));
if ~op.pairflag
	fprintf('\tFermi Operator wedge^%d H -> wedge^%d H (orbs == %d)\n\n',...
		sum(op.pFrom),sum(op.pTo),sum(op.orbs));
else
	fprintf('\tFermi Operator wedge^%d H -> wedge^%d H (orbs == %d, pair structure)\n\n',...
		sum(op.pFrom),sum(op.pTo),sum(op.orbs));
end
% output matrix
if any(op.data(:)~=0)
	fprintf('Matrix representation w.r.t. ordered Slater basis\n%s -> %s:\n\n',...
		print_base(op.orbs,op.pFrom,op.orbnames,op.pairflag),...
		print_base(op.orbs,op.pTo,op.orbnames,op.pairflag));
	disp(op.data);
else
	fprintf('\t<zero matrix>\n\n');
end


function s = print_base(orbs,N,orbnames,pairflag)

% determine maximum orbital name length
lg_max = 0;
for j=1:length(orbnames)
	lg = length(orbnames{j});
	if lg > lg_max, lg_max = lg; end;
end
% choose delimiter accordingly
if lg_max>1, delim = ' ';
else delim = ''; end

if sum(N)==0
	s = '|vac>';
elseif ~pairflag
	X = fermi2coords(orbs,N);

	if size(X,2) < 5
		s = '(';
		for j=1:size(X,2)
			s = [s,'|'];
			for k=1:size(X,1)
				s = [s,orbnames{X(k,j)},delim];
			end
			s = [s(1:end-length(delim)),'>, '];
		end
		% remove trailing ', '
		s = [s(1:end-2),')'];
	else
		s = '(|';
		for k=1:size(X,1)
			s = [s,orbnames{X(k,1)},delim];
		end
		% remove trailing delimiter
		s = [s(1:end-length(delim)),'>, |'];
		for k=1:size(X,1)
			s = [s,orbnames{X(k,2)},delim];
		end
		s = [s(1:end-length(delim)),'>, ... |'];
		for k=1:size(X,1)
			s = [s,orbnames{X(k,end)},delim];
		end
		s = [s(1:end-length(delim)),'>)'];
	end
else	% pairflag
	X = fermi2coords(orbs/2,N/2);

	if size(X,2) < 5
		s = '(';
		for j=1:size(X,2)
			s = [s,'|'];
			for k=1:size(X,1)
				s = [s,orbnames{2*X(k,j)-1},':',orbnames{2*X(k,j)},' '];
			end
			s = [s(1:end-1),'>, '];
		end
		% remove trailing ', '
		s = [s(1:end-2),')'];
	else
		s = '(|';
		for k=1:size(X,1)
			s = [s,orbnames{2*X(k,1)-1},':',orbnames{2*X(k,1)},' '];
		end
		% remove trailing delimiter
		s = [s(1:end-1),'>, |'];
		for k=1:size(X,1)
			s = [s,orbnames{2*X(k,2)-1},':',orbnames{2*X(k,2)},' '];
		end
		s = [s(1:end-1),'>, ... |'];
		for k=1:size(X,1)
			s = [s,orbnames{2*X(k,end)-1},':',orbnames{2*X(k,end)},' '];
		end
		s = [s(1:end-1),'>)'];
	end
end
