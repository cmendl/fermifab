function s = char(psi)
%CHAR - Convert Fermi state to a string
%
%   s = CHAR(psi)
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

if ~psi.pairflag
	orbseff = psi.orbs;
	Neff = psi.N;
else
	orbseff = psi.orbs/2;
	Neff = psi.N/2;
end

% ignore near-zero entries
psi = elim_zeros(psi);

X = fermi2coords(orbseff,Neff);

% determine maximum orbital name length
lg_max = 0;
for j=1:length(psi.orbnames)
	lg = length(psi.orbnames{j});
	if lg > lg_max, lg_max = lg; end;
end
% choose delimiter accordingly
if lg_max>1 && sum(Neff)>0, delim = ' ';
else delim = ''; end

s = [];
for j=find(psi.data ~= 0).'
	if psi.data(j) == 1
		s = [s,' + |'];
	elseif isa(psi.data,'sym')
		s = sprintf('%s + (%s)|',s,char(psi.data(j)));
	else
		% 'num2str': handle complex data correctly
		s = sprintf('%s + (%s)|',s,num2str(psi.data(j)));
	end
	if ~psi.pairflag
		for k=1:sum(Neff)
			s = [s,psi.orbnames{X(k,j)},delim];
		end
		s = [s(1:end-length(delim)),'>'];
	else
		for k=1:sum(Neff), s = [s,psi.orbnames{2*X(k,j)-1},':',psi.orbnames{2*X(k,j)},' ']; end;
		s = [s(1:end-1),'>'];
	end
end

if isempty(s)
	s = '0';
else
	% remove leading ' + '
	s = s(4:end);
end

s = strrep(s,'|>','|vac>');
