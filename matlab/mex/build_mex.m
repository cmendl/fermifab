function build_mex()
%BUILD_MEX - Compile the C source code of the FermiFab toolbox (Matlab MEX interface)
%
%BUILD_MEX()
%
%    Copyright (c) 2008-2015, Christian B. Mendl
%    All rights reserved.

clear gen_rdm gen_rdm_pairs gen_annihil gen_annihil_pairs calc_tensor_op...
	coords2fermi fermi2coords coords2boson boson2coords map_pairs...
	gen_comprise gen_spintrace_sing gen_spintrace_coul gen_config;

fprintf('Compiling mex back-end files...\n');
cd([fileparts(mfilename('fullpath')),'/src']);
v = ver('Matlab');
% mex options
opts = '-I../include/ -O CFLAGS="$CFLAGS -std=c99" ';
if str2double(v.Version) > 7.2
	opts = [opts,'-largeArrayDims '];
else
	opts = [opts,'-DDEFINE_INDEX '];
end

eval(['mex ',opts,'gen_rdm.c bitfield.c fermi_map.c util.c mextools.c']);
eval(['mex ',opts,'gen_rdm_pairs.c bitfield.c fermi_map.c util.c mextools.c']);
eval(['mex ',opts,'gen_annihil.c bitfield.c fermi_map.c util.c mextools.c']);
eval(['mex ',opts,'gen_annihil_pairs.c bitfield.c fermi_map.c util.c mextools.c']);
eval(['mex ',opts,'calc_tensor_op.c bitfield.c fermi_map.c util.c']);
eval(['mex ',opts,'coords2fermi.c bitfield.c fermi_map.c util.c']);
eval(['mex ',opts,'fermi2coords.c bitfield.c fermi_map.c util.c']);
eval(['mex ',opts,'coords2boson.c bitfield.c boson_map.c fermi_map.c util.c']);
eval(['mex ',opts,'boson2coords.c bitfield.c boson_map.c fermi_map.c util.c']);
eval(['mex ',opts,'map_pairs.c bitfield.c fermi_map.c util.c mextools.c']);
eval(['mex ',opts,'gen_comprise.c bitfield.c fermi_map.c util.c']);
eval(['mex ',opts,'gen_spintrace_sing.c bitfield.c fermi_map.c util.c mextools.c']);
eval(['mex ',opts,'gen_spintrace_coul.c bitfield.c fermi_map.c boson_map.c util.c mextools.c']);
eval(['mex ',opts,'gen_config.c util.c']);

fprintf('Done.\n');
fprintf('Moving files to the ''fermifab'' folder...');
movefile *.mex* ../../fermifab/
cd ..
fprintf(' Done.\n');
