%COORDS2FERMI - Convert N-particle orbital "coordinates" to a Slater basis index and sign
%
%   n = COORDS2FERMI(x,orbs)
%   [n,sign] = COORDS2FERMI(x,orbs)
%
%   Examples:
%       COORDS2FERMI([1,3,4,6],6) -> [8,1]
%       COORDS2FERMI([3,1,4,6],6) -> [8,-1]     (sign has changed due to flipped coordinates)
%       COORDS2FERMI([1,3,4,4],6) -> [-1,0]     (invalid coordinates, '4' is occupied twice)
%
%   C mex function
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.
