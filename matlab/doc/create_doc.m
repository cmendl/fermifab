function create_doc()
%CREATE_DOC - Create documentation by 'publishing' Matlab files
%
%   Copyright (c) 2008-2015, Christian B. Mendl
%   All rights reserved.

%%
options.format = 'html';
options.maxOutputLines = 16;
options.stylesheet = 'mxdom2mathjax.xsl';	% use a custom stylesheet to enable MathJax
publish('usrguide.m',options);

%%
web(['file:///',pwd,'/html/usrguide.html']);
