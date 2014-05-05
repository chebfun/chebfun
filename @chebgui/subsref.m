function varargout = subsref(cg,index)
% SUBSREF   Obtain info from chebgui

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Allow calls on the form guifile.options.plotting
idx = index(1).subs;
if size(index,2) > 1 && strcmp(idx,'options')
    idx = index(2).subs;
end

switch index(1).type
    case '.'
        varargout = { get(cg,idx) };
    otherwise
        error('CHEBGUI:subsref:indexType',['Unexpected index.type of ' index(1).type]);
end