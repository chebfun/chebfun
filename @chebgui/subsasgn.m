function varargout = subsasgn(cg, index, varargin)
% SUBSASGN   Modify a CHEBGUI object.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Allow calls on the form guifile.options.plotting
idx = index(1).subs;
if ( (size(index, 2) > 1) && strcmp(idx, 'options') )
    idx = index(2).subs;
end

vin = varargin{:};
switch ( index(1).type )
    case '.'
        varargout = {set(cg, idx, vin)};
        
    otherwise
        error('CHEBFUN:CHEBGUI:subsasgn:indexType', ...
            ['Unexpected index.type of ' index(1).type]);
        
end

end
