function out = get(f, prop)
%GET   GET method for the CLASSICFUN class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the CLASSICFUN F.  The string PROP may be the name of a CLASSICFUN property
%   (see the CLASSICFUN classdef file for a list) or any of the following
%   strings:
%       'EXPONENTS'        - Exponents of F.ONEFUN or [0 0] if none.
%       'DELTAS'           - CLASSICFUNs have no delta functions.  Returns [].
%   If PROP is a string other than those specified above, GET(F, PROP) returns
%   the result of GET(F.ONEFUN, PROP).
%
% See also CLASSICFUN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case 'domain'
        out = f.domain;
    case 'onefun'
        out = f.onefun;
    case 'mapping'
        out = f.mapping;
    case 'exponents'
        if ( issing(f) )
            out = get(f.onefun, prop);
        else
            out = [0 0];
        end
    case 'deltas'
        out = [];
    otherwise
        out = get(f.onefun, prop);
end

end
