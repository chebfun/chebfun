function out = get(f, prop)
%GET   GET method for the SINGFUN class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the SINGFUN F.  The string PROP may be the name of a SINGFUN property (see
%   the SINGFUN classdef file for a list) or any of the following strings:
%       'VALUES'          - Values of F at its underlying gridpoints.
%       'LVAL'            - Value of F at -1.
%       'RVAL'            - Value of F at +1.
%   If PROP is a string other than those specified above, GET(F, PROP) returns
%   the result of GET(F.SMOOTHPART, PROP).
%
% See also SINGFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

switch prop      
    case 'values'
        % TODO:  This breaks encapsulation: not all techs have "points".
        pts = get(f.smoothPart, 'points');
        out = feval(f, pts);
    case 'lval'
        out = feval(f, -1);
    case 'rval'
        out = feval(f, 1);
    case 'smoothPart'
        out = f.smoothPart;
    case 'exponents'
        out = f.exponents;
    otherwise
        out = get(f.smoothPart, prop);
end

end
