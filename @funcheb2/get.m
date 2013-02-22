function out = get(f, prop)
%GET    GET method for the FUNCHEB2 class
%   P = GET(F,PROP) returns the property P specified in the string PROP from
%   the fun F. Valid entries for the string PROP are:
%       'VALUES' - Values of F at Chebyshev points.
%       'COEFFS' - Chebyshev coefficients of F.
%       'VSCALE' - Vertical scale of F.
%       'EPSLEVEL' - Happiness level of F.
%       'POINTS' - 2nd-kind Chebyshev grid corresponding to F.
%       'LVAL'
%       'RVAL'

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        % Allow access to any of F's properties via GET.
        out = f.(prop);
    case 'points'
        % GET the Chebyshev grid corresponding to the VALUES.
        out = funcheb2.chebpts(size(f.values,1));
    case 'lval'
        % The value at -1.
        out = f.values(1,:);
    case 'rval'
        % The value at +1.
        out = f.values(end,:);
    otherwise
        error('CHEBFUN:FUNCHEB2:GET:proname', ...
            'Unknown property name ''%s'' for object of type funcheb2.', prop);
end