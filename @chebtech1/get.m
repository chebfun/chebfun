function out = get(f, prop)
%GET   GET method for the CHEBTECH1 class
%   P = GET(F,PROP) returns the property P specified in the string PROP from
%   the fun F. Valid entries for the string PROP are:
%       'VALUES' - Values of F at Chebyshev points.
%       'COEFFS' - Chebyshev coefficients of F.
%       'VSCALE' - Vertical scale of F.
%       'EPSLEVEL' - Happiness level of F.
%       'POINTS' - 1st-kind Chebyshev grid corresponding to F.
%       'LVAL' - Value of F at -1.
%       'RVAL' - Value of F at +1.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        % Allow access to any of the properties of F via GET:
        out = f.(prop);
    case 'points'
        % Get the Chebyshev grid corresponding to the VALUES:
        out = chebtech1.chebpts(size(f.values,1));
    case 'lval'
        % The value at -1:
        out = feval(f, -1);
    case 'rval'
        % The value at +1:
        out = feval(f, 1);
    case 'values' 
        out = f.coeffs2vals(f.coeffs);
    otherwise
        error('CHEBFUN:CHEBTECH1:GET:proname', ...
            'Unknown property name ''%s'' for object of type chebtech1.', prop);
end

end
