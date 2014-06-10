function out = get(f, prop)
%GET   GET method for the CHEBTECH2 class
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the CHEBTECH2 F.  The string PROP may be the name of a CHEBTECH2 property
%   (see the CHEBTECH and CHEBTECH2 classdef files for a list) or any of the
%   following strings:
%       'POINTS'          - 1st-kind Chebyshev grid corresponding to F.
%       'VALUES'          - Values of F at Chebyshev points.
%       'LVAL'            - Value of F at -1.
%       'RVAL'            - Value of F at +1.
%       'TECHCONSTRUCTOR' - Handle to the CHEBTECH1 constructor.
%
% See also CHEBTECH, CHEBTECH2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        out = f.(prop);
    case 'points'
        out = f.points();
    case 'lval'
        out = feval(f, -1);
    case 'rval'
        out = feval(f, 1);
    case 'values'
        out = f.coeffs2vals(f.coeffs);
    case 'techConstructor'
        out = @chebtech2;
    otherwise
        error('CHEBFUN:CHEBTECH2:GET:proname', ...
            'Unknown property name ''%s'' for object of type chebtech2.', prop);
end

end
