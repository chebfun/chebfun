function out = get(f, prop)
%GET   GET method for the FOURTECH class.
%   P = GET(F,PROP) returns the property P specified in the string PROP from
%   the fun F. Valid entries for the string PROP are:
%       'VALUES' - Values of F at Chebyshev points.
%       'COEFFS' - Fourier coefficients of F.
%       'VSCALE' - Vertical scale of F.
%       'EPSLEVEL' - Happiness level of F.
%       'POINTS' - Equally spaced points where F is sampled.
%       'LVAL' - Value of F at -1.
%       'RVAL' - Value of F at +1.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        % Allow access to any of the properties of F via GET:
        out = f.(prop);
    case 'points'
        % Get the Fourier grid corresponding to the VALUES:
        out = f.points();
    case 'lval'
        % The value at -1:
        out = f.values(1,:);
    case 'rval'
        % The value at 1:
        out = f.values(end,:);
    otherwise
        error('FOURTECH:GET:proname', ...
            'Unknown property name ''%s'' for object of type fourtech.', prop);
end

end
