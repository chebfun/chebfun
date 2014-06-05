function out = get(f, prop)
%GET   GET method for the DELTAFUN class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the DELTAFUN object F. Valid entries for the string PROP are:
%
%       'LOCATION'     - Location of the delta functions 
%       'DELTAMAG'     - Magnitude of the delta functions
%       'FUNPART'      - The smooth function contained in DELTAFUN.
%       'LVAL', 'RVAL' - Evaluate a DELTAFUN at an end point of its domain.
%                        Whether there is a delta function at the left or the 
%                        right end point, this is equivalent to just evaluating
%                        the funPart at the left or right end. 
%                        See DELTAFUN/FEVAL for further details.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        % Allow access to any of the properties of F via GET:
        out = f.(prop);                
    otherwise
        % Delegate to the get method of funPart. All error handling will also be
        % done here:
        out = get(f.funPart, prop);                    
end

end