function out = get(f, prop)
%GET   GET method for the DELTAFUN class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the DELTAFUN object F. Valid entries for the string PROP are:
%
%       'LOCATION'     - Location of the delta functions 
%       'DELTAMAG'     - Magnitude of the delta functions
%       'FUNPART'      - The smooth function contained in DELTAFUN.
%       'LVAL', 'RVAL' - Evaluate a DELTAFUN at an end point of the domain.
%                        If there is no delta function at the left or the right 
%                        end point, this is equivalent to evaluating the funPart 
%                        at the left or right end, otherwise, appropriately, 
%                        a  signed infinitey or a NaN is returned. 
%                        See DELTAFUN/FEVAL for fruther help on this.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        % Allow access to any of the properties of F via GET:
        out = f.(prop);
        
    case 'lval'        
        % Evaluate at the lefthand end of the domain:
        dom = f.domain;
        out = feval(f, dom(1));
        
    case 'rval'
        % Evaluate at the righthand end of the domain:
        dom = f.domain;
        out = feval(f, dom(end));
        
    otherwise
        % Delegate to the get method of funPart. All error handling will also be
        % done here:
        out = get(f.funPart, prop);                    
end

end