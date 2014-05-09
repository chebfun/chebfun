function out = get(f, prop)
%GET   GET method for the DELTAFUN class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the DELTAFUN object F. Valid entries for the string PROP are:
%
%       'LOCATION'    - Location of the delta functions 
%       'DELTAMAG'    - Magnitude of the delta functions
%       'FUNPART'     - The smooth function contained in DELTAFUN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        % Allow access to any of the properties of F via GET:
        out = f.(prop);
        
    case {'lval'}        
        % If there is a delta functions at the left end point, we weil either
        % get an Inf or a NaN:
        dom = f.domain;
        out = feval(f, dom(1));
        if ( ~isinf(out) && ~isnan(out) )
            % If there is no delta function at the left end, call the property
            % of the funPart:
            out = get(f.funPart, prop);
        end
    case {'rval'}
        % If there is a delta functions at the right end point, we weil either
        % get an Inf or a NaN:
        dom = f.domain;
        out = feval(f, dom(2));
        if ( ~isinf(out) && ~isnan(out) )
            % If there is no delta function at the right end, call the property
            % of the funPart:
            out = get(f.funPart, prop);
        end                       
    otherwise
        out = get(f.funPart, prop);                    
end

end
