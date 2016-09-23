function out = get(f, prop)
%GET   GET method for the DELTAFUN class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the DELTAFUN object F.  The string PROP may be the name of a DELTAFUN
%   property (see the DELTAFUN classdef file for a list) or any of the
%   following strings:
%       'DELTAS'       - Return the locations and magnitude of delta functions
%                        as a single matrix with the first row corresponding to
%                        the locations of the delta functions, and the
%                        magnitudes appended below the first row.
%   If PROP is a string other than those specified above, GET(F, PROP) returns
%   the result of GET(F.FUNPART, PROP).
%
% See also DELTAFUN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        out = f.(prop);
    case 'deltas'
        out = [f.deltaLoc ; f.deltaMag];
    otherwise
        out = get(f.funPart, prop);
end

end
