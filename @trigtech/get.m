function out = get(f, prop)
%GET   GET method for the TRIGTECH class
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the TRIGTECH F.  The string PROP may be the name of a TRIGTECH property (see
%   the TRIGTECH classdef files for lists) or any of the following strings:
%       'POINTS'          - Equally spaced points where F is sampled.
%       'VALUES'          - Values of F at trig points.
%       'LVAL'            - Value of F at -1.
%       'RVAL'            - Value of F at +1.
%       'TECH'            - Handle to the TRIGTECH constructor. *
%
% See also CHEBTECH, CHEBTECH2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%  * Currently get(f, 'tech') returns a function handle to the tech constructor.
%    This may change in future to return instead an empty instance of the tech.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch prop
    case fieldnames(f)
        % Allow access to any of the properties of F via GET:
        out = f.(prop);
    case 'points'
        % Get the trigier grid corresponding to the VALUES:
        out = f.points();
    case 'lval'
        % The value at -1:
        out = feval(f, -1);
    case 'rval'
        % The value at 1:
        out = feval(f, 1);
    case 'tech'
        % TODO: Return function handle, or empty instance of the tech?
        out = @trigtech;        
    otherwise
        error('CHEBFUN:TRIGTECH:GET:proname', ...
            'Unknown property name ''%s'' for object of type trigtech.', prop);
end

end
