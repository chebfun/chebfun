function f = rdivide(f, g, varargin)
%./   Right array divide for a CLASSICFUN.
%   F ./ C divides the CLASSICFUN F by a double C. If F is an array-valued CLASSICFUN with M
%   columns, then C must be either a scalar or a 1xM array.
%
%   Alternatively if C is a CLASSICFUN with the same number of columns as F or F is a
%   scalar, then the resulting division is returned if C is found to have no
%   roots in its domain. The division is performed column-wise. Note that F and
%   C are assumed to have the same domain. The method gives no warning if their
%   domains don't agree, but the output of the method will be meaningless.
%
% See also MRDIVIDE, TIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% CLASSICFUN ./ [] = []:
if ( isempty(f) || isempty(g) )
    f = [];
    return
end

% Cast ONEFUN of G to SINGFUN if G has vanishing values at the endpoints:
if ( ~isnumeric(g) )
    g = extractBoundaryRoots(g);
end

% Look at different cases:

if ( isa(g, 'double') )         % CLASSICFUN ./ DOUBLE
    
    % Divide the ONEFUN:
    f.onefun = rdivide(f.onefun, g, varargin{:});
    
elseif ( isa(f, 'double')  )    % DOUBLE ./ CLASSICFUN
    
    % If g has roots in [a, b], we're going to get an error. However, don't
    % worry about catching errors here; rather, just allow the error generated
    % by the ONEFUN method to be thrown.
    g.onefun = rdivide(f, g.onefun, varargin{:});
    
    % Assign g to be the output argument
    f = g;
    
else                            % CLASSICFUN ./ CLASSICFUN
    
    % Both f and g are FUN objects. Call RDIVIDE() on their onefuns. Similar to
    % above, don't worry about error catching.
    f.onefun = rdivide(f.onefun, g.onefun, varargin{:});
    
end

end
