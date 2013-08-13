function f = times(f, g, varargin)
%.*   FUN multiplication.
%   F.*G multiplies FUN objects F and G or a FUN by a scalar if either F or G is
%   a scalar.
%
%   If F is a array-valued FUN, then F.*C is supported if C is a row vector of
%   doubles with the same number of columns as F. Similarly if F if a row
%   vector of doubles and C is an array-valued FUN
%
%   If F and G are both FUN objects, they are assumed to have the same domain
%   and, if they are array-valued, the same number of columns. The method gives
%   no warning if their domains don't agree, but the output of the method will
%   be meaningless.
%
% See also MTIMES, RDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % FUN * [] = []
    
    % Return empty:
    f = []; 
    
elseif ( ~isa(f, 'fun') )       % double * FUN
    
    % If f is not a FUN, g must be. Call ONEFUN/TIMES():
    g.onefun = times(f, g.onefun, varargin{:});
    % Swap arguments for output variable:
    f = g;
    
elseif ( isa(g, 'fun'))         % FUN * FUN
    
    % Multiply the ONEFUN objects of f and g together.
    f.onefun = times(f.onefun, g.onefun, varargin{:});
    
elseif ( isa(g, 'double') )     % FUN * double
    
    % Multiply the ONEFUN of f with g.
    f.onefun = times(f.onefun, g, varargin{:});
    
end

end
