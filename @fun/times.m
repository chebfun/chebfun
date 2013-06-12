function f = times(f, g, varargin)
%.*   FUN multiplication.
%   F.*G multiplies FUN objects F and G or a FUN by a scalar if either
%   F or G is a scalar.
%
%   If F is a array-valued FUN, then F.*C is supported if C is a row
%   vector of doubles with the same number of columns as F.
%
%   If F and G are both FUN objects, they are assumed to have the same domain.
%   The method gives no warning if their domains don't agree, but the output of
%   the method will be gibberish.
%
% See also MTIMES, RDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% FUN * [] = []:
if ( isempty(f) || isempty(g) )
    f = []; 
    return
end

if ( ~isa(f, 'fun') )
    % If f is not a FUN, g must be
    g.onefun = times(f, g.onefun, varargin{:});
    % Swap arguments for output variable
    f = g;
    return
elseif ( isa(g, 'fun'))      % FUN * FUN
    % Multiply the onefun objects of f and g together.
    f.onefun = times(f.onefun, g.onefun, varargin{:});
elseif ( isa(g, 'double') )     % FUN * double
    % Multiply the onefun of f with g.
    f.onefun = times(f.onefun, g, varargin{:});
end


end
