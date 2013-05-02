function f = times(f, g, varargin)
%.*   BNDFUN multiplication.
%   F.*G multiplies BNDFUN objects F and G or a BNDFUN by a scalar if either
%   F or G is a scalar.
%
%   If F is a array-valued BNDFUN, then F.*C is supported if C is a row
%   vector of doubles with the same number of columns as F.
%
% See also MTIMES, RDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% BNDFUN * [] = []:
if ( isempty(f) || isempty(g) )
    f = []; 
    return
end

if ( ~isa(f, 'bndfun') )        % Ensure F is a BNDFUN
    f = times(g, f, varargin{:});
    return
elseif ( isa(g, 'bndfun'))      % BNDFUN * BNDFUN
    % Multiply the onefun objects of f and g together.
    f.onefun = times(f.onefun, g.onefun, varargin{:});
elseif ( isa(g, 'double') )     % BNDFUN * double
    % Multiply the onefun of f with g.
    f.onefun = times(f.onefun, g, varargin{:});
end


end
