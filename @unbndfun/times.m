function f = times(f, g, varargin)
%.*   UNBNDFUN multiplication.
%   F.*G multiplies UNBNDFUN objects F and G or an UNBNDFUN by a scalar if
%   either F or G is a scalar.
%
%   If F is an array-valued UNBNDFUN, then F.*C is supported if C is a row
%   vector of doubles with the same number of columns as F.
%
% See also MTIMES, RDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% UNBNDFUN * [] = []:
if ( isempty(f) || isempty(g) )
    f = [];
    return
end

if ( ~isa(f, 'unbndfun') )        % Ensure F is an UNBNDFUN
    f = times(g, f, varargin{:});
    return
elseif ( isa(g, 'unbndfun'))      % UNBNDFUN * UNBNDFUN
    
    if checkDomain(f ,g)
        % Multiply the onefun objects of f and g together.
        f.onefun = times(f.onefun, g.onefun, varargin{:});
    else
        error('CHEBFUN:UNBNDFUN:times:domainMismatch',...
            'The domain of unbndfun f and g do not match.');
    end
    
elseif ( isa(g, 'double') )     % UNBNDFUN * double
    % Multiply the onefun of f with g.
    f.onefun = times(f.onefun, g, varargin{:});
end


end
