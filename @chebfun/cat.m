function out = cat(dim, varargin)
%CAT   Concatenation of CHEBFUN objects.
%   CAT of a CHEBFUN is not yet supported.

if ( dim == 1 )
    out = chebmatrix(varargin);
elseif ( dim == 2 )
    out = chebmatrix(varargin.');
else
    error('CHEBFUN:cat:noSupport', 'CAT of a CHEBFUN is not yet supported.');
end

% [TODO]: Implement this.

end
