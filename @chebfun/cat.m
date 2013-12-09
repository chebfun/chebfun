function out = cat(varargin)
%CAT   Concatenation of CHEBFUN objects.
%   CAT of a CHEBFUN is not yet supported.

% [TODO]: Implement this.

if ( dim == 1 )
    out = vertcat(varargin);
elseif ( dim == 2 )
    out = horzcat(varargin);
else
    error('CHEBFUN:cat:moo', 'Derp?');
end

end
