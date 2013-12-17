function out = cat(varargin)
%CAT   Concatenation of CHEBFUN objects.
%   CAT of a CHEBFUN is not yet supported.

% [TODO]: Implement this.

if ( isa(dim, 'chebfun') )
    error('CHEBFUN:cat:dim',['First input to  CHEBFUN/CAT must be 1 or 2.\n',...
    '\n',...
    ' )\\_/(    ChebMeow\n',...
    ' (o.o)\n',...
    ' (_|_)_/']);
elseif ( dim == 1 )
    out = vertcat(varargin);
elseif ( dim == 2 )
    out = horzcat(varargin);
else
    error('CHEBFUN:cat:invalidDim', 'CHEBFUN dimension must be 1 or 2.');
end

end
