function out = cat(dim, varargin)
%CAT   Concatenation of CHEBFUN objects.
%    CAT(2, F, G) is the same as [F, G].
%    CAT(1, F, G) is the same as [F ; G].
% 
%    G = CAT(DIM, F1, F2, F3, F4,...) concatenates the input CHEBFUN objects F1,
%    F2, etc. along the dimension DIM.
%
% See also HORZCAT, VERTCAT, NUM2CELL.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(dim, 'chebfun') )
    error('CHEBFUN:CHEBFUN:cat:dim', ...
        ['First input to CHEBFUN/CAT must be 1 or 2.\n',...
            '\n',...
            ' )\\_/(    ChebMeow\n',...
            ' (o.o)\n',...
            ' (_|_)_/']);

elseif ( dim == 1 )
    out = vertcat(varargin{:});

elseif ( dim == 2 )
    out = horzcat(varargin{:});

else
    error('CHEBFUN:CHEBFUN:cat:invalidDim', ...
        'CHEBFUN dimension must be 1 or 2.');

end

end
