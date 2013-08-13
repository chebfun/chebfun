function f = subsasgn(f, index, val)
%SUBSREF   Chebfun subsref.
% ( )
%   F(X) = VAL assigns the values of vector VAL at locations specified in vector
%   X in the chebfun F. size(X, 1) should be equal to length(VALS) and size(X,
%   2) should be the number of columns in F. SUBSASGN introduces new breakpoints
%   in F at points in X that were not originally in F.DOMAIN. See DEFINEPOINT
%   for further details.
%
% .
%   CHEBFUN properties are restricted, so F.PROP = VAL has no effect.
%
% {}
%   F{A, B} = G redefines the CHEBFUN F in the interval [A, B] using G. See
%   CHEBFUN/DEFINEINTERVAL for further details.
%
% See also SUBSREF, DEFINEPOINT, DEFINEINTERVAL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

idx = index(1).subs;

switch index(1).type
    case '.'

        % [TODO]: Restrict access to this.
        f = builtin('subsasgn', f, index, val);

    case '()'

        % Define a point value:
        f = definePoint(f, idx{:}, val);

    case '{}'
        
        % Define an interval:
        f = defineInterval(f, [idx{:}], val);

    otherwise
        error('CHEBFUN:UnexpectedType',...
            ['??? Unexpected index.type of ' index(1).type]);
end

end
