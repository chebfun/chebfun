function f = subsasgn(f, index, val)
%SUBSASGN   Chebfun SUBSASGN.
% ( )
%   F(X) = VAL assigns the values VAL at locations specified by X to the 
%   CHEBFUN F. SIZE(X, 1) should be equal to LENGTH(VAL) and SIZE(X, 2) should 
%   be the number of columns in F. SUBSASGN introduces new breakpoints
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
% See http://www.chebfun.org for Chebfun information.

idx = index(1).subs;

switch index(1).type
    case '.'

        % [TODO]: Restrict access to this.
        f = builtin('subsasgn', f, index, val);

    case '()'

        % [TODO]: Check inputs more carefully. (In particular for ROW chebfuns)
        if ( ischar(idx{1}) && strcmp(idx{1}, ':') )
            % Assign a column:
            f = assignColumns(f, idx{2}, val);
        else 
            % Define a point value:
            f = definePoint(f, idx{1}, val);
        end

    case '{}'
        
        % Define an interval:
        f = defineInterval(f, [idx{:}], val);

    otherwise
        error('CHEBFUN:UnexpectedType',...
            ['??? Unexpected index.type of ' index(1).type]);
end

end
