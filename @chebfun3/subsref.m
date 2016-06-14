function varargout = subsref(f, index)
%SUBSREF   CHEBFUN3 subsref.
%   F(X, Y, Z) returns the values of the CHEBFUN3 object F evaluated at the
%   point (X, Y, Z).
% 
%   F.PROP returns the property of F specified in PROP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Implement vector inputs like f(0,0,0:1).  This kind of input
% is permitted in Chebfun2.

idx = index(1).subs;
switch index(1).type
    case '()'
        % FEVAL / COMPOSE
        if numel(idx) == 3
            % Find where to evaluate:
            x = idx{1}; 
            y = idx{2}; 
            z = idx{3};
            out = feval(f, x, y, z); 
            varargout = {out};
        else
            error('CHEBFUN:CHEBFUN3:subsref:inputs', ...
                'Can only evaluate at triples (X,Y,Z)')
        end
    case '.'
        % Call GET() for .PROP access.
        out = get(f, idx);
        if ( numel(index) > 1 )
            % Recurse on SUBSREF():
            index(1) = [];
            out = subsref(out, index);
        end
        varargout = {out};

    case '{}'
        % RESTRICT
        if ( length(idx) == 6 ) 
            out = restrict(f, [ idx{:} ]);
            varargout = {out};
        else
            error('CHEBFUN:CHEBFUN3:subsref:dimensions', ...
                'Index exceeds chebfun3 dimensions.')
        end

end

end