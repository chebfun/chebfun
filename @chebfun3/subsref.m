function varargout = subsref(f, index)
%SUBSREF   CHEBFUN3 subsref.
%   F(X, Y, Z) returns the values of the CHEBFUN3 object F evaluated at the
%   point (X, Y, Z).
%
%   F(G) where G is a CHEBFUN3V returns the CHEBFUN3 representing the
%   composition F(G).  If G is a CHEBFUN2V, then F(G) is a CHEBFUN2.  If G is a
%   CHEBFUN with three columns, then F(G) is a CHEBFUN.  If G is a SPHEREFUNV,
%   then F(G) is a SPHEREFUN.
%
%   F(X, Y, Z) where X, Y, Z are CHEBFUN3 objects is a CHEBFUN3 representing the
%   composition.  Similarly if X, Y, Z are CHEBFUN or CHEBFUN2 objects.
%
%   F.PROP returns the property of F specified in PROP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Implement vector inputs like f(0,0,0:1).  This kind of input
% is permitted in Chebfun2.

idx = index(1).subs;
switch index(1).type
    case '()'
        % FEVAL / COMPOSE
        if ( numel(idx) == 3 )
            % Find where to evaluate:
            x = idx{1};
            y = idx{2};
            z = idx{3};
            % If x, y, z are numeric or ':' call feval().  If x, y, z are
            % CHEBFUN, CHEBFUN2 or CHEBFUN3, concatenate (also checks that
            % domains are compatible) and call compose().
            if ( ( isnumeric(x) || strcmpi(x, ':') ) && ...
                    ( isnumeric(y) || strcmpi(y, ':') ) && ...
                    ( isnumeric(z) || strcmpi(z, ':') ) )
                out = feval(f, x, y, z);
            elseif ( isa(x, 'chebfun') && isa(y, 'chebfun') && isa(z, 'chebfun') )
                out = compose([x, y, z], f);
            elseif ( isa(x, 'chebfun2') && isa(y, 'chebfun2') && isa(z, 'chebfun2') )
                out = compose([x; y; z], f);
            elseif ( isa(x, 'chebfun3') && isa(y, 'chebfun3') && isa(z, 'chebfun3') )
                out = compose([x; y; z], f);
            else
                % Don't know what to do.
                error('CHEBFUN:CHEBFUN3:subsref:inputs3', ...
                    'Unrecognized inputs.')
            end
            varargout = {out};
            
        elseif ( isa(idx{1}, 'chebfun') || isa(idx{1}, 'chebfun2v') || ...
                isa(idx{1}, 'chebfun3v') || isa(idx{1}, 'spherefunv') )
            % Composition with a CHEBFUN, CHEBFUN2V, CHEBFUN3V, or SPHEREFUN:
            out = compose(idx{1}, f);
            varargout = {out};
            
        else
            error('CHEBFUN:CHEBFUN3:subsref:inputs', ...
                'Can only evaluate at triples (X,Y,Z), a CHEBFUN with 3 columns, a CHEBFUN2V or a CHEBFUN3V.')
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