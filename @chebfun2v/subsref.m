function varargout = subsref(F, ref)
%SUBSREF   CHEBFUN2V subsref.
%
% ( )
%   F(X,Y) returns the values of the CHEBFUN2 F evaluated on the array (X,Y).
%
%   F(k) returns the first component of F if k=1, the second if k=2, and
%   the third if k=3.
%
%   F(G) computes the composition with a CHEBFUN with two columns, a CHEBFUN2V
%   or a CHEBFUN3V with two components, or a DISKFUNV.  If G is a CHEBFUN with
%   one column, a CHEBFUN2, a CHEBFUN3, a DISKFUN or a SPHEREFUN, F(G) is
%   interpreted as F(real(G), imag(G)), regardless whether G is real or complex.
%
%   F(X, Y) where X and Y are CHEBFUN2 objects returns the CHEBFUN2V
%   representing the composition.  If X and Y are CHEBFUN3 objects, then
%   F(X, Y) is a CHEBFUN3V.  If X and Y are CHEBFUNs, then F(X, Y) is a
%   CHEBFUN.
%
%  .
%   F.PROP returns the property PROP of F as defined by GET(F,'PROP').
%
% { }
%    Throws an error.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% check for empty CHEBFUN2V object.
if ( isempty(F) )
    varargout = {[]};
    return
end

% Recursive going through the ref structure:
indx = ref(1).subs;

switch ( ref(1).type )
    
    case '.'
        if ( numel(ref) == 1 )
            % This is a get call to get a property.
            varargout = {get(F, indx)};
        else
            %
            t2 = ref(2).type;
            if ( strcmp(t2,'.') )
                out = get(F, indx, ref(2).subs{:});
            else
                out = get(F, indx);
                out = out( ref(2).subs{:} );
            end
            if ( numel(ref) > 2 )
                varargout = {subsref(out, ref(3:end))};
            else
                varargout = {out};
            end
        end
        
    case '()'
        if ( length(indx) > 1 )
            % where to evaluate
            x = indx{1};
            y = indx{2};
            % If x, y are numeric or ':' call feval().  If x, y are CHEBFUN,
            % CHEBFUN2 or CHEBFUN3, concatenate (also checks that domains are
            % compatible) and call compose().
            if ( ( isnumeric(x) || strcmpi(x, ':') ) && ...
                    ( isnumeric(y) || strcmpi(y, ':') ) )
                out = feval(F, x, y);
            elseif ( isa(x, 'chebfun') && isa(y, 'chebfun') )
                out = compose([ x, y ], F);
            elseif ( isa(x, 'chebfun2') && isa(y, 'chebfun2') )
                out = compose([ x; y ], F);
            elseif ( isa(x, 'chebfun3') && isa(y, 'chebfun3') )
                out = compose([ x; y ], F);
            else
                % Don't know what to do.
                error('CHEBFUN:CHEBFUN2V:subsref:input', ...
                    'Unrecognized input.')
            end
            varargout = {out};
            
        elseif ( isa(indx{1}, 'chebfun') )
            x = indx{1};
            % Composition F(X(t)):
            if ( size(x, 2) == 1 )
                % Interpret as F(real(X(t)), iag(X(t))):
                out = compose([ real(x), imag(x) ], F);
            elseif ( size(x, 2) == 2 )
                % Composition F([X(t), Y(t)]):
                out = compose([ x(:,1), x(:,2) ], F);
            else
                error('CHEBFUN:CHEBFUN2:subsref:ChebfunSize', ...
                    'Can compose only with a CHEBFUN with values in C or R^2.')
            end
            varargout = {out};
            
        elseif ( isa(indx{1}, 'chebfun2') || isa(indx{1}, 'chebfun3') || ...
            isa(indx{1}, 'diskfun') || isa(indx{1}, 'spherefun') )
            % Composition F(G) where G is a CHEBFUN2, CHEBFUN3, DISKFUN or
            % SPHEREFUN, interpreted as F(real(G), imag(G)), regardless whether
            % G is real or complex.
            out = compose(indx{1}, F);
            varargout = {out};
            
        elseif ( isa(indx{1}, 'chebfun2v') || isa(indx{1}, 'chebfun3v') || ...
                isa(indx{1}, 'diskfunv') )
            % Composition F(CHEBFUN2V), F(CHEBFUN3V) or F(DISKFUNV):
            out = compose(indx{1}, F);
            varargout = {out};
            
        else
            if ( isa(indx{1}, 'double') )
                if all( indx{1} == 1  )
                    varargout = F.components(1);
                elseif ( all( indx{1} == 2 ) )
                    varargout = F.components(2);
                elseif ( ( all(indx{1} == 3) ) && ( ~isempty(F.components(3)) ) )
                    varargout = F.components(3);
                else
                    error('CHEBFUN:CHEBFUN2V:subsref:index', ...
                        'CHEBFUN2V only contains two/three components');
                end
            end
        end
        
    otherwise
        error('CHEBFUN:CHEBFUN2V:subsref:unexpectedType', ...
            ['??? Unexpected index.type of ' index(1).type]);
        
end

% Recurse down:
if ( numel(ref) > 1 )
    ref(1) = [];
    varargout = { subsref( varargout{ : }, ref ) };
end

end
