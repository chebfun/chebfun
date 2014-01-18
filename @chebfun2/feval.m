function out = feval(f, x, y)
%FEVAL  evaluate a chebfun2 at one or more points.
%
%  FEVAL(F,X,Y) evaluates the chebfun2 F and the point(s) in (X,Y), where
%    X and Y are doubles.
%
%  FEVAL(F,X) evaluates the chebfun2 F along the complex valued chebfun X
%    and returns  g(t) = F(real(X(t)),imag(X(t)))
%
%  FEVAL(F,X,Y) returns g(t) = F(X(t),Y(t)), where X and Y are real valued
%  chebfuns with the same domain.
%
% See also SUBSREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    varargout = {[]};
    return
end

% Get the low rank representation for f.
cols = f.cols;
rows = f.rows;
piv = f.pivotValues;
d = 1./piv;
d(d==inf) = 0;  % set infinite values to zero.


if ( strcmpi(x, ':') && strcmpi(y, ':') )    % f(:, :)
    out = f;
elseif ( strcmpi(x, ':') && isnumeric( y ) ) % f(:, y)
    % Make evaluation points a vector.
    y = y(:);
    % Evaluate (returns a column chebfun):
    out = feval( cols, y ) * diag( d ) * rows';
    % Simplify:
    out = simplify( out );
elseif ( isnumeric( x ) && strcmpi(y, ':') ) % f(x, :)
    % Make evaluation points a vector.
    x = x( : );
    % Evaluate (returns a row chebfun):
    out = Cols * diag( 1./pivotValues ) * feval( Rows, x )';
    % Simplify:
    out = simplify( out );
elseif ( isnumeric( x ) && isnumeric( y ) )    % f(x, y)
    
    sx = size(x);
    sy = size(y);
    
    if ( min(sx) > 1 && all(sx == sy) )
        if ( rank(x) == 1 && rank(y) == 1 )
            x = x(1,:);
            y = y(:,1);
        end
    end
    
    % Evaluate:
    out = feval( cols, y(:) ) * diag( d ) * feval( rows, x(:)) .';
elseif ( isa(x,'chebfun') )
    if ( ~isreal(x) ) % complex valued chebfun.
        % Extract chebfun along the path
        out = chebfun(@(t) feval(f, real(x(t))', imag(x(t))'), x.domain, 'vectorize' );  % F(real(X(t)),imag(X(t)))
    elseif ( isa(y, 'chebfun') )
        if ( isreal( y ) ) % both x and y are real valued.
            % check domains of x and y match
            if( domainCheck(x, y) )
                out = chebfun( @(t) feval(f, x(t), y(t)), x.domain, 'vectorize');
            else
                error('CHEBFUN2:feval:path','Chebfun path has domain inconsistency');
            end
        else
            error('CHEBFUN2:feval:complex','Cannot evaluate along complex-valued chebfun');
        end
    end
else
    error('CHEBFUN2:FEVAL:INPUTS','Unrecognized arguments for evaluation');
end

end