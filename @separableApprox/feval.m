function out = feval(f, x, y)
%FEVAL  Evaluate a SEPARABLEAPPROX at one or more points.
%   FEVAL(F,X,Y) evaluates the SEPARABLEAPPROX F and the point(s) in (X,Y), where X and
%   Y are doubles.
%
%   FEVAL(F,X) evaluates the SEPARABLEAPPROX F along the complex valued CHEBFUN X and
%   returns  g(t) = F(real(X(t)),imag(X(t)))
%
%   FEVAL(F,X,Y) returns g(t) = F(X(t),Y(t)), where X and Y are real valued
%   CHEBFUN objects with the same domain.
%
% See also SUBSREF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    out = [];
    return
end

% Get the low rank representation for f.
[cols, D, rows] = cdr(f);

if ( strcmpi(x, ':') && strcmpi(y, ':') )    % f(:, :)
    % Simply return the SEPARABLEAPPROX:
    out = f;
    
elseif ( strcmpi(x, ':') && isnumeric( y ) ) % f(:, y)
    % Make evaluation points a vector.
    y = y(:);
    % Evaluate (returns a column chebfun):
    out = feval( cols, y ) * D * rows.';
    % Simplify:
    out = simplify( out, [], 'globaltol' );
    
elseif ( isnumeric( x ) && strcmpi(y, ':') ) % f(x, :)
    % Make evaluation points a vector.
    x = x( : );
    % Evaluate (returns a row chebfun):
    out = cols * D * feval( rows, x ).';
    % Simplify:
    out = simplify( out, [], 'globaltol' );
    
elseif ( isnumeric(x) && isnumeric(y) )  % f(x, y)    
    if ( ndims(x) >= 3 && isequal(size(x), size(y)) )
        % x and y are tensors. Call from CHEBFUN3.
        sizeX = size(x);
        x = x(:);
        y = y(:);
        out = feval(f, x, y);
        %% RESHAPE FOR OUTPUT:
        out = reshape(out, sizeX);
        return
    end
   
    takeTranspose = 0;
    
    % If the evaluation points are derived from meshgrid, then there is a
    % fast way to evaluate a separableApprox. Check for this property.
    if ( min(size(x)) > 1 && all(size(x) == size(y)) && numel(size(x)) == 2 )
        % Check to see if the input is a meshgrid:
        if ( max(max(abs(bsxfun(@minus, x, x(1,:))))) <= 10*eps  && ...
                max(max(abs(bsxfun(@minus, y, y(:,1))))) <= 10*eps )
            % This allows someone to do:
            % [xx,yy] = meshgrid(linspace(-1,1));
            % f(xx,yy)
            
            x = x(1,:);
            y = y(:,1);
            
        elseif ( max(max(abs(bsxfun(@minus, y, y(1,:))))) <= 10*eps && ...
                max(max(abs(bsxfun(@minus, x, x(:,1))))) <= 10*eps )
            % This allows someone to do:
            % [yy,xx] = meshgrid(linspace(-1,1));
            % f(xx,yy)
            
            x = x(:,1);
            y = y(1,:);
            takeTranspose = 1;
        else
            % Evaluate at matrices, but they're not from meshgrid:
            [m,n] = size( x );
            out = zeros( m, n );
            % Unroll the loop that is the longest
            if m > n
                for ii = 1:n
                    out(:,ii) = dot(feval(cols, ...
                        y(:,ii))*D,feval(rows, x(:,ii)),2);
                end
            else
                for jj = 1:m
                    out(jj,:) = dot(feval(cols, ...
                        y(jj,:).')*D,feval(rows, x(jj,:).'),2);
                end
            end
            return
        end
    else
    end
    
    % Evaluate:
    if ( isvector(x) && ~isscalar(x) && all(size(x) == size(y)) )
        % Determine whether inputs are pure vectors.
        out = feval(cols, y(:)) .* feval(rows, x(:)) * diag(D);
    else
        out = feval(cols, y(:)) * D * feval(rows, x(:)).';
    end
    
    % Take transpose:
    if ( takeTranspose )
        out = transpose( out );
    end
    
elseif ( isa(x, 'chebfun') )
    if ( min( size( x ) ) > 1 )
        error('CHEBFUN:SEPARABLEAPPROX:feval:arrayValued', ...
            'Cannot evaluate a SEPARABLEAPPROX at an array-valued CHEBFUN.');
    end
    
    if ( ~isreal(x) )      % Complex valued chebfun.
        % Extract chebfun along the path,  F(real(X(t)),imag(X(t)))
        out = chebfun(@(t) feval(f, real(x(t))', imag(x(t))'), x.domain, 'vectorize' );
    elseif ( isa(y, 'chebfun') )
        if ( isreal( y ) ) % Both x and y are real valued.
            % Check domains of x and y match
            if( domainCheck(x, y) )
                out = chebfun( @(t) feval(f, x(t), y(t)), x.domain, 'vectorize');
            else
                error('CHEBFUN:SEPARABLEAPPROX:feval:path', ...
                    'CHEBFUN path has domain inconsistency.');
            end
        else
            error('CHEBFUN:SEPARABLEAPPROX:feval:complex', ...
                'Cannot evaluate along complex-valued CHEBFUN.');
        end
    end
elseif ( isa(x, 'chebfun2v') )
    components = x.components;
    % [TODO]: Check domain and range are compatible?
    domain = components{1}.domain;
    out = f.constructor(@(s,t) feval(f, feval(components{1},s,t),...
        feval(components{2},s,t)), domain);
else
    error('CHEBFUN:SEPARABLEAPPROX:feval:inputs', ...
        'Unrecognized arguments for evaluation.');
    
end

end