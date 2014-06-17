function out = feval(f, x, y)
%FEVAL  Evaluate a CHEBFUN2 at one or more points.
%   FEVAL(F,X,Y) evaluates the CHEBFUN2 F and the point(s) in (X,Y), where X and
%   Y are doubles.
%
%   FEVAL(F,X) evaluates the CHEBFUN2 F along the complex valued chebfun X and
%   returns  g(t) = F(real(X(t)),imag(X(t)))
%
%   FEVAL(F,X,Y) returns g(t) = F(X(t),Y(t)), where X and Y are real valued
%   CHEBFUN objects with the same domain.
%
% See also SUBSREF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    out = [];
    return
end

% Get the low rank representation for f.
[cols, D, rows] = cdr(f);

if ( strcmpi(x, ':') && strcmpi(y, ':') )    % f(:, :)
    % Simply return the CHEBFUN2:
    out = f;
    
elseif ( strcmpi(x, ':') && isnumeric( y ) ) % f(:, y)
    % Make evaluation points a vector.
    y = y(:);
    % Evaluate (returns a column chebfun):
    out = feval( cols, y ) * D * rows.'; 
    % Simplify:
    out = simplify( out );
    
elseif ( isnumeric( x ) && strcmpi(y, ':') ) % f(x, :)
    % Make evaluation points a vector.
    x = x( : );
    % Evaluate (returns a row chebfun):
    out = cols * D * feval( rows, x ).'; 
    % Simplify:
    out = simplify( out );
    
elseif ( isnumeric( x ) && isnumeric( y ) )  % f(x, y)
    
    takeTranspose = 0;
    
    % If the evaluation points are derived from meshgrid, then there is a
    % fast way to evaluate a chebfun2. Check for this property. 
    if ( min(size(x)) > 1 && all(size(x) == size(y)) )
        % Check to see if the input is a meshgrid:
        if ( max(max(abs(bsxfun(@minus, x, x(1,:))))) == 0  && ... 
                max(max(abs(bsxfun(@minus, y, y(:,1)))) == 0 ))
            % This allows someone to do: 
            % [xx,yy] = meshgrid(linspace(-1,1));
            % f(xx,yy)
            
            x = x(1,:);
            y = y(:,1);
            takeDiag = 0;
            
        elseif ( max(max(abs(bsxfun(@minus, y, y(1,:))))) == 0 && ... 
                max(max(abs(bsxfun(@minus, x, x(:,1)))) == 0 ) )
            % This allows someone to do: 
            % [yy,xx] = meshgrid(linspace(-1,1));
            % f(xx,yy)
            
            x = x(:,1); 
            y = y(1,:);
            takeTranspose = 1;
            takeDiag = 0;
        else
            % Evaluate at matrices, but they're not from meshgrid: 
            out = zeros( size( x ) ); 
            for jj = 1:size(out,1)
                for kk = 1: size(out,2)
                    out(jj, kk) = feval( f, x(jj,kk), y(jj,kk) ); 
                end
            end
            return
        end
    else 
        takeDiag = 1; 
    end
    
    % Evaluate:
    out = feval(cols, y(:)) * D * feval(rows, x(:)).';
    
    % Take transpose: 
    if ( takeTranspose ) 
        out = transpose( out ); 
    end 
    
    % Take diagonal:
    if ( takeDiag ) 
        out = diag( out ); 
    end 
    
elseif ( isa(x, 'chebfun') )
    if ( min( size( x ) ) > 1 )
        error('CHEBFUN:CHEBFUN2:feval:arrayValued', ...
            'Cannot evaluate a CHEBFUN2 at an array-valued CHEBFUN.');
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
                error('CHEBFUN:CHEBFUN2:feval:path', ...
                    'CHEBFUN path has domain inconsistency.');
            end
        else
            error('CHEBFUN:CHEBFUN2:feval:complex', ...
                'Cannot evaluate along complex-valued CHEBFUN.');
        end
    end
elseif ( isa(x, 'chebfun2v') ) 
    components = x.components; 
    % [TODO]: Check domain and range are compatible? 
    domain = components{1}.domain;
    out = chebfun2(@(s,t) feval(f, feval(components{1},s,t),...
                                        feval(components{2},s,t)), domain);
else
    error('CHEBFUN:CHEBFUN2:feval:inputs', ...
        'Unrecognized arguments for evaluation.');
    
end

end
