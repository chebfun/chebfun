function u = mldivide(N, f, varargin )
% Solve partial differential equation using a tensor products.
% The only ACA stuff going on here is in the determining the discretisation
% size required to approximate the solution.

% Try and make a chebfun2 out of right hand side.
if isa(f,'chebfun')
    f = @(x,y) f(x) + 0*y;
    warning('Chebop2:mldivide:rhs','Univariate righthand side.');
end
if isa(f,'double')
    f = @(x,y) f + 0*x;
end
rect = N.domain; f = chebfun2(f, rect);

prefs = chebfunpref(); 
tol = max(prefs.cheb2Prefs.eps,1e-14); 
maxDiscretise = 2*prefs.cheb2Prefs.maxRank;
minsample = 9; 

if ( nargin == 3 && isa(varargin{1},'double') )
    % Nonadaptive solve with an m-by-m discretization
    m = varargin{1};
    n = m;
    if ( abs( round(n) - n ) > 0 )
        error('CHEBOP2:MLDIVIDE:DISC','Discretization size should be an integer');
    end
    X = chebop2.denseSolve(N, f, m, n);
elseif ( nargin == 4 && isa(varargin{1},'double') && ~isinf(varargin{1}) ...
        && isa(varargin{2},'double') && ~isinf(varargin{2}) )
    % Nonadaptive solve with an m-by-n discretization
    m = varargin{1};
    n = varargin{2};
    if ( abs( round(n) - n ) > 0 || abs( round(m) - m ) > 0)
        error('CHEBOP2:MLDIVIDE:DISC','Discretization size should be an integer');
    end
    X = chebop2.denseSolve(N, f, m, n);
elseif ( nargin == 4 && isa(varargin{1},'double') && ~isinf(varargin{1})...
        && isinf(varargin{2}) )
    % Adaptive solve in the horizontal variable, nonadaptive in the other.
    % (An m-by-inf discretization.)
    m = varargin{1};
    n = minsample;
    
    Resolved = 0;
    while ( ( ~Resolved ) && ( max(n) < maxDiscretise ))
        X = chebop2.denseSolve(N, f, m, n);               % Solve PDE on a m x n grid
        
        if ( ~Resolved )                       % Resolved in vertical direction?
            clcfs = max(abs(X(:,end:-1:end-8)));
            Resolved = all( ( clcfs < 20*n*tol ) );
            n = 2^(floor(log2(n))+1)+1;
        end
        if ( max(m,n) > 250 ), tol = max(tol,1e-10); end % Increase tolerance on large grids.
        
        if any(isnan(X))
            error('Nonunique solution to PDE');
        end
        
    end
    
    if max(n) >= maxDiscretise
        warning('Maximum discretization size reached. Solution may not be resolved.')
    end
    
%     %TODO:
%     % cut trailing coefficients
%     idx = find( max(abs(X))/max(abs(X(:))) > tol, 1, 'last');
%     if ~isempty(idx)
%         X = X(:,1:idx);
%     end
    
elseif ( nargin == 4 && isa(varargin{2},'double') && ~isinf(varargin{2})...
        && isinf(varargin{1}) )
    % Adaptive solve in the vertical variable, nonadaptive in the other.
    % (An inf-by-n discretization.)
    n = varargin{2};
    m = minsample;
    
    Resolved = 0;
    while ( ( ~Resolved ) && ( max(m) < maxDiscretise ))
        X = chebop2.denseSolve(N, f, m, n);               % Solve PDE on a m x n grid
        
        if ( ~Resolved )                             % Resolved in x-direction
            rwcfs = max(abs(X(end:-1:end-8,:)));
            Resolved = all( ( rwcfs < 20*m*tol ) );
            m = 2^(floor(log2(m))+1)+1;
        end
        if ( max(m,n) > 250 ), tol = max(tol,1e-10); end % Increase tolerance on large grids.
        
        if any(isnan(X))
            error('Nonunique solution to PDE');
        end
        
    end
    
    if max(m) >= maxDiscretise
        warning('Maximum discretization size reached. Solution may not be resolved.')
    end
    
    %TODO:
%     % cut trailing coefficients
%     idy = find( max(abs(X),[],2)/max(abs(X(:))) > tol, 1, 'last');
%     if ~isempty(idy)
%         X = X(1:idy,:);
%     end
    
elseif ( nargin == 2 || (nargin == 4 && isinf(varargin{2}) && isinf(varargin{1})))
    n = minsample;
    m = minsample;
    
    Resolved_x = 0; Resolved_y = 0; Resolved = Resolved_x & Resolved_y;
    while ( ( ~Resolved ) && ( max(m, n) < maxDiscretise ))
        X = chebop2.denseSolve(N, f, m, n);               % Solve PDE on a m x n grid
        
        if ( ~Resolved_y )                             % Resolved in y-direction
            clcfs = max(abs(X(:,end:-1:end-8)));
            Resolved_y = all( ( clcfs < 20*m*tol ) );
            m = 2^(floor(log2(m))+1)+1;
        end
        if ( ~Resolved_x )                             % Resolved in x-direction
            rwcfs = max(abs(X(end:-1:end-8,:)));
            Resolved_x = all( ( rwcfs < 20*n*tol ) );
            n = 2^(floor(log2(n))+1)+1;
        end
        if ( max(m,n) > 250 ), tol = max(tol,1e-10); end % Increase tolerance on large grids.
        Resolved = Resolved_x & Resolved_y;
        
        if any(isnan(X))
            error('Nonunique solution to PDE');
        end
        
    end
    
    if max(m,n) >= maxDiscretise
        warning('Maximum discretization size reached. Solution may not accurate.')
    end
     
    %TODO:
%     % cut trailing coefficients
%     idy = find( max(abs(X))/max(abs(X(:))) > eps, 1, 'last');
%     idx = find( max(abs(X),[],2)/max(abs(X(:))) > eps, 1, 'last');
%     if isempty(idx) && ~isempty(idy)
%         X = X(1:idy,:);
%     elseif ~isempty(idx) && isempty(idy)
%         X = X(:,1:idx);
%     else
%         X = X(1:idy,1:idx);
%     end
    
else
    error('CHEBOP2:MLDIVIDE:SYNTAX','Unrecognized input syntax')
end

% Form a chebfun2 object.
u = chebfun2( polyval( rot90( X, 2 ) ), rect);
end


function X = ImposeBoundaryConditions(X, bb, gg, Px, Py, m, n)

% Recombine in the boundary conditions.
cs = size(bb{3},2) + size(bb{4},2);
rs = size(bb{1},2) + size(bb{2},2);

By = [bb{3}.' ; bb{4}.'].'; Gy = [gg{3}.' ; gg{4}.'].';
if ( ~isempty(By) )
    By = Py.' * By; % Gy = Gy * Py;
    X12 = By(1:cs,:).' \ (Gy(rs+1:size(X,2)+rs,:).' - By(cs+1:size(X,1)+cs,:).'*X);
    X = [ X12; X ];
end

Bx = [bb{1}.' ; bb{2}.'].'; Gx = [gg{1}.' ; gg{2}.'].';
if ( ~isempty(Bx) )
    Bx = Px.' * Bx; %Gx = Px.' * Gx;
    X2 = ( Bx(1:rs,:).' \ ( Gx(1:size(X,1),:).' - Bx(rs+1:size(X,2)+rs,:).'*X.') ).' ;
    X = [X2 X];
end

if ( size(X,1) < m )
    X(size(X,1)+1:m, :) = 0;
end

if ( size(X,2) < n )
    X(:, size(X,2)+1:n) = 0;
end

if ( ~isempty(Px) )
    X = X * Px.';
end

if ( ~isempty(Py) )
    X = Py * X;
end

end


function c = polyval(c)
% Convert bivariate chebyshev coeffs to values on a Chebyshev grid.
[sy, sx] = size(c);
for j = 1:sx
    c(:,j) = chebtech2.coeffs2vals(c(:,j));
end
for k = 1:sy
    c(k,:) = chebtech2.coeffs2vals(c(k,:).').';
end
end


function [Resolved, newDisc] = ResolveCheck( coeffs, oldDisc, tol )
% Basic Resolution check: 

Resolved = all( ( coeffs < 20*oldDisc*tol ) );
newDisc = 2^(floor(log2(oldDisc))+1)+1;

end
