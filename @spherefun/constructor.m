function g = constructor( g, op, varargin )
%CONSTRUCTOR   The main SPHEREFUN constructor.
%
% The algorithm for constructing a SPHEREFUN comes in two phases:
%
% PHASE 1: The first phase attempts to determine the numerical rank of the
% function by performing Gaussian elimination with special 2x2 pivoting matrices 
% on a tensor grid of sample values. GE is perform until the sample matrix
% is approximated.  At the end of this stage we have candidate pivot locations
% and pivot elements.
%
% PHASE 2: The second phase attempts to resolve the corresponding column and row
% slices by sampling along the slices and performing GE (pivoting at 2x2 matrices) 
% on the skeleton.   Sampling along each slice is increased until the Fourier 
% coefficients of the slice fall below machine precision.

if ( nargin == 0 )          % SPHEREFUN( )
    return
end

if ( isa(op, 'spherefun') )  % SPHEREFUN( SPHEREFUN )
    g = op;
    return
end

tol = 50*eps;       % Tolerance
max_rank = 4000; 

% If f is defined in terms of x,y,z; then convert: 
h = redefine_function_handle( op );

% PHASE ONE  
% Sample at square grids, determine the numerical rank of the
% function.
n = 4;
happy_rank = 0;     % Happy with phase one? 
while ( ~happy_rank )
    n = 2*n;
    F = sample(h, n, n);
    [ PivotLocations, PivotMatrices, happy_rank ] = PhaseOne( F, tol );
    if ( n >= max_rank  )
        warning('SPHEREFUN:CONSTRUCTOR:MAXRANK', ... 
                                'Unresolved with maximum rank.')
    end
end

% PHASE TWO 
% Find the appropriate discretizations in the columns and rows. 
[Cols, BlockDiag, Rows] = PhaseTwo(h, PivotLocations, PivotMatrices, n, tol );

% Make a spherefun, we are done. 
g.Cols = trigtech( Cols );
g.Rows = trigtech( Rows.' );
g.BlockDiag = BlockDiag;
g.PivotLocations = PivotLocations;

end

function [PivotLocations, PivotMatrices, happy] = PhaseOne( F, tol )

% Phase 1: Go find rank, plus pivot locations, ignore cols and rows.
alpha = 50; 
[m, n] = size( F );
PivotLocations = []; PivotMatrices = [];
vscl = norm( F( : ), inf);
rank_count = 0;    % keep track of the rank of the approximation. 

while ( norm( F( : ), inf ) > tol*vscl )
    % Find pivot:
    % Calculate the maximum 1st singular value of all the special 2x2
    % submatrices.
    S1 = max(  abs(F - flipud(F)), abs(F + flipud(F)) );
    S1 = S1(1:m/2, 1:n/2);
    [ignored, idx] = max( S1(:) );
    [j, k] = ind2sub( size( S1 ), idx );
    
    % Max determinant submatrix is:
    M = [ F(j,k) F(j,k+n/2) ; F(j,k+n/2) F(j,k)];
    s = sort( abs([ diff(M(1,:)) ; sum(M(1,:)) ]), 1, 'descend' );  % equivalent to svd( M )
    
    PivotLocations = [ PivotLocations ; j k m-j+1 k+n/2];
    PivotMatrices = [PivotMatrices ; M ];
    
    if ( s(1) <= alpha*s(2) )  % Theoretically, should be s1 <= 2s2.
        % Calculate inverse of pivot matrix:
        F = F - F(:, [k k+n/2] ) * ( M \ F( [j m-j+1],: ) );
        rank_count = rank_count + 2; 
    else
        % Calculate pseudoinverse of pivot matrix, there is
        % no full rank pivot matrix:
        [U, S, V] = svd( M );
        S(1,1) = 1./S(1,1); S(2,2) = 0;
        invM = U * S * V';
        F = F - F(:, [k k+n/2] ) * ( invM *  F( [j m-j+1],: ) );
        rank_count = rank_count + 1; 
    end
    
end

% If the rank of the matrix is less than 1/4 its size. We are happy:
if ( rank_count < min(size(F))/8 )
    happy = 1;
else
    happy = 0;
end

end



function [Cols, BlockDiag, Rows] = PhaseTwo( h, PivotLocations, PivotMatrices, n, tol)

alpha = 50; 
happy_columns = 0;   % Not happy, until proven otherwise.
happy_rows = 0;
m = n;
[x, y] = getPoints( m, n );
rk = size(PivotLocations, 1);
BlockDiag = spalloc( 2*rk, 2*rk, 6*rk-2 );
id = PivotLocations'; id = id(:);
id_cols = id(2:2:end); id_rows = id(1:2:end);
col_pivots = x(id_cols);
row_pivots = y(id_rows);

% Phase 2: Calculate decomposition on sphere.

while ( ~happy_columns || ~happy_rows )
    
    [x, y] = getPoints( m, n );
    [xx, yy] = meshgrid( col_pivots, y);
    Cols_new = h( xx ,yy );
    [xx, yy] = meshgrid( x, row_pivots);
    Rows_new = h( xx, yy );
    Cols = zeros( 2*m, 2*rk );
    Rows = zeros( 2*rk, 2*n );
    for ii = 1:rk
        
        M = PivotMatrices( 2*ii-1:2*ii, :);
        s = sort( abs([ diff(M(1,:)) ; sum(M(1,:)) ]), 1, 'descend' );  % equivalent to svd( M )
        
        Cols(:,2*ii-1:2*ii) = Cols_new(:,2*ii-1:2*ii); %F(:, [k k+n] );
        Rows(2*ii-1:2*ii,:) = Rows_new(2*ii-1:2*ii,:); %F([j 2*m-j+1],: );
        if ( s(1) <= alpha*s(2) )
            % Calculate inverse of pivot matrix:
            row_correction = Cols_new(id_rows,2*ii-1:2*ii) * ( M \ Rows_new(2*ii-1:2*ii,:) );
            Cols_new = Cols_new - Cols_new(:,2*ii-1:2*ii) * ( M \ Rows_new(2*ii-1:2*ii,id_cols) );
            Rows_new = Rows_new - row_correction;
            BlockDiag(2*ii-1:2*ii, 2*ii-1:2*ii) = inv( M );
        else
            % Calculate pseudoinverse of pivot matrix, there is
            % no full rank pivot matrix:
            [U, S, V] = svd( M );
            S(1,1) = 1./S(1,1); S(2,2) = 0;
%             invM = U * S * V';
            invM = V * S * U';
            row_correction = Cols_new(id_rows,2*ii-1:2*ii) * ( invM * Rows_new(2*ii-1:2*ii,:) );
            Cols_new = Cols_new - Cols_new(:,2*ii-1:2*ii) * ( invM * Rows_new(2*ii-1:2*ii,id_cols) );
            Rows_new = Rows_new - row_correction;
            BlockDiag(2*ii-1:2*ii, 2*ii-1:2*ii) = invM;
        end
        
    end
    
    
    % Happiness check for columns:
    col_coeffs = trigtech.vals2coeffs( Cols ); 
    % GBW: This only tests the coefficients in the last column.
    % if ( all( abs( col_coeffs(end:-1:end-8) ) <= 1e2*tol*norm(Cols,inf) ) )
    % Below tests all the coefficients
    if ( all( abs( col_coeffs(end:-1:end-8,:) ) <= 1e2*tol*norm(Cols,inf) ) )
        happy_columns = 1;
    end
    
    % Happiness check for rows:
    row_coeffs = trigtech.vals2coeffs( Rows.' ); 
    % GBW: This only tests the coefficients in the last row.    
    % if ( all(abs( row_coeffs(end:-1:end-8) ) <= 1e2*tol*norm(Rows,inf)) )
    if ( all(abs( row_coeffs(end:-1:end-8,:) ) <= 1e2*tol*norm(Rows,inf)) )
        happy_rows = 1; 
    end
    
    % Adaptive:
    if( ~happy_columns )
        m = 2*m;
        ii = [1:2:m-1 m+2:2:2*m]; 
        id_rows = ii(id_rows); 
    end
    
    if ( ~happy_rows )       
        n = 2*n; 
        id_cols = 2*id_cols - 1;
    end
    
    if ( max(m, n) >= 4000 ) 
        error
    end 
    
end


end

function [x, y] = getPoints( m, n )

% x = linspace(-pi, pi, 2*n+1)';  x( end ) = [ ];
% GBW: You can't just remove the pole and keep everything else equally
% spaced between -pi/2 and 3*pi/2.  The issue is that you can't keep
% both point -pi/ and 3*pi/2.
% y = linspace(-pi/2, 3*pi/2, 2*m+1)'; y( m+1 ) = [ ];

x = trigpts(2*n,[-pi pi]);
% GBW: I believe we have to sample at equally spaced points shifted by h/2
% to not sample the poles and keep an even total number of points.
y = trigpts(2*m,[-pi/2 3*pi/2]);
y = y+0.5*pi/m; % Shift y by h/2 to avoid the poles
end

function F = sample( h, m, n )
[x, y] = getPoints( m, n );
[L2, T2] = meshgrid(x, y);
F = h(L2, T2);
end


function f = redefine_function_handle( f )
% nargin( f ) = 2, then we are already on the sphere, if nargin( f ) = 3,
% then do change of variables:

if ( nargin( f ) == 3 )
    % Wrap f so it can be evaluated in spherical coordinates
    f = @(lam, th) sphf2cartf(f,lam,th);
    % Double g up.
    f = @(lam, th) sph2torus(f,lam,th);
end

end


function fdf = sph2torus(f,lam,th)

fdf = real(f(lam,th));

id = th-pi/2 > 100*eps;

if ~isempty(id) && any(id(:))
    fdf(id) = f(lam(id)-pi,pi-th(id));
end

end

function fdf = sphf2cartf(f,lam,th)

x = cos(lam).*cos(th);
y = sin(lam).*cos(th);
z = sin(th);

fdf = f(x,y,z);

end