function spherefun( )
tol = 10*eps;
clc
n = 501;
m = 501;

% Define f in terms of x,y,z to make things easier.
% f = @(x,y,z) exp(-cos(pi*(x+y)));
f = @(x,y,z) exp(-cos(pi*(x+y+z)));
%f = @(x,y,z) exp(x + y + z);

h = redefine_function_handle( f );

[x, y] = getPoints( m, n );
[L2, T2] = meshgrid(x, y);
F = h(L2, T2);

[ PivotLocations, PivotMatrices ] = PhaseOne( F, tol );


[Cols, BlockDiag, Rows] = PhaseTwo(h, PivotLocations, PivotMatrices, tol );
%% TEST
x = linspace(-pi, pi, 2*n+1);  x( end ) = [ ];
y = linspace(-pi/2, 3*pi/2, 2*m+1); y( m+1 ) = [ ];
[L2, T2] = meshgrid(x, y);
F = h(L2, T2);
norm( F - Cols * BlockDiag * Rows, inf )

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


function [PivotLocations, PivotMatrices] = PhaseOne( F, tol )

%F = [A A ; flipud( A ) flipud( A ) ]
% Phase 1: Go find rank, plus pivot locations, ignore cols and rows.
[m, n] = size( F ); m = m/2; n = n/2;
PivotLocations = []; PivotMatrices = [];
vscl = norm( F( : ), inf);

while ( norm( F( : ), inf ) > tol*vscl )
    % Find pivot:
    % This find the maximum determinant of 2x2 submatrices of A:
    DET = F.^2 - flipud(F.^2);
    DET = DET(1:m,1:n);
    [mx, idx] = max( abs( DET(:) ) );
    [j, k] = ind2sub( size( DET ), idx );
    
    if ( abs(mx) < 1e8*tol )
        % Redo if all dets are zero:
        Fsub = F(1:m, 1:n);
        [ignored, idx] = max( abs( Fsub(:) ) );
        [j, k] = ind2sub( size( Fsub ), idx);
    end
    
    % Max determinant submatrix is:
    M = [ F(j,k) F(j,k+n) ; F(j,k+n) F(j,k)];
    mx = det( M );
    PivotLocations = [ PivotLocations ; j k 2*m-j+1 k+n];
    PivotMatrices = [PivotMatrices ; M ];
    
    
    if ( abs(mx) > 1e8*tol )
        % Calculate inverse of pivot matrix:
        F = F - F(:, [k k+n] ) * ( M \ F( [j 2*m-j+1],: ) );
    else
        % Calculate pseudoinverse of pivot matrix, there is
        % no full rank pivot matrix:
        [U, S, V] = svd( M );
        S(1,1) = 1./S(1,1); S(2,2) = 0;
        invM = U * S * V';
        F = F - F(:, [k k+n] ) * ( invM *  F( [j 2*m-j+1],: ) );
    end
end

end



function [Cols, BlockDiag, Rows] = PhaseTwo( h, PivotLocations, PivotMatrices, tol)
n = 501;
m = 501;
[x, y] = getPoints( m, n );
% Phase 2: Calculate decomposition on sphere.
rk = size(PivotLocations, 1);
BlockDiag = zeros( 2*rk );
id = PivotLocations'; id = id(:);
id_cols = id(2:2:end); id_rows = id(1:2:end);
[xx, yy] = meshgrid( x(id_cols), y);
Cols_new = h( xx ,yy );
[xx, yy] = meshgrid( x, y(id_rows));
Rows_new = h( xx, yy );
Cols = zeros( 2*m, 2*rk );
Rows = zeros( 2*rk, 2*n );
for ii = 1:rk
    
    M = PivotMatrices( 2*ii-1:2*ii, :);
    mx = det( M );
    
    Cols(:,2*ii-1:2*ii) = Cols_new(:,2*ii-1:2*ii); %F(:, [k k+n] );
    Rows(2*ii-1:2*ii,:) = Rows_new(2*ii-1:2*ii,:); %F([j 2*m-j+1],: );
    if ( abs(mx) > 1e8*tol )
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
        invM = U * S * V';
        row_correction = Cols_new(id_rows,2*ii-1:2*ii) * ( invM * Rows_new(2*ii-1:2*ii,:) );
        Cols_new = Cols_new - Cols_new(:,2*ii-1:2*ii) * ( invM * Rows_new(2*ii-1:2*ii,id_cols) );
        Rows_new = Rows_new - row_correction;
        BlockDiag(2*ii-1:2*ii, 2*ii-1:2*ii) = invM;
    end
end

end

function [x, y] = getPoints( m, n )

    x = linspace(-pi, pi, 2*n+1);  x( end ) = [ ];
    y = linspace(-pi/2, 3*pi/2, 2*m+1); y( m+1 ) = [ ];
end