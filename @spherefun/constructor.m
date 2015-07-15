function g = constructor( g, op, dom, varargin )
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

if ( nargin == 0 )          % SPHEREFUN( )
    return
end

if ( isa(op, 'spherefun') )  % SPHEREFUN( SPHEREFUN )
    g = op;
    return
end

% If domain is empty take it to be co-latitude (doubled up).
if ( nargin < 3 || isempty(dom) )
    dom = [-pi pi 0 pi];
end

if ~isempty(varargin)
    alpha = varargin{1};
else
    alpha = 2;
end

% TODO: Should we allow any other domains that latitude and co-latitude?

if ( isa(op, 'spherefun') )  % SPHEREFUN( SPHEREFUN )
    g = op;
    return
end

% TODO: 
% 1. Need to allow for different domains.
% 2. Add support for preferences
% 3. Add non-adaptive construction
% 4. Add fixed-rank.
% 5. Add tensor-product.

maxRank = 8192; 
maxSample = 8192;
pseudoLevel = eps;

% If f is defined in terms of x,y,z; then convert it to
% (longitude,latitude).
h = redefine_function_handle( op );

% PHASE ONE  
% Sample at square grids, determine the numerical rank of the
% function.
n = 4;
happy_rank = 0;     % Happy with phase one? 
failure = false;
pivotIndices = [];
pivotMatrices = [];
while ( ~happy_rank && ~failure )
    n = 2*n;
        
    % Sample function on a tensor product grid.
    % TODO: Add a more sophisticated evaluate function that does
    % vectorization like chebfun2.
    F = evaluate(h, n, n, dom);
    
    tol = GetTol(F, pi/n, pi/n, dom, pseudoLevel);
    
    pivotIndices2 = pivotIndices;
    pivotMatrices2 = pivotMatrices;
    [ pivotIndices, pivotMatrices, happy_rank, removePoles ] = PhaseOne( F, tol, alpha );
    if size(pivotIndices,1) > size(pivotIndices2,1)
        happy_rank = 0;
    elseif happy_rank
        pivotIndices = pivotIndices2;
        pivotMatrices = pivotMatrices2;
        n = n/2;
    end
    
    if ( n >= maxRank  )
        warning('SPHEREFUN:CONSTRUCTOR:MAXRANK', ... 
                                'Unresolved with maximum rank.');
        failure = true;
    end
end

% PHASE TWO 
% Find the appropriate discretizations in the columns and rows. 
[cols, pivots, rows, pivotLocations, idxPlus, idxMinus] = PhaseTwo( h, pivotIndices, pivotMatrices, n, dom, tol, maxSample, removePoles, alpha );

g.cols = chebfun( cols, dom(3:4)-[pi 0], 'trig');
g.rows = chebfun( rows, dom(1:2), 'trig');
g.pivotValues = pivots;
g.pivotIndices = pivotIndices;
g.domain = dom;
g.idxPlus = idxPlus;
g.idxMinus = idxMinus;

% Adjust the pivot locations so that they correspond to 
% -pi < lam < pi and 0 < th < pi or -pi/2 < th < pi/2
if iscolat(g)
    pivotLocations(:,2) = -pivotLocations(:,2);
else
    pivotLocations(:,2) = -(pivotLocations(:,2)+pi);
end
pivotLocations(:,1) = pivotLocations(:,1) + pi;
g.pivotLocations = pivotLocations;

% Sort according to the maginuted of the pivots using the partition and
% combine functions.
% [gp,gm] = partition(g);
% g = combine(gp,gm);

end

function [pivotIndices, pivotMatrices, happy, removePole] = PhaseOne( F, tol, alpha )

% Phase 1: Go find rank, plus pivot locations, ignore cols and rows.
% alpha = spherefun.alpha; % get growth rate factor.
[m, n] = size( F );
pivotIndices = []; pivotMatrices = [];
vscl = norm( F( : ), inf);
rank_count = 0;    % keep track of the rank of the approximation.

B = F(:,1:n/2);    %% (1,1) Block of F.
C = F(:,n/2+1:n);  % (1,2) block of F.
Fp = 0.5*(B + C);
Fm = 0.5*(B - C);

%
% Deal with the poles by removing them from Fp.
%
pole1 = mean(Fp(1,:));     % Take the value at the poles to be 
pole2 = mean(Fp(m,:));     % the mean of all the samples at the poles.

% If the the values at both poles are not zero then we need to add zero
% them out before removing these entries from F.
removePole = false;
if abs(pole1) > vscl*tol || abs(pole2) > vscl*tol
    % Determine the column with maximum inf-norm
    [ignored, poleCol] = max(max(abs(Fp),[],1));
    % Zero out the pole using the poleCol.
    Fp = Fp - Fp(:, poleCol)*ones(1,n/2);
    % Do we need to keep track of the pivot locations?
    removePole = true;
    % Update the rank count
    rank_count = rank_count + 1;
end

% Remove the rows corresponding to the poles before determining the
% rank.  We do this because then F is an even BMC matrix, which is what the
% code below requires.
Fp = Fp( 2:m-1, : );
Fm = Fm( 2:m-1, : );

% while ( norm( F( : ), inf ) > tol )
% while ( max(norm(B(:),inf),norm(C(:),inf) ) > tol )
while ( max(norm(Fp(:),inf),norm(Fm(:),inf) ) > tol )    
    % Find pivots:
    S1 = max( abs(Fp), abs(Fm) );
    [ignored, idx] = max( S1(:) );
    [j, k] = myind2sub( size( S1 ), idx );
    
    % Use maximum of the Fp and Fm matrices for pivots
    ev = [Fp(j,k) Fm(j,k)];
    evp = ev(1);
    evm = ev(2);
    
    % Singular-values sorted by magnitude
    sv = [max(abs(ev)) min(abs(ev))];
        
    pivotIndices = [ pivotIndices ; j k];

    if ( sv(1) <= alpha*sv(2) )  % Theoretically, should be s1 <= 2*s2.
        % Calculate inverse of pivot matrix:
        Fp = Fp - Fp(:,k)*(Fp(j,:)*(1/evp));
        Fm = Fm - Fm(:,k)*(Fm(j,:)*(1/evm));
        rank_count = rank_count + 2; 
        pivotMatrices = [pivotMatrices ; ev ];
    else
        % Calculate pseudoinverse of pivot matrix, there is
        % no full rank pivot matrix:
        if abs(evp) > abs(evm)
            Fp = Fp - Fp(:,k)*(Fp(j,:)*(1/evp));
            ev(2) = 0;
        else
            Fm = Fm - Fm(:,k)*(Fm(j,:)*(1/evm));
            ev(1) = 0;
        end
        pivotMatrices = [pivotMatrices ; ev ];
        rank_count = rank_count + 1; 
    end
end

% Adjust the pivot locations so that they now correspond to F having
% the poles.
if ~isempty( pivotIndices )
    pivotIndices(:,1) = pivotIndices(:,1) + 1;
end

% Put the poles at the begining the pivot locations array and also include
% the pivot matrix.
if removePole
    pivotIndices = [ 1 poleCol; pivotIndices ];
    ev = [1 0];
    pivotMatrices = [ev ; pivotMatrices];
end

% If the rank of the matrix is less than 1/8 its size. We are happy:
if ( rank_count < min(size(F))/8 )
    happy = 1;
else
    happy = 0;
end

end

function [cols, pivots, rows, pivotLocations, idxPlus, idxMinus] = PhaseTwo( h, pivotIndices, pivotMatrices, n, dom, tol, maxSample, removePoles, alpha )

% alpha = spherefun.alpha; % get growth rate factor.
happy_columns = 0;   % Not happy, until proven otherwise.
happy_rows = 0;
m = n;

[x, y] = getPoints( m, n, dom );

rk = size( pivotIndices, 1);
id_rows = pivotIndices(:,1); id_cols = pivotIndices(:,2);

row_pivots = y(id_rows);
col_pivots = x(id_cols);

% % Need to also include id_cols+n to account for the entries in the C
% % block
% id_cols = reshape([id_cols id_cols+n].',[],1);
% 
% col_pivots = x(id_cols);
% row_pivots = y(id_rows);

numPosPivots = sum( abs( pivotMatrices(:,1) ) > 0 );
numMinusPivots = sum( abs( pivotMatrices(:,2) ) > 0 );
totalPivots = numPosPivots + numMinusPivots;

pivotPlus = zeros(numPosPivots,1);
pivotMinus = zeros(numMinusPivots,1);
pivots = zeros(totalPivots,1);
idxPlus = zeros(numPosPivots,1);
idxMinus = zeros(numMinusPivots,1);

% Phase 2: Calculate decomposition on sphere.
failure = false;
while ( ~(happy_columns && happy_rows) && ~failure)
    
    [x, y] = getPoints( m, n, dom );
    [xx, yy] = meshgrid( col_pivots, y);
    newCols = h( xx, yy ); temp = h( xx + pi, yy );
    newColsPlus = 0.5*(newCols + temp);
    newColsMinus = 0.5*(newCols - temp);
    
    [xx, yy] = meshgrid( x, row_pivots );
    newRows = h( xx, yy );

    % This code will be unnecessary once ticket #1532 is addressed on the
    % chebfun tracker.  Don't forget to remove it.
    if numel(row_pivots) == 1
        newRows = newRows(:).';
    end
    
    newRowsPlus = 0.5*( newRows(:,1:n) + newRows(:,n+1:2*n) );
    newRowsMinus = 0.5*( newRows(:,1:n) - newRows(:,n+1:2*n) );
    
    
    colsPlus = zeros( m+1, numPosPivots );
    colsMinus = zeros( m+1, numMinusPivots );
    rowsPlus = zeros( numPosPivots, n );
    rowsMinus = zeros( numMinusPivots, n );
    plusCount = 1;
    minusCount = 1;
    pivotCount = 1;

    % Need to remove pole, which means we use the column with the largest
    % max norm (repeated) with rows of all ones in the elimination
    % algorithm.
    if removePoles
        newRowsPlus(1,:) = 1;
    end
    
    
    for ii = 1:rk
        
        % Get the eigenvalues of the pivot matrix M
        ev = pivotMatrices( ii, : );
        s = [max(abs(ev)) min(abs(ev))];
                                
        % Do GE step on both matrices
        if ( s(1) <= alpha*s(2) )
            colPlus = newColsPlus(:,ii);
            rowPlus = newRowsPlus(ii,:);
            
            colMinus = newColsMinus(:,ii);
            rowMinus = newRowsMinus(ii,:);

            % Store the columns and rows.
            colsPlus(:,plusCount) = colPlus;
            rowsPlus(plusCount,:) = rowPlus;
            pivotPlus(plusCount) = ev(1);
            idxPlus(plusCount) = pivotCount;
            pivots(pivotCount) = ev(1);
            
            plusCount = plusCount + 1;
            pivotCount = pivotCount + 1;

            colsMinus(:,minusCount) = colMinus;
            rowsMinus(minusCount,:) = rowMinus;
            pivotMinus(minusCount) = ev(2);
            idxMinus(minusCount) = pivotCount;
            pivots(pivotCount) = ev(2);

            minusCount = minusCount + 1;
            pivotCount = pivotCount + 1;

            newColsPlus = newColsPlus - ...
                    colPlus*(rowPlus(id_cols)*(1/ev(1)));
            newRowsPlus = newRowsPlus - ...
                    ((1/ev(1))*colPlus(id_rows))*rowPlus;
            newColsMinus = newColsMinus - ...
                    colMinus*(rowMinus(id_cols)*(1/ev(2)));
            newRowsMinus = newRowsMinus - ...
                    ((1/ev(2))*colMinus(id_rows))*rowMinus;
        else
            % Use the pseudoinverse of the pivot matrix, there is
            % no full rank pivot matrix:
            if abs(ev(1)) > abs(ev(2))                
                colPlus = newColsPlus(:,ii);
                rowPlus = newRowsPlus(ii,:);

                % Store the columns and rows.
                colsPlus(:,plusCount) = colPlus;
                rowsPlus(plusCount,:) = rowPlus;
                pivotPlus(plusCount) = ev(1);
                idxPlus(plusCount) = pivotCount;
                pivots(pivotCount) = ev(1);

                plusCount = plusCount + 1;
                pivotCount = pivotCount + 1;
                
                newColsPlus = newColsPlus - ...
                        colPlus*(rowPlus(id_cols)*(1/ev(1)));
                newRowsPlus = newRowsPlus - ...
                        ((1/ev(1))*colPlus(id_rows))*rowPlus;
            else
                colMinus = newColsMinus(:,ii);
                rowMinus = newRowsMinus(ii,:);

                % Store the columns and rows.
                colsMinus(:,minusCount) = colMinus;
                rowsMinus(minusCount,:) = rowMinus;
                pivotMinus(minusCount) = ev(2);
                idxMinus(minusCount) = pivotCount;
                pivots(pivotCount) = ev(2);

                minusCount = minusCount + 1;
                pivotCount = pivotCount + 1;

                newColsMinus = newColsMinus - ...
                        colMinus*(rowMinus(id_cols)*(1/ev(2)));
                newRowsMinus = newRowsMinus - ...
                        ((1/ev(2))*colMinus(id_rows))*rowMinus;
            end
        end
    end    
    % Happiness check for columns:
    % TODO: Make this more similar to hapiness check in trigtech.

    % Double up the columns.
    temp1 = sum([colsPlus colsMinus],2); temp2 = sum([colsPlus -colsMinus],2);
    col_coeffs = trigtech.vals2coeffs( [temp1;temp2(m:-1:2)] );

    % Length of tail to test.
    testLength = min(m, max(3, round((m-1)/8)));
    tail = col_coeffs(1:testLength);

    if ( all( abs( tail ) <= 1e1*tol ) )
        happy_columns = 1;
    end
    
    % Happiness check for rows:
    % TODO: Make this more similar to hapiness check in trigtech.

    % Double up the rows.
    temp1 = sum([rowsPlus; rowsMinus],1); temp2 = sum([rowsPlus; -rowsMinus],1);
    row_coeffs = trigtech.vals2coeffs( [temp1 temp2].' );

    % Length of tail to test.
    testLength = min(n, max(3, round((n-1)/8)));
    tail = row_coeffs(1:testLength);
    if ( all( abs( tail ) <= 1e1*tol ) )
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
    
    if ( max(m, n) >= maxSample ) 
        warning('SPHEREFUN:constructor:notResolved', ...
        'Unresolved with maximum length: %u.', maxSample);
        failure = true;
    end 
end

% Combine the types of pivots and set-up indices to track them
cols = zeros( 2*size(colsPlus,1)-2, totalPivots );
cols(:,idxPlus) = [ colsPlus; flipud(colsPlus(2:end-1,:)) ];
cols(:,idxMinus) = [ colsMinus; -flipud(colsMinus(2:end-1,:)) ];

rows = zeros( 2*size(rowsPlus,2), totalPivots );
rows(:,idxPlus) = [rowsPlus rowsPlus].';
rows(:,idxMinus) = [rowsMinus -rowsMinus].';

pivotLocations = [col_pivots row_pivots];

end

function F = evaluate( h, m, n, dom )
% Evaluate h on a m-by-n tensor product grid.

[x, y] = getPoints( m, n, dom );

% Tensor product grid
[xx,yy] = meshgrid(x, y);

% Evaluate h on the non-doubled up grid
F = h( xx, yy );

end


function [x, y] = getPoints( m, n, dom )

colat = [-pi pi 0 pi]; % Colatitude (doubled up)
lat = [-pi pi -pi/2 pi/2]; % Latitude (doubled up)

% Sample at an even number of points so that the poles are included.
if all( (dom-colat) == 0 )
    x = trigpts( 2*n, [-pi,pi] );   % azimuthal angle, lambda
    y = linspace( -pi, 0, m+1 ).';  % elevation angle, theta
%     y = linspace( 0, pi, m+1 ).';  % elevation angle, theta
elseif all( (dom-lat) == 0 )
    x = trigpts( 2*n, [-pi,pi] );           % azimuthal angle, lambda
    y = linspace( -3*pi/2, -pi/2, m+1 ).';  % elevation angle, theta
%     y = linspace( -pi/2, pi/2, m+1 ).';  % elevation angle, theta
else
    error('SPHEREFUN:constructor:points2D:unkownDomain', ...
        'Unrecognized domain.');
end

end

% function pinvM = getPseudoInv( M )
% lam1 = M(1,1)+M(1,2);  % Eigenvalues of M (which is symmetric)
% lam2 = M(1,1)-M(1,2);
% if abs(lam1) > abs(lam2)
%     pinvM = ones(2)/(2*lam1);
% else
%     pinvM = [[1 -1];[-1 1]]/(2*lam2);
% end
% end

function [row, col] = myind2sub(siz, ndx)
% Alex's version of ind2sub. In2sub is slow because it has a varargout. Since this
% is at the very inner part of the constructor and slowing things down we will
% make our own. This version is about 1000 times faster than MATLAB ind2sub.

vi = rem( ndx - 1, siz(1) ) + 1 ;
col = ( ndx - vi ) / siz(1) + 1;
row = ( vi - 1 ) + 1;

end


function f = redefine_function_handle( f )
% nargin( f ) = 2, then we are already on the sphere, if nargin( f ) = 3,
% then do change of variables:

if ( nargin( f ) == 3 )
    % Wrap f so it can be evaluated in spherical coordinates
    f = @(lam, th) spherefun.sphf2cartf(f,lam,th,0);
%     % Double g up.
%     f = @(lam, th) sph2torus(f,lam,th);
end

end

function tol = GetTol(F, hx, hy, dom, pseudoLevel)
% GETTOL     Calculate a tolerance for the spherefun constructor.
%
%  This is the 2D analogue of the tolerance employed in the trigtech
%  constructors. It is based on a finite difference approximation to the
%  gradient, the size of the approximation domain, the internal working
%  tolerance, and an arbitrary (2/3) exponent. 

[m, n] = size( F ); 
grid = max( m, n );

% Remove some edge values so that df_dx and df_dy have the same size. 
dfdx = diff(F(1:m-1,:),1,2) / hx; % xx diffs column-wise.
dfdy = diff(F(:,1:n-1),1,1) / hy; % yy diffs row-wise.
% An approximation for the norm of the gradient over the whole domain.
Jac_norm = max( max( abs(dfdx(:)), abs(dfdy(:)) ) );
vscale = max( abs( F(:) ) );
tol = grid.^(2/3) * max( abs( dom(:) ) ) * max( Jac_norm, vscale) * pseudoLevel;

end

% function f = redefine_function_handle_pole( f, poleColPivot )
% % Set f to f - f(poleColPivot,theta) where poleColPivot is the value of
% % lambda (the column) used to zero out the poles of f.
% 
% f = @(lam, th) f(lam,th) - f(poleColPivot,th);
% 
% end

