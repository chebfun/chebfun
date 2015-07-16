function g = constructorPlusMinus( g, op, dom, varargin )
%CONSTRUCTOR   The main SPHEREFUN constructor.
%
% A smooth function on the doubled-up sphere has the splitting:
% [B    C]  = 0.5*[ (B+C)   (B+C)]  + 0.5*[  (B-C)    -(B-C)]
% [JC  JB]        [J(B+C)  J(B+C)]        [-J(B-C)    J(B-C)] 
%
% This version of the constructor applies Gaussian elimination with
% complete pivoting separately on the two functions on the right. The code
% exploits symmetries in the problem and the fact that the rank of the
% doubled up function is equal to rank(B+C) + rank(B-C) so that it is more
% efficient than directly calling the chebfun2 constructor on these
% functions.
%

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

% TODO: Should we allow any other domains than latitude and co-latitude?

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
happyRankPlus = 0;     % Happy with phase one? 
happyRankMinus = 0;
happyRank = 0;
failure = false;
pivotIndicesPlus = []; pivotsPlus = [];
pivotIndicesMinus = []; pivotsMinus = [];
while ( ~happyRank && ~failure )
    n = 2*n;
        
    % Sample function on a tensor product grid.
    % TODO: Add a more sophisticated evaluate function that does
    % vectorization like chebfun2.
    F = evaluate(h, n, n, dom);
    
    % Split into even and odd parts:
    B = F(:,1:n);
    C = F(:,n+1:2*n);
    Fplus = 0.5*(B + C);
    Fminus = 0.5*(B - C);

    tolplus = GetTol(Fplus, pi/n, pi/n, [-pi 0 0 pi], pseudoLevel);
    tolminus = GetTol(Fminus, pi/n, pi/n, [-pi 0 0 pi], pseudoLevel);
    tol = max(tolplus,tolminus);
  
    if ~happyRankPlus
        pivotIndicesPlus2 = pivotIndicesPlus;
        pivotsPlus2 = pivotsPlus;
        [ pivotIndicesPlus, pivotsPlus, happyRankPlus, removePoles ] = PhaseOnePlus( Fplus, tol );
        if size(pivotsPlus,1) > size(pivotsPlus2,1)
            happyRankPlus = 0;
            nplus = n;
        elseif happyRankPlus
            pivotIndicesPlus = pivotIndicesPlus2;
            pivotsPlus = pivotsPlus2;
            nplus = n/2;
        else
            nplus = n;
        end
    end
    
    if ~happyRankMinus
        pivotIndicesMinus2 = pivotIndicesMinus;
        pivotsMinus2 = pivotsMinus;
        
        [ pivotIndicesMinus, pivotsMinus, happyRankMinus] = PhaseOne( Fminus, tol, 0 );
        if size(pivotsMinus,1) > size(pivotsMinus2,1)
            happyRankMinus = 0;
            nminus = n;
        elseif happyRankMinus
            pivotIndicesMinus = pivotIndicesMinus2;
            pivotsMinus = pivotsMinus2;
            nminus = n/2;
        else
            nminus = n;
        end
    end
    
    happyRank = happyRankPlus && happyRankMinus;
    
    if ( nplus >= maxRank || nminus >= maxRank )
        warning('SPHEREFUN:CONSTRUCTOR:MAXRANK', ... 
                                'Unresolved with maximum rank.');
        failure = true;
    end
end

% PHASE TWO 
% Find the appropriate discretizations in the columns and rows.
[colsPlus, pivotsPlus, rowsPlus, pivotLocationsPlus ] = PhaseTwo( h, pivotIndicesPlus, pivotsPlus, nplus, dom, tol, maxSample, 1, removePoles );
[colsMinus, pivotsMinus, rowsMinus, pivotLocationsMinus ] = PhaseTwo( h, pivotIndicesMinus, pivotsMinus, nminus, dom, tol, maxSample, -1, 0);

g.cols = [chebfun( colsPlus, dom(3:4)-[pi 0], 'trig') chebfun( colsMinus, dom(3:4)-[pi 0], 'trig')];
g.rows = [chebfun( rowsPlus.', dom(1:2), 'trig') chebfun( rowsMinus.', dom(1:2), 'trig')];
g.pivotValues = [pivotsPlus;pivotsMinus];
g.pivotIndices = [pivotIndicesPlus; pivotIndicesMinus];
g.idxPlus = 1:length(pivotsPlus);
g.idxMinus = length(pivotsPlus)+1:length(g.pivotValues);
g.domain = dom;

% Adjust the pivot locations so that they correspond to 
% -pi < lam < pi and 0 < th < pi or -pi/2 < th < pi/2
pivotLocations = [pivotLocationsPlus;pivotLocationsMinus];
if iscolat(g)
    pivotLocations(:,2) = -pivotLocations(:,2);
else
    pivotLocations(:,2) = -(pivotLocations(:,2)+pi);
end
pivotLocations(:,1) = pivotLocations(:,1) + pi;
g.pivotLocations = pivotLocations;

% Sort according to the maginuted of the pivots using the partition and
% combine functions.
[gp,gm] = partition(g);
g = combine(gp,gm);

end

function [pivotIndices, pivots, happy, poleRemoved] = PhaseOnePlus( F, tol )

[m, n] = size( F );

vscl = norm( F( : ), inf);

%
% Deal with the poles by removing them from F.
%
pole1 = mean(F(1,:));     % Take the value at the poles to be 
pole2 = mean(F(m,:));     % the mean of all the samples at the poles.

% If the the values at both poles are not zero then we need to add zero
% them out before removing these entries from F.
poleRemoved = 0;
if abs(pole1) > tol || abs(pole2) > tol
    % Determine the column with maximum inf-norm
    [ignored, poleCol] = max(max(abs(F),[],1));
    % Zero out the pole using the poleCol.
    F = F - F(:, poleCol)*ones(1,n);
    % Do we need to keep track of the pivot locations?
    poleRemoved = 1;
end

% Remove the rows corresponding to the poles in F before determining the
% rank.  We do this because then F is an even BMC matrix, which is what the
% code below requires.
F = F( 2:m-1, : );

[pivotIndices, pivots, happy] = PhaseOne( F, tol, poleRemoved );

% Adjust the pivot locations so that they now correspond to F having
% the poles.
if ~isempty( pivotIndices )
    pivotIndices(:,1) = pivotIndices(:,1) + 1;
end

% Put the poles at the begining the pivot locations array and also include
% the pivot matrix.
if poleRemoved
    pivotIndices = [ 1 poleCol; pivotIndices ];
    pivots = [1 ; pivots];
end

end

function [pivotIndices, pivots, happy] = PhaseOne( F, tol, rank_count )

% Phase 1: Go find rank, plus pivot locations, ignore cols and rows.
pivotIndices = []; pivots = [];

% Do GE with complete pivoting

while ( norm( F( : ), inf ) > tol )
    % Find pivot:
    [ignored, idx] = max( abs( F(:) ) );
    [j, k] = myind2sub( size( F ), idx );        
    pivotIndices = [ pivotIndices ; j k];
    pivot = F(j,k);
    
    F = F - F(:,k)*F(j,:)/pivot;
    
    pivots = [pivots; pivot];

    rank_count = rank_count + 1; 
end

% If the rank of the matrix is less than 1/8 its size. We are happy:
if ( rank_count < min(size(F))/8 )
    happy = 1;
else
    happy = 0;
end

end

function [cols, pivots, rows, pivotLocations] = PhaseTwo( h, pivotIndices, pivots, n, dom, tol, maxSample, sgn, removePoles)

if isempty( pivotIndices )
    cols = [];
    pivots = [];
    rows = [];
    pivotLocations = [];
    return;
end

happy_columns = 0;   % Not happy, until proven otherwise.
happy_rows = 0;
m = n;

[x, y] = getPoints( m, n, dom );

rk = size( pivotIndices, 1);
id_rows = pivotIndices(:,1); id_cols = pivotIndices(:,2);

% % Need to also include id_cols+n to account for the entries in the C
% % block
% id_cols = reshape([id_cols id_cols+n].',[],1);

col_pivots = x(id_cols);
row_pivots = y(id_rows);

numPivots = sum( abs( pivots ) > 0 );

% Phase 2: Calculate decomposition on sphere.
failure = false;
while ( ~(happy_columns && happy_rows) && ~failure )
    
    [x, y] = getPoints( m, n, dom );
    [xx, yy] = meshgrid( col_pivots, y);
    newCols = 0.5*( h( xx, yy ) + sgn*h( xx+pi, yy ) );
    
    
    [xx, yy] = meshgrid( x, row_pivots );
    newRows = h( xx, yy );
    % This code will be unnecessary once ticket #1532 is addressed on the
    % chebfun tracker.  Don't forget to remove it.
    if numel(row_pivots) == 1
        newRows = newRows(:).';
    end
    
    newRows = 0.5*( newRows(:,1:n) + sgn*newRows(:,n+1:2*n) );

    cols = zeros( m+1, numPivots );
    rows = zeros( numPivots, n );
    count = 1;

    % Need to remove pole, which means we use the column with the largest
    % max norm (repeated) with rows of all ones in the elimination
    % algorithm.
    if removePoles
        newRows(1,:) = 1;
    end
    
    
    for ii = 1:rk
        
        % Get the pivot
        pivot = pivots( ii, : );
        col = newCols( :, ii );
        row = newRows(ii,:);

        % Store the columns and rows.
        cols(:,count) = col;
        rows(count,:) = row;
        count = count + 1;
                
        newCols = newCols - (col*row(id_cols))/pivot;
        newRows = newRows - (col(id_rows)*row)/pivot;
    end    
    % Happiness check for columns:
    % TODO: Make this more similar to hapiness check in trigtech.

    % Double up the columns
    colsSum = sum([ cols ; sgn*flipud(cols(2:m,:)) ],2);
    col_coeffs = trigtech.vals2coeffs( colsSum ); 
    % Length of tail to test.
    testLength = min(m, max(3, round((m-1)/8)));
    tail = col_coeffs(1:testLength);

    if ( all( abs( tail ) <= 1e1*tol ) )
        happy_columns = 1;
    end
    
    % Happiness check for rows:
    % TODO: Make this more similar to hapiness check in trigtech.
    rowsSum = sum([ rows sgn*rows ],1).';
    row_coeffs = trigtech.vals2coeffs( rowsSum ); 
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
    
%     if happy_rows && n < 256 
%         happy_rows = false;
%         n = 2*n; 
%         id_cols = 2*id_cols - 1;
%     end
%     
%     if happy_columns && m < 256 
%         happy_columns = false;
%         m = 2*m; 
%         ii = [1:2:m-1 m+2:2:2*m]; 
%         id_rows = ii(id_rows); 
%     end
        
    if ( max(m, n) >= maxSample ) 
        warning('SPHEREFUN:constructor:notResolved', ...
        'Unresolved with maximum length: %u.', maxSample);
        failure = true;
    end 
end
cols = [ cols ; sgn*flipud(cols(2:end-1,:)) ];
rows = [ rows sgn*rows ];

% Adjust pivots to track indicies.
pivotLocations = [col_pivots row_pivots];

end

function [cols, pivots, rows, pivotLocations, idxPlus, idxMinus] = PhaseTwoOld( h, pivotIndices, pivotMatrices, n, dom, tol, maxSample, removePoles)

alpha = spherefun.alpha; % get growth rate factor.
happy_columns = 0;   % Not happy, until proven otherwise.
happy_rows = 0;
m = n;

[x, y] = getPoints( m, n, dom );

rk = size( pivotIndices, 1);
id_rows = pivotIndices(:,1); id_cols = pivotIndices(:,2);
% Need to also include id_cols+n to account for the entries in the C
% block
id_cols = reshape([id_cols id_cols+n].',[],1);

col_pivots = x(id_cols);
row_pivots = y(id_rows);

numPosPivots = sum( abs( pivotMatrices(:,1) ) > 0 );
numMinusPivots = sum( abs( pivotMatrices(:,2) ) > 0 );

pivotPlus = zeros(numPosPivots,1);
pivotMinus = zeros(numMinusPivots,1);

% Phase 2: Calculate decomposition on sphere.
failure = false;
while ( ~(happy_columns && happy_rows) && ~failure)
    
    [x, y] = getPoints( m, n, dom );
    [xx, yy] = meshgrid( col_pivots, y);
    newCols = h( xx, yy );
    
    [xx, yy] = meshgrid( x, row_pivots );
    newRows = h( xx, yy );

    % This code will be unnecessary once ticket #1532 is addressed on the
    % chebfun tracker.  Don't forget to remove it.
    if numel(row_pivots) == 1
        newRows = newRows(:).';
    end
    
    colsPlus = zeros( m+1, numPosPivots );
    colsMinus = zeros( m+1, numMinusPivots );
    rowsPlus = zeros( numPosPivots, 2*n );
    rowsMinus = zeros( numMinusPivots, 2*n );
    plusCount = 1;
    minusCount = 1;

    % Need to remove pole, which means we use the column with the largest
    % max norm (repeated) with rows of all ones in the elimination
    % algorithm.
    if removePoles
        newRows(1,:) = 1;
    end
    
    
    for ii = 1:rk
        
        % Get the eigenvalues of the pivot matrix M
        ev = pivotMatrices( ii, : );
        s = [max(abs(ev)) min(abs(ev))];
                                
        if ( s(1) <= alpha*s(2) )
            % Calculate inverse of pivot matrix:
            colPlus = newCols(:,2*ii-1) + newCols(:,2*ii);
            temp = newRows(ii,1:n) + newRows(ii,n+1:2*n);
            rowPlus = [temp temp];
            
            colMinus = newCols(:,2*ii-1) - newCols(:,2*ii);
            temp = newRows(ii,1:n) - newRows(ii,n+1:2*n);
            rowMinus = [temp -temp];
            
            % Store the columns and rows.
            colsPlus(:,plusCount) = colPlus;
            rowsPlus(plusCount,:) = rowPlus;
            pivotPlus(plusCount) = 2*ev(1);
            plusCount = plusCount + 1;

            colsMinus(:,minusCount) = colMinus;
            rowsMinus(minusCount,:) = rowMinus;
            pivotMinus(minusCount) = 2*ev(2);
            minusCount = minusCount + 1;

            newCols = newCols - ...
                    1/(2*ev(1))*(colPlus*rowPlus(id_cols)) - ...
                    1/(2*ev(2))*(colMinus*rowMinus(id_cols));                
            newRows = newRows - ...
                    1/(2*ev(1))*(colPlus(id_rows)*rowPlus) - ...
                    1/(2*ev(2))*(colMinus(id_rows)*rowMinus);
        else
            % Use the pseudoinverse of the pivot matrix, there is
            % no full rank pivot matrix:
            if abs(ev(1)) > abs(ev(2))                
                colPlus = newCols(:,2*ii-1) + newCols(:,2*ii);
                temp = newRows(ii,1:n) + newRows(ii,n+1:2*n);
                rowPlus = [temp temp];

                % Store the columns and rows.
                colsPlus(:,plusCount) = colPlus;
                rowsPlus(plusCount,:) = rowPlus;
                pivotPlus(plusCount) = 2*ev(1);
                plusCount = plusCount + 1;
                
                newCols = newCols - ...
                        1/(2*ev(1))*(colPlus*rowPlus(id_cols));
                newRows = newRows - ...
                        1/(2*ev(1))*(colPlus(id_rows)*rowPlus);
            else
                colMinus = newCols(:,2*ii-1) - newCols(:,2*ii);
                temp = newRows(ii,1:n) - newRows(ii,n+1:2*n);
                rowMinus = [temp -temp];

                % Store the columns and rows.
                colsMinus(:,minusCount) = colMinus;
                rowsMinus(minusCount,:) = rowMinus;
                pivotMinus(minusCount) = 2*ev(2);
                minusCount = minusCount + 1;
                
                newCols = newCols - ...
                        1/(2*ev(2))*(colMinus*rowMinus(id_cols));                
                newRows = newRows - ...
                        1/(2*ev(2))*(colMinus(id_rows)*rowMinus);
            end
        end
    end    
    % Happiness check for columns:
    % TODO: Make this more similar to hapiness check in trigtech.

    % Double up the columns
    cols = [ [colsPlus colsMinus] ; [flipud(colsPlus(2:m,:)) -flipud(colsMinus(2:m,:))] ];
    rows = [ rowsPlus ;  rowsMinus ].';
    col_coeffs = trigtech.vals2coeffs( cols ); 
    % Length of tail to test.
    testLength = min(m, max(3, round((m-1)/8)));
    tail = col_coeffs(1:testLength,:);

    if ( all( abs( tail ) <= 1e2*tol*norm(cols,inf) ) )
        happy_columns = 1;
    end
    
    % Happiness check for rows:
    % TODO: Make this more similar to hapiness check in trigtech.
    row_coeffs = trigtech.vals2coeffs( rows ); 
    % Length of tail to test.
    testLength = min(n, max(3, round((n-1)/8)));
    tail = row_coeffs(1:testLength,:);
    if ( all(abs( tail ) <= 1e2*tol*norm(rows,inf)) )
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
pivots = [pivotPlus;pivotMinus];
idxPlus = 1:numPosPivots;
idxMinus = (numPosPivots+1):(numPosPivots+numMinusPivots);

pivotLocations = [col_pivots(1:length(col_pivots)/2) row_pivots];

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

