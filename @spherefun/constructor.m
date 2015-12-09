function g = constructor( g, op, dom, varargin )
%CONSTRUCTOR   The main SPHEREFUN constructor.
%
% This code is when functions on the surface of the sphere are represented
% as SPHEREFUN objects. A SPHEREFUN object is a low rank representation and
% expresses a function as a sum of rank-0 or 1 outerproduct of univariate
% functions in spherical coordinates.
%
% The algorithm for constructing a SPHEREFUN comes in two phases:
%
% PHASE 1: The first phase attempts to determine the numerical rank of the
% function by performing Gaussian elimination with special 2x2 pivoting
% matrices on a tensor grid of sample values. GE is perform until the
% sample matrix is approximated.  At the end of this stage we have
% candidate pivot locations and pivot elements.
%
% PHASE 2: The second phase attempts to resolve the corresponding column
% and row slices by sampling along the slices and performing GE (pivoting
% at 2x2 matrices) on the skeleton.   Sampling along each slice is
% increased until the Fourier coefficients of the slice fall below machine
% precision.
%
% The algorithm is fully described in:
%  A. Townsend, H. Wilber, and G. Wright, Computing with function on
%  spherical and polar geometries I. The sphere, submitted, 2015. 
%
% See also SPHEREFUN.

if ( nargin == 0 )          % SPHEREFUN( )
    return
end

if ( isa(op, 'spherefun') )  % SPHEREFUN( SPHEREFUN )
    g = op;
    return
end

% If domain is empty take it to be co-latitude.
if ( nargin < 3 || isempty( dom ) )
     dom = [-pi pi 0 pi]; 
elseif ( numel( dom ) ~= 4 )
    error('CHEBFUN:SPHEREFUN:CONSTRUCTOR:domain',... 
          ['A domain is rarely given for spherefun, ' ... 
          'but it needs to be given by four corner values' 
          'in intrinstic coordinates.'])
elseif ( numel( dom ) == 4 && norm( dom(:)' - [-pi pi 0 pi] ) >0 )
    error('CHEBFUN:SPHEREFUN:CONSTRUCTOR:domain',...
        'The domain of a spherefun is always [-pi pi]x[0 pi] in intrinstic coordinates');
else
    dom = [-pi pi 0 pi];
end

% Default value for coupling parameter
alpha = 100;

% TODO: Should we allow any other domains than latitude and co-latitude?

% TODO: 
% 1. Need to allow for different domains.
% 2. Add support for preferences
% 3. Add non-adaptive construction
% 4. Add fixed-rank.
% 5. Add tensor-product.

maxRank = 8192; 
maxSample = 8192;
pseudoLevel = eps;

if ( isa(op, 'double') )    % SPHEREFUN( DOUBLE )
    % Should we allow coefficients to be passed in?
    
    % Only do Phase I on the values.
    F = op;
    [n, m] = size(F);
    
    if mod(m,2) ~= 0
        error('SPHEREFUN:CONSTRUCTOR:VALUES', ... 
         'When constructing from values the number of columns must be even.');
    end
    
    % Flip F arround since Phase I operates on the doubled-up portion of
    % the sphere [-pi pi] x [-pi, 0] or [-pi pi] x [-3*pi/2 -pi/2]
    F = [F(n:-1:1,m/2+1:m) F(n:-1:1,1:m/2)];
    
    % TODO: Add a way to loosen tolerances for this type of construction.
    tol = GetTol(F, 2*pi/m, pi/(n-1), dom, pseudoLevel);
    [pivotIndices, pivotArray, removePoles, happyRank, cols, pivots, ...
        rows, idxPlus, idxMinus ] = PhaseOne( F, tol, alpha, 0 );
    [x, y] = getPoints( n, m, dom );
    pivotLocations = [x(pivotIndices(:,2)) y(pivotIndices(:,1))];
    
else  % SPHEREFUN( FUNCTION )
    % If f is defined in terms of x,y,z; then convert it to
    % (longitude,latitude).
    h = redefine_function_handle( op );

    % PHASE ONE  
    % Sample at square grids, determine the numerical rank of the
    % function.
    n = 4;
    happyRank = 0;     % Happy with phase one? 
    failure = 0;
    pivotIndices = [];
    pivotArray = [];
    while ( ~happyRank && ~failure )
        n = 2*n;

        % Sample function on a tensor product grid.
        % TODO: Add a more sophisticated evaluate function that does
        % vectorization like chebfun2.
        F = evaluate(h, n, n, dom);

        tol = GetTol(F, pi/n, pi/n, dom, pseudoLevel);

        pivotIndices2 = pivotIndices;
        pivotArray2 = pivotArray;
        [ pivotIndices, pivotArray, removePoles, happyRank ] = ...
            PhaseOne( F, tol, alpha, 8 );

        if size(pivotIndices,1) > size(pivotIndices2,1)
            happyRank = 0;
        else
            % If the norm of the pivots for the n/2 case is within 1/10 of
            % the norm for the n case then just keep the smaller n as this
            % will be faster.
            if norm(pivotArray2,inf) > 0.1*norm(pivotArray,inf) && happyRank
                pivotIndices = pivotIndices2;
                pivotArray = pivotArray2;
                n = n/2;
            end
        end

        if ( n >= maxRank  )
            warning('SPHEREFUN:CONSTRUCTOR:MAXRANK', ... 
                                    'Unresolved with maximum rank.');
            failure = 1;
        end
    end

    % PHASE TWO 
    % Find the appropriate discretizations in the columns and rows. 
    [cols, pivots, rows, pivotLocations, idxPlus, idxMinus] = ...
        PhaseTwo( h, pivotIndices, pivotArray, n, dom, tol, maxSample, removePoles );
end

g.cols = chebfun( cols, dom(3:4)-[pi 0], 'trig');
% c = g.cols.funs{1}.onefun.coeffs;
% offset = real(sum(c));
% n = size(c,1);
% c(n/2+1,:) = c(n/2+1,:)-offset;
% g.cols.funs{1}.onefun.coeffs = c;
g.rows = chebfun( rows, dom(1:2), 'trig');
if all(pivots) == 0
    pivots = inf;
end
g.pivotValues = pivots;
g.domain = dom;
g.idxPlus = idxPlus;
g.idxMinus = idxMinus;
g.nonZeroPoles = removePoles;
g.pivotLocations = adjustPivotLocations(pivotLocations, pivotArray, iscolat(g) ); 

% Simplifying rows and columns after they are happy.
g = simplify( g );

end

function [pivotIndices, pivotArray, removePole, ihappy, cols, pivots, ...
        rows, idxPlus, idxMinus ] = PhaseOne( F, tol, alpha, factor )

% Phase 1: Go find rank and pivot locations

% Setup
[m, n] = size( F );
minSize = min(m,n);
width = minSize/factor;
pivotIndices = []; pivotArray = [];
ihappy = 0;  % Assume we are not happy

% If given a 1xn matrix, then this only gives us a function samples at the 
% two poles, which is simple to deal with.
if ( m <= 1 ) 
    error('CHEBFUN:SPHEREFUN:constructor:poleSamples',...
        ['Matrix of function samples contains < 2 rows. ',...
        'This is not enough information to reconstruct the function. Please increase samples in the latitudinal direction.'])
end

% Only information at the poles is given.
if ( m == 2) 
    cols = F(:,1);
    rows = F(1,:).';
    idxPlus = 1;
    idxMinus = [];
    pivotArray = [1 0];
    pivotIndices = [1 1];
    removePole = false;
    pivots = 1;
    ihappy = 1;
    return;
end

B = F(:,1:n/2);    %% (1,1) Block of F.
C = F(:,n/2+1:n);  % (1,2) block of F.
Fp = 0.5*(B + C);
Fm = 0.5*(B - C);

%
% Deal with the poles by removing them from Fp.
%

% Check if the poles are numerically constant and get the value.
[pole1,constValue1] = checkPole(Fp(1,:),tol);
[pole2,constValue2] = checkPole(Fp(m,:),tol);

if ~(constValue1 || constValue1)
    warning('CHEBFUN:SPHEREFUN:constructor:constPoles',...
        ['Results may be inaccurate as the function may not be constant '...
         'at either the north or south poles.']);
end

colsPlus = []; rowsPlus = []; kplus = 0;  idxPlus = [];
colsMinus = []; rowsMinus = []; kminus = 0; idxMinus = [];

rankCount = 0;    % keep track of the rank of the approximation.

% If the the values at both poles are not zero then we need to add zero
% them out before removing these entries from F.
removePole = false;
if abs(pole1) > tol || abs(pole2) > tol
    % Determine the column with maximum inf-norm
    [rowVal, poleCol] = max(max(abs(Fp),[],1));
    % Zero out the pole using the poleCol.
    rowPole = rowVal*ones(1,n/2);
    colPole = Fp(:, poleCol);
    Fp = Fp - colPole*(rowPole/rowVal);
    % Do we need to keep track of the pivot locations?
    removePole = true;
    % Update the rank count
    rankCount = rankCount + 1;
end

% If given a 1xn matrix, then this only gives us a function samples at the 
% two poles, which is simple to deal with.
if ( m <= 1 ) 
    error('CHEBFUN:SPHEREFUN:constructor:poleSamples',...
        ['Matrix of function samples contains < 2 rows. ',...
        'This is not enough information to reconstruct the function. Please increase samples in the latitudinal direction.'])
end

% Remove the rows corresponding to the poles before determining the
% rank.  We do this because then F is an even BMC matrix, which is what the
% code below requires.
Fp = Fp( 2:m-1, : );
Fm = Fm( 2:m-1, : );

[maxp,idxp] = max(abs(Fp(:)));
[maxm,idxm] = max(abs(Fm(:)));

% Zero function
if ( maxp == 0 ) && ( maxm == 0 ) && ~( removePole )
    m = 3; n = 3;
    cols = zeros( 2*m-2, 1 );
    rows = zeros( n, 1 );
    idxPlus = 1;
    idxMinus = [];
    pivotArray = [0 0];
    pivotIndices = [1 1];
    pivots = inf;
    ihappy = 1;
    return;
end

while ( ( max( maxp, maxm ) > tol ) && ( rankCount < width ) && ...
    ( rankCount < minSize ) )    
    % Find pivots:
    if maxp >= maxm
        idx = idxp;
    else
        idx = idxm;
    end
    [j, k] = myind2sub( [m-2 n/2], idx );
    
    % Use maximum of the Fp and Fm matrices for pivots
    evp = Fp(j,k); absevp = abs(evp);
    evm = Fm(j,k); absevm = abs(evm);
            
    pivotIndices = [ pivotIndices ; j k];
    
    % Smallest pivots is within an acceptable multiple of larger pivot so
    % do a rank 2 update.
    if ( max( absevp, absevm ) <= alpha*min( absevp, absevm ) )
        kplus = kplus + 1;
        colsPlus(:,kplus) = Fp(:,k);
        rowsPlus(kplus,:) = Fp(j,:);
        Fp = Fp - colsPlus(:,kplus)*(rowsPlus(kplus,:)*(1/evp));
        
        kminus = kminus + 1;
        colsMinus(:,kminus) = Fm(:,k);
        rowsMinus(kminus,:) = Fm(j,:);
        Fm = Fm - colsMinus(:,kminus)*(rowsMinus(kminus,:)*(1/evm));        
        
        rankCount = rankCount + 1;
        if absevp >= absevm
            idxPlus(kplus) = rankCount;
            rankCount = rankCount + 1;
            idxMinus(kminus) = rankCount;
        else
            idxMinus(kminus) = rankCount;
            rankCount = rankCount + 1;
            idxPlus(kplus) = rankCount;
        end            
        pivotArray = [pivotArray ; [evp evm] ];
        [maxp,idxp] = max(abs(Fp(:)));
        [maxm,idxm] = max(abs(Fm(:)));
    else
        % Positive pivot dominates
        if absevp > absevm
            kplus = kplus + 1;
            rankCount = rankCount + 1;

            colsPlus(:,kplus) = Fp(:,k);
            rowsPlus(kplus,:) = Fp(j,:);
            Fp = Fp - colsPlus(:,kplus)*(rowsPlus(kplus,:)*(1/evp));
            idxPlus(kplus) = rankCount;
            
            % Minus pivot is zero
            evm = 0;
            [maxp,idxp] = max(abs(Fp(:)));
        % Negative pivot dominates
        else
            kminus = kminus + 1;
            rankCount = rankCount + 1;

            colsMinus(:,kminus) = Fm(:,k);
            rowsMinus(kminus,:) = Fm(j,:);
            Fm = Fm - colsMinus(:,kminus)*(rowsMinus(kminus,:)*(1/evm));
            idxMinus(kminus) = rankCount;

            % Plus pivot is zero
            evp = 0;
            [maxm,idxm] = max(abs(Fm(:)));
        end
        pivotArray = [pivotArray ; [evp evm] ];
    end
end

if ( max( maxp, maxm ) <= tol )
    ihappy = 1;                               % We are happy
end
if ( rankCount >= width )
    ihappy = 0;                               % We are not happy
end

% No sense in giving row and column values if they are not wanted.
if ( nargout > 4 )
    % Combine the types of pivots and set-up indices to track them
    cols = zeros( 2*m-2, rankCount );
    rows = zeros( n, rankCount );
    pivots = zeros( rankCount, 1);
    if kplus ~= 0
        cols(2:m-1,idxPlus) = colsPlus;
        cols(m+1:2*m-2,idxPlus) = flipud(colsPlus);
        rows(:,idxPlus) = [rowsPlus rowsPlus].';
        pivotPlus = pivotArray(pivotArray(:,1) ~= 0,1);
        pivots(idxPlus) = pivotPlus;
    end
    
    if kminus ~= 0
        cols(2:m-1,idxMinus) = colsMinus;
        cols(m+1:2*m-2,idxMinus) = -flipud(colsMinus);
        rows(:,idxMinus) = [rowsMinus -rowsMinus].';
        pivotMinus = pivotArray(pivotArray(:,2) ~= 0,2);
        pivots(idxMinus) = pivotMinus;
    end
    
%     pivots = reshape(pivotArray.',[],1);
%     pivots = pivots(pivots ~= 0 );
%     pivots = pivots([idxPlus idxMinus]);

    if removePole
        cols(:,1) = [colPole;flipud(colPole(2:m-1))];
        rows(:,1) = [rowPole rowPole];
        pivots(1) = rowVal;
    end

end

% Adjust the pivot locations so that they now correspond to F having
% the poles.
if ~isempty( pivotIndices )
    pivotIndices(:,1) = pivotIndices(:,1) + 1;
end

% Put the poles at the begining of the pivot locations array and also include
% the pivot matrix.
if removePole
    pivotIndices = [ 1 poleCol; pivotIndices ];
    pivotArray = [[rowVal 0] ; pivotArray];
    idxPlus = [1 idxPlus];
end

end

function [cols, pivots, rows, pivotLocations, idxPlus, idxMinus] = PhaseTwo( h, pivotIndices, pivotArray, n, dom, tol, maxSample, removePoles )

% alpha = spherefun.alpha; % get growth rate factor.
happy_columns = 0;   % Not happy, until proven otherwise.
happy_rows = 0;
m = n;

[x, y] = getPoints( m, n, dom );

rk = size( pivotIndices, 1);
id_rows = pivotIndices(:,1); id_cols = pivotIndices(:,2);

row_pivots = y(id_rows);
col_pivots = x(id_cols);

numPosPivots = sum( abs( pivotArray(:,1) ) > 0 );
numMinusPivots = sum( abs( pivotArray(:,2) ) > 0 );
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
        newRowsPlus(1,:) = pivotArray(1,1);
    end
    
    
    for ii = 1:rk
        
        % Get the pivots
        evp = pivotArray( ii, 1 );
        evm = pivotArray( ii, 2 );
                                        
        % Do GE step on both matrices
        if evp ~=0 && evm ~= 0
            colPlus = newColsPlus(:,ii);
            rowPlus = newRowsPlus(ii,:);
            
            colMinus = newColsMinus(:,ii);
            rowMinus = newRowsMinus(ii,:);

            % Store the columns and rows.
            colsPlus(:,plusCount) = colPlus;
            rowsPlus(plusCount,:) = rowPlus;
            pivotPlus(plusCount) = evp;
            
            colsMinus(:,minusCount) = colMinus;
            rowsMinus(minusCount,:) = rowMinus;
            pivotMinus(minusCount) = evm;

            newColsPlus = newColsPlus - ...
                    colPlus*(rowPlus(id_cols)*(1/evp));
            newRowsPlus = newRowsPlus - ...
                    ((1/evp)*colPlus(id_rows))*rowPlus;
            newColsMinus = newColsMinus - ...
                    colMinus*(rowMinus(id_cols)*(1/evm));
            newRowsMinus = newRowsMinus - ...
                    ((1/evm)*colMinus(id_rows))*rowMinus;
                
            if abs(evp) >= abs(evm)
                idxPlus(plusCount) = pivotCount;
                idxMinus(minusCount) = pivotCount+1;
                pivots(pivotCount) = evp;
                pivots(pivotCount+1) = evm;
            else
                idxMinus(minusCount) = pivotCount;
                idxPlus(plusCount) = pivotCount+1;
                pivots(pivotCount) = evm;
                pivots(pivotCount+1) = evp;
            end 
            plusCount = plusCount + 1;
            minusCount = minusCount + 1;
            pivotCount = pivotCount + 2;
        else
            if ( evp ~= 0 )
                colPlus = newColsPlus(:,ii);
                rowPlus = newRowsPlus(ii,:);

                % Store the columns and rows.
                colsPlus(:,plusCount) = colPlus;
                rowsPlus(plusCount,:) = rowPlus;
                pivotPlus(plusCount) = evp;
                idxPlus(plusCount) = pivotCount;
                pivots(pivotCount) = evp;

                plusCount = plusCount + 1;
                pivotCount = pivotCount + 1;
                
                newColsPlus = newColsPlus - ...
                        colPlus*(rowPlus(id_cols)*(1/evp));
                newRowsPlus = newRowsPlus - ...
                        ((1/evp)*colPlus(id_rows))*rowPlus;
            else
                colMinus = newColsMinus(:,ii);
                rowMinus = newRowsMinus(ii,:);

                % Store the columns and rows.
                colsMinus(:,minusCount) = colMinus;
                rowsMinus(minusCount,:) = rowMinus;
                pivotMinus(minusCount) = evm;
                idxMinus(minusCount) = pivotCount;
                pivots(pivotCount) = evm;

                minusCount = minusCount + 1;
                pivotCount = pivotCount + 1;

                newColsMinus = newColsMinus - ...
                        colMinus*(rowMinus(id_cols)*(1/evm));
                newRowsMinus = newRowsMinus - ...
                        ((1/evm)*colMinus(id_rows))*rowMinus;
            end
        end
    end    
    % Happiness check for columns:
    % TODO: Make this more similar to hapiness check in trigtech.
    
    if removePoles
        colsPlus(1,2:end) = 0;
        colsPlus(end,2:end) = 0;
    elseif ~isempty(colsPlus)
        colsPlus(1,:) = 0;
        colsPlus(end,:) = 0;
    end
    
    if ~isempty(colsMinus)
        colsMinus(1,:) = 0;
        colsMinus(end,:) = 0;
    end

    % Double up the columns.
    temp1 = sum([colsPlus colsMinus],2); temp2 = sum([colsPlus -colsMinus],2);
    col_coeffs = trigtech.vals2coeffs( [temp1;temp2(m:-1:2)] );

%     colData.hscale = norm(dom(3:4), inf);
%     colValues = [temp1;temp2(m:-1:2)];
%     colTrigtech = trigtech.make(colValues, colData);
%     resolvedCol = happinessCheck(colTrigtech, [], colValues, colData);
    
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

%     rowValues = [temp1 temp2].';
%     rowData.hscale = norm(dom(1:2), inf);
%     rowTrigtech = trigtech.make(rowValues, rowData);
%     resolvedRows = happinessCheck(rowTrigtech, [], rowValues, rowData);
    
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

function pivLocNew = adjustPivotLocations(pivLoc, pivArray, colat)

% Adjust the pivot locations so that they correspond to 
% -pi < lam < pi and 0 < th < pi or -pi/2 < th < pi/2
if colat
    pivLoc(:,2) = -pivLoc(:,2);
else
    pivLoc(:,2) = -(pivLoc(:,2)+pi);
end
pivLoc(:,1) = pivLoc(:,1) + pi;

% We will store the pivotLocations for both the plus and minus
% pieces, which could result in duplicate values being stored.  This
% happens whenever a rank-2 deflation step occurs. Really only one set of
% pivotLocations need to be stored in this case, but we will double up the
% information as it makes the combine and partition methods much easier to
% write.  If this is changed then look the combine and partition methods
% need to also be changed.

pivLocNew = zeros(sum(sum(pivArray ~= 0)),2);
count = 1;
for j=1:size(pivLoc,1)
    if pivArray(j,1) ~= 0 && pivArray(j,2) ~= 0
        pivLocNew(count,:) = pivLoc(j,:);
        pivLocNew(count+1,:) = pivLoc(j,:);
        count = count + 2;
    else
        pivLocNew(count,:) = pivLoc(j,:);
        count = count + 1;
    end
end

end

% Check that the values at the pole are constant.
function [pole,constVal] = checkPole(val,tol)

% Take the mean of the values that are at the pole.
pole = mean(val);
% Compute there standard deviation
stddev = std(val);

% If the standard deviation does not exceed the 100*tolearnce then the pole
% is "constant".
% TODO: Get a better feel for the tolerance check.
if stddev <= 1e8*tol || stddev < eps
    constVal = 1;
else
    constVal = 0;
end

end

% function f = redefine_function_handle_pole( f, poleColPivot )
% % Set f to f - f(poleColPivot,theta) where poleColPivot is the value of
% % lambda (the column) used to zero out the poles of f.
% 
% f = @(lam, th) f(lam,th) - f(poleColPivot,th);
% 
% end


