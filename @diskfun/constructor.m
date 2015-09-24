function g = constructor( g, op, coords, dom, varargin )
%CONSTRUCTOR   The main DISKFUN constructor.
%
% The algorithm for constructing a DISKFUN comes in two phases:
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

if ( nargin == 0 )          % DISKFUN( )
    return
end

if ( isa(op, 'diskfun') )  % DISKFUN( DISKFUN )
    g = op;
    return
end

if ( isa(op, 'chebfun2') )  % DISKFUN( CHEBFUN2 ) 
    g = diskfun; 
    g.cols = op.cols; 
    g.rows = op.rows; 
    g.domain = op.domain; 
    g.pivotLocations = op.pivotLocations;
    g.pivotValues = op.pivotValues; 
    return
end

% If domain is empty, assign.
if ( nargin < 4 || isempty(dom) )
    dom = [-pi pi 0 1];
end

if ( nargin < 3 || isempty(coords) ) %for now assume polar if not specified
    coords = 1;
end

if strcmpi(coords,'cart') 
    coords = 0;
end

%coords=1 -> polar
%coords=0 -> cartesian

if ~isempty(varargin)
    alpha = varargin{1};
else
    alpha = 2;
end



maxRank = 8192; 
maxSample = 8192;
pseudoLevel = eps;

if ( isa(op, 'double') )    % DISKFUN( DOUBLE )
    % Should we allow coefficients to be passed in?
    
    % Only do Phase I on the values.
    F = op;
    [n, m] = size(F);
    
    % Flip F around since Phase I operates on the doubled-up portion of
    % the disk [-pi pi] x [-1, 0] 
    F = [F(n:-1:1,m/2+1:m) F(n:-1:1,1:m/2)];
    
    % TODO: Add a way to loosen tolerances for this type of construction.
    tol = GetTol(F, 2*pi/m, pi/(n-1), dom, 50*pseudoLevel);
    [pivotIndices, pivotArray, removePoles, happyRank, cols, pivots, ...
        rows, idxPlus, idxMinus ] = PhaseOne( F, tol, alpha, 0 );
    [x, y] = getPoints( n, m);
    pivotLocations = [x(pivotIndices(:,2)) y(pivotIndices(:,1))];
    
else  % DISKFUN( FUNCTION )
    % If f is defined in terms of x,y; then convert it to
    % (longitude,latitude).
    h = redefine_function_handle( op, coords );

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
        F = evaluate(h, n, n);

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
            warning('DISKFUN:CONSTRUCTOR:MAXRANK', ... 
                                    'Unresolved with maximum rank.');
            failure = 1;
        end
    end

    % PHASE TWO 
    % Find the appropriate discretizations in the columns and rows. 
    [cols, pivots, rows, pivotLocations, idxPlus, idxMinus, removePoles] = ...
        PhaseTwo( h, pivotIndices, pivotArray, n, dom, tol, maxSample, removePoles );
end

g.cols = chebfun( cols, dom(3:4)-[1 0]);
g.rows = chebfun( rows, dom(1:2), 'trig');
g.pivotValues = pivots;
g.pivotIndices = pivotIndices;
g.domain = dom;
g.idxPlus = idxPlus;
g.idxMinus = idxMinus;
g.nonZeroPoles = removePoles;

% Adjust the pivot locations so that they correspond to 
% -pi < th < pi and 0 < r <1

pivotLocations(:,2) = -pivotLocations(:,2); %vertical in [0 1]

pivotLocations(:,1) = pivotLocations(:,1) + pi; %adjust horz using BMC sym.
g.pivotLocations = pivotLocations;

% Sort according to the maginuted of the pivots using the partition and
% combine functions. 
% [gp,gm] = partition(g);
% g = combine(gp,gm);

end

function [pivotIndices, pivotArray, removePoles, ihappy, cols, pivots, ...
        rows, idxPlus, idxMinus ] = PhaseOne( F, tol, alpha, factor )

% Phase 1: Go find rank and pivot locations

% Setup
[m, n] = size( F );
minSize = min(m,n);
width = minSize/factor;
pivotIndices = []; pivotArray = [];
ihappy = 0;  % Assume we are not happy

B = F(:,1:n/2);    %% (1,1) Block of F.
C = F(:,n/2+1:n);  % (1,2) block of F.
Fp = 0.5*(B + C);
Fm = 0.5*(B - C);

%
% Deal with the pole by removing them from Fp.
%
                           % Take the value at the pole to be 
pole = mean(Fp(m,:));     % the mean of all the samples at the poles.

colsPlus = []; rowsPlus = []; kplus = 0;  idxPlus = [];
colsMinus = []; rowsMinus = []; kminus = 0; idxMinus = [];

rankCount = 0;    % keep track of the rank of the approximation.

% If the the value at pole is not zero then we need to zero
% out before removing these entries from F.
removePoles = false;
if  abs(pole) > tol
    % Determine the column with maximum inf-norm
    [ignored, poleCol] = max(max(abs(Fp),[],1));
    % Zero out the pole using the poleCol.
    rowPole = ones(1,n/2);
    colPole = Fp(:, poleCol);
    Fp = Fp - colPole*rowPole;
%     kplus = kplus + 1;
    % Do we need to keep track of the pivot locations?
    removePoles = true;
    % Update the rank count
    rankCount = rankCount + 1;
%     idxPlus(kplus) = rankCount;
end

% Remove the rows corresponding to the poles before determining the
% rank.  We do this because then F is an even BMC matrix, which is what the
% code below requires.
Fp = Fp( 1:m-1, : );
Fm = Fm( 1:m-1, : );

[maxp,idxp] = max(abs(Fp(:)));
[maxm,idxm] = max(abs(Fm(:)));

% Zero function
if ( maxp == 0 ) && ( maxm == 0 ) && ~( removePoles )
    colsPlus = 0;
    rowsPlus = 0;
    pivotArray = [0 0];
    pivotIndices = [1 1];
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
    [j, k] = myind2sub( [m-1 n/2], idx );
    
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
    cols = zeros( 2*m-1, rankCount );
    rows = zeros( n, rankCount );
    pivots = zeros( rankCount, 1);
    if kplus ~= 0
        cols(1:m-1,idxPlus) = colsPlus;
        cols(m+1:2*m-1,idxPlus) = flipud(colsPlus);
        rows(:,idxPlus) = [rowsPlus rowsPlus].';
        pivotPlus = pivotArray(pivotArray(:,1) ~= 0,1);
        pivots(idxPlus) = pivotPlus;
    end
    
    if kminus ~= 0
        cols(1:m-1,idxMinus) = colsMinus;
        cols(m+1:2*m-1,idxMinus) = -flipud(colsMinus);
        rows(:,idxMinus) = [rowsMinus -rowsMinus].';
        pivotMinus = pivotArray(pivotArray(:,2) ~= 0,2);
        pivots(idxMinus) = pivotMinus;
    end
    
%     pivots = reshape(pivotArray.',[],1);
%     pivots = pivots(pivots ~= 0 );
%     pivots = pivots([idxPlus idxMinus]);

    if removePoles
        cols(:,1) = [colPole;flipud(colPole(1:m-1))];
        rows(:,1) = [rowPole rowPole];
        pivots(1) = 1;
    end

end

% Adjust the pivot locations so that they now correspond to F having
% the poles.


% Put the poles at the begining of the pivot locations array and also include
% the pivot matrix.
if removePoles
    pivotIndices = [ 1 poleCol; pivotIndices ];
    pivotArray = [[1 0] ; pivotArray];
    idxPlus = [1 idxPlus];
end

end

function [cols, pivots, rows, pivotLocations, idxPlus, idxMinus, removePoles] = PhaseTwo( h, pivotIndices, pivotArray, n, dom, tol, maxSample, removePoles )

% alpha = diskfun.alpha; % get growth rate factor.
happy_columns = 0;   % Not happy, until proven otherwise.
happy_rows = 0;
m = n;

[x, y] = getPoints( m, n );

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

% Phase 2: Calculate decomposition on disk.
failure = false;
while ( ~(happy_columns && happy_rows) && ~failure)
    
    [x, y] = getPoints( m, n);
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

    % Double up the columns.
    temp1 = sum([colsPlus colsMinus],2); temp2 = sum([colsPlus -colsMinus],2);
    col_coeffs = chebtech2.vals2coeffs( [temp1;temp2(m:-1:1)] );

    % Length of tail to test.
    testLength = min(m, max(3, round((m-1)/8)));
    %tail = col_coeffs(1:testLength);
    tail = col_coeffs(end-testLength+1:end,:);

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
        warning('DISKFUN:constructor:notResolved', ...
        'Unresolved with maximum length: %u.', maxSample);
        failure = true;
    end 
end

% Combine the types of pivots and set-up indices to track them
cols = zeros( 2*size(colsPlus,1)-1, totalPivots );
cols(:,idxPlus) = [ colsPlus; flipud(colsPlus(1:end-1,:)) ];
cols(:,idxMinus) = [ colsMinus; -flipud(colsMinus(1:end-1,:)) ];

rows = zeros( 2*size(rowsPlus,2), totalPivots );
rows(:,idxPlus) = [rowsPlus rowsPlus].';
rows(:,idxMinus) = [rowsMinus -rowsMinus].';

pivotLocations = [col_pivots row_pivots];

end

function F = evaluate( h, m, n)
% Evaluate h on a m-by-n tensor product grid.

[x, y] = getPoints( m, n);

% Tensor product grid
[xx,yy] = meshgrid(x, y);

% Evaluate h on the non-doubled up grid
F = h( xx, yy );

end


function [x, y] = getPoints( m, n)

x = trigpts(2*n, [-pi, pi]);
y =  chebpts(2*m+1, [-1, 1]); %doubled;
y = y(1:(2*m)/2+1); %this includes pole (r=0), which needs evaluated. 


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


function f = redefine_function_handle( f, coords )

    % Wrap f so it can be evaluated in polar coordinates
    if ~(coords==1)
    f = @(th, r) diskfun.pol2cartf(f,th, r);
    end

end

function tol = GetTol(F, hx, hy, dom, pseudoLevel)
% GETTOL     Calculate a tolerance for the diskfun constructor.
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



