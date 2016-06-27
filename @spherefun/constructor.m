function g = constructor(g, op, varargin)
%CONSTRUCTOR   The main SPHEREFUN constructor.
%
% This code is when functions on the surface of the sphere are represented
% as SPHEREFUN objects. A SPHEREFUN object is a real-valued function as a 
% sum of rank 1 outerproducts of univariate functions in spherical 
% coordinates.
%
% The algorithm for constructing a SPHEREFUN comes in two phases:
%
% PHASE 1: The first phase attempts to determine the numerical rank of the
% function by performing Gaussian elimination with special 2x2 pivoting
% matrices on a tensor grid of sample values. GE is performed until the
% sample matrix is approximated. At the end of this stage we have candidate
% pivot locations and pivot elements.
%
% PHASE 2: The second phase attempts to resolve the corresponding column
% and row slices by sampling along the slices and performing GE (pivoting
% at 2x2 matrices) on the skeleton. Sampling along each slice is increased 
% until the Fourier coefficients of the slice fall below machine precision.
%
% The algorithm is fully described in:
%  A. Townsend, H. Wilber, and G. Wright, Computing with function on
%  spherical and polar geometries I: The sphere, SIAM J. Sci. Comput., 
%  Accepted, 2016. 
%
% See also SPHEREFUN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 0 )          % SPHEREFUN( )
    return
end

% Parse the inputs:
[op, dom, pref, fixedRank, vectorize] = parseInputs(op, varargin{:});

% Return op if construction is from coefficients which is handled in
% parseInputs.
if ( isa(op, 'spherefun') )  
    g = op;
    % Fix the rank:
    g = fixTheRank(g, fixedRank);
    return
end

% Preferences:
tech        = trigtech;
tpref       = tech.techPref;
minSample   = 4;
maxSample   = tpref.maxLength;
cheb2Prefs  = pref.cheb2Prefs;
sampleTest  = cheb2Prefs.sampleTest;
maxRank     = cheb2Prefs.maxRank;
pseudoLevel = cheb2Prefs.chebfun2eps;

alpha = 100; % Default value for coupling parameter

if ( isa(op, 'char') )     % SPHEREFUN( CHAR )
    op = str2op( op );
end

% TODO: 
% 1. Need to allow for different domains.
% 2. Add non-adaptive construction
% 3. Add tensor-product.

% Deal with constructions from numeric data:
if ( isa(op, 'double') )    % SPHEREFUN( DOUBLE )
    g = constructFromDouble(op, dom, alpha, pref);
    % Fix the rank:
    g = fixTheRank(g, fixedRank);
    return
end

%
% Construction is from a function handle.
%

% Check for op = @(lam,th) constant function
[ll, tt] = meshgrid(dom(1:2), dom(3:4));
if ( numel(op(ll,tt)) == 1 )
    op1 = op;
    op = @(ll, tt) op1(ll, tt) + 0*ll;
end

factor  = 8; % Ratio between size of matrix and no. pivots.
isHappy = 0; % If we are currently unresolved.
failure = 0; % Reached max discretization size without being happy.

while ( ~isHappy && ~failure )
    %
    % Setup Phase I: GE with block 2-by-2 pivoting to determine the
    % numerical rank and pivot locations.  Sampling is done at equally
    % spaced square grids.
    %
    grid = minSample;          
    happyRank = 0;             % Happy with phase one? 
    strike = 1;
    while ( ~happyRank && ~failure && strike < 3)
        grid = 2*grid;

        % Sample function on a tensor product grid.
        [x, y] = getPoints(grid, grid, dom);
        [xx, yy] = meshgrid(x, y);
        F = evaluate(op, xx, yy, vectorize);
        
        if ( ~isreal( F ) ) 
            warning('SPHEREFUN:CONSTRUCTOR:COMPLEX', ...
                    ['Only real-valued spherefuns are supported. The '...
                     'imaginary part is being set to zero now.'])
             F = real( F );   
        end

        [tol, vscale] = getTol(F, pi/grid, pi/grid, dom, pseudoLevel);
        pref.chebfuneps = tol;
        
        % Does the function blow up or evaluate to NaN?:
        if ( isinf(vscale) )
            error('CHEBFUN:SPHEREFUN:constructor:inf', ...
                'Function returned INF when evaluated');
        elseif ( any(isnan(F(:)) ) )
            error('CHEBFUN:SPHEREFUN:constructor:nan', ...
                'Function returned NaN when evaluated');
        end
        
        % Do GE
        [pivotIndices, pivotArray, removePoles, happyRank] = ...
            PhaseOne(F, tol, alpha, factor);

        if ( grid > factor*(maxRank-1) )
            warning('SPHEREFUN:CONSTRUCTOR:MAXRANK', ... 
                                    'Unresolved with maximum rank.');
            failure = 1;
        end
        
        % If the function is 0+noise then stop after three strikes.
        if ( max(abs(pivotArray(1,:))) < 1e4*tol )
            strike = strike + 1;
        end
    end

    % Do Phase 2: resolve along the column and row slices.
    [cols, pivots, rows, pivotLocations, idxPlus, idxMinus, isHappy, failure] = ...
        PhaseTwo(op, pivotIndices, pivotArray, grid, dom, vscale, ...
        maxSample, removePoles, vectorize, pref);
    
    g.cols = chebfun(cols, dom(3:4)-[pi 0], 'trig', pref);
    g.rows = chebfun(rows, dom(1:2), 'trig', pref);
    if ( all(pivots) == 0 )
        pivots = inf;
    end
    g.pivotValues = pivots;
    g.domain = dom;
    g.idxPlus = idxPlus;
    g.idxMinus = idxMinus;
    g.nonZeroPoles = removePoles;
    g.pivotLocations = adjustPivotLocations(pivotLocations, pivotArray); 

    % Sample Test:
    if ( sampleTest )
        % wrap the op with evaluate in case the 'vectorize' flag is on:
        sampleOP = @(lam,th) evaluate(op, lam, th, vectorize);
        
        % Evaluate at points in the domain:
        pass = g.sampleTest(sampleOP, tol, vectorize);
        if ( ~pass )
            % Increase minSamples and try again.
            minSample = 2*minSample;
            isHappy = 0;
        end
    end
end

% Simplifying rows and columns after they are happy.
g = simplify( g, pref.chebfuneps );

% Fix the rank, if in nonadaptive mode.
g = fixTheRank( g , fixedRank );

% Project onto BMC-I symmetry so the function is smooth on the sphere.
g = projectOntoBMCI( g );

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = constructFromDouble(op, dom, alpha, pref)

g = spherefun();

if ( ~isreal( op ) ) 
    warning('SPHEREFUN:CONSTRUCTOR:COMPLEX', ...
            ['Only real-valued spherefuns are supported. The '...
             'imaginary part is being set to zero now.'])
    op = real( op );
end

% If single numerical value given
if ( (numel( op ) == 1) )
    g = constructor(g, @(lam,th) op + 0*lam, dom);
    return
end

% Only do Phase I on the values.
F = op;
[n, m] = size(F);

if ( mod(m,2) ~= 0 )
    error('SPHEREFUN:CONSTRUCTOR:VALUES', ... 
     'When constructing from values the number of columns must be even.');
end

% TODO: Add a way to loosen tolerances for this type of construction.
tol = getTol(F, 2*pi/m, pi/(n-1), dom, pref.cheb2Prefs.chebfun2eps);
pref.chebfuneps = tol;

% Perform GE with complete pivoting
[pivotIndices, pivotArray, removePoles, unused, cols, pivots, ...
    rows, idxPlus, idxMinus] = PhaseOne(F, tol, alpha, 0);

[x, y] = getPoints(n, m, dom);
pivotLocations = [ x(pivotIndices(:, 2)) y(pivotIndices(:, 1)) ];

g.cols = chebfun(cols, dom(3:4)-[pi 0], 'trig', pref);
g.rows = chebfun(rows, dom(1:2), 'trig', pref);
if ( all(pivots) == 0 )
    pivots = inf;
end
g.pivotValues = pivots;
g.domain = dom;
g.idxPlus = idxPlus;
g.idxMinus = idxMinus;
g.nonZeroPoles = removePoles;
g.pivotLocations = adjustPivotLocations(pivotLocations, pivotArray); 

g = projectOntoBMCI(g);

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pivotIndices, pivotArray, removePole, isHappy, cols, pivots, ...
        rows, idxPlus, idxMinus] = PhaseOne(F, tol, alpha, factor)

% Phase 1: Go find rank and pivot locations

% Setup
[m, n] = size(F);
minSize = min(m, n);
width = minSize/factor;
pivotIndices = []; 
pivotArray = [];
isHappy = 0;  % Assume we are not happy

% If given a 1xn matrix, then this only gives us a function samples at the 
% two poles, which is simple to deal with.
if ( m <= 1 ) 
    error('CHEBFUN:SPHEREFUN:constructor:poleSamples',...
        ['Matrix of function samples contains < 2 rows. ',...
        'This is not enough information to reconstruct the function. ',...
        'Please increase samples in the latitudinal direction.'])
end

% Only information at the poles is given.
if ( m == 2 ) 
    cols = F(:, 1);
    rows = F(1, :).';
    idxPlus = 1;
    idxMinus = [];
    pivotArray = [1 0];
    pivotIndices = [1 1];
    removePole = true;
    pivots = 1;
    isHappy = 1;
    return;
end

C = F(:, 1:n/2);    % (2,2) Block of F.
B = F(:, n/2+1:n);  % (1,2) block of F.
Fp = 0.5*(B + C);
Fm = 0.5*(B - C);

%
% Deal with the poles by removing them from Fp.
%

% Check if the poles are numerically constant and get the value.
pole1 = checkPole(Fp(1, :), tol);
pole2 = checkPole(Fp(m, :), tol);

% TODO: Figure out if we really need to warn the user about their function
% % Check if the poles are numerically constant and get the value.
% [pole1, constValue1] = checkPole(Fp(1, :), tol);
% [pole2, constValue2] = checkPole(Fp(m, :), tol);
% not being constant along the poles.
% if ~(constValue1 || constValue1)
%     warning('CHEBFUN:SPHEREFUN:constructor:constPoles',...
%         ['Results may be inaccurate as the function may not be constant '...
%          'at either the north or south poles.']);
% end

colsPlus = []; 
rowsPlus = []; 
kplus = 0;  
idxPlus = [];
colsMinus = []; 
rowsMinus = []; 
kminus = 0; 
idxMinus = [];
rankCount = 0;    % keep track of the rank of the approximation.

% If the the values at both poles are not zero then we need to add zero
% them out before removing these entries from F.
removePole = false;
if ( (abs(pole1) > tol) || (abs(pole2) > tol) )
    
    % Determine the column with maximum inf-norm
    [rowVal, poleCol] = max(max(abs(Fp), [], 1));
    
    % Zero out the pole using the poleCol.
    rowPole = rowVal*ones(1, n/2);
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
        'This is not enough information to reconstruct the function. '...
        'Please increase samples in the latitudinal direction.'])
end

% Remove the rows corresponding to the poles before determining the
% rank.  We do this because then F is an even BMC matrix, which is what the
% code below requires.
Fp = Fp(2:m-1, :);
Fm = Fm(2:m-1, :);

[maxp, idxp] = max(abs(Fp(:)));
[maxm, idxm] = max(abs(Fm(:)));

% Zero function
if ( (maxp == 0) && (maxm == 0) && ~removePole )
    % Pass back a zero matrix that is the same size as F. 
    % This ensures that spherefun( zeros(5) ) has a 5x5 (zero) coefficient 
    % matrix.      
    cols = zeros(m, 1);
    rows = zeros(n, 1);
    idxPlus = 1;
    idxMinus = [];
    pivotArray = [0 0];
    pivotIndices = [1 1];
    pivots = inf;
    isHappy = 1;
    return;
end

while ( (max(maxp, maxm) > tol) && (rankCount < width) && ...
    (rankCount < minSize) )
    % Find pivots:
    if ( maxp >= maxm )
        idx = idxp;
    else
        idx = idxm;
    end
    [j, k] = myind2sub([m-2 n/2], idx);
    
    % Use maximum of the Fp and Fm matrices for pivots
    evp = Fp(j, k); 
    absevp = abs(evp);
    evm = Fm(j,k); 
    absevm = abs(evm);
    pivotIndices = [ pivotIndices; j k];
    
    % Smallest pivots is within an acceptable multiple of larger pivot so
    % do a rank 2 update.
    if ( max(absevp, absevm) <= alpha*min(absevp, absevm) )
        kplus = kplus + 1;
        colsPlus(:, kplus) = Fp(:, k);
        rowsPlus(kplus, :) = Fp(j, :);
        Fp = Fp - colsPlus(:, kplus)*(rowsPlus(kplus, :)*(1/evp));
        
        kminus = kminus + 1;
        colsMinus(:, kminus) = Fm(:, k);
        rowsMinus(kminus, :) = Fm(j, :);
        Fm = Fm - colsMinus(:, kminus)*(rowsMinus(kminus, :)*(1/evm));
        
        rankCount = rankCount + 1;
        if ( absevp >= absevm )
            idxPlus(kplus) = rankCount;
            rankCount = rankCount + 1;
            idxMinus(kminus) = rankCount;
        else
            idxMinus(kminus) = rankCount;
            rankCount = rankCount + 1;
            idxPlus(kplus) = rankCount;
        end            
        pivotArray = [ pivotArray; [evp evm] ];
        [maxp, idxp] = max(abs(Fp(:)));
        [maxm, idxm] = max(abs(Fm(:)));
    else
        % Positive pivot dominates
        if ( absevp > absevm )
            kplus = kplus + 1;
            rankCount = rankCount + 1;

            colsPlus(:, kplus) = Fp(:, k);
            rowsPlus(kplus, :) = Fp(j, :);
            Fp = Fp - colsPlus(:, kplus)*(rowsPlus(kplus, :)*(1/evp));
            idxPlus(kplus) = rankCount;
            
            % Minus pivot is zero
            evm = 0;
            [maxp, idxp] = max(abs(Fp(:)));
            
        % Negative pivot dominates
        else
            kminus = kminus + 1;
            rankCount = rankCount + 1;

            colsMinus(:, kminus) = Fm(:, k);
            rowsMinus(kminus, :) = Fm(j, :);
            Fm = Fm - colsMinus(:, kminus)*(rowsMinus(kminus, :)*(1/evm));
            idxMinus(kminus) = rankCount;

            % Plus pivot is zero
            evp = 0;
            [maxm, idxm] = max(abs(Fm(:)));
        end
        pivotArray = [ pivotArray; [evp evm] ];
    end
end

if ( max(maxp, maxm) <= tol )
    isHappy = 1;                               % We are happy
end

if ( rankCount >= width )
    isHappy = 0;                               % We are not happy
end

% No sense in giving row and column values if they are not wanted.
if ( nargout > 4 )
    % Combine the types of pivots and set-up indices to track them
    cols = zeros(2*m-2, rankCount);
    rows = zeros(n, rankCount);
    pivots = zeros(rankCount, 1);
    if ( kplus ~= 0 )
        cols(m+1:2*m-2, idxPlus) = colsPlus;
        cols(2:m-1, idxPlus) = flipud(colsPlus);
        rows(:, idxPlus) = [ rowsPlus rowsPlus ].';
        pivotPlus = pivotArray(pivotArray(:,1) ~= 0,1);
        pivots(idxPlus) = pivotPlus;
    end
    
    if ( kminus ~= 0 )
        cols(m+1:2*m-2, idxMinus) = colsMinus;
        cols(2:m-1, idxMinus) = -flipud(colsMinus);
        rows(:, idxMinus) = [ -rowsMinus rowsMinus ].';
        pivotMinus = pivotArray(pivotArray(:,2) ~= 0,2);
        pivots(idxMinus) = pivotMinus;
    end
    
    if removePole
        cols(:, 1) = [ flipud(colPole); colPole(2:m-1)];
        rows(:, 1) = [ rowPole rowPole ];
        pivots(1) = rowVal;
    end

end

% Adjust the pivot locations so that they now correspond to F having the 
% poles.
if ( ~isempty(pivotIndices) )
    pivotIndices(:, 1) = pivotIndices(:, 1) + 1;
end

% Put the poles at the begining of the pivot locations array and also 
% include the pivot matrix.
if ( removePole )
    pivotIndices = [ 1 poleCol; pivotIndices ];
    pivotArray = [ [ rowVal 0 ]; pivotArray ];
    idxPlus = [ 1 idxPlus ];
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cols, pivots, rows, pivotLocations, idxPlus, idxMinus, isHappy, failure] = ...
    PhaseTwo(h, pivotIndices, pivotArray, n, dom, vscale, maxSample, ...
    removePoles, vectorize, pref)

happy_columns = 0;   % Not happy, until proven otherwise.
happy_rows = 0;
m = n;

% Set up the structs for the column and row trigtechs.
colData.hscale = norm(dom(3:4), inf);
colData.vscale = vscale;
rowData.hscale = norm(dom(1:2), inf);
rowData.vscale = vscale;

[x, y] = getPoints(m, n, dom);

rk = size(pivotIndices, 1);
id_rows = pivotIndices(:, 1); 
id_cols = pivotIndices(:, 2);

row_pivots = y(id_rows);
col_pivots = x(id_cols);

numPosPivots = sum(abs(pivotArray(:,1)) > 0);
numMinusPivots = sum(abs(pivotArray(:,2)) > 0);
totalPivots = numPosPivots + numMinusPivots;

pivotPlus = zeros(numPosPivots, 1);
pivotMinus = zeros(numMinusPivots, 1);
pivots = zeros(totalPivots, 1);
idxPlus = zeros(numPosPivots, 1);
idxMinus = zeros(numMinusPivots, 1);

% Phase 2: Calculate decomposition on sphere.
failure = false;
while ( ~(happy_columns && happy_rows) && ~failure )
    
    [x, y] = getPoints(m, n, dom);
    [xx, yy] = meshgrid(col_pivots, y);
    newCols = real(evaluate(h, xx + pi, yy, vectorize)); 
    temp = real(evaluate(h, xx, yy, vectorize));
    newColsPlus = 0.5*(newCols + temp);
    newColsMinus = 0.5*(newCols - temp);
    
    [xx, yy] = meshgrid(x, row_pivots);
    newRows = real(evaluate(h, xx, yy, vectorize));

    % This code will be unnecessary once ticket #1532 is addressed on the
    % chebfun tracker.  Don't forget to remove it.
    if ( numel(row_pivots) == 1 )
        newRows = newRows(:).';
    end
    
    newRowsPlus = 0.5*(newRows(:, n+1:2*n) + newRows(:, 1:n));
    newRowsMinus = 0.5*(newRows(:, n+1:2*n) - newRows(:, 1:n));
    
    colsPlus = zeros(m+1, numPosPivots);
    colsMinus = zeros(m+1, numMinusPivots);
    rowsPlus = zeros(numPosPivots, n);
    rowsMinus = zeros(numMinusPivots, n);
    plusCount = 1;
    minusCount = 1;
    pivotCount = 1;

    % Need to remove pole, which means we use the column with the largest
    % max norm (repeated) with rows of all ones in the elimination
    % algorithm.
    if removePoles
        newRowsPlus(1, :) = pivotArray(1, 1);
    end
    
    for ii = 1:rk
        
        % Get the pivots
        evp = pivotArray(ii, 1);
        evm = pivotArray(ii, 2);
                                        
        % Do GE step on both matrices
        if ( (evp ~=0) && (evm ~= 0) )
            colPlus = newColsPlus(:, ii);
            rowPlus = newRowsPlus(ii, :);
            
            colMinus = newColsMinus(:, ii);
            rowMinus = newRowsMinus(ii, :);
            
            % Store the columns and rows
            colsPlus(:, plusCount) = colPlus;
            rowsPlus(plusCount, :) = rowPlus;
            pivotPlus(plusCount) = evp;
            
            colsMinus(:, minusCount) = colMinus;
            rowsMinus(minusCount, :) = rowMinus;
            pivotMinus(minusCount) = evm;

            newColsPlus = newColsPlus - ...
                    colPlus*(rowPlus(id_cols)*(1/evp));
            newRowsPlus = newRowsPlus - ...
                    ((1/evp)*colPlus(id_rows))*rowPlus;
            newColsMinus = newColsMinus - ...
                    colMinus*(rowMinus(id_cols)*(1/evm));
            newRowsMinus = newRowsMinus - ...
                    ((1/evm)*colMinus(id_rows))*rowMinus;
                
            if ( abs(evp) >= abs(evm) )
                idxPlus(plusCount) = pivotCount;
                idxMinus(minusCount) = pivotCount + 1;
                pivots(pivotCount) = evp;
                pivots(pivotCount+1) = evm;
            else
                idxMinus(minusCount) = pivotCount;
                idxPlus(plusCount) = pivotCount + 1;
                pivots(pivotCount) = evm;
                pivots(pivotCount+1) = evp;
            end 
            plusCount = plusCount + 1;
            minusCount = minusCount + 1;
            pivotCount = pivotCount + 2;
        else
            if ( evp ~= 0 )
                colPlus = newColsPlus(:, ii);
                rowPlus = newRowsPlus(ii, :);

                % Store the columns and rows
                colsPlus(:, plusCount) = colPlus;
                rowsPlus(plusCount, :) = rowPlus;
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
                colMinus = newColsMinus(:, ii);
                rowMinus = newRowsMinus(ii, :);

                if any(isnan(newRowsMinus(:))) || any(isinf(newRowsMinus(:)))
                    fprintf('Shit\n');
                end
                
                % Store the columns and rows
                colsMinus(:, minusCount) = colMinus;
                rowsMinus(minusCount, :) = rowMinus;
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
    
    if ( removePoles )
        colsPlus(1, 2:end) = 0;
        colsPlus(end, 2:end) = 0;
    elseif ( ~isempty(colsPlus) )
        colsPlus(1,:) = 0;
        colsPlus(end,:) = 0;
    end
    
    if ~isempty(colsMinus)
        colsMinus(1, :) = 0;
        colsMinus(end, :) = 0;
    end
    
    % Happiness check for columns:
    % Double up the columns
    temp1 = sum([ colsPlus colsMinus ], 2); 
    temp2 = sum([ colsPlus -colsMinus ], 2);

    colValues = [ flipud(temp2); temp1(2:m) ];
    colTrigtech = trigtech.make(colValues, colData);
    happy_columns = happinessCheck(colTrigtech, [], colValues, colData, pref);
    
    % Happiness check for rows:
    % Double up the rows
    temp1 = sum([ rowsPlus; rowsMinus ], 1); 
    temp2 = sum([ rowsPlus; -rowsMinus ], 1);

    rowValues = [ temp1 temp2 ].';
    rowTrigtech = trigtech.make(rowValues, rowData);
    happy_rows = happinessCheck(rowTrigtech, [], rowValues, rowData, pref);
    
    % Adaptive:
    if( ~happy_columns )
        m = 2*m;
        ii = [ 1:2:m-1 m+2:2:2*m ];
        id_rows = ii(id_rows);
    end
    
    if ( ~happy_rows )
        n = 2*n;
        id_cols = 2*id_cols - 1;
    end
    
    if ( max(m, n) >= maxSample ) 
        warning('SPHEREFUN:constructor:notResolved', ...
        'Unresolved with maximum length: %u.', ...
        maxSample + mod(maxSample+1,2) );
        % Note the reason for adding one to maxSample when it's even: Since
        % we call simplify after construction the columns will always be an
        % odd length, but maxSample should always be even since we resample
        % in powers of 2. The warning message needs to have the right value
        % in the case that maxSample is even.
        failure = true;
    end 
end

% Combine the types of pivots and set-up indices to track them
cols = zeros(2*size(colsPlus, 1)-2, totalPivots);
cols(:, idxPlus) = [ flipud(colsPlus); colsPlus(2:end-1, :)  ];
cols(:, idxMinus) = [ -flipud(colsMinus); colsMinus(2:end-1, :)  ];

rows = zeros(2*size(rowsPlus, 2), totalPivots);
rows(:, idxPlus) = [ rowsPlus rowsPlus ].';
rows(:, idxMinus) = [ -rowsMinus rowsMinus ].';

pivotLocations = [ col_pivots row_pivots ];

isHappy = happy_rows & happy_columns;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vals = evaluate(h, xx, yy, vectorize)
% Evaluate h on an m-by-n tensor product grid.

% Evaluate h on the non-doubled up grid
if ( vectorize )
    vals = zeros(size( yy, 1), size( xx, 2));
    for jj = 1:size(yy, 1)
        for kk = 1:size(xx, 2)
            vals(jj, kk) = feval(h, xx(1, kk) , yy(jj, 1) );
        end
    end
else
    vals = feval(h, xx, yy );
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, y] = getPoints(m, n, dom)

colat = [ -pi pi 0 pi ];   % Colatitude (doubled up)
lat = [ -pi pi -pi/2 pi/2 ]; % Latitude (doubled up)

% Sample at an even number of points so that the poles are included.
if ( all((dom - colat) == 0) )
    x = trigpts(2*n, [-pi, pi]);   % azimuthal angle, lambda
    y = linspace(0, pi, m+1).';   % elevation angle, theta
elseif ( all((dom - lat) == 0) )
    x = trigpts(2*n, [-pi, pi]);          % azimuthal angle, lambda
    y = linspace(-pi/2,pi/2, m+1).';
else
    error('SPHEREFUN:constructor:points2D:unkownDomain', ...
        'Unrecognized domain.');
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [row, col] = myind2sub(siz, ndx)
% Alex's version of ind2sub. In2sub is slow because it has a varargout. 
% Since this is at the very inner part of the constructor and slowing 
% things down we will make our own. This version is about 1000 times faster
% than MATLAB ind2sub.

vi = rem(ndx - 1, siz(1)) + 1 ;
col = (ndx - vi ) / siz(1) + 1;
row = (vi - 1) + 1;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tol, vscale] = getTol(F, hx, hy, dom, pseudoLevel)
% GETTOL     Calculate a tolerance for the spherefun constructor.
%
%  This is the 2D analogue of the tolerance employed in the trigtech
%  constructors. It is based on a finite difference approximation to the
%  gradient, the size of the approximation domain, the internal working
%  tolerance, and an arbitrary (2/3) exponent. 

[m, n] = size(F);
grid = max(m, n);

% Remove some edge values so that df_dx and df_dy have the same size. 
dfdx = diff(F(1:m-1, :), 1, 2) / hx; % xx diffs column-wise.
dfdy = diff(F(:, 1:n-1), 1, 1) / hy; % yy diffs row-wise.

% An approximation for the norm of the gradient over the whole domain.
Jac_norm = max(max(abs(dfdx(:)), abs(dfdy(:))) );
vscale = max(abs(F(:)));
tol = grid.^(2/3) * max(abs(dom(:))) * max(Jac_norm, vscale) * pseudoLevel;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pivLocNew = adjustPivotLocations(pivLoc, pivArray)

% We will store the pivotLocations for both the plus and minus pieces, 
% which could result in duplicate values being stored.  This happens 
% whenever a rank-2 deflation step occurs. Really only one set of
% pivotLocations need to be stored in this case, but we will double up the
% information as it makes the combine and partition methods much easier to
% write. If this is changed then look the combine and partition methods
% need to also be changed.

pivLocNew = zeros(sum(sum(pivArray ~= 0)), 2);
count = 1;
for j = 1:size(pivLoc, 1)
    if ( (pivArray(j, 1) ~= 0) && (pivArray(j,2) ~= 0) )
        pivLocNew(count, :) = pivLoc(j, :);
        pivLocNew(count+1, :) = pivLoc(j, :);
        count = count + 2;
    else
        pivLocNew(count, :) = pivLoc(j, :);
        count = count + 1;
    end
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pole, constVal] = checkPole(val, tol)
% Check that the values at the pole are constant.

% Take the mean of the values that are at the pole.
pole = mean(val);

% Compute their standard deviation
stddev = std(val);

% If the standard deviation does not exceed the 100*tolearnce then the pole
% is "constant".
% TODO: Get a better feel for the tolerance check.
if ( (stddev <= 1e8*tol) || (stddev < eps) )
    constVal = 1;
else
    constVal = 0;
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [op, dom, pref, fixedRank, vectorize] = parseInputs(op, varargin)

if ( isa(op, 'char') )     % SPHEREFUN( CHAR )
    op = str2op(op);
end

% If the operator has one argument, then throw an error
if ( isa(op, 'function_handle') )
    % If the operator has one argument, then throw an error
    if ( nargin(op) <= 1 )
        error('CHEBFUN:SPHEREFUN:CONSTRUCTOR:toFewInputArgs',...
            'The function %s must accept 2 or 3 input arguments.',op);
    % If f is defined in terms of x,y,z; then convert it to
    % (longitude,latitude).
    elseif ( nargin(op) == 3 )
        % Wrap op so it can be evaluated in spherical coordinates
        op = @(lam, th) spherefun.sphf2cartf(op, lam, th, 0);
    end
end

% Get the domain: (Always first if given)
% The only domain presently supported is [-pi pi 0 pi], which corresponds
% to co-latitude in spherical coordinates.
% If domain is empty take it to be co-latitude.
dom = [-pi, pi, 0, pi]; 
fixedRank = NaN;
fixedLength = 0;

while ( numel(varargin) > 0 && isnumeric(varargin{1}) )
    d = varargin{1};
    varargin(1) = [];
    
    if ( numel(d) == 4 )                 % SPHEREFUN(OP, [A B C D])
        dom = d;
        if ( norm(dom(:)' - [-pi pi 0 pi]) > 0 )
            error('CHEBFUN:SPHEREFUN:CONSTRUCTOR:domain',...
                ['The only domain presently supported in spherefun is [-pi pi]x[0 pi] in '...
                'intrinstic (spherical) coordinates, which corresponds to colatitude.']);
        end
    elseif ( numel(d) == 2 )             % SPHEREFUN(OP, [M N])
        % Interpret this as the user wants a degree (M,N)
        % spherefun
        fixedLength = 1;
        m = d(1); 
        n = d(2);        
    elseif ( numel(d) == 1 )             % SPHEREFUN(OP, K)
        fixedRank = d;
    else
        error('CHEBFUN:SPHEREFUN:CONSTRUCTOR:domain',... 
              ['A domain is rarely given for spherefun, ', ... 
              'but it needs to be given by four corner values',... 
              'in intrinstic coordinates.'])
    end
end

if ( fixedLength )  % Check that m and n are positive integers
    if ( ( m <= 0 ) || ( n <= 0 ) || ( abs(round(m)-m)  > eps ) || ...
            ( abs(round(n)-n) > eps ) )
        error('CHEBFUN:SPHEREFUN:constructor:parseInputs:domain2', ...
            ['When constructing with fixed lengths, the values '...
             'for the lengths must be positive integers.']);
    end
end

if ( ( fixedRank < 0 ) || ( abs(round(fixedRank)-fixedRank) > eps ) )
        error('CHEBFUN:SPHEREFUN:constructor:parseInputs:domain3', ...
            ['When constructing with a fixed rank, the value must '...
             'be a positive integer.']);
end    

% Preferences structure given?
isPref = find(cellfun(@(p) isa(p, 'chebfunpref'), varargin));
if ( any(isPref) )
    pref = varargin{isPref};
    varargin(isPref) = [];
else
    pref = chebfunpref();
end

isEpsGiven = find(cellfun(@(p) strcmpi(p, 'eps'), varargin));
if ( isEpsGiven )
    pseudoLevel = varargin{isEpsGiven+1};
    varargin(isEpsGiven+(0:1)) = [];
else
    pseudoLevel = 0;
end
pref.cheb2Prefs.chebfun2eps = max(pref.cheb2Prefs.chebfun2eps, pseudoLevel);

% Look for vectorize flag:
vectorize = find(cellfun(@(p) strncmpi(p, 'vectori', 7), varargin));
if ( vectorize )
    varargin(vectorize) = [];
    vectorize = true;
else
    vectorize = false;
end

% If the vectorize flag is off, do we need to give user a warning?
if ( ~vectorize && ~isnumeric(op) ) % another check
    [vectorize, op] = vectorCheck(op, dom, pref.chebfuneps);
end

isCoeffs = find(cellfun(@(p) strcmpi(p, 'coeffs'), varargin));
if ( isCoeffs )
    varargin(isCoeffs) = [];
    op = spherefun.coeffs2spherefun(op);
end

% Deal with SPHEREFUN(OP, [M N]) now that all the other things are set.
if ( fixedLength )
    [x, y] = getPoints(m, n, dom);
    [xx, yy] = meshgrid(x, y);
    % Handle the special case of the input being a spherefun.  We can't call
    % evaluate here because, we have to use feval(op,xx,yy).
    if ( isa(op,'spherefun') )
        op = feval(op, xx, yy);
    else
        op = evaluate(op, xx, yy, vectorize);
    end    
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = fixTheRank( g , fixedRank )
% Fix the rank of a SPHEREFUN. Used for nonadaptive calls to the constructor.

if ( fixedRank < 0 )
    error('CHEBFUN:SPHEREFUN:constructor:fixTheRank:negative', ...
        'Nonadaptive rank should be positive.')
elseif ( fixedRank > 0 )
    if ( length(g.pivotValues) > fixedRank )
        % Truncate things:
        g.cols = g.cols(:,1:fixedRank);
        g.rows = g.rows(:,1:fixedRank);
        g.pivotValues = g.pivotValues(1:fixedRank);
        g.idxPlus = g.idxPlus( g.idxPlus <= fixedRank );
        g.idxMinus = g.idxMinus( g.idxMinus <= fixedRank );
    elseif ( length(g.pivotValues) < fixedRank )
        % Pad things with zero columns:
        zcols = chebfun(0, g.cols.domain, 'trig');
        zrows = chebfun(0, g.rows.domain, 'trig');
        for jj = length(g.pivotValues) : fixedRank - 1
            g.cols = [g.cols zcols];
            g.rows = [g.rows zrows];
            g.pivotValues = [g.pivotValues ; 0];
        end
    end
elseif ( fixedRank == 0 )
    g.cols = chebfun(0, g.cols.domain, 'trig');
    g.rows = chebfun(0, g.rows.domain, 'trig'); 
    g.pivotValues = Inf;
    g.idxPlus = [];
    g.idxMinus = 1;
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vectorize, op] = vectorCheck(op, dom, pseudoLevel)
% Check for cases: @(x,y) x*y, and @(x,y) x*y'

vectorize = false;

if isa(op,'spherefun')
    return;
end

% Evaluate at a 2-by-2 grid on the interior of the domain.
[xx, yy] = meshgrid( dom(1:2)/3 + diff(dom(1:2))/3,...
                     dom(3:4)/2 + diff(dom(3:4))/3);
try
    A = op(xx, yy);
catch
    throwVectorWarning();
    vectorize = true;
    return
end
if ( isscalar(A) )
    op = @(x,y) op(x,y) + 0*x + 0*y;
    A = op(xx, yy);
end
B = zeros(2);
for j = 1:2
    for k = 1:2
        B(j,k) = op(xx(j,k), yy(j,k));
    end
end
if ( any(any( abs(A - B) > min( 1000*pseudoLevel, 1e-4 ) ) ) )
    % Function handle probably needs vectorizing.
    % Give user a warning and then vectorize.
    throwVectorWarning();
    vectorize = true;
end
    function throwVectorWarning()
        warning('CHEBFUN:SPHEREFUN:constructor:vectorize',...
            ['Function did not correctly evaluate on an array.\n', ...
            'Turning on the ''vectorize'' flag. Did you intend this?\n', ...
            'Use the ''vectorize'' flag in the SPHEREFUN constructor\n', ...
            'call to avoid this warning message.']);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function op = str2op( op )
% OP = STR2OP(OP), finds the dependent variables in a string and returns an op
% handle than can be evaluated.

depvar = symvar( op );
if ( numel(depvar) > 3)
    error('CHEBFUN:SPHEREFUN:constructor:str2op:depvars', ...
        'Too many dependent variables in string input.');
elseif ( numel(depvar) == 1 )
    % Treat as a complex variable:
    op = eval(['@(' real(depvar{1}) + 1i*imag(depvar{1}) ')' op]);
elseif ( numel(depvar) == 2 )
    op = eval(['@(' depvar{1} ',' depvar{2} ')' op]);
else
    op = eval(['@(' depvar{1} ',' depvar{2} ',' depvar{3} ')' op]);
end

end
