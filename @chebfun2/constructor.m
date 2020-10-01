function g = constructor(g, op, varargin)
%CONSTRUCTOR   The main CHEBFUN2 constructor.
%
% This code is when functions of two variables are represented as CHEBFUN2
% objects. A CHEBFUN2 object is a low rank representation and expresses a
% function as a sum of rank-0 or 1 outer products of univariate functions.
%
% The algorithm for constructing a CHEBFUN2 comes in two phases:
%
% PHASE 1: The first phase attempts to determine the numerical rank of the
% function by performing Gaussian elimination with complete pivoting on a tensor
% grid of sample values. GE is perform until the pivoting elements fall below
% machine precision.  At the end of this stage we have candidate pivot locations
% and pivot elements.
%
% PHASE 2: The second phase attempts to resolve the corresponding column and row
% slices by sampling along the slices and performing GE on the skeleton.
% Sampling along each slice is increased until the Chebyshev coefficients of the
% slice fall below machine precision.
%
% The algorithm is fully described in:
%  A. Townsend and L. N. Trefethen, An extension of Chebfun to two dimensions,
%  SISC, 35 (2013), C495-C518.
%
% See also CHEBFUN2.

% Copyright 2017 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

% Parse the inputs:
[op, dom, prefx, prefy, isEqui, isTrig, fixedRank, vectorize] = parseInputs(op, varargin{:});

% Preferences:
techx       = prefx.tech();
techy       = prefy.tech();
tprefx      = techx.techPref;
tprefy      = techy.techPref;
minSample   = [tprefx.minSamples, tprefy.minSamples];
maxSample   = [tprefx.maxLength, tprefy.maxLength];
sampleTest  = prefx.cheb2Prefs.sampleTest | prefy.cheb2Prefs.sampleTest;
maxRank     = [prefx.cheb2Prefs.maxRank, prefy.cheb2Prefs.maxRank];
pseudoLevel = min(prefx.cheb2Prefs.chebfun2eps, prefy.cheb2Prefs.chebfun2eps);

% minSample needs to be a power of 2 when building periodic CHEBFUN2 objects or
% a ones plus power of 2 otherwise.  See #1771.
minSample( isTrig) = 2.^(floor(log2(minSample(isTrig))));
minSample(~isTrig) = 2.^(floor(log2(minSample(~isTrig) - 1))) + 1;

factor  = 4; % Ratio between size of matrix and no. pivots.
isHappy = 0; % If we are currently unresolved.
failure = 0; % Reached max discretization size without being happy.

if ( isa(op, 'chebfun2') )  % CHEBFUN2( CHEBFUN2 )
    g = fixTheRank(op, fixedRank);
    return
end

% The 'equi' flag can be used only with numeric data:
if ( any(isEqui) && ~isa(op, 'double') )
    error('CHEBFUN:CHEBFUN2:constructor:equi', ...
        'The EQUI flag is valid only when constructing from numeric data');
end
% Deal with constructions from numeric data:
if ( isa(op, 'double') )    % CHEBFUN2( DOUBLE )
    g = constructFromDouble(op, dom, prefx, prefy, isEqui, isTrig);
    % Fix the rank:
    g = fixTheRank(g, fixedRank);
    return
end

while ( ~isHappy && ~failure )
    grid = minSample;
    
    % Sample function on a Chebyshev tensor grid:
    [xx, yy] = points2D(grid(1), grid(2), dom, prefx, prefy);
    vals = evaluate(op, xx, yy, vectorize);
    
    % Does the function blow up or evaluate to NaN?:
    vscale = max(abs(vals(:)));
    if ( isinf(vscale) )
        error('CHEBFUN:CHEBFUN2:constructor:inf', ...
            'Function returned INF when evaluated');
    elseif ( any(isnan(vals(:)) ) )
        error('CHEBFUN:CHEBFUN2:constructor:nan', ...
            'Function returned NaN when evaluated');
    end
    
    % Two-dimensional version of CHEBFUN's tolerance:
    [relTol, absTol] = getTol(xx, yy, vals, dom, pseudoLevel);
    prefx.chebfuneps = relTol;
    prefy.chebfuneps = relTol;
    
    %% %%% PHASE 1: %%%
    % Do GE with complete pivoting:
    [pivotVal, pivotPos, rowVals, colVals, iFail] = completeACA(vals, absTol, factor);
    
    strike = 1;
    % grid <= 4*(maxRank-1)+1, see Chebfun2 paper.
    while ( iFail && all(grid <= factor*(maxRank-1)+1) && strike < 3)
        % Refine sampling on tensor grid:
        grid(1) = gridRefine(grid(1), prefx);
        grid(2) = gridRefine(grid(2), prefy);
        [xx, yy] = points2D(grid(1), grid(2), dom, prefx, prefy);
        vals = evaluate(op, xx, yy, vectorize); % resample
        vscale = max(abs(vals(:)));
        % New tolerance:
        [relTol, absTol] = getTol(xx, yy, vals, dom, pseudoLevel);
        prefx.chebfuneps = relTol;
        prefy.chebfuneps = relTol;
        % New GE:
        [pivotVal, pivotPos, rowVals, colVals, iFail] = ...
                                 completeACA(vals, absTol, factor);
        % If the function is 0+noise then stop after three strikes.
        if ( abs(pivotVal(1)) < 1e4*vscale*relTol )
            strike = strike + 1;
        end
    end
    
    % If the rank of the function is above maxRank then stop.
    if ( any(grid > factor*(maxRank-1) + 1) )
        warning('CHEBFUN:CHEBFUN2:constructor:rank', ...
            'Not a low-rank function.');
        failure = 1;
    end
    
    % Check if the column and row slices are resolved.
    colData.hscale = norm(dom(3:4), inf);
    colData.vscale = vscale;
    colTech = techy.make(sum(colVals,2), colData);
    resolvedCols = happinessCheck(colTech, [], sum(colVals, 2), colData, prefy);
    rowData.hscale = norm(dom(1:2), inf);
    rowData.vscale = vscale;
    rowTech = techx.make(sum(rowVals.',2), rowData);
    resolvedRows = happinessCheck(rowTech, [], sum(rowVals.', 2), rowData, prefx);
    isHappy = resolvedRows & resolvedCols;
    
    % If the function is zero, set midpoint of domain as pivot location.
    if ( length(pivotVal) == 1 && pivotVal == 0 )
        pivPos = [0, 0];
        isHappy = 1;
    else
        pivPos = [xx(1, pivotPos(:, 2)); yy(pivotPos(:, 1), 1).'].';
        PP = pivotPos;
    end
    
    %% %%% PHASE 2: %%%
    % Now resolve along the column and row slices:
    n = grid(2);  m = grid(1);
    while ( ~isHappy && ~failure  )
        if ( ~resolvedCols )
            % Double sampling along columns
            [n, nesting] = gridRefine( n , prefy );
            [xx, yy] = meshgrid(pivPos(:, 1), myPoints(n, dom(3:4), prefy));
            colVals = evaluate(op, xx, yy, vectorize);
            % Find location of pivots on new grid (using nesting property).
            PP(:, 1) = nesting(PP(:, 1));
        else
            [xx, yy] = meshgrid(pivPos(:, 1), myPoints(n, dom(3:4), prefy));
            colVals = evaluate(op, xx, yy, vectorize);
        end
        if ( ~resolvedRows )
            [m, nesting] = gridRefine( m , prefx );
            [xx, yy] = meshgrid(myPoints(m, dom(1:2), prefx), pivPos(:, 2));
            rowVals = evaluate(op, xx, yy, vectorize);
            % find location of pivots on new grid  (using nesting property).
            PP(:, 2) = nesting(PP(:, 2));
        else
            [xx, yy] = meshgrid(myPoints(m, dom(1:2), prefx), pivPos(:, 2));
            rowVals = evaluate(op, xx, yy, vectorize);
        end
        
        % Do GE on the skeleton to update slices:
        nn = numel(pivotVal);
        for kk = 1:nn-1
            colVals(:, kk+1:end) = colVals(:, kk+1:end) -...
                colVals(:, kk)*(rowVals(kk, PP(kk+1:nn, 2))./pivotVal(kk));
            rowVals(kk+1:end, :) = rowVals(kk+1:end, :) -...
                colVals(PP(kk+1:nn, 1), kk)*(rowVals(kk, :)./pivotVal(kk));
        end
        
        % If function is of rank-1 then make rowValues a row vector:
        if ( nn == 1 )
            rowVals = rowVals(:).';
        end
        
        % Are the columns and rows resolved now?
        if ( ~resolvedCols )
            colTech = techy.make(sum(colVals,2));
            resolvedCols = happinessCheck(colTech,[],sum(colVals,2), colData, prefy);
        end
        if ( ~resolvedRows )
            rowTech = techx.make(sum(rowVals.',2));
            resolvedRows = happinessCheck(rowTech,[],sum(rowVals.',2), rowData, prefx);
        end
        isHappy = resolvedRows & resolvedCols;

        % STOP if degree is over maxLength:
        sampleCheck = ( [m n] < maxSample );
        if ( ~all(sampleCheck) )
            k = find( ~sampleCheck );
            warning('CHEBFUN:CHEBFUN2:constructor:notResolved', ...
                'Unresolved with maximum CHEBFUN length: %u.', maxSample(k(1)));
            failure = 1;
        end
        
    end
    
    % For some reason, on some computers simplify is giving back a scalar zero.
    % In which case the function is numerically zero. Artificially set the
    % columns and rows to zero.
    if ( (norm(colVals) == 0) || (norm(rowVals) == 0) )
        colVals = 0;
        rowVals = 0;
        pivotVal = Inf;
        pivPos = [0, 0];
        isHappy = 1;
    end
    
    % Construct a CHEBFUN2:
    g.pivotValues = pivotVal;
    g.cols = chebfun(colVals,   dom(3:4), prefy);
    g.rows = chebfun(rowVals.', dom(1:2), prefx);
    g.pivotLocations = pivPos;
    g.domain = dom;
    
    % Sample Test:
    if ( sampleTest )
        % wrap the op with evaluate in case the 'vectorize' flag is on:
        sampleOP = @(x,y) evaluate(op, x, y, vectorize);
        
        % Evaluate at points in the domain:
        pass = g.sampleTest(sampleOP, absTol, vectorize);
        if ( ~pass )
            % Increase minSamples and try again.
            minSample = [gridRefine(minSample(1), prefx), ...
                         gridRefine(minSample(2), prefy)];
            isHappy = 0;
        end
    end
    
end

% Simplifying rows and columns after they are happy.
mineps = min(prefx.chebfuneps, prefy.chebfuneps);
g = simplify( g, mineps );

% Fix the rank, if in nonadaptive mode.
g = fixTheRank( g , fixedRank );

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = constructFromDouble(op, dom, prefx, prefy, isEqui, isTrig)

g = chebfun2();

if ( ~any(isEqui) && (numel( op ) == 1) )
    % LNT wants this:
    g = constructor(g, @(x,y) op + 0*x, dom);
    return
end
    
% Calculate a tolerance and find numerical rank to this tolerance:
% The tolerance assumes the samples are from a function. It depends
% on the size of the sample matrix, hscale of domain, vscale of
% the samples, condition number of the function, and the accuracy
% target in chebfun2 preferences.
if ( ~isEqui(1) || isTrig(1) )
    x = myPoints(size(op,2), dom(1:2), prefx);
else
    x = linspace(dom(1), dom(2), size(op,2));
end
if ( ~isEqui(2) || isTrig(2) )
    y = myPoints(size(op,1), dom(3:4), prefy);
else
    y = linspace(dom(3), dom(4), size(op,1));
end
[xx, yy] = meshgrid(x, y);

mineps = min(prefx.cheb2Prefs.chebfun2eps, prefy.cheb2Prefs.chebfun2eps);
[relTol, absTol] = getTol(xx, yy, op, dom, mineps);
prefx.chebfuneps = relTol;
prefy.chebfuneps = relTol;

% Perform GE with complete pivoting:
[pivotValue, ~, rowValues, colValues] = completeACA(op, absTol, 0);

% Construct a CHEBFUN2:
g.pivotValues = pivotValue;
if ( ~isEqui(1) || isTrig(1) )
    g.rows = chebfun(rowValues.', dom(1:2), prefx);
else
    g.rows = chebfun(rowValues.', dom(1:2), 'equi', prefx);
end
if ( ~isEqui(2) || isTrig(2) )
    g.cols = chebfun(colValues,   dom(3:4), prefy);
else
    g.cols = chebfun(colValues,   dom(3:4), 'equi', prefy);
end

g.domain = dom;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pivotValue, pivotElement, rows, cols, ifail] = ...
    completeACA(A, absTol, factor)
% Adaptive Cross Approximation with complete pivoting. This command is
% the continuous analogue of Gaussian elimination with complete pivoting.
% Here, we attempt to adaptively find the numerical rank of the function.

% Set up output variables.
[nx, ny] = size(A);
width = min(nx, ny);        % Use to tell us how many pivots we can take.
pivotValue = zeros(1);      % Store an unknown number of Pivot values.
pivotElement = zeros(1, 2); % Store (j,k) entries of pivot location.
ifail = 1;                  % Assume we fail.

% Main algorithm
zRows = 0;                  % count number of zero cols/rows.
[infNorm, ind] = max( abs ( reshape(A,numel(A),1) ) );
[row, col] = myind2sub( size(A) , ind);

% Bias toward diagonal for square matrices (see reasoning below):
if ( ( nx == ny ) && ( max( abs( diag( A ) ) ) - infNorm ) > -absTol )
    [infNorm, ind] = max( abs ( diag( A ) ) );
    row = ind;
    col = ind;
end

scl = infNorm;

% The function is the zero function.
if ( scl == 0 )
    % Let's pass back the zero matrix that is the same size as A. 
    % This ensures that chebfun2( zeros(5) ) has a 5x5 (zero) coefficient 
    % matrix.  
    pivotValue = 0;
    rows = zeros(1, size(A,2));
    cols = zeros(size(A,1), 1);
    ifail = 0;
else
    rows(1,:) = zeros(1, size(A, 2));
    cols(:,1) = zeros(size(A, 1), 1);
end

while ( ( infNorm > absTol ) && ( zRows < width / factor) ...
        && ( zRows < min(nx, ny) ) )
    rows(zRows+1,:) = A(row,:);
    cols(:,zRows+1) = A(:,col);              % Extract the columns.
    PivVal = A(row,col);
    A = A - cols(:,zRows+1)*(rows(zRows+1,:)./PivVal); % One step of GE.
    
    % Keep track of progress.
    zRows = zRows + 1;                       % One more row is zero.
    pivotValue(zRows) = PivVal;              % Store pivot value.
    pivotElement(zRows,:)=[row col];         % Store pivot location.
    
    % Next pivot.
    [ infNorm , ind ] = max( abs ( A(:) ) ); % Slightly faster.
    [ row , col ] = myind2sub( size(A) , ind );
    
    % Have a bias towards the diagonal of A, so that it can be used as a test
    % for nonnegative definite functions. (Complete GE and Cholesky are the
    % same as nonnegative definite functions have an absolute maximum on the
    % diagonal, except there is the possibility of a tie with an off-diagonal
    % absolute maximum. Bias toward diagonal maxima to prevent this.)
    if ( ( nx == ny ) && ( max( abs( diag( A ) ) ) - infNorm ) > -absTol )
        [infNorm, ind] = max( abs ( diag( A ) ) );
        row = ind;
        col = ind;
    end
end

if ( infNorm <= absTol )
    ifail = 0;                               % We didn't fail.
end
if ( zRows >= (width/factor) )
    ifail = 1;                               % We did fail.
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [row, col] = myind2sub(siz, ndx)
% My version of ind2sub. In2sub is slow because it has a varargout. Since this
% is at the very inner part of the constructor and slowing things down we will
% make our own. This version is about 1000 times faster than MATLAB ind2sub.

vi = rem( ndx - 1, siz(1) ) + 1 ;
col = ( ndx - vi ) / siz(1) + 1;
row = ( vi - 1 ) + 1;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vals = evaluate( op, xx, yy, flag )
% EVALUATE  Wrap the function handle in a FOR loop if the vectorize flag is
% turned on.

if ( flag )
    vals = zeros( size( yy, 1), size( xx, 2) );
    for jj = 1 : size( yy, 1)
        for kk = 1 : size( xx , 2 )
            vals(jj, kk) = op( xx( 1, kk) , yy( jj, 1 ) );
        end
    end
else
    vals = op( xx, yy );  % Matrix of values at cheb2 pts.
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function op = str2op( op )
% OP = STR2OP(OP), finds the dependent variables in a string and returns an op
% handle than can be evaluated.

depvar = symvar( op );
if ( numel(depvar) > 2)
    error('CHEBFUN:CHEBFUN2:constructor:str2op:depvars', ...
        'Too many dependent variables in string input.');
elseif ( numel(depvar) == 1 )
    % Treat as a complex variable:
    op = eval(['@(' real(depvar{1}) + 1i*imag(depvar{1}) ')' op]);

elseif ( numel(depvar) == 2 ) 

    op = eval(['@(' depvar{1} ',' depvar{2} ')' op]);

elseif ( isempty( depvar ) ) 
    
    op = @(x,y) str2double(op) + 0*x;
    
else 
    error('CHEBFUN2:STR2OP','Function as string has too many independent variables.');
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = fixTheRank( g , fixedRank )
% Fix the rank of a CHEBFUN2. Used for nonadaptive calls to the constructor.

if ( fixedRank < 0 )
    error('CHEBFUN:CHEBFUN2:constructor:fixTheRank:negative', ...
        'Nonadaptive rank should be positive.')
elseif ( fixedRank )
    if ( length(g.pivotValues) > fixedRank )
        % Truncate things:
        g.cols = g.cols(:,1:fixedRank);
        g.rows = g.rows(:,1:fixedRank);
        g.pivotValues = g.pivotValues(1:fixedRank);
    elseif ( length(g.pivotValues) < fixedRank )
        % Pad things with zero columns:
        zcols = chebfun(0, g.cols.domain);
        zrows = chebfun(0, g.rows.domain);
        for jj = length(g.pivotValues) : fixedRank - 1
            g.cols = [g.cols zcols];
            g.rows = [g.rows zrows];
            g.pivotValues = [g.pivotValues 0];
        end
    end
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xx, yy] = points2D(m, n, dom, prefx, prefy)
% Get the sample points that correspond to the right grid for a particular
% technology.

if ( nargin == 4 )
    prefy = prefx;
end

% What tech am I based on?:
techx = prefx.tech();
techy = prefy.tech();

if ( isa(techx, 'chebtech2') )
    x = chebpts( m, dom(1:2), 2 );
elseif ( isa(techx, 'chebtech1') )
    x = chebpts( m, dom(1:2), 1 );
elseif ( isa(techx, 'trigtech') )
    x = trigpts( m, dom(1:2) );
else
    error('CHEBFUN:CHEBFUN2:constructor:points2D:tecType', ...
        'Unrecognized technology');
end

if ( isa(techy, 'chebtech2') )
    y = chebpts( n, dom(3:4), 2 );
elseif ( isa(techy, 'chebtech1') )
    y = chebpts( n, dom(3:4), 1 );
elseif ( isa(techy, 'trigtech') )
    y = trigpts( n, dom(3:4) );
else
    error('CHEBFUN:CHEBFUN2:constructor:points2D:tecType', ...
        'Unrecognized technology');
end

[xx, yy] = meshgrid( x, y );

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = myPoints(n, dom, pref)
% Get the sample points that correspond to the right grid for a particular
% technology.

% What tech am I based on?:
tech = pref.tech();

if ( isa(tech, 'chebtech2') )
    x = chebpts( n, dom, 2 );   % x grid.
elseif ( isa(tech, 'chebtech1') )
    x = chebpts( n, dom, 1 );   % x grid.
elseif ( isa(tech, 'trigtech') )
    x = trigpts( n, dom );   % x grid.
else
    error('CHEBFUN:CHEBFUN2:constructor:mypoints:techType', ...
        'Unrecognized technology');
end

end

function [grid, nesting] = gridRefine( grid, pref )
% Hard code grid refinement strategy for tech.

% What tech am I based on?:
tech = pref.tech();

% What is the next grid size?
if ( isa(tech, 'chebtech2') )
    % Double sampling on tensor grid:
    grid = 2^( floor( log2( grid ) ) + 1) + 1;
    nesting = 1:2:grid;
elseif ( isa(tech, 'trigtech') )
    % Double sampling on tensor grid:
    grid = 2^( floor( log2( grid ) + 1 ));
    nesting = 1:2:grid;
elseif ( isa(tech, 'chebtech1' ) )
    grid = 3 * grid;
    nesting = 2:3:grid;
else
    error('CHEBFUN:CHEBFUN2:constructor:gridRefine:techType', ...
        'Technology is unrecognized.');
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [relTol, absTol] = getTol(xx, yy, vals, dom, pseudoLevel)
% GETTOL     Calculate a tolerance for the Chebfun2 constructor.
%
%  This is the 2D analogue of the tolerance employed in the chebtech
%  constructors. It is based on a finite difference approximation to the
%  gradient, the size of the approximation domain, the internal working
%  tolerance, and an arbitrary (2/3) exponent.

[m, n] = size( vals );
grid = max( m, n );
dfdx = 0;
dfdy = 0;
if ( m > 1 && n > 1 )
    % Remove some edge values so that df_dx and df_dy have the same size.
    dfdx = diff(vals(1:m-1,:),1,2) ./ diff(xx(1:m-1,:),1,2); % xx diffs column-wise.
    dfdy = diff(vals(:,1:n-1),1,1) ./ diff(yy(:,1:n-1),1,1); % yy diffs row-wise.
elseif ( m > 1 && n == 1 )
    % Constant in x-direction
    dfdy = diff(vals,1,1) ./ diff(yy,1,1);
elseif ( m == 1 && n > 1 )
    % Constant in y-direction
    dfdx = diff(vals,1,2) ./ diff(xx,1,2);
end
% An approximation for the norm of the gradient over the whole domain.
Jac_norm = max( max( abs(dfdx(:)), abs(dfdy(:)) ) );
vscale = max( abs( vals(:) ) );
relTol = grid.^(2/3) * pseudoLevel; % this should be vscale and hscale invariant
absTol = max( abs(dom(:) ) ) * max( Jac_norm, vscale) * relTol;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [op, dom, prefx, prefy, isEqui, isTrig, fixedRank, vectorize] = parseInputs(op, varargin)

if ( isa(op, 'char') )     % CHEBFUN2( CHAR )
    op = str2op(op);
end

% If the operator has one argument, then make it complex.
if ( isa(op, 'function_handle') && (nargin(op) == 1) )
    op = @(x, y) op(x + 1i*y);
end

% Determine the domain, fixed rank, and fixed length
dom = [-1, 1, -1, 1];
fixedRank = 0;
fixedLength = 0;
while ( numel(varargin) > 0 && isnumeric(varargin{1}) )
    d = varargin{1};
    varargin(1) = [];
    
    if ( numel(d) == 4 )                 % CHEBFUN2(OP, [A B C D])
        dom = d;
        
    elseif ( numel(d) == 2 )
        if ( (numel(varargin) > 0) && isa(varargin{1}, 'double') )
            ends = varargin{1};
            if ( numel( ends ) == 2 )    % CHEBFUN2(OP, [A B], [C D])
                dom = [d(:) ; ends(:)].';
                varargin(1) = [];                
            elseif ( numel( ends ) == 1 || numel( ends ) == 4) 
                % CHEBFUN2(OP, [M N], K) or CHEBFUN2(OP, [M N], [A B C D])
                % Just get fixed length [M N].
                fixedLength = 1;
                m = d(1);
                n = d(2);
            else
                error('CHEBFUN:CHEBFUN2:constructor:parseInputs:domain', ...
                    'Domain not valid or fully determined.');
            end
        else                             % CHEBFUN2(OP, [M N])
            % Just get fixed length [M N].
            fixedLength = 1;
            m = d(1);
            n = d(2);
       end
    elseif ( numel(d) == 1 )             % CHEBFUN2(OP, K)
        fixedRank = d;        
    else
        error('CHEBFUN:CHEBFUN2:constructor:parseInputs:unknown', ...
            'Unknown optional arguments for constructor.');
    end
end

if ( fixedLength )  % Check that m and n are positive integers
    if ( ( m <= 0 ) || ( n <= 0 ) || ( abs(round(m))-m  > eps ) || ...
            ( abs(round(n))-n > eps ) )
        error('CHEBFUN:CHEBFUN2:constructor:parseInputs:domain2', ...
            ['When constructing with fixed lengths, the values '...
             'for the lengths must be positive integers.']);
    end
end

if ( ( fixedRank < 0 ) || ( abs(round(fixedRank))-fixedRank > eps ) )
        error('CHEBFUN:CHEBFUN2:constructor:parseInputs:domain3', ...
            ['When constructing with a fixed rank, the value must '...
             'be a positive integer.']);
end    
    
% Check for infinite domains:
if ( any(isinf(dom) ) )
    error('CHEBFUN:CHEBFUN2:constructor:parseInputs:infDomain', ...
        'Chebfun2 cannot approximation functions on infinite domains.');
end

% Preferences structure given?
isPref = cellfun(@(p) (iscell(p) && numel(p)<3 && any(cellfun(@(q) isa(q, 'chebfunpref'), p))) || ...
    isa(p, 'chebfunpref'), varargin);
prefx = chebfunpref();
prefy = chebfunpref();
if ( any(isPref) )
    args = varargin(isPref);
    varargin(isPref) = [];
    if ( iscell(args{end}) )
        if ( isa(args{end}{1}, 'chebfunpref') ), prefx = args{end}{1}; end
        if ( isa(args{end}{2}, 'chebfunpref') ), prefy = args{end}{2}; end
    else
        prefx = args{end};
        prefy = args{end};
    end
end

isEqui = [false false];
match = cellfun(@(p) any(strcmpi(p, {'equi', 'equix', 'equiy'})), varargin);
if ( any(match) )
    args = varargin(match);
    varargin(match) = [];
    if ( strcmpi(args{end}, 'equix') )
        isEqui = [true false];
    elseif ( strcmpi(args{end}, 'equiy') )
        isEqui = [false true];
    else
        isEqui = [true true];
    end
end

match = cellfun(@(p) any(strcmpi(p, {'trig', 'trigx', 'trigy', ...
    'periodic', 'periodicx', 'periodicy'})), varargin);
if ( any(match) )
    args = varargin(match);
    varargin(match) = [];
    if ( any(strcmpi(args{end}, {'trigx', 'periodicx'})) )
        isTrig = [true false];
    elseif ( any(strcmpi(args{end}, {'trigy', 'periodicy'})) )
        isTrig = [false true];
    else
        isTrig = [true true];
    end
    if ( isTrig(1) ), prefx.tech = @trigtech; end
    if ( isTrig(2) ), prefy.tech = @trigtech; end
else
    % Even if the user didn't supply the 'trig' flag, we could still be doing a
    % periodic construction if the tech preference is 'trigtech'.
    %
    % TODO:  The only reason this is necessary is because of the adjustments we
    % have to make to minSample in the main construction routine above.  Can we
    % avoid this?
    isTrig = [isa(prefx.tech(), 'trigtech'), ...
              isa(prefy.tech(), 'trigtech')];
end

isEpsGiven = find(cellfun(@(p) strcmpi(p, 'eps'), varargin));
if ( isEpsGiven )
    pseudoLevel = varargin{isEpsGiven+1};
    varargin(isEpsGiven+(0:1)) = [];
else
    pseudoLevel = 0;
end
prefx.cheb2Prefs.chebfun2eps = max(prefx.cheb2Prefs.chebfun2eps, pseudoLevel);
prefy.cheb2Prefs.chebfun2eps = max(prefy.cheb2Prefs.chebfun2eps, pseudoLevel);

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
    mineps = min(prefx.chebfuneps, prefy.chebfuneps);
    [vectorize, op] = vectorCheck(op, dom, mineps);
end

% Deal with fixed length construction CHEBFUN(OP,[M N])
if ( fixedLength && ~isnumeric(op) )
    x = myPoints(m, dom(1:2), prefx);
    y = myPoints(n, dom(3:4), prefy);
    [xx, yy] = meshgrid(x,y);
    % Handle the special case of the input being a chebfun2.  We can't call
    % evaluate here because, we have to use feval(op,xx,yy).
    if ( isa(op,'chebfun2') )
        op = feval(op, xx, yy);
    else
        op = evaluate(op, xx, yy, vectorize);
    end
end

isPadua = find(cellfun(@(p) strcmpi(p, 'padua'), varargin));
if ( isPadua )
    varargin(isPadua) = [];
    op = chebfun2.paduaVals2coeffs(op);
    op = chebfun2.coeffs2vals(op);
end

match = cellfun(@(p) any(strcmpi(p, {'coeffs', 'coeffsx', 'coeffsy'})), varargin);
if ( any(match) )
    args = varargin(match);
    varargin(match) = [];
    if ( strcmpi(args{end}, 'coeffsx') )
        isCoeffs = [true false];
    elseif ( strcmpi(args{end}, 'coeffsy') )
        isCoeffs = [false true];
    else
        isCoeffs = [true true];
    end

    % Get the coeffs2vals transform corresponding to each tech
    if ( ~( isa(prefx.tech(), 'chebtech2') || isa(prefy.tech(), 'chebtech2') || ...
            isa(prefx.tech(), 'chebtech1') || isa(prefy.tech(), 'chebtech1') || ...
            isa(prefx.tech(), 'trigtech')  || isa(prefy.tech(), 'trigtech') ) )
        error('CHEBFUN:CHEBFUN2:constructor:parseInputs:techType', ...
            'Unrecognized technology');
    end
    transformX = @(x) x;
    transformY = @(y) y;
    if ( isCoeffs(1) )
        transformX = str2func([func2str(prefx.tech) '.coeffs2vals']);
    end
    if ( isCoeffs(2) )
        transformY = str2func([func2str(prefy.tech) '.coeffs2vals']);
    end
    op = transformX( transformY( op ).' ).';
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vectorize, op] = vectorCheck(op, dom, pseudoLevel)
% Check for cases: @(x,y) x*y, and @(x,y) x*y'

vectorize = false;
[xx, yy] = meshgrid( dom(1:2), dom(3:4));

if ( isa(op,'chebfun2') )
    vectorize = false;
    return;
end

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
        B(j,k) = op(dom(j), dom(2+k));
    end
end
if ( any(any( abs(A - B.') > min( 1000*pseudoLevel, 1e-4 ) ) ) )
    % Function handle probably needs vectorizing.
    % Give user a warning and then vectorize.
    throwVectorWarning();
    vectorize = true;
end
    function throwVectorWarning()
        warning('CHEBFUN:CHEBFUN2:constructor:vectorize',...
            ['Function did not correctly evaluate on an array.\n', ...
            'Turning on the ''vectorize'' flag. Did you intend this?\n', ...
            'Use the ''vectorize'' flag in the CHEBFUN2 constructor\n', ...
            'call to avoid this warning message.']);
    end
end
