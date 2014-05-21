function g = constructor(g, op, domain, varargin)
%CONSTRUCTOR   The main CHEBFUN2 constructor.
%
% This code is when functions of two variables are represented as CHEBFUN2
% objects. A CHEBFUN2 object is a low rank representation and expresses a
% function as a sum of rank-0 or 1 outerproduct of univariate functions.
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

if ( nargin == 0 )          % CHEBFUN2( )
    return
end

if ( isa(op, 'chebfun2') )  % CHEBFUN2( CHEBFUN2 )
    g = op;
    return
end

% If domain is empty take [-1 1 -1 1]:
if ( nargin < 3 || isempty(domain) )
    domain = [-1 1 -1 1];
end

if ( nargin > 3 && isa(varargin{1}, 'chebfunpref') )
    defaults = chebfunpref();
    pref = chebfunpref.mergePrefs(defaults, varargin{1});
else
    pref = chebfunpref();
end

if ( isa(op, 'double') )    % CHEBFUN2( DOUBLE )
    if ( numel( op ) == 1 )
        % LNT wants this:
        g = constructor(g, @(x,y) op + 0*x, domain);
        % Look for coeffs flag:
    elseif ( any(strcmpi(domain, 'coeffs')) )
        op = chebfun2.coeffs2vals( op );
        g = chebfun2( op, varargin{:} );
        return
    elseif ( any(strcmpi(domain, 'padua')) )
        [ignored, op] = chebfun2.paduaVals2coeffs( op );
        g = chebfun2( op );
        return        
    elseif ( (nargin > 3) && (any(strcmpi(varargin{1}, 'coeffs'))) )
        op = chebfun2.coeffs2vals( op );
        g = chebfun2( op, domain );
        return
    elseif ( (nargin > 3) && (any(strcmpi(varargin{1}, 'padua'))) )
        [ignored, op] = chebfun2.paduaVals2coeffs( op, domain );
        g = chebfun2( op, domain );
        return        
    else
        % If CHEBFUN2(f, rk), then nonadaptive call:
        if ( numel(domain) == 1 )
            fixedRank = domain;
            domain = [-1 1 -1 1];
        else
            % Otherwise its an adaptive call:
            fixedRank = 0;
        end
        % Perform GE with complete pivoting:
        [pivotValue, ignored, rowValues, colValues] = CompleteACA(op, 0);
        % Construct a CHEBFUN2:
        g.pivotValues = pivotValue;
        g.cols = chebfun(colValues, domain(3:4) );
        g.rows = chebfun(rowValues.', domain(1:2) );
        g.domain = domain;
        % Did we have a nonadaptive construction?:
        g = fixTheRank(g, fixedRank);
    end
    return
end

if ( isa(op, 'char') )     % CHEBFUN2( CHAR )
    op = str2op( op );
end

% If the operator has one argument, then make it complex.
if ( nargin(op) == 1 )
    op = @(x, y) op( x + 1i*y );
end

% Look for vectorize flag:
vectorize = 0;
if (any(strcmpi(domain, 'vectorize')) || any(strcmpi(domain, 'vectorise')))
    vectorize = 1;
    domain = [-1 1 -1 1];
elseif ( (nargin > 3) && (any(strcmpi(varargin{1}, 'vectorize')) ||...
        any(strcmpi(varargin{1}, 'vectorise'))))
    vectorize = 1;
end

fixedRank = 0;
% If the domain isn't of length 4, search for the other 2 endpoints: For
% instance, allow CHEBFUN2( OP, [-1 1], [-1 1]).
if ( numel(domain) == 2 )
    if ( ( nargin > 3) && isa(varargin{1}, 'double') )
        ends = varargin{1};
        if ( numel( ends ) == 2 )
            domain = [domain(:) ; ends(:)].';
        else
            error('CHEBFUN2:CONSTRUCTOR:DOMAIN', 'Domain not fully determined.');
        end
    else
        error('CHEBFUN2:CONSTRUCTOR:DOMAIN', 'Domain not fully determined.');
    end
elseif ( numel(domain) == 1 )
    fixedRank = domain;
    domain = [-1 1 -1 1];
elseif ( numel(domain) ~= 4 )
    error('CHEBFUN2:CONSTRUCTOR:DOMAIN', 'Domain not fully determined.');
end

% Get default preferences from chebfunpref:
prefStruct = pref.cheb2Prefs;
maxRank = prefStruct.maxRank;
maxLength = prefStruct.maxLength;
pseudoLevel = prefStruct.eps;
sampleTest = prefStruct.sampleTest;
minsample = 17;   % minsample

% If the vectorize flag is off, do we need to give user a warning?
if ( vectorize == 0 ) % another check
    % Check for cases: @(x,y) x*y, and @(x,y) x*y'
    [xx, yy] = meshgrid( domain(1:2), domain(3:4));
    A = op(xx, yy);
    B = zeros(2);
    for j = 1:2
        for k = 1:2
            B(j,k) = op(domain(j), domain(2+k));
        end
    end
    if ( any(any( abs(A - B.') > min( 1000*pseudoLevel, 1e-4 ) ) ) )
        % Function handle probably needs vectorizing, give user a warning and
        % then vectorize.
        warning('CHEBFUN2:CTOR:VECTORIZE','Function did not correctly evaluate on an array. Turning on the ''vectorize'' flag. Did you intend this? Use the ''vectorize'' flag in the chebfun2 constructor call to avoid this warning message.');
        g = chebfun2(op, domain, 'vectorize');
        return
    end
end

isHappy = 0;
while ( ~isHappy )
    grid = minsample; 
    
    % Sample function on a Chebyshev tensor grid:
    [xx, yy] = chebfun2.chebpts2(grid, grid, domain);
    vals = evaluate(op, xx, yy, vectorize);
    
    % Does the function blow up or evaluate to NaN?:
    vscale = max(abs(vals(:)));
    if ( isinf(vscale) )
        error('CHEBFUN2:CTOR', 'Function returned INF when evaluated');
    elseif ( any(isnan(vals(:)) ) )
        error('CHEBFUN2:CTOR', 'Function returned NaN when evaluated');
    end
    
    % Two-dimensional version of CHEBFUN's tolerance:
    tol = grid.^(2/3) * max( max( abs(domain(:))), 1) * vscale * pseudoLevel;
    
    %%% PHASE 1: %%%
    % Do GE with complete pivoting:
    [pivotValue, pivotPosition, rowValues, colValues, iFail] = CompleteACA(vals, tol);
    
    strike = 1;
    % grid <= 4*(maxRank-1)+1, see Chebfun2 paper. 
    while ( iFail && grid <= 4*(maxRank-1)+1 && strike < 3)
        % Double sampling on tensor grid:
        grid = 2^( floor( log2( grid ) ) + 1) + 1;
        [xx, yy] = chebfun2.chebpts2(grid, grid, domain);
        vals = evaluate(op, xx, yy, vectorize); % resample
        vscale = max(abs(vals(:)));
        % New tolerance:
        tol = grid.^(2/3) * max( max( abs(domain(:))), 1) * vscale * pseudoLevel;
        % New GE:
        [pivotValue, pivotPosition, rowValues, colValues, iFail] = CompleteACA(vals, tol);
        % If the function is 0+noise then stop after three strikes.
        if ( abs(pivotValue(1))<1e4*vscale*tol )
            strike = strike + 1;
        end
    end
    
    % If the rank of the function is above maxRank then stop.
    if ( grid > 4*(maxRank-1)+1 )
        error('CHEBFUN2:CTOR', 'Not a low-rank function.');
    end
    
    % Check if the column and row slices are resolved.
    colChebtech = chebtech2(sum(colValues,2), domain(3:4) );
    resolvedCols = happinessCheck(colChebtech,[],sum(colValues,2));
    rowChebtech = chebtech2(sum(rowValues.',2), domain(1:2) );
    resolvedRows = happinessCheck(rowChebtech,[],sum(rowValues.',2));
    isHappy = resolvedRows & resolvedCols;
    
    % If the function is zero, set midpoint of domain as pivot location.
    if ( length(pivotValue) == 1 && pivotValue == 0 )
        PivPos = [0, 0];
        isHappy = 1;
    else
        PivPos = [xx(1, pivotPosition(:, 2)); yy(pivotPosition(:, 1), 1).'].';
        PP = pivotPosition;
    end
    
    %%% PHASE 2: %%%
    % Now resolve along the column and row slices:
    n = grid;  m = grid;
    while ( ~isHappy )
        if ( ~resolvedCols )
            % Double sampling along columns
            n = 2^( floor( log2( n ) ) + 1) + 1;
            [xx, yy] = meshgrid(PivPos(:, 1), chebpts(n, domain(3:4)));
            colValues = evaluate(op, xx, yy, vectorize);
            % Find location of pivots on new grid (using nesting property).
            oddn = 1:2:n;
            PP(:, 1) = oddn(PP(:, 1));
        else
            [xx, yy] = meshgrid(PivPos(:, 1), chebpts(n, domain(3:4)));
            colValues = evaluate(op, xx, yy, vectorize);
        end
        if ( ~resolvedRows )
            % Double sampling along rows
            m = 2^( floor( log2( m ) ) + 1 ) + 1;
            [xx, yy] = meshgrid(chebpts(m, domain(1:2)), PivPos(:, 2));
            rowValues = evaluate(op, xx, yy, vectorize);
            % find location of pivots on new grid  (using nesting property).
            oddm = 1:2:m;
            PP(:, 2) = oddm(PP(:, 2));
        else
            [xx, yy] = meshgrid(chebpts(m, domain(1:2)), PivPos(:, 2));
            rowValues = evaluate(op, xx, yy, vectorize);
        end
        
        % Do GE on the skeleton to update slices:
        nn = numel(pivotValue);
        for kk = 1:nn-1
            colValues(:, kk+1:end) = colValues(:, kk+1:end) -...
                colValues(:, kk)*(rowValues(kk, PP(kk+1:nn, 2))./pivotValue(kk));
            rowValues(kk+1:end, :) = rowValues(kk+1:end, :) -...
                colValues(PP(kk+1:nn, 1), kk)*(rowValues(kk, :)./pivotValue(kk));
        end
        
        % If function is on rank-1 then make rowValues a row vector:
        if ( nn == 1 )
            rowValues = rowValues(:).';
        end
        
        % Are the columns and rows resolved now?
        if ( ~resolvedCols )
            colChebtech = chebtech2(sum(colValues,2));
            resolvedCols = happinessCheck(colChebtech,[],sum(colValues,2));
        end
        if ( ~resolvedRows )
            rowChebtech = chebtech2(sum(rowValues.',2));
            resolvedRows = happinessCheck(rowChebtech,[],sum(rowValues.',2));
        end
        isHappy = resolvedRows & resolvedCols;
        
        % STOP if degree is over maxLength:
        if ( max(m, n) >= maxLength )
            error('CHEBFUN2:CTOR', 'Unresolved with maximum CHEBFUN length: %u.', maxLength);
        end
        
    end
    
    % For some reason, on some computers simplify is giving back a scalar zero.
    % In which case the function is numerically zero. Artifically set the
    % columns and rows to zero.
    if ( norm(colValues) == 0 || norm(rowValues) == 0)
        colValues = 0;
        rowValues = 0;
        pivotValue = 0;
        PivPos = [0, 0];
        isHappy = 1;
    end
    
    % Construct a CHEBFUN2:
    g.pivotValues = pivotValue;
    g.cols = chebfun(colValues, domain(3:4) );
    g.rows = chebfun(rowValues.', domain(1:2) );
    g.pivotLocations = PivPos;
    g.domain = domain;
    
    % Sample Test:
    if ( sampleTest )
        % Evaluate at arbitrary point in domain:
        r = 0.029220277562146;
        s = 0.237283579771521;
        r = (domain(2)+domain(1))/2 + r*(domain(2)-domain(1));
        s = (domain(4)+domain(3))/2 + s*(domain(4)-domain(3));
        if ( abs( op(r,s) - feval(g, r, s) ) > 1e5 * tol )
           % out of lives.
           warning('CHEBFUN2:SAMPLETEST:FAILURE',...
               'Function may not be resolved. Is the function discontinuous?');
        end
    end
    
end

% Fix the rank, if in nonadaptive mode.
g = fixTheRank( g , fixedRank );

end

function [pivotValue, pivotElement, rows, cols, ifail] = CompleteACA(A, tol)
% Adaptive Cross Approximation with complete pivoting. This command is
% the continuous analogue of Gaussian elimination with complete pivoting.
% Here, we attempt to adaptively find the numerical rank of the function.

% Set up output variables.
[nx, ny] = size(A);
width = min(nx, ny);        % Use to tell us how many pivots we can take.
pivotValue = zeros(1);      % Store an unknown number of Pivot values.
pivotElement = zeros(1, 2); % Store (j,k) entries of pivot location.
ifail = 1;                  % Assume we fail.
factor = 4*(tol > 0);       % ratio between size of matrix and no. pivots. If tol = 0, then do full no. of steps.

% Main algorithm
zRows = 0;                  % count number of zero cols/rows.
[infNorm, ind] = max( abs ( reshape(A,numel(A),1) ) );
[row, col] = myind2sub( size(A) , ind);

% Bias toward diagonal for square matrices (see reasoning below):
if ( ( nx == ny ) && ( max( abs( diag( A ) ) ) - infNorm ) > -tol )
    [infNorm, ind] = max( abs ( diag( A ) ) );
    row = ind; 
    col = ind;
end

scl = infNorm;

% If the function is the zero function.
if ( scl == 0 )
    pivotValue = 0;
    rows = 0;
    cols = 0;
    ifail = 0;
else
    rows(1,:) = zeros(1, size(A, 2));
    cols(:,1) = zeros(size(A, 1), 1);
end

while ( ( infNorm > tol ) && ( zRows < width / factor)...
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
    if ( ( nx == ny ) && ( max( abs( diag( A ) ) ) - infNorm ) > -tol )
        [infNorm, ind] = max( abs ( diag( A ) ) );
        row = ind; 
        col = ind;
    end
end

if ( infNorm <= tol )
    ifail = 0;                               % We didn't fail.
end
if (zRows >= width / factor)
    ifail = 1;                               % We did fail.
end

end



function [row, col] = myind2sub(siz, ndx)
% My version of ind2sub. In2sub is slow because it has a varargout. Since this
% is at the very inner part of the constructor and slowing things down we will
% make our own. This version is about 1000 times faster than MATLAB ind2sub.

vi = rem( ndx - 1, siz(1) ) + 1 ;
col = ( ndx - vi ) / siz(1) + 1;
row = ( vi - 1 ) + 1;

end


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


function op = str2op( op )
% OP = STR2OP(OP), finds the dependent variables in a string and returns an op
% handle than can be evaluated.

depvar = symvar( op );
if ( numel(depvar) > 2 )
    error('CHEBFUN2:fun2:depvars', 'Too many dependent variables in string input.');
end
op = eval(['@(' depvar{1} ',' depvar{2} ')' op]);

end

function g = fixTheRank( g , fixedRank )
% Fix the rank of a CHEBFUN2. Used for nonadaptive calls to the constructor.

if ( fixedRank < 0 )
    error('CHEBFUN2:CONSTRUCTOR','Nonadaptive rank should be positive.')
end

if ( fixedRank > 0 )
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

