function g = constructor(g, op, dom, varargin)
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

% Copyright 2014 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

% TODO: The input parsing of CHEBFUN2 needs to be improved.

if ( nargin == 0 )          % CHEBFUN2( )
    return
end

if ( isa(op, 'chebfun2') )  % CHEBFUN2( CHEBFUN2 )
    g = op;
    return
end

% If domain is empty take [-1 1 -1 1]:
if ( nargin < 3 || isempty(dom) )
    dom = [-1 1 -1 1];
end
% Get preferences:
if ( nargin > 3 && isa(varargin{1}, 'chebfunpref') )
    pref = chebfunpref(varargin{1});
else
    pref = chebfunpref();
end

if ( isa(dom, 'chebfunpref') )
    % Support chebfun2(op, pref);
    pref = dom;
    dom = [-1 1 -1 1];
end
    
% Get default preferences from chebfunpref:
prefStruct = pref.cheb2Prefs;
sampleTest = prefStruct.sampleTest;
maxRank = prefStruct.maxRank;

% Get default preferences from the techPref:
tech = pref.tech();
tpref = chebfunpref.mergeTechPrefs(pref, tech.techPref);
minSample = tpref.minSamples; 
maxSample = tpref.maxLength;
pseudoLevel = tpref.eps;

% Deal with periodic functions: 
if ( any(strcmpi(dom, 'periodic')) )
        % If periodic flag, then map chebfun2 with TRIGTECHs. 
        pref.tech = @trigtech;
        tpref = chebfunpref.mergeTechPrefs(pref, tech.techPref);
        minSample = tpref.minSamples;
        maxSample = tpref.maxLength;
        pseudoLevel = tpref.eps;
        dom = [-1 1 -1 1];
elseif ( (nargin > 3) && (any(strcmpi(varargin{1}, 'periodic'))) )
        % If periodic flag, then map chebfun2 with TRIGTECHs
        pref.tech = @trigtech; 
        tpref = chebfunpref.mergeTechPrefs(pref, tech.techPref);
        minSample = tpref.minSamples;
        maxSample = tpref.maxLength;
        pseudoLevel = tpref.eps;
end

% Deal with constructions from equally spaced data:
if ( any(strcmpi(dom, 'equi')) || ((nargin > 3) && (any(strcmpi(varargin{1}, 'equi')))) )
        % Equally spaced data: 
        if ( any(strcmpi(dom, 'equi')) ) 
            dom = [-1 1 -1 1];
        end
        % Calculate a tolerance and find numerical rank to this tolerance: 
        % The tolerance assumes the samples are from a function. It depends
        % on the size of the sample matrix, hscale of domain, vscale of
        % the samples, and the accuracy target in chebfun2 preferences. 
        grid = max( size( op ) ); 
        vscale = max( abs(op(:)) ); 
        tol = grid.^(2/3) * max( max( abs(dom(:))), 1) * vscale * pseudoLevel;
        [pivotValue, ignored, rowValues, colValues] = CompleteACA(op, tol, 0); % Do ACA on matrices
        
        % Make a chebfun2: 
        g.pivotValues = pivotValue;
        g.cols = chebfun(colValues, dom(3:4), 'equi' );
        g.rows = chebfun(rowValues.', dom(1:2), 'equi'  );
        g.domain = dom;
        return
end

if ( isa(op, 'double') )    % CHEBFUN2( DOUBLE )
    if ( numel( op ) == 1 )
        % LNT wants this:
        g = constructor(g, @(x,y) op + 0*x, dom);
        
    elseif ( any(strcmpi(dom, 'coeffs')) )
        % Look for coeffs flag:
        op = chebfun2.coeffs2vals( op );
        g = chebfun2( op, varargin{:} );
        return
    elseif ( any(strcmpi(dom, 'padua')) )
        op = chebfun2.paduaVals2coeffs( op );
        op = chebfun2.coeffs2vals( op );
        g = chebfun2( op, 'coeffs' );
        return
    elseif ( (nargin > 3) && (any(strcmpi(varargin{1}, 'coeffs'))) )
        op = chebfun2.coeffs2vals( op );
        g = chebfun2( op, dom );
        return
    elseif ( (nargin > 3) && (any(strcmpi(varargin{1}, 'padua'))) )
        op = chebfun2.paduaVals2coeffs( op, dom );
        op = chebfun2.coeffs2vals( op );
        g = chebfun2( op, dom );
        return
    else
        % If CHEBFUN2(f, rk), then nonadaptive call:
        if ( numel(dom) == 1 )
            fixedRank = dom;
            dom = [-1 1 -1 1];
        else
            % Otherwise its an adaptive call:
            fixedRank = 0;
        end

        % Calculate a tolerance and find numerical rank to this tolerance: 
        % The tolerance assumes the samples are from a function. It depends
        % on the size of the sample matrix, hscale of domain, vscale of
        % the samples, and the accuracy target in chebfun2 preferences. 
        grid = max( size( op ) ); 
        vscale = max( abs(op(:)) ); 
        tol = grid.^(2/3) * max( max( abs(dom(:))), 1) * vscale * pseudoLevel;
        
        % Perform GE with complete pivoting:
        [pivotValue, ignored, rowValues, colValues] = CompleteACA(op, tol, 0);
        
        % Construct a CHEBFUN2:
        g.pivotValues = pivotValue;
        g.cols = chebfun(colValues, dom(3:4), pref );
        g.rows = chebfun(rowValues.', dom(1:2), pref );
        g.domain = dom;
        
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
if (any(strcmpi(dom, 'vectorize')) || any(strcmpi(dom, 'vectorise')))
    vectorize = 1;
    dom = [-1 1 -1 1];
elseif ( (nargin > 3) && (any(strcmpi(varargin{1}, 'vectorize')) ||...
        any(strcmpi(varargin{1}, 'vectorise'))))
    vectorize = 1;
end

fixedRank = 0;
% If the domain isn't of length 4, search for the other 2 endpoints: For
% instance, allow CHEBFUN2( OP, [-1 1], [-1 1]).
if ( numel(dom) == 2 )
    if ( ( nargin > 3) && isa(varargin{1}, 'double') )
        ends = varargin{1};
        if ( numel( ends ) == 2 )
            dom = [dom(:) ; ends(:)].';
        elseif ( numel(ends) == 4 ) 
            % Interpret this as the user wants a degree (dom(1),dom(2)) 
            % chebfun2 on the domain [ends]. 
            [xx, yy] = chebfun2.chebpts2(dom(1), dom(2), ends);
            g = chebfun2( op(xx, yy), varargin{:} ); 
            return
        else 
            error('CHEBFUN:CHEBFUN2:constructor:domain1', ...
                'Domain not valid or fully determined.');
        end
    else
        % The domain is not given, but perhaps the user 
        % wants a degree (dom(1),dom(2)) representation.
        if ( dom(2) - dom(1) > 0 && dom(1)>0 &&...   % A valid bivariate degree? 
                abs(round(dom(1)) - dom(1))< eps &&...
                abs(round(dom(2)) - dom(2))< eps) 
            [xx, yy] = chebfun2.chebpts2(dom(1), dom(2));
            g = chebfun2( op(xx, yy), varargin ); 
            return
        else
            error('CHEBFUN:CHEBFUN2:constructor:domain2', ...
                'Domain not valid or fully determined.');
        end
    end
elseif ( numel(dom) == 1 )
    fixedRank = dom;
    dom = [-1 1 -1 1];
elseif ( numel(dom) ~= 4 )
    error('CHEBFUN:CHEBFUN2:constructor:DOMAIN', ...
        'Domain not fully determined.');
end

% If the vectorize flag is off, do we need to give user a warning?
if ( vectorize == 0 ) % another check
    % Check for cases: @(x,y) x*y, and @(x,y) x*y'
    [xx, yy] = meshgrid( dom(1:2), dom(3:4));
    A = op(xx, yy);
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
        % Function handle probably needs vectorizing, give user a warning and
        % then vectorize.
        
        warning('CHEBFUN:CHEBFUN2:constructor:vectorize',...
            ['Function did not correctly evaluate on an array.\n', ...
             'Turning on the ''vectorize'' flag. Did you intend this?\n', ...
             'Use the ''vectorize'' flag in the CHEBFUN2 constructor\n', ...
             'call to avoid this warning message.']);
        g = chebfun2(op, dom, 'vectorize', pref);
        return
    end
end

factor = 4;  % ratio between size of matrix and no. pivots. 
isHappy = 0; % If we are currently unresolved. 
failure = 0; % Reached max discretization size without being happy. 
while ( ~isHappy && ~failure )
    grid = minSample; 
    
    % Sample function on a Chebyshev tensor grid:
    [xx, yy] = points2D(grid, grid, dom, pref);
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
    tol = grid.^(2/3) * max( max( abs(dom(:))), 1) * vscale * pseudoLevel;
    
    %%% PHASE 1: %%%
    % Do GE with complete pivoting:
    [pivotValue, pivotPosition, rowValues, colValues, iFail] = CompleteACA(vals, tol, factor);
    
    strike = 1;
    % grid <= 4*(maxRank-1)+1, see Chebfun2 paper. 
    while ( iFail && grid <= factor*(maxRank-1)+1 && strike < 3)
        % Refine sampling on tensor grid:
        grid = gridRefine( grid , pref);
        [xx, yy] = points2D(grid, grid, dom, pref);
        vals = evaluate(op, xx, yy, vectorize); % resample
        vscale = max(abs(vals(:)));
        % New tolerance:
        tol = grid.^(2/3) * max( max( abs(dom(:))), 1) * vscale * pseudoLevel;
        % New GE:
        [pivotValue, pivotPosition, rowValues, colValues, iFail] = CompleteACA(vals, tol, factor);
        % If the function is 0+noise then stop after three strikes.
        if ( abs(pivotValue(1))<1e4*vscale*tol )
            strike = strike + 1;
        end
    end
    
    % If the rank of the function is above maxRank then stop.
    if ( grid > factor*(maxRank-1)+1 )
        warning('CHEBFUN:CHEBFUN2:constructor:rank', ...
            'Not a low-rank function.');
        failure = 1; 
    end
    
    % Check if the column and row slices are resolved.
    colData.vscale = dom(3:4);
    tech = pref.tech(); 
    colChebtech = tech.make(sum(colValues,2), colData);
    resolvedCols = happinessCheck(colChebtech,[],sum(colValues,2));
    rowData.vscale = dom(1:2);
    rowChebtech = tech.make(sum(rowValues.',2), rowData);
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
            [n, nesting] = gridRefine( n , pref );
            [xx, yy] = meshgrid(PivPos(:, 1), mypoints(n, dom(3:4), pref));
            colValues = evaluate(op, xx, yy, vectorize);
            % Find location of pivots on new grid (using nesting property).
            PP(:, 1) = nesting(PP(:, 1));
        else
            [xx, yy] = meshgrid(PivPos(:, 1), mypoints(n, dom(3:4), pref));
            colValues = evaluate(op, xx, yy, vectorize);
        end
        if ( ~resolvedRows )
            [m, nesting] = gridRefine( m , pref ); 
            [xx, yy] = meshgrid(mypoints(m, dom(1:2), pref), PivPos(:, 2));
            rowValues = evaluate(op, xx, yy, vectorize);
            % find location of pivots on new grid  (using nesting property).
            PP(:, 2) = nesting(PP(:, 2));
        else
            [xx, yy] = meshgrid(mypoints(m, dom(1:2), pref), PivPos(:, 2));
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
            colChebtech = tech.make(sum(colValues,2));
            resolvedCols = happinessCheck(colChebtech,[],sum(colValues,2));
        end
        if ( ~resolvedRows )
            rowChebtech = tech.make(sum(rowValues.',2));
            resolvedRows = happinessCheck(rowChebtech,[],sum(rowValues.',2));
        end
        isHappy = resolvedRows & resolvedCols;
        
        % STOP if degree is over maxLength:
        if ( max(m, n) >= maxSample )
            warning('CHEBFUN:CHEBFUN2:constructor:notResolved', ...
                'Unresolved with maximum CHEBFUN length: %u.', maxSample);
            failure = 1;
        end
        
    end
    
    % For some reason, on some computers simplify is giving back a scalar zero.
    % In which case the function is numerically zero. Artificially set the
    % columns and rows to zero.
    if ( (norm(colValues) == 0) || (norm(rowValues) == 0) )
        colValues = 0;
        rowValues = 0;
        pivotValue = Inf;
        PivPos = [0, 0];
        isHappy = 1;
    end
    
    % Construct a CHEBFUN2:
    g.pivotValues = pivotValue;
    g.cols = chebfun(colValues, dom(3:4), pref);
    g.rows = chebfun(rowValues.', dom(1:2), pref );
    g.pivotLocations = PivPos;
    g.domain = dom;
    
    % Sample Test:
    if ( sampleTest )
        % wrap the op with evaluate in case the 'vectorize' flag is on: 
        sampleOP = @(x,y) evaluate( op, x, y, vectorize);
        
        % Evaluate at points in the domain:
        pass = g.sampleTest( sampleOP, tol, vectorize);
        if ( ~pass )
            % Increase minSamples and try again.
            minSample = gridRefine( minSample, pref );
            isHappy = 0;
        end
    end
    
end

% Fix the rank, if in nonadaptive mode.
g = fixTheRank( g , fixedRank );

end

function [pivotValue, pivotElement, rows, cols, ifail] = ...
    CompleteACA(A, tol, factor)
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

while ( ( infNorm > tol ) && ( zRows < width / factor) ...
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
if ( zRows >= (width/factor) )
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
    error('CHEBFUN:CHEBFUN2:constructor:str2op:depvars', ...
        'Too many dependent variables in string input.');
end
op = eval(['@(' depvar{1} ',' depvar{2} ')' op]);

end


function g = fixTheRank( g , fixedRank )
% Fix the rank of a CHEBFUN2. Used for nonadaptive calls to the constructor.

if ( fixedRank < 0 )
    error('CHEBFUN:CHEBFUN2:constructor:fixTheRank:negative', ...
        'Nonadaptive rank should be positive.')
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


function [xx, yy] = points2D(m, n, dom, pref)
% Get the sample points that correspond to the right grid for a particular
% technology.

% What tech am I based on?:
tech = pref.tech();

if ( isa(tech, 'chebtech2') )
    x = chebpts( m, dom(1:2), 2 );   % x grid.
    y = chebpts( n, dom(3:4), 2 ); 
    [xx, yy] = meshgrid( x, y ); 
elseif ( isa(tech, 'chebtech1') )
    x = chebpts( m, dom(1:2), 1 );   % x grid.
    y = chebpts( n, dom(3:4), 1 ); 
    [xx, yy] = meshgrid( x, y ); 
elseif ( isa(tech, 'trigtech') )
    x = trigpts( m, dom(1:2) );   % x grid.
    y = trigpts( n, dom(3:4) );
    [xx, yy] = meshgrid( x, y );
else
    error('CHEBFUN:CHEBFUN2:constructor:points2D:tecType', ...
        'Unrecognized technology');
end

end


function x = mypoints(n, dom, pref)
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
