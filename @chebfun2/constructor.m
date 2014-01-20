function g = constructor(g, op, domain, varargin)
% The main Chebfun2 constructor.


% Remove two trivial cases:
if ( nargin == 0 )   % chebfun2();
    return
end

if ( isa(op, 'chebfun2') )  % chebfun2(f), f = chebfun2
    g = op;
    return
end

if ( nargin < 3 || isempty(domain) )
    domain = [-1 1 -1 1];
end

if ( nargin > 3)
    if ( isa(varargin{1}, 'chebpref') )
        defaults = chebpref;
        pref = mergePrefs(defaults, varargin{1});
    end
else
    pref = chebpref;
end

if ( isa(op, 'double') )                        % chebfun2( double )
    if ( numel( op ) == 1 )
        g = constructor(g, @(x,y) op + 0*x, domain);
    else
        op = flipud( op );
        [pivotValue, ignored, rowValues, colValues] = CompleteACA(op, 0);
        g.pivotValues = pivotValue;
        g.cols = chebfun(colValues, domain(3:4) );
        g.rows = chebfun(rowValues.', domain(1:2) );
        g.domain = domain;
    end
    return
end

if ( isa(op, 'char') )                          % chebfun2('fh')
    op = str2op( op );
end

% Check the operator has one argument, then make it complex.
if ( nargin(op) == 1 )
    op = @(x, y) op( x + 1i*y );
end

% Look for vectorize and coeffs flag:
vectorize = 0;
if (any(strcmpi(domain,'vectorize')) || any(strcmpi(domain,'vectorise')))
    vectorize = 1;
    domain = [-1 1 -1 1];
end
if ( (nargin > 3) && (any(strcmpi(varargin{1},'vectorize')) || any(strcmpi(varargin{1},'vectorise'))))
    vectorize = 1;
end

if (any(strcmpi(domain,'coeffs')) || any(strcmpi(domain,'coeffs')) )
    op = coeffs2vals( op );
    g = chebfun2( op );
    return
end
if (( nargin > 3 ) && ( any(strcmpi(varargin{1},'coeffs')) || any(strcmpi(varargin{1},'coeffs'))))
    op = coeffs2vals( op );
    g = chebfun2( op, domain );
    return
end

% If the domain isn't of length 4, search for the other 2 endpoints:
if ( numel(domain) == 2 )
    if ( ( nargin > 3) && isa(varargin{1}, 'double') )
        ends = varargin{1};
        if ( numel( ends ) == 2 )
            domain = [domain(:);ends(:)]';
        else
            error('CHEBFUN2:CONSTRUCTOR:DOMAIN','Domain not fully determined.');
        end
    else
        error('CHEBFUN2:CONSTRUCTOR:DOMAIN','Domain not fully determined.');
    end
elseif ( numel(domain) ~= 4 )
    error('CHEBFUN2:CONSTRUCTOR:DOMAIN','Domain not fully determined.');
end

% Get default preferences from chebPref:
prefs = chebpref;
prefStruct = prefs.cheb2Prefs;
% maxRank = prefStruct.maxRank;
maxRank = 1025;
maxLength = prefStruct.maxLength;
pseudoLevel = prefStruct.eps;
exactLength = prefStruct.exactLength;
sampleTest = prefStruct.sampleTest;
grid = 9;   % minsample

% Check if we need to turn on vectorize flag:
% m1 = mean( domain(1:2) );
% m2 = mean( domain(3:4) );
% 
% E = ones(2,1);
% if ( (vectorize == 0) && all( size(op(m1*E,m2*E)) == [1 1]) )   % scalar check
%     % sizes are not going to match so let's try with the vectorizeflag on.
%     if ~( numel(op(m1*E,m2*E))==1 && norm(op(m1*E,m2*E))==0 )
%         warning('CHEBFUN2:CTOR:VECTORIZE','Function did not correctly evaluate on an array. Turning on the ''vectorize'' flag. Did you intend this? Use the ''vectorize'' flag in the chebfun2 constructor call to avoid this warning message.');
%         g = chebfun2( op, domain, 'vectorize' );
%         return
%     end
% else
if ( vectorize == 0 )                                       % another check
    % check for cases: @(x,y) x*y, and @(x,y) x*y'
    [xx, yy] = meshgrid( domain(1:2), domain(3:4));
    A = op(xx, yy);
    B = zeros(2);
    for j = 1:2
        for k = 1:2
            B(j,k) = op(domain(j), domain(2+k));
        end
    end
    if ( any(any( abs(A - B.') > min( 1000*pseudoLevel, 1e-4 ) ) ) )
        warning('CHEBFUN2:CTOR:VECTORIZE','Function did not correctly evaluate on an array. Turning on the ''vectorize'' flag. Did you intend this? Use the ''vectorize'' flag in the chebfun2 constructor call to avoid this warning message.');
        g = chebfun2(op, domain, 'vectorize');
        return
    end
end

isHappy = 0;  % If unhappy, selected pivots were not good enough.
while ( ~isHappy )
    [xx, yy] = chebfun2.chebpts2(grid, grid, domain);
    vals = evaluate(op, xx, yy, vectorize);             % Matrix of values at cheb2 pts.
    
    vscale = max(abs(vals(:)));
    if ( isinf(vscale) )
        error('FUN2:CTOR', 'Function returned INF when evaluated');
    end
    if ( any(isnan(vals(:)) ) )
        error('FUN2:CTOR', 'Function returned NaN when evaluated');
    end
    tol = grid.^(4/3) * max( max( abs(domain(:))), 1) * vscale * pseudoLevel;
    
    %% FIND NUMERICAL RANK:
    [pivotValue, pivotPosition, rowValues, colValues, iFail] = CompleteACA(vals, tol);
    % Use chebtech's happiness check
    strike = 1;
    while ( iFail && grid <= maxRank && strike < 3)
        grid = 2^(floor(log2(grid)) + 1) + 1;                % Double the sampling
        [xx, yy] = chebfun2.chebpts2(grid, grid, domain);
        vals = evaluate(op, xx, yy, vectorize);                        % Resample on denser grid.
        vscale = max(abs(vals(:)));
        tol = grid.^(4/3) * max( max( abs(domain(:))), 1) * vscale * pseudoLevel;
        [pivotValue, pivotPosition, rowValues, colValues, iFail] = CompleteACA(vals, tol);
        if ( abs(pivotValue(1))<1e4*vscale*tol )
            % If the function is 0+noise then stop after three strikes.
            strike = strike + 1;
        end
    end
    
    if ( grid > maxRank )
        error('FUN2:CTOR', 'Not a low-rank function.');
    end
    
    colChebtech = chebtech2(sum(colValues,2), domain(3:4) );
    resolvedCols = happinessCheck(colChebtech);
    rowChebtech = chebtech2(sum(rowValues.',2), domain(1:2) );
    resolvedRows = happinessCheck(rowChebtech);
    
    isHappy = resolvedRows & resolvedCols;
    
    if ( length(pivotValue) == 1 && pivotValue == 0 )
        PivPos = [0, 0];
        isHappy = 1;
    else
        PivPos = [xx(1, pivotPosition(:, 2)); yy(pivotPosition(:, 1), 1).'].';
        PP = pivotPosition;
    end
    
    n = grid;  m = grid;
    
    % If unresolved then perform ACA on selected slices.
    while ( ~isHappy )
        if ( ~resolvedCols )
            n = 2^(floor(log2(n))+1) + 1;
            [xx, yy] = meshgrid(PivPos(:, 1), chebpts(n, domain(3:4)));
            colValues = evaluate(op, xx, yy, vectorize);
            % Find location of pivots on new grid.
            oddn = 1:2:n;
            PP(:, 1) = oddn(PP(:, 1));
        else
            [xx, yy] = meshgrid(PivPos(:, 1), chebpts(n, domain(3:4)));
            colValues = evaluate(op, xx, yy, vectorize);
        end
        if ( ~resolvedRows )
            m = 2^(floor(log2(m))+1) + 1;
            [xx, yy] = meshgrid(chebpts(m, domain(1:2)), PivPos(:, 2));
            rowValues = evaluate(op, xx, yy, vectorize);
            oddm = 1:2:m; PP(:, 2) = oddm(PP(:, 2)); % find location of pivots on new grid.
        else
            [xx, yy] = meshgrid(chebpts(m, domain(1:2)), PivPos(:, 2));
            rowValues = evaluate(op, xx, yy, vectorize);
        end
        
        nn = numel(pivotValue);
        
        % ACA on selected Pivots.
        for kk = 1:nn-1
            colValues(:, kk+1:end) = colValues(:, kk+1:end) - colValues(:, kk)*(rowValues(kk, PP(kk+1:nn, 2))./pivotValue(kk));
            rowValues(kk+1:end, :) = rowValues(kk+1:end, :) - colValues(PP(kk+1:nn, 1), kk)*(rowValues(kk, :)./pivotValue(kk));
        end
        
        if ( nn == 1 )
            rowValues = rowValues(:).';  % make a row vector.
        end
        
        % Are the columns and rows resolved now?
        if ( ~resolvedCols )
            colChebtech = chebtech2(sum(colValues,2));
            resolvedCols = happinessCheck(colChebtech);
        end
        if ( ~resolvedRows )
            rowChebtech = chebtech2(sum(rowValues.',2));
            resolvedRows = happinessCheck(rowChebtech);
        end
        
        isHappy = resolvedRows & resolvedCols;
        if ( max(m, n) >= maxLength )  % max number of degrees allows.
            error('FUN2:CTOR', 'Unresolved with maximum Chebfun length: %u.', maxLength);
        end
        
    end
    
    % For some reason, on some computers simplify is giving back a
    % scalar zero.  In which case the function is numerically zero.
    % Artifically set the columns and rows to zero.
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
        r = 0.029220277562146; s = 0.237283579771521;
        r = (domain(2)+domain(1))/2 + r*(domain(2)-domain(1));
        s = (domain(4)+domain(3))/2 + s*(domain(4)-domain(3));
        if ( abs( op(r,s) - feval(g, r, s) ) > 1e5 * tol )
            isHappy = 0;
        end
    end
    
end

end

function [PivotValue, PivotElement, Rows, Cols, ifail] = CompleteACA(A, tol)
% Adaptive Cross Approximation with complete pivoting. This command is
% the continuous analogue of Gaussian elimination with complete pivoting.
% Here, we attempt to adaptively find the numerical rank of the function.

% Set up output variables.
[nx,ny]=size(A);
width = min(nx,ny);         % Use to tell us how many pivots we can take.
PivotValue = zeros(1);      % Store an unknown number of Pivot values.
PivotElement = zeros(1,2);  % Store (j,k) entries of pivot location.
ifail = 1;                  % Assume we fail.
factor = 4*(tol>0);         % ratio between size of matrix and no. pivots. If tol = 0, then do full no. of steps.

% Main algorithm
zrows = 0;                  % count number of zero cols/rows.
[ infnorm , ind ]=max( abs ( reshape(A,numel(A),1) ) );
[ row , col ]=myind2sub( size(A) , ind);
scl = infnorm;

% If the function is the zero function.
if scl == 0
    PivotValue=0;
    Rows = 0; Cols = 0;
    ifail = 0;
else
    Rows(1, :) = zeros(1, size(A, 2));
    Cols(:, 1) = zeros(size(A, 1), 1);
end
while ( ( infnorm > tol ) && ( zrows < width / factor) )
    Rows(zrows+1,:) = A( row , : ) ;
    Cols(:,zrows+1) = A( : , col ) ;    % Extract the columns out
    PivVal = A(row,col);
    A = A - Cols(:,zrows+1)*(Rows(zrows+1,:)./PivVal); % One step of GE.
    
    % Keep track of progress.
    zrows = zrows + 1;                  % One more row is zero.
    PivotValue(zrows) = PivVal;         % Store pivot value.
    PivotElement(zrows,:)=[row col];    % Store pivot location.
    
    %Next pivot.
    [ infnorm , ind ]=max( abs ( A(:) ) );  % slightly faster.
    [ row , col ] = myind2sub( size(A) , ind );
end

if infnorm <= tol, ifail = 0; end  % We didn't fail.
if (zrows >= width / factor), ifail = 1; end  % We did fail.

end



function [row, col] = myind2sub(siz, ndx)
% My version of ind2sub. In2sub is slow because it has a varargout. Since
% this is at the very inner part of the constructor and slowing things down
% we will make our own.
% This version is about 1000 times faster than MATLAB ind2sub.

vi = rem( ndx - 1, siz(1) ) + 1 ;
col = ( ndx - vi ) / siz(1) + 1;
row = ( vi - 1 ) + 1;

end


function vals = evaluate( op, xx, yy, flag )
if ( flag )
    vals = zeros( size( yy, 1), size( xx, 2) );
    for jj = 1 : size( yy, 1)
        for kk = 1 : size( xx , 2 )
            vals(jj, kk) = op( xx( 1, kk) , yy( jj, 1 ) );
        end
    end
else
    vals = op( xx, yy );              % Matrix of values at cheb2 pts.
end
end


function op = str2op( op )
% OP = STR2OP(OP), finds the dependent variables in a string and returns
% an op handle than can be evaluated.
depvar = symvar( op );
if ( numel(depvar) > 2 )
    error('FUN2:fun2:depvars',...
        'Too many dependent variables in string input.');
end
if ( numel(depvar) == 1 )
    warning('CHEBFUN2:chebfun2:depvars',...
        'Not a bivariate function handle.');  % exclude the case @(x) for now..
    
    % Not sure if this should be a warning or not.
    
end
op = eval(['@(' depvar{1} ',' depvar{2} ')' op]);
end