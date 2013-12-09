function g = constructor(g, op, domain, varargin)
% The main Chebfun2 constructor. 

if ( nargin < 3 || isempty(domain) )
    domain = [-1 1 -1 1];
end

if ( isa(op, 'double') )
   g = constructor(g, @(x,y) op + 0*x, domain, varargin); 
   return
end
    
% Check that the operator then make it complex.
if ( nargin(op) == 1 )
    op = @(x, y) op(x+1i*y);
end

% [TODO]: Get preferences:
% if ( nargin > 4 )
%     pref = chebfun2.pref(varargin{:});
% else
%     pref = chebfun2.pref;
% end
tol = eps;
maxRank = 500;
maxDegree = 2^16;
r = 9;
isHappy = 0;

% Initialise:
vectorize = 0;
vscale = 1;
hscale = 1;

notHappy = 1;  % If unhappy, selected pivots were not good enough.
while ( ~isHappy )
    [xx, yy] = chebpts2(r, r, domain);
    vals = evaluate(op, xx, yy, vectorize);             % Matrix of values at cheb2 pts.
    
    vscale = max(abs(vals(:)));
    % scale has been overloaded.
    if ( ~isempty(vscale) )
        vscale = max(vscale, 1);
    end
    if ( isinf(vscale) )
        error('FUN2:CTOR', 'Function returned INF when evaluated');
    end
    if ( any(isnan(vals(:)) ) )
        error('FUN2:CTOR', 'Function returned NaN when evaluated');
    end
    
    %% FIND NUMERICAL RANK:
    [pivotValue, pivotPosition, rowValues, colValues, iFail] = CompleteACA(vals, tol);
    % Use chebtech's happiness check
    strike = 1;
    while ( iFail && r <= maxRank && strike < 3 && r < 65)
        r = 2^(floor(log2(r)) + 1) + 1;                % Double the sampling
        [xx, yy] = chebpts2(r, r, domain);
        vals = evaluate(op, xx, yy, vectorize);                        % Resample on denser grid.
        [pivotValue, pivotPosition, rowValues, colValues, iFail] = CompleteACA(vals, tol);
        if ( abs(pivotValue(1))<1e4*vscale*tol )
            % If the function is 0+noise then stop after three strikes.
            strike = strike + 1;
        end
    end
    
    if ( r >= maxRank )
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

    n = r;  m = r;
    
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
        
        % Are the columns and rows resolved now?
        if ( ~resolvedCols )
            colChebtech = chebtech2(sum(colValues,2), domain(3:4) );
            resolvedCols = happinessCheck(colChebtech);    
        end
        if ( ~resolvedRows )
            rowChebtech = chebtech2(sum(rowValues.',2), domain(1:2) );
            resolvedRows = happinessCheck(rowChebtech);    
        end
        
        isHappy = resolvedRows & resolvedCols;
        if ( max(m, n) >= maxDegree )  % max number of degrees allows.
            error('FUN2:CTOR', 'Unresolved with maximum Chebfun length: %u.', maxDegree);
        end
        
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
g.cols = simplify(chebfun(colValues, domain(3:4) ));
g.rows = simplify(chebfun(rowValues.', domain(1:2) ));
g.domain = domain;

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
% [xx,yy]=cheb2pts(nx,ny,g.map);  % points sampling from.
[ infnorm , ind ]=max( abs ( reshape(A,numel(A),1) ) );
[ row , col ]=myind2sub( size(A) , ind);
% [row,col,infnorm]=rook_pivot(A);
scl = infnorm;

% If the function is the zero function.
if scl == 0
    PivotValue=0;
    Rows = 0; Cols = 0;
    ifail = 0;
end

while ( ( infnorm > 10*tol*scl ) && ( zrows < width / factor) )
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

if infnorm <= 10*tol*scl, ifail = 0; end  % We didn't fail.
if (zrows >= width / factor), ifail = 1; end  % We did fail.

end



function [row, col] = myind2sub(siz,ndx)
% My version of ind2sub. In2sub is slow because it has a varargout. Since
% this is at the very inner part of the constructor and slowing things down
% we will make our own.
% This version is about 1000 times faster than MATLAB ind2sub.

vi = rem(ndx-1,siz(1)) + 1 ;
col = (ndx - vi)/siz(1) + 1;
row = (vi-1) + 1;

end

function vals = evaluate(op,xx,yy,flag)
if flag
    vals = zeros(size(yy,1),size(xx,2));
    for jj = 1:size(yy,1)
        for kk = 1:size(xx,2)
            vals(jj,kk) = op(xx(1,kk),yy(jj,1));
        end
    end
else
    vals = op(xx,yy);              % Matrix of values at cheb2 pts.
end
end