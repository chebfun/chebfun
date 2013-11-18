function g = constructor(g, op, domain, varargin)

if ( nargin < 3 || isempty(domain) )
    domain = [-1 1 -1 1];
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
approxRank = 9;

% Initialise:
vectorize = 0;
vscale = 1;
hscale = 1;

notHappy = 1;  % If unhappy, selected pivots were not good enough.
while ( notHappy )
    [xx, yy] = chebpts2(approxRank, approxRank, domain);
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
    [PivotValue, PivotPos, rowValues, colValues, ifail] = CompleteACA(vals, tol);
    strike = 1;
    while ( ifail && approxRank <= maxRank && strike < 3 )
        [xx, yy] = chebpts2(approxRank, approxRank, domain);
        vals = evaluate(op, xx, yy, vectorize);                        % Resample on denser grid.
        [PivotValue, PivotPos, rowValues, colValues, ifail] = CompleteACA(vals, tol);
        if ( abs(PivotValue(1))<1e4*vscale*tol )
            % If the function is 0+noise then stop after three strikes.
            strike = strike + 1;
        end
        approxRank = 2^(floor(log2(approxRank)) + 1) + 1;                % Double the sampling
    end
    
    if ( approxRank >= maxRank )
        error('FUN2:CTOR', 'Not a low-rank function.');
    end
    
    %% See if the slices are resolved:
    colChebtech = chebtech2(colValues);
    ResolvedCols = happinessCheck(colChebtech);
    rowChebtech = chebtech2(rowValues.');
    ResolvedRows = happinessCheck(rowChebtech);   
    
    ResolvedSlices = ResolvedRows & ResolvedCols;
    if ( length(PivotValue) == 1 && PivotValue == 0 )
        PivPos = [0, 0]; 
        ResolvedSlices = 1;
    else
        PivPos = [xx(1, PivotPos(:, 2)); yy(PivotPos(:, 1), 1).'].'; 
        PP = PivotPos;
    end
    
    n = approxRank; 
    m = approxRank;
    % If unresolved then perform ACA on selected slices.
    while ( ~ResolvedSlices )
        if ( ~ResolvedCols )
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
        if ( ~ResolvedRows )
            m = 2^(floor(log2(m))+1) + 1;
            [xx, yy] = meshgrid(chebpts(m, domain(1:2)), PivPos(:, 2));
            rowValues = evaluate(op, xx, yy, vectorize);
            oddm = 1:2:m; PP(:, 2) = oddm(PP(:, 2)); % find location of pivots on new grid.
        else
            [xx, yy] = meshgrid(chebpts(m, domain(1:2)), PivPos(:, 2));
            rowValues = evaluate(op, xx, yy, vectorize);
        end
        
        nn = numel(PivotValue);
        
%         % ACA on selected Pivots.
%         for kk = 1:nn-1
%             colValues(:, kk+1:end) = colValues(:, kk+1:end) - colValues(:, kk)*(rowValues(kk, PP(kk+1:nn, 2))./PivotValue(kk));
%             rowValues(kk+1:end, :) = rowValues(kk+1:end, :) - colValues(PP(kk+1:nn, 1), kk)*(rowValues(kk, :)./PivotValue(kk));
%         end
        
        for kk=1:nn-1
            selx = PP(kk+1:nn,1); sely = PP(kk+1:nn,2);
            colValues(:,kk+1:end) = colValues(:,kk+1:end) - colValues(:,kk)*(rowValues(kk,sely)./PivotValue(kk));
            rowValues(kk+1:end,:) = rowValues(kk+1:end,:) - colValues(selx,kk)*(rowValues(kk,:)./PivotValue(kk));
        end

        % Are the columns and rows resolved now?
        if ( ~ResolvedCols )
            colChebtech = chebtech2(colValues);
            ResolvedCols = happinessCheck(colChebtech);    
            ResolvedCols = true;
        end
        if ( ~ResolvedRows )
            rowChebtech = chebtech2(rowValues.');
            ResolvedRows = happinessCheck(rowChebtech);    
            ResolvedRows = true;
        end
        ResolvedSlices = ResolvedRows & ResolvedCols;
        if ( max(m, n) >= maxDegree )  % max number of degrees allows.
            error('FUN2:CTOR', 'Unresolved with maximum Chebfun length: %u.', maxDegree);
        end
        
    end
    
    % Pivots locations are probably ok..
    if ( ResolvedSlices )
        notHappy = 0; 
    end 
    
    % check the first three pivots are at least okay. Though a rank
    % 1 function doesn't need this check (almost any nonzero pivot location
    % will do).
    if ( length( PivotValue ) > 1 )
        for jj = 1:min(3, length(PivotValue))
            if max(max(abs(colValues(:, jj))), max(abs(rowValues(jj, :)))) - abs(PivotValue(jj)) > hscale*vscale*1e-2;
                notHappy = 1; break
            end
        end
    end
    
end

% For some reason, on some computers simplify is giving back a
% scalar zero.  In which case the function is numerically zero.
% Artifically set the columns and rows to zero.
if ( norm(colValues) == 0 || norm(rowValues) == 0)
    colValues = 0; 
    rowValues = 0; 
    PivotValue = 0;
    PivPos = [0, 0]; 
    ResolvedSlices = 1;
end

% Construct a CHEBFUN2:
g.pivotValues = PivotValue;
g.cols = colChebtech;
g.rows = rowChebtech;

end

function [PivotValue, PivotElement, Rows, Cols, ifail] = CompleteACA(A, tol)
% Adaptive Cross Approximation with complete pivoting. This command is
% completely analogous to Gaussian elimination with complete pivoting. We
% adaptively find the rank of the approximant in this command.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

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
    A = A - Cols(:,zrows+1)*(Rows(zrows+1,:)./PivVal);    % Rank one update.
    
    % Keep track of progress.
    zrows = zrows + 1;                  % One more row is zero.
    PivotValue(zrows) = PivVal;         % Store pivot value.
    PivotElement(zrows,:)=[row col];    % Store pivot location.
    
    %Next pivot.
    %     [row,col,infnorm]=rook_pivot(A);
    %     [ infnorm , ind ]=max( abs ( reshape(A,numel(A),1) ) );
    [ infnorm , ind ]=max( abs ( A(:) ) );  % slightly faster.
    [ row , col ] = myind2sub( size(A) , ind );
end

if infnorm <= 10*tol*scl, ifail = 0; end  % We didn't fail.
if (zrows >= width / factor), ifail = 1; end  % We did fail.



if ifail == 0, return; end

% Toy Plateau detection.
if length(PivotValue )  > 1
    if all(abs(PivotValue(end-2:end))<100*tol*scl)
        % perhaps we did plateau, but there was a bit of noise.
        if all(diff(abs(PivotValue(end-2:end)))<10*tol*scl)
            % If we think we plateau, we probably didn't fail.
            ifail = 0;
        end
    end
end
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