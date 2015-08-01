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

if ( nargin == 0 )          % SPHEREFUN( )
    return
end

if ( nargin == 0 )          % SPHEREFUN( )
    return
end

if ( isa(op, 'diskfun') )  % SPHEREFUN( SPHEREFUN )
    g = op;
    return
end


if ( nargin < 4 || isempty(dom) )
   dom = [-pi pi 0 1];
end

if ( nargin < 3 || isempty(coords) ) %for now assume polar if not specified
    coords = 1
end

%coords=1 -> polar
%coords=0 -> cartesian


if ( isa(op, 'diskfun') )  % SPHEREFUN( SPHEREFUN )
    g = op;
    return
end

% TODO: 
% 1. Need to allow for different domains.
% 2. Add support for preferences
% 3. Add non-adaptive construction
% 4. Add fixed-rank.
% 5. Add tensor-product.

maxRank = 4000; 
maxSample = 4000;
pseudoLevel = eps;

% If f is defined in terms of x,y,z; then convert it to
% (longitude,latitude).

%if numel(varargin) > 0
               % if ischar(varargin{1})
                  %  if strcmpi(varargin{1},'p')
                    %    polarcoords=1;
                    %    varargin{1} = [];
                  %  end
               % end
           % end


h = redefine_function_handle( op, coords );

% PHASE ONE  
% Sample at square grids, determine the numerical rank of the
% function.


n = 4; %sample is 2n+1 points across -1 to 1, n+1 points from 0 to 1; 2n points from -pi to pi in theta
happy_rank = 0;     % Happy with phase one? 
failure = false;
while ( ~happy_rank && ~failure )
   n=2*n;
    % Sample function on a tensor product grid.
    % TODO: Add a more sophisticated evaluate function that does
    % vectorization like chebfun2.
    F = evaluate(h, n, n);
    
    tol = GetTol(F, pi/(n-1), pi/n, dom, pseudoLevel); 
     %tol=3e-14;
    [ pivotIndices, pivotMatrices, happy_rank, removePoles ] = PhaseOne( F, tol );
    if ( n >= maxRank  )
        warning('SPHEREFUN:CONSTRUCTOR:MAXRANK', ... 
                                'Unresolved with maximum rank.');
        failure = true;
    end
end

% PHASE TWO 
% Find the appropriate discretizations in the columns and rows. 
[cols, pivots, rows, pivotLocations, idxPlus, idxMinus] = PhaseTwo( h, pivotIndices, pivotMatrices, n, dom, tol, maxSample, removePoles );

g.cols = chebfun( cols, dom(3:4)-[1 0]);
g.rows = chebfun( rows, dom(1:2), 'trig');
g.pivotValues = pivots;
g.pivotIndices = pivotIndices;
g.pivotLocations = pivotLocations;
g.idxPlus = idxPlus;
g.idxMinus = idxMinus;
g.domain = dom;

end

function [pivotIndices, pivotMatrices, happy, removePole] = PhaseOne( F, tol )

% Phase 1: Go find rank, plus pivot locations, ignore cols and rows.
alpha = 2;%spherefun.alpha; % get growth rate factor.
[m, n] = size( F );
pivotIndices = []; pivotMatrices = [];
vscl = norm( F( : ), inf);
rank_count = 0;    % keep track of the rank of the approximation.

%
% Deal with the poles by removing them from F.
%
%pole1 = mean(F(1,:));     % Take the value at the poles to be 
pole2 = mean(F(m,:));     % the mean of all the samples at the poles.

% If the the values at both poles are not zero then we need to add zero
% them out before removing these entries from F.
removePole = false;
if abs(pole2) > vscl*tol 
    % Determine the column with maximum inf-norm
    [ignored, poleCol] = max(max(abs(F(:,1:n/2)),[],1));
    % Zero out the pole using the poleCol.
    F = F - 0.25*F(:, [poleCol poleCol+n/2])*(ones(2)*ones(2,n)); 
    % Do we need to keep track of the pivot locations?
    removePole = true;
end

% Remove the rows corresponding to the poles in F before determining the
% rank.  We do this because then F is an even BMC matrix, which is what the
% code below requires.
F = F( 1:m-1, : );

% Update the number of rows F now contains.
m = m-1; 


while ( norm( F( : ), inf ) > tol )
   
    % Find pivot:
    % Calculate the maximum 1st singular value of all the special 2x2
    % submatrices.
    B = F(:,1:n/2);    %% (1,1) Block of F.
    C = F(:,n/2+1:n);  % (1,2) block of F.
    Fp = B + C;
    Fm = B - C;
    S1 = max( abs(Fp), abs(Fm) );
    [ignored, idx] = max( S1(:) );
    [j, k] = myind2sub( size( S1 ), idx );
    
    % Calculate the eigenvalues of the pivot matrix:
    % This is what we really care about in the algorithm.
    ev = [ Fp(j,k) , Fm(j,k) ];
    % Singular-values sorted by magnitude
    sv = [max(abs(ev)) min(abs(ev))];
        
    pivotIndices = [ pivotIndices ; j k];

    if ( sv(1) <= alpha*sv(2) )  % Theoretically, should be s1 <= 2*s2.
        % Calculate inverse of pivot matrix:
        plusBlk = 1/(2*ev(1))*(Fp(:,k)*Fp(j,:));
        minusBlk = 1/(2*ev(2))*(Fm(:,k)*Fm(j,:));
        F =  F - [plusBlk plusBlk] - [minusBlk -minusBlk];                  
        rank_count = rank_count + 2; 
        pivotMatrices = [pivotMatrices ; ev ];
    else
        % Calculate pseudoinverse of pivot matrix, there is
        % no full rank pivot matrix:
        if abs(ev(1)) > abs(ev(2))
            plusBlk = 1/(2*ev(1))*(Fp(:,k)*Fp(j,:));
            F =  F - [plusBlk plusBlk];
            ev(2) = 0;
        else
            minusBlk = 1/(2*ev(2))*(Fm(:,k)*Fm(j,:));
            F =  F - [minusBlk -minusBlk];                  
            ev(1) = 0;
        end
        pivotMatrices = [pivotMatrices ; ev ];
        rank_count = rank_count + 1; 
    end
     
end



% Put the poles at the begining the pivot locations array and also include
% the pivot matrix.
if removePole
    pivotIndices = [ 1 poleCol; pivotIndices ];
    ev = [2 0];
    pivotMatrices = [ev ; pivotMatrices];
end

% If the rank of the matrix is less than 1/4 its size. We are happy:
if ( rank_count < min(size(F))/4 )
    happy = 1;
else
    happy = 0;
end

end

function [cols, pivots, rows, pivotLocations, idxPlus, idxMinus] = PhaseTwo( h, pivotIndices, pivotMatrices, n, dom, tol, maxSample, removePoles)

alpha = 2;%spherefun.alpha; % get growth rate factor.
happy_columns = 0;   % Not happy, until proven otherwise.
happy_rows = 0;
m = n;

[x, y] = getPoints( m, n);

rk = size( pivotIndices, 1);
id = pivotIndices'; id = id(:);
id_rows = id(1:2:end); id_cols = id(2:2:end); 
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
    
    [x, y] = getPoints( m, n);
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
    cols = [ [colsPlus colsMinus] ; [flipud(colsPlus(1:m,:)) -flipud(colsMinus(1:m,:))] ];
    rows = [ rowsPlus ;  rowsMinus ].';
    
    
    col_coeffs = chebtech2.vals2coeffs( cols ); 
    %Length of tail to test.
    testLength = min(m, max(3, round((m-1)/8)));
    tail = col_coeffs(end-testLength+1:end,:);
    

    if ( all( abs( tail ) <= 1e2*tol*norm(cols,inf) ) )
       happy_columns = 1;
    end
    
    %colData.vscale = dom(3:4);
    %tech = pref.tech(); 
   % colChebtech = chebtech1.make(sum(cols,2), [-1 1]);
    %resolvedCols = happinessCheck(colChebtech,[],sum(cols,2));
   % happy_columns=resolvedCols;
   
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
        warning('DISKFUN:constructor:notResolved', ...
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





    % Wrap f so it can be evaluated in spherical coordinates
    if ~(coords==1)
    f = @(th, r) diskfun.pol2cartf(f,th, r);
    end
%     % Double g up.
%     f = @(lam, th) sph2torus(f,lam,th);
end

%end

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
tol = grid.^(2/3) * max( abs(dom(:) ) ) * max( Jac_norm, vscale) * pseudoLevel;

end

% function f = redefine_function_handle_pole( f, poleColPivot )
% % Set f to f - f(poleColPivot,theta) where poleColPivot is the value of
% % lambda (the column) used to zero out the poles of f.
% 
% f = @(lam, th) f(lam,th) - f(poleColPivot,th);
% 
%end

