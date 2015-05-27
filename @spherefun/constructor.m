function g = constructor( g, op, varargin )
%CONSTRUCTOR   The main SPHEREFUN constructor.
%
% The algorithm for constructing a SPHEREFUN comes in three phases:
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
%
% PHASE 3: The decomposition is reduced to CDR form such that the BMC
% structure is maintained.

if ( nargin == 0 )          % SPHEREFUN( )
    return
end

if ( isa(op, 'spherefun') )  % SPHEREFUN( SPHEREFUN )
    g = op;
    return
end

% TODO: 
% 1. Need to allow for different domains.
% 2. Add support for preferences
% 3. Add non-adaptive construction
% 4. Add fixed-rank.
% 5. Add tensor-product.

tol = 50*eps;       % Tolerance
maxRank = 4000; 
maxSample = 4000;

% If f is defined in terms of x,y,z; then convert: 
h = redefine_function_handle( op );

% PHASE ONE  
% Sample at square grids, determine the numerical rank of the
% function.
n = 4;
happy_rank = 0;     % Happy with phase one? 
failure = false;
while ( ~happy_rank && ~failure )
    n = 2*n;
    F = sample(h, n, n);
    [ pivotLocations, pivotMatrices, happy_rank, removePoles ] = PhaseOne( F, tol );
    if ( n >= maxRank  )
        warning('SPHEREFUN:CONSTRUCTOR:MAXRANK', ... 
                                'Unresolved with maximum rank.');
        failure = true;
    end
end

% PHASE TWO 
% Find the appropriate discretizations in the columns and rows. 
[cols, blockDiag, rows] = PhaseTwo( h, pivotLocations, pivotMatrices, n, tol, maxSample, removePoles );

% PHASE THREE
% Put in CDR form
[cols, rows, pivots, idxPlus, idxMinus] = PhaseThree( cols, blockDiag, rows.' );

% Make a spherefun, we are done. 
g.cols = trigtech( cols );
g.rows = trigtech( rows );
g.blockDiag = blockDiag;
g.pivotValues = pivots;
g.pivotLocations = pivotLocations;
g.idxPlus = idxPlus;
g.idxMinus = idxMinus;

end

function [pivotLocations, pivotMatrices, happy, removePole] = PhaseOne( F, tol )

% Phase 1: Go find rank, plus pivot locations, ignore cols and rows.
alpha = spherefun.alpha; % get growth rate factor.
[m, n] = size( F );
pivotLocations = []; pivotMatrices = [];
vscl = norm( F( : ), inf);
rank_count = 0;    % keep track of the rank of the approximation.

%
% Deal with the poles by removing them from F.
%
pole1 = mean(F(1,:));       % Take the value at the poles to be 
pole2 = mean(F(m/2+1,:));   % the mean of all the samples at the poles.

% If the the values at both poles are not zero then we need to add zero
% them out before removing these entries from F.
removePole = false;
if abs(pole1) > vscl*10*eps || abs(pole2) > vscl*10*eps
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
F = F([2:m/2 m/2+2:m],:);

% Update the number of rows F now contains.
m = m-2;

while ( norm( F( : ), inf ) > tol*vscl )
    % Find pivot:
    % Calculate the maximum 1st singular value of all the special 2x2
    % submatrices.
    S1 = max(  abs(F - flipud(F)), abs(F + flipud(F)) );
    S1 = S1(1:m/2, 1:n/2);
    [ignored, idx] = max( S1(:) );
    [j, k] = myind2sub( size( S1 ), idx );
    
    % Max singular value submatrix is:
    M = [ F(j,k) F(j,k+n/2) ; F(j,k+n/2) F(j,k)];
    s = sort( abs([ diff(M(1,:)) ; sum(M(1,:)) ]), 1, 'descend' );  % equivalent to svd( M )
    
    pivotLocations = [ pivotLocations ; j k m-j+1 k+n/2];
    pivotMatrices = [pivotMatrices ; M ];
    
    if ( s(1) <= alpha*s(2) )  % Theoretically, should be s1 <= 2s2.
        % Calculate inverse of pivot matrix:
        invM = [M(1,1) -M(1,2);-M(1,2) M(1,1)]./(M(1,1)^2 - M(1,2)^2);
        F = F - F(:, [k k+n/2] ) * ( invM * F( [j m-j+1],: ) );
        rank_count = rank_count + 2; 
    else
        % Calculate pseudoinverse of pivot matrix, there is
        % no full rank pivot matrix:
        pinvM = getPseudoInv( M );
        F = F - F(:, [k k+n/2] ) * ( pinvM *  F( [j m-j+1],: ) );
        rank_count = rank_count + 1; 
    end
    
end

% Adjust the pivot locations so that they now correspond to F having
% the poles.
if ~isempty( pivotLocations )
    pivotLocations(:,1) = pivotLocations(:,1) + 1;
    pivotLocations(:,3) = pivotLocations(:,3) + 2;
end

% Put the poles at the begining the pivot locations array and also include
% the pivot matrix.
if removePole
    pivotLocations = [ 1 poleCol m/2+1 poleCol+n/2 ; pivotLocations ];
%     M = diag([1/pole1 1/pole2]);
    M = ones(2);
    pivotMatrices = [M ; pivotMatrices];
end

% If the rank of the matrix is less than 1/4 its size. We are happy:
if ( rank_count < min(size(F))/8 )
    happy = 1;
else
    happy = 0;
end

end



function [cols, blockDiag, rows] = PhaseTwo( h, pivotLocations, pivotMatrices, n, tol, maxSample, removePoles)

alpha = spherefun.alpha; % get growth rate factor.
happy_columns = 0;   % Not happy, until proven otherwise.
happy_rows = 0;
m = n;
[x, y] = getPoints( m, n );

rk = size( pivotLocations, 1);
blockDiag = spalloc( 2*rk, 2*rk, 6*rk-2 );
id = pivotLocations'; id = id(:);
id_cols = id(2:2:end); id_rows = id(1:2:end);
col_pivots = x(id_cols);
row_pivots = y(id_rows);

% Phase 2: Calculate decomposition on sphere.
failure = false;
while ( ~(happy_columns || happy_rows) && ~failure)
    
    [x, y] = getPoints( m, n );
    [xx, yy] = meshgrid( col_pivots, y);
    Cols_new = h( xx ,yy );
    [xx, yy] = meshgrid( x, row_pivots);
    Rows_new = h( xx, yy );
    cols = zeros( 2*m, 2*rk );
    rows = zeros( 2*rk, 2*n );
    
    % Need to remove pole, which means we use the column with the largest
    % max norm (repeated) with rows of all ones in the elimination
    % algorithm.
    if removePoles
        Rows_new(1:2,:) = 1 + 0*Rows_new(1:2,:);
%         poleCol = PivotLocations( 1, 2 );
%         poleColPivot = x(poleCol);
%         h = redefine_function_handle_pole(h,poleColPivot);
%         M = PivotMatrices(1:2,:);
%         cols(:,1:2) = Cols_new(:,1:2);
%         % rows will be all ones.
%         rows(1:2,:) = ones( 2, 2*n ); 
%         % Remove the pivots corresponding to the poles.
    end
    
    
    for ii = 1:rk
        
        M = pivotMatrices( 2*ii-1:2*ii, :);
        s = sort( abs([ diff(M(1,:)) ; sum(M(1,:)) ]), 1, 'descend' );  % equivalent to svd( M )
        
        cols(:,2*ii-1:2*ii) = Cols_new(:,2*ii-1:2*ii); %F(:, [k k+n] );
        rows(2*ii-1:2*ii,:) = Rows_new(2*ii-1:2*ii,:); %F([j 2*m-j+1],: );
        if ( s(1) <= alpha*s(2) )
            % Calculate inverse of pivot matrix:
            invM = [M(1,1) -M(1,2);-M(1,2) M(1,1)]./(M(1,1)^2 - M(1,2)^2);
            row_correction = Cols_new(id_rows,2*ii-1:2*ii) * ( invM * Rows_new(2*ii-1:2*ii,:) );
            Cols_new = Cols_new - Cols_new(:,2*ii-1:2*ii) * ( invM * Rows_new(2*ii-1:2*ii,id_cols) );
            Rows_new = Rows_new - row_correction;
            blockDiag(2*ii-1:2*ii, 2*ii-1:2*ii) = invM;
        else
            % Calculate pseudoinverse of the pivot matrix, there is
            % no full rank pivot matrix:
            pinvM = getPseudoInv( M );
            row_correction = Cols_new(id_rows,2*ii-1:2*ii) * ( pinvM * Rows_new(2*ii-1:2*ii,:) );
            Cols_new = Cols_new - Cols_new(:,2*ii-1:2*ii) * ( pinvM * Rows_new(2*ii-1:2*ii,id_cols) );
            Rows_new = Rows_new - row_correction;
            blockDiag(2*ii-1:2*ii, 2*ii-1:2*ii) = pinvM;
        end
        
    end
    
    % Happiness check for columns:
    col_coeffs = trigtech.vals2coeffs( cols ); 
    if ( all( abs( col_coeffs(end:-1:end-8,:) ) <= 1e2*tol*norm(cols,inf) ) )
        happy_columns = 1;
    end
    
    % Happiness check for rows:
    row_coeffs = trigtech.vals2coeffs( rows.' ); 
    if ( all(abs( row_coeffs(end:-1:end-8,:) ) <= 1e2*tol*norm(rows,inf)) )
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

end

% Switch from BMC to CDR format
function [cols, rows, pivots, idxPlus, idxMinus ] = PhaseThree( cols, blockDiag, rows )

tol = 50*eps;

C = cols;
R = rows;
d = zeros(size(blockDiag,1),1);

for j = 1 : size(blockDiag, 1)/2
    
    % Look at 2x2 blocks at one time. 2x2 blocks separate. 
    ii = 2*j-1:2*j;
    
    % Follows equation (4) of Grady's memo.  Basic idea is to compute
    % eigenvalue decomposition of each 2-by-2 matrix on the block diagonal
    % M = U*S*U' 
    %   = 1/sqrt(2)[1 1;1 -1]*[M(1,1)+M(1,2) 0;0 M(1,1)-M(1,2)]*1/sqrt(2)*[1 1;1 -1]
    % Then Compute C*U and U'*R.  This will preserve the BMC structure
    
    % TODO: Perhaps we should switch to the SVD of M so all pivots are
    % positive?
    
    C(:,ii) = [C(:,ii(1))+C(:,ii(2)) C(:,ii(1))-C(:,ii(2))]/sqrt(2);
    R(:,ii) = [R(:,ii(1))+R(:,ii(2)) R(:,ii(1))-R(:,ii(2))]/sqrt(2);    
    d(ii) = full([blockDiag(ii(1),ii(1))+blockDiag(ii(1),ii(2));blockDiag(ii(1),ii(1))-blockDiag(ii(1),ii(2))]);
end

% Extract the plus/minus terms from the columns and rows
Cplus = C(:,1:2:end);  
Cminus = C(:,2:2:end);
Rplus = R(:,1:2:end);
Rminus = R(:,2:2:end);
% Split the pivots as well.
dplus = d(1:2:end);
dminus = d(2:2:end);

% Compress the results by removing the terms corresponding to a zero
% diagonal.
idx = find( abs(dplus) > tol );
dplus = dplus(idx);
Cplus = Cplus(:,idx);
Rplus = Rplus(:,idx);

idx = find( abs(dminus) > tol );
dminus = dminus(idx);
Cminus = Cminus(:,idx);
Rminus = Rminus(:,idx);

% Put the plus/minus terms back together and sort them according to the
% the pivot size (which is 1/d in this setting).  We will keep track
% of these terms in case they are useful later on.
cols = [Cplus Cminus];
rows = [Rplus Rminus];
pivots = 1./[dplus;dminus];
idx = 1:size(pivots,1);

[ignore, perm] = sort( abs( pivots ), 1, 'descend');
pivots = pivots(perm);
cols = cols(:,perm);
rows = rows(:,perm);

% Get the indices of the plus and minus terms so they can be recovered
% later.
idx = idx(perm);
plusFlag = idx <= numel(dplus);
idxPlus = find( plusFlag );
idxMinus = find( ~plusFlag );

% % Now make a new spherefun:
% fp = f;
% fp.cols = trigtech( Cplus );
% fp.rows = trigtech( Rplus );
% fp.BlockDiag = spdiags(dplus,0,size(dplus,1),size(dplus,1));
% 
% fm = f;
% fm.cols = trigtech( Cminus );
% fm.rows = trigtech( Rminus );
% fm.BlockDiag = spdiags(dminus,0,size(dminus,1),size(dminus,1));

end


function [x, y] = getPoints( m, n )

% x = linspace(-pi, pi, 2*n+1)';  x( end ) = [ ];
% GBW: You can't just remove the pole and keep everything else equally
% spaced between -pi/2 and 3*pi/2.  The issue is that you can't keep
% both point -pi/ and 3*pi/2.
% y = linspace(-pi/2, 3*pi/2, 2*m+1)'; y( m+1 ) = [ ];

% x = trigpts(2*n,[-pi pi]);
% % GBW: I believe we have to sample at equally spaced points shifted by h/2
% % to not sample the poles and keep an even total number of points.
% y = trigpts(2*m,[-pi/2 3*pi/2]);
% y = y+0.5*pi/m; % Shift y by h/2 to avoid the poles

x = trigpts(2*n,[-pi pi]);
% Sample at an even number of points so that the poles are included.
% TODO: Change domain.
y = trigpts(2*m,[-pi/2 3*pi/2]);

end

function F = sample( h, m, n )
[x, y] = getPoints( m, n );
[L2, T2] = meshgrid(x, y);
F = h(L2, T2);
end

function pinvM = getPseudoInv( M )
lam1 = M(1,1)+M(1,2);  % Eigenvalues of M (which is symmetric)
lam2 = M(1,1)-M(1,2);
if abs(lam1) > abs(lam2)
    pinvM = ones(2)/(2*lam1);
else
    pinvM = [[1 -1];[-1 1]]/(2*lam2);
end
end

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
    f = @(lam, th) sphf2cartf(f,lam,th);
    % Double g up.
    f = @(lam, th) sph2torus(f,lam,th);
end

end

function f = redefine_function_handle_pole( f, poleColPivot )
% Set f to f - f(poleColPivot,theta) where poleColPivot is the value of
% lambda (the column) used to zero out the poles of f.

f = @(lam, th) f(lam,th) - f(poleColPivot,th);

end


function fdf = sph2torus(f,lam,th)

fdf = real(f(lam,th));

id = th-pi/2 > 100*eps;

if ~isempty(id) && any(id(:))
    fdf(id) = f(lam(id)-pi,pi-th(id));
end

end

function fdf = sphf2cartf(f,lam,th)

x = cos(lam).*cos(th);
y = sin(lam).*cos(th);
z = sin(th);

fdf = f(x,y,z);

end