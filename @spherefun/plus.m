function h = plus(f, g)
%+   Plus for SPHEREFUN objects.
%
% F + G adds F and G. F and G can be scalars or SPHEREFUN objects.

if ( ~isa(f, 'spherefun') ) % ??? + SPHEREFUN
    
    h = plus(g, f);
    
elseif ( isempty(g) ) % SPHEREFUN + []
    
    h = f; 
    
elseif ( isempty(f) ) % [] + SPHEREFUN
    
    h = g; 
    
elseif ( isa( g, 'double' ) )           % SPHEREFUN + DOUBLE
    
    g = compose( 0*f,@plus, g);   % promote double to object class.  
    h = plus(f, g); 
    
elseif ( ~isa(g, 'spherefun') )          % SPHEREFUN + ???
    
    error( 'SPHEREFUN:plus:unknown', ...
        ['Undefined function ''plus'' for input arguments of type %s ' ...
        'and %s.'], class( f ), class( g ));
        
else                                     % SPHEREFUN + SPHEREFUN
    
    % Domain Check:
    if ( ~domainCheck(f, g) )
        error('SPHEREFUN:plus:domain', 'Inconsistent domains.');
    end
    
    % Check for zero SPHEREFUN objects:
    if ( iszero(f) )
        h = g;
    elseif ( iszero(g) )
        h = f;
    elseif ( isequal(f,-g) )
        h = 0*f;
    else
        % Add together two nonzero SPHEREFUN objects:
        % The algorithm is as follows: Split f and g into their plus/minus
        % components.  Do the compression_plus algorithm described in
        % @separableApprox/compression plus on each pair of plus and minus
        % components.
        
        % Check if the pole is non-zero.  If it is then we need to strip
        % out the column and row that deal with this in the plus piece
        % before doing compression plus.  The reason is that what we feed
        % compression plus needs to be zero at the poles if what is to be
        % returned is zero at the poles.  Including one column that is
        % non-zero at the poles can screw everything up.
        [f,fPole] = extractPole(f);
        [g,gPole] = extractPole(g);

        [fp,fm] = partition(f);
        [gp,gm] = partition(g);
        
        hp = plus@separableApprox(fp,gp);
        r = size(hp.cols,2);
        hp.idxPlus = 1:r;
        % Indices or locations of the pivots do not make sense after a
        % compression plus, so we set them to NaN.
        hp.pivotLocations = nan(r,2);  % This should be done at the separableApprox level.
        
        hm = plus@separableApprox(fm,gm);
        r = size(hm.cols,2);
        hm.idxMinus = 1:r;
        hm.pivotLocations = nan(r,2);
                
        if ( ~isempty( fPole ) ) || ( ~isempty( gPole ) )
            % Set tolerance for determining if fPole+gPole=0.
            tol = eps*max(vscale(f),vscale(g)); 
            g = addPoles(fPole,gPole,tol);
            % Handle the rare case that g is zero and hp is not empty
            if ( g.pivotValues == 0 ) && ( ~isempty( hp ) )
                % Set g to empty spherefun
                g = spherefun([]);
            end
            hp.cols = [g.cols hp.cols];
            hp.rows = [g.rows hp.rows];
            hp.pivotValues = [g.pivotValues;hp.pivotValues];
            hp.pivotLocations = [g.pivotLocations;hp.pivotLocations];
            hp.idxPlus = 1:size(hp.cols,2);
            hp.nonZeroPoles = ~isempty(g);
        end
        
        hp = projectOntoEvenBMCI( hp );
        hm = projectOntoOddBMCI( hm );
        
        % Put pieces back together.
        h = combine(hp,hm);

    end 
    
end

end

function f = addPoles(f,g,tol)

if isempty(g)
    return;
elseif isempty(f)
    f = g;
    return;
end

fmean = mean(f.rows);
gmean = mean(g.rows);

cols = (fmean/f.pivotValues)*f.cols + (gmean/g.pivotValues)*g.cols;

% If cols is numerically zero then return an empty result
if norm(cols) <= tol
    fmean = 0;
    gmean = 0;
    pivot = 0;
    nonZeroPoles = 0;
else
    pivot = 1;
    nonZeroPoles = 1;
end

f.cols = (fmean/f.pivotValues)*f.cols + (gmean/g.pivotValues)*g.cols;
f.rows = chebfun('1',f.domain(1:2),'trig');
f.pivotValues = pivot;
% No idea what indices or locations should be after plus
f.pivotLocations = [nan nan];
f.nonZeroPoles = nonZeroPoles;

end

function f = projectOntoEvenBMCI( f )
% Project a spherefun to have even BMC-I symmetry, i.e., a spherefun that
% is pi-periodic in lambda and even in theta. The projection is orthogonal,
% i.e., the correction matrix to fix up the structure has the smallest
% possible Frobenius norm.

% Nothing to project
if isempty( f )
    return;
end

% Operate on the column coefficients first to project them onto even
% functions.
X = f.cols.funs{1}.onefun.coeffs;

% Get size: 
[m, n] = size(X); 

isevenM = false;
if mod(m,2) == 0
    X(1,:) = 0.5*X(1,:);
    X = [X;X(1,:)];
    m = m+1;
    isevenM = true;
end

% Only project the nonzero Fourier modes:
waveNumbers = -(m-1)/2:(m-1)/2;

evenModes = 1:n;
A = [];
if f.nonZeroPoles
    zeroMode = 1;
    % Need to handle the zero mode in lambda separately
    % Enforce the expansion is even in theta
    I = eye(m); A = I - fliplr(I); A = A(1:(m-1)/2,:);
    % Solution to underdetermined system A*(X + Y) = 0 with smallest Frobenius
    % norm: 
    C = A\(A*X(:,zeroMode));

    % Update coeff matrix: 
    X(:,zeroMode) = X(:,zeroMode) - C; 

    % The result of the code now needs to operate on the remaining even,
    % non-zero modes.
    evenModes = 2:n;
end

% Second do the even, non-zero modes in lambda
% Enforce these are zero at the poles and that the expansion is even in 
% theta
A = [[ones(1,m); (-1).^waveNumbers];A];

% Solution to underdetermined system A*(X + Y) = 0 with smallest Frobenius
% norm: 
C = A\(A*X(:,evenModes));

% Update coeff matrix: 
X(:,evenModes) = X(:,evenModes) - C; 

% If m is even we need to remove the mode that was appended 
if ( isevenM )
    X(1,:) = (X(1,:)+X(end,:));
    X = X(1:m-1,:);
end

ctechs = real(trigtech({'',X}));
f.cols.funs{1}.onefun = ctechs;

% Now operate on the rows. The coefficients for the rows of an even BMCI
% function should only contain even wave numbers. The projection is to
% simply zero out the odd wave numbers.
X = f.rows.funs{1}.onefun.coeffs;
n = size(X,1); 
zeroMode = floor(n/2)+1;
oddModes = [fliplr(zeroMode-1:-2:1) zeroMode+1:2:n];
X(oddModes,:) = 0;
rtechs = real(trigtech({'',X}));
f.rows.funs{1}.onefun = rtechs;

% Weird feval behavior in chebfun requires this
f.cols.pointValues = feval(ctechs,[-1;1]);
f.rows.pointValues = feval(rtechs,[-1;1]); 

end

function f = projectOntoOddBMCI( f )
% Project a spherefun to have odd BMC-I symmetry, i.e., a spherefun that is
% pi-anti-periodic in lambda and even in theta. The projection is
% orthogonal, i.e., the correction matrix to fix up the structure has the
% smallest possible Frobenius norm.

% Nothing to project
if isempty( f )
    return;
end

% Operate on the column coefficients first to project them onto odd
% functions.
X = f.cols.funs{1}.onefun.coeffs;

% Get size: 
m = size(X,1);

isevenM = false;
if mod(m,2) == 0
    X(1,:) = 0.5*X(1,:);
    X = [X;X(1,:)];
    m = m+1;
    isevenM = true;
end

I = eye(m); A = I + fliplr(I); 
A = A(1:(m-1)/2+1,:); A((m-1)/2+1,(m-1)/2+1) = 1;

% Solution to underdetermined system A*(X + Y) = 0 with smallest Frobenius
% norm: 
C = A\(A*X);
% Update coeff matrix: 
X = X - C; 

% If m is even we need to remove the mode that was appended 
if ( isevenM )
    X(1,:) = (X(1,:)+X(end,:));
    X = X(1:m-1,:);
end

ctechs = real(trigtech({'',X}));
f.cols.funs{1}.onefun = ctechs;

% Now operate on the rows. The coefficients for the rows of an odd BMCI
% function should only contain odd wave numbers. The projection is to
% simply zero out the even wave numbers.
X = f.rows.funs{1}.onefun.coeffs;
n = size(X,1); 
zeroMode = floor(n/2)+1;
evenModes = [fliplr(zeroMode-2:-2:1) zeroMode:2:n];
X(evenModes,:) = 0;

rtechs = real(trigtech({'',X}));
f.rows.funs{1}.onefun = rtechs;

% Weird feval behavior in chebfun requires this
f.cols.pointValues = feval(ctechs,[-1;1]);
f.rows.pointValues = feval(rtechs,[-1;1]); 

end

