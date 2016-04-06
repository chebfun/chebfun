function f = projectOntoBMCI(f)
% PROJECTONTOBMCI  Projection onto BMC-I symmetry.
%
% f = projectOntoBMCI(f) is the orthogonal projection of f onto BMC-I 
% symmetry, i.e., a function that is
% 1. even in theta for every even wave number in lambda;
% 2. odd in theta for every odd wave number in lambda;
% Additionally, for all but k=0 wavenumber lambda the resulting projection
% enforces the spherefun is zero at the poles. 
%
% The projection is orthogonal, i.e., the correction matrix to fix up the
% structure has the smallest possible Frobenius norm.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Even part
feven = f;
feven.cols = feven.cols(:, feven.idxPlus);
feven.rows = feven.rows(:, feven.idxPlus);
feven = projectOntoEvenBMCI(feven);

% Odd part
fodd = f;
fodd.cols = fodd.cols(:, fodd.idxMinus);
fodd.rows = fodd.rows(:, fodd.idxMinus);
fodd = projectOntoOddBMCI(fodd);

% Put pieces back together.
f.cols(:, f.idxPlus) = feven.cols;
f.rows(:, f.idxPlus) = feven.rows;

f.cols(:, f.idxMinus) = fodd.cols;
f.rows(:, f.idxMinus) = fodd.rows;

end

function f = projectOntoEvenBMCI(f)
% Project a spherefun to have even BMC-I symmetry, i.e., a spherefun that
% is pi-periodic in lambda and even in theta. The projection is orthogonal,
% i.e., the correction matrix to fix up the structure has the smallest
% possible Frobenius norm.

% Nothing to project
if isempty(f)
    return;
end

% Operate on the column coefficients first to project them onto even
% functions.
X = f.cols.funs{1}.onefun.coeffs;

% Get size: 
[m, n] = size(X); 

isevenM = false;
if ( mod(m, 2) == 0 )
    X(1, :) = 0.5*X(1, :);
    X = [X; X(1,:)];
    m = m+1;
    isevenM = true;
end

% Goal is to enforce the expansion is even in theta

evenModes = 1:n;
if ( f.nonZeroPoles )
    zeroMode = 1;
    % Need to handle the zero mode in lambda separately as the expansion
    % only has to be forced to be even in theta.

    % Let 
    % I = eye(m); A = I - fliplr(I); A = A(1:(m-1)/2, :);
    % then we want to find C with smallest two-norm such that
    % A*(X + C) = 0
    % The solution is 
    % C = A'*((A*A')\(A*X))
    % However, because of the form of the matrix equation the solution can
    % be worked out to just be
    C = 0.5*(X(1:m, zeroMode) - X(m:-1:1, zeroMode));

    % Update coeff matrix: 
    X(:, zeroMode) = X(:, zeroMode) - C; 

    % The result of the code now needs to operate on the remaining even,
    % non-zero modes.
    evenModes = 2:n;
end

if ( ~isempty(evenModes) )
    % Second do the even, non-zero modes in lambda

    % First enforce these are even as above.
    C = 0.5*(X(1:m, evenModes) - X(m:-1:1, evenModes));
    X(:, evenModes) = X(:, evenModes) - C;

    % Now enforce these are zero at the poles (evenness is preserved)
    % Letting 
    % A = [[ones(1,m); (-1).^waveNumbers];
    % We want to find the C with smallest two-norm such that 
    % A*(X + C) = 0
    % The solution is 
    % C = A'*((A*A')\(A*X))
    % However, we can again work out the solution in close form because of
    % the special structure of the matrix equation.
    X(1:2:m, evenModes) = bsxfun(@minus,...
        X(1:2:m, evenModes), (2/(m+1))*sum(X(1:2:m, evenModes), 1));
    X(2:2:m, evenModes) = bsxfun(@minus,...
        X(2:2:m, evenModes), (2/(m-1))*sum(X(2:2:m, evenModes), 1));    
end

% If m is even we need to remove the mode that was appended
if ( isevenM )
    X(1, :) = (X(1, :) + X(end, :));
    X = X(1:m-1, :);
end

ctechs = real(trigtech({'', X}));
f.cols.funs{1}.onefun = ctechs;

% Now operate on the rows. The coefficients for the rows of an even BMCI
% function should only contain even wave numbers. The projection is to
% simply zero out the odd wave numbers.
X = f.rows.funs{1}.onefun.coeffs;
n = size(X, 1); 
zeroMode = floor(n/2) + 1;
oddModes = [fliplr(zeroMode-1:-2:1) zeroMode+1:2:n];
X(oddModes, :) = 0;
rtechs = real(trigtech({'', X}));
f.rows.funs{1}.onefun = rtechs;

% Weird feval behavior in chebfun requires this
f.cols.pointValues = feval(ctechs, [-1; 1]);
f.rows.pointValues = feval(rtechs, [-1; 1]); 

end

function f = projectOntoOddBMCI(f)
% Project a spherefun to have odd BMC-I symmetry, i.e., a spherefun that is
% pi-anti-periodic in lambda and even in theta. The projection is
% orthogonal, i.e., the correction matrix to fix up the structure has the
% smallest possible Frobenius norm.

% Nothing to project
if isempty(f)
    return;
end

% Operate on the column coefficients first to project them onto odd
% functions.
X = f.cols.funs{1}.onefun.coeffs;

% Get size: 
m = size(X, 1);

isevenM = false;
if ( mod(m, 2) == 0 )
    X(1, :) = 0.5*X(1, :);
    X = [X; X(1,:)];
    m = m + 1;
    isevenM = true;
end

% I = eye(m); A = I + fliplr(I); 
% A = A(1:(m-1)/2+1, :); A((m-1)/2+1, (m-1)/2+1) = 1;
% 
% % Solution to underdetermined system A*(X + Y) = 0 with smallest Frobenius
% % norm: 
% C = A\(A*X);
C = 0.5*(X(1:m, :)+X(m:-1:1, :));

% Update coeff matrix: 
X = X - C;

% If m is even we need to remove the mode that was appended 
if ( isevenM )
    X(1, :) = (X(1, :)+X(end, :));
    X = X(1:m-1, :);
end

ctechs = real(trigtech({'', X}));
f.cols.funs{1}.onefun = ctechs;

% Now operate on the rows. The coefficients for the rows of an odd BMCI
% function should only contain odd wave numbers. The projection is to
% simply zero out the even wave numbers.
X = f.rows.funs{1}.onefun.coeffs;
n = size(X, 1); 
zeroMode = floor(n/2) + 1;
evenModes = [fliplr(zeroMode-2:-2:1) zeroMode:2:n];
X(evenModes, :) = 0;

rtechs = real(trigtech({'', X}));
f.rows.funs{1}.onefun = rtechs;

% Weird feval behavior in chebfun requires this
f.cols.pointValues = feval(ctechs, [-1; 1]);
f.rows.pointValues = feval(rtechs, [-1; 1]); 

end