function f = projectOntoBMCII(f)
%PROJECTONTOBMCII   Projection onto BMC-II symmetry.
%   F = projectOntoBMCII(F) is the orthogonal projection of f onto BMC-II 
%   symmetry, i.e., a function that is
%   1. even in r for every even wave number in theta.
%   2. odd in r for every odd wave number in theta.
%   Additionally, for all but the k=0 mode in theta, the resulting 
%   projection enforces the diskfun is zero at the origin.
%
%   The projection is orthogonal, i.e., the correction matrix to fix up the
%   structure has the smallest possible Frobenius norm.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Even part
feven = f;
feven.cols = feven.cols(:, feven.idxPlus);
feven.rows = feven.rows(:, feven.idxPlus);
feven = projectOntoEvenBMCII(feven);

% Odd part
fodd = f;
fodd.cols = fodd.cols(:, fodd.idxMinus);
fodd.rows = fodd.rows(:, fodd.idxMinus);
fodd = projectOntoOddBMCII(fodd);

% Put pieces back together.
f.cols(:, f.idxPlus) = feven.cols;
f.rows(:, f.idxPlus) = feven.rows;

f.cols(:, f.idxMinus) = fodd.cols;
f.rows(:, f.idxMinus) = fodd.rows;

end

function f = projectOntoEvenBMCII(f)
% Project a diskfun to have even BMC-II symmetry, i.e., a diskfun that
% is pi-periodic in theta and even in r. The projection is orthogonal,
% i.e., the correction matrix to fix up the structure has the smallest
% possible Frobenius norm.

% Nothing to project
if isempty(f)
    return;
end

% Operate on the column coefficients first.
X = f.cols.funs{1}.onefun.coeffs;

% Get size: 
[m, n] = size(X); 

 % First we enforce that the every function in r is even: 
 X(2:2:end, :) = 0; 
 
% The rest of the code now needs to operate on the remaining even,
% non-zero modes in theta.
evenModes = 2:n;
    

if ( ~isempty(evenModes) )
    % to enforce that f(r) = 0 for each column function, 
    % we want C with smallest 2-norm such that  A*(X(1:2:end, 2:end)+C) = 0, 
    % where A =[1 -1 1 -1..] . 
    % The solution is 
    % C = A'*((A*A')\(A*X)). 
    %odd modes in r are zero, they won't contribute.
    even = 1:2:m;
    Xe = X(even, evenModes); 
    factor = 1/length(even)*(sum(Xe(1:2:end, :),1)-sum(Xe(2:2:end, :),1));
    C = ((-1*ones(length(even), 1)).^((2:length(even)+1)'))*factor; 
    %now add C to X
    X(even, evenModes) = Xe+C; 
     
end

ctechs = chebtech2({'',X}); 
f.cols.funs{1}.onefun = ctechs;

% Now operate on the rows. The coefficients for the rows of an even BMC-II
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

function f = projectOntoOddBMCII(f)
% Project a diskfun to have odd BMC-II symmetry, i.e., a diskfun that is
% pi-anti-periodic in theta and even in r. The projection is orthogonal, 
% i.e., the correction matrix to fix up the structure has the smallest 
% possible Frobenius norm.

% Nothing to project
if isempty(f)
    return;
end

% to enforce odd, simply zero out even coeffs: 

% Operate on the column coefficients first to project them onto odd
% functions.
X = f.cols.funs{1}.onefun.coeffs;
X(1:2:end, :) = 0; 

ctechs = chebtech2({'',X}); 
f.cols.funs{1}.onefun = ctechs;

% Now operate on the rows. The coefficients for the rows of an odd BMC-II
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