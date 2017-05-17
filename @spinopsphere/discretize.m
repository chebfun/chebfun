function [L, Nc] = discretize(S, N)
%DISCRETIZE   Discretize a SPINOPSPHERE.
%   [L, NC] = DISCRETIZE(S, N) uses a Fourier spectral method in coefficient 
%   space to discretize the SPINOPSPHERE S with N grid points in each direction. 
%   L is the linear part, a N^2xN^2 matrix, and NC is the differentiation term 
%   of the nonlinear part (and hence is linear).
%
% Remark: DISCRETIZE will fail to discretize SPINOPSPHERE objects which are not
%         of the right form. See HELP/SPINOPSPHERE.
%
% See also SPINOPSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Set-up:
 
% Get the domain DOM, the linear part LFUN, the nonlinear part NFUN, and the 
% number of variables NVARS from S:
funcL = S.lin;
nVars = nargin(funcL);

% Get the variables of the workspace:
func = functions(funcL);
wrk = func.workspace{1};
names = fieldnames(wrk);
if ( isempty(names) == 0 )
    lengthNames = size(names, 1);
    for k = 1:lengthNames
        eval(sprintf('%s = wrk.(names{k});', names{k}));
    end
end
 
%% Discretize the linear part:

% Construct Laplacian matrix (multiplied by Tsin2):
m = N; 
n = N;
Dm = spdiags(1i*[0,-m/2+1:m/2-1]', 0, m, n);
D2m = spdiags(-(-m/2:m/2-1).^2', 0, m, n);
D2n = spdiags(-(-n/2:n/2-1).^2', 0, m, n);
Im = speye(m); In = speye(n);
P = speye(m+1); P = P(:, 1:m); P(1,1) = .5; P(m+1,1) = .5;
Q = speye(m+1+4); Q = Q(3:m+2,:); Q(1,3) = 1; Q(1,m+3) = 1;
Msin2 = toeplitz([1/2, 0, -1/4, zeros(1, m+2)]);
Msin2 = sparse(Msin2(:, 3:m+3));
Tsin2 = round(Q*Msin2*P, 15);
Mcossin = toeplitz([0, 0, 1i/4, zeros(1, m+2)]);
Mcossin = sparse(Mcossin(:, 3:m+3));
Tcossin = round(Q*Mcossin*P, 15);
lapmat = kron(In, Tsin2*D2m + Tcossin*Dm) + kron(D2n, Im);
 
% Convert to a string and initialize L:
strL = func2str(funcL);

% Get the constants A in front of the Laplacians:
str = strrep(strL, 'laplacian', '');
str = strrep(str, 'lap', '');
func = eval(str);
inputs = cell(1, nVars);
for k = 1:nVars
   inputs{k} = 1; 
end
A = feval(func, inputs{:}); 

% Compute L:
L = cell(nVars, 1);
for k = 1:nVars
    L{k} = A(k)*lapmat;
end
L = blkdiag(L{:});

%% Disretize the differentiation term of the nonlinear part:

% We only support no differentiation, i.e., NC = 1.
% [TODO]: Support diff order > 1.
Nc = 1;

end