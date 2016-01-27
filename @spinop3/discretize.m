function [L, Nc] = discretize(S, N)
%DISCRETIZE   Discretize a SPINOP3.
%   [L, NC] = DISCRETIZE(S, N) uses a Fourier spectral method in coefficient 
%   space to discretize the SPINOP3 S with N grid points in each direction. L is 
%   the linear part, a N^3xN^3 diagonal matrix stored as a NxNxN tensor, and NC 
%   is the differentiation term of the nonlinear part (and hence is linear).
%
% Remark: DISCRETIZE will fail to discretize SPINOP3 objects which are not of 
%         the right form. See HELP/SPINOP3.
%
% See also SPINOP3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Set-up:
 
% Get the domain DOM, the linear part LFUN, the nonlinear part NFUN, and the 
% number of variables NVARS from S:
dom = S.domain;
funcL = S.linearPart;
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
 
% Create a CHEBOPPREF object with TRIGSPEC discretization:
pref = cheboppref();
pref.discretization = @trigspec;

%% Discretize the linear part:

% Second order Fourier differentiation matrix with TRIGSPEC (sparse diagonal 
% matrix):
D2 = trigspec.diffmat(N,2)*(2*pi/(dom(2) - dom(1)))^2;
if ( mod(N,2) == 0 )
    D2 = fftshift(D2);
else
    D2 = ifftshift(D2);
end

% Look for 'laplacian' and 'biharmonic':
strL = func2str(funcL);
isLap = ~isempty(strfind(strL, 'laplacian'));
isBih = ~isempty(strfind(strL, 'biharmonic'));

% NxN identity matrix for the Kronecker products:
I = eye(N);

% The linear part has a A*laplacian(u) term:
if ( isLap == 1 )
    
    % Compute the N^3xN^3 Laplacian with KRON:
    lapmat = kron(kron(I, I), D2) + kron(kron(I, D2), I) + kron(kron(D2, I), I);
    
    % Create a NxNxN tensor with the diagonal of the N^3xN^3 Laplacian:
    lapmat = reshape(full(diag(lapmat)), N, N, N);
    
else
    lapmat = 0;
end

% The linear part has a B*biharmonic(u) term:
if ( isBih == 1 )
    
    % Fourth order Fourier differentiation matrix:
    D4 = trigspec.diffmat(N,4)*(2*pi/(dom(2) - dom(1)))^4;
    if ( mod(N,2) == 0 )
        D4 = fftshift(D4);
    else
        D4 = ifftshift(D4);
    end
    
    % Compute the N^3xN^3 biharmonic operator with KRON:
    bihmat = kron(kron(I, I), D4) + kron(kron(I, D4), I) + ...
        kron(kron(D4, I), I) + 2*kron(kron(I, I), D2)*kron(kron(I, D2), I) + ...
        + 2*kron(kron(I, I), D2)*kron(kron(D2, I), I) + ...
        + 2*kron(kron(D2, I), I)*kron(kron(I, D2), I);
    
    % Create a NxNxN tensor with the diagonal of the N^3xN^3 biharmonic
    % operator:
    bihmat = reshape(full(diag(bihmat)), N, N, N);
    
else
    bihmat = 0;
end

% Convert to a string and initialize L:
strL = func2str(funcL);
L = [];

% Get the constants A in front of the Laplacians:
str = strrep(strL, 'laplacian', '');
str = strrep(str, 'biharmonic', '0*');
func = eval(str);
inputs = cell(1, nVars);
for k = 1:nVars
   inputs{k} = 1; 
end
A = feval(func, inputs{:}); 

% Get the constants B in front of the biharmonic operators:
str = strrep(strL, 'laplacian', '0*');
str = strrep(str, 'biharmonic', '');
func = eval(str);
B = feval(func, inputs{:}); 

% Compute L:
for k = 1:nVars
    L = [L; A(k)*lapmat + B(k)*bihmat]; %#ok<*AGROW>
end

%% Disretize the differentiation term of the nonlinear part:

% We only support no differentiation, i.e., NC = 1.
% [TODO]: Support diff order > 1.
Nc = 1;

end