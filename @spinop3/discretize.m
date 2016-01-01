function [L, Nc] = discretize(S, N)
%DISCRETIZE   Discretize a SPINOP3.
%   [L, NC] = DISCRETIZE(S, N) uses a Fourier spectral method in coefficient 
%   space to discretize the SPINOP3 S with N grid points in each direction. L is 
%   the linear part, a N^3xN^3 diagonal matrix stored as a NxNxN tensor, and NC 
%   is the diffenriation term of the nonlinear part (and hence is linear).

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

% Compute the N^3xN^3 Laplacian with KRON:
I = eye(N);
lapmat = kron(kron(D2, I), I)  + kron(kron(I, D2), I) + kron(kron(I, I), D2);

% Create a NxNxN tensor with the diagonal of the N^3xN^3 Laplacian:
lapmat = reshape(full(diag(lapmat)), N, N, N);

% Get the constants in front of the Laplacians:
strL = func2str(funcL);
strL = strrep(strL, 'laplacian', '');
funcL = eval(strL);
inputs = cell(1, nVars);
for k = 1:nVars
   inputs{k} = 1; 
end
constants = feval(funcL, inputs{:}); 
L = [];
for k = 1:nVars
    L = [L; constants(k)*lapmat]; %#ok<*AGROW>
end

%% Disretize the differentiation term of the nonlinear part:

% We only support no differentiation, i.e., NC = 1.
% [TODO]: Support diff order > 1.
Nc = 1;

end