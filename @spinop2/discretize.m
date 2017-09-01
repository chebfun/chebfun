function [L, Nc] = discretize(S, N)
%DISCRETIZE   Discretize a SPINOP2.
%   [L, NC] = DISCRETIZE(S, N) uses a Fourier spectral method in coefficient 
%   space to discretize the SPINOP2 S with N grid points in each direction. L is 
%   the linear part, a N^2xN^2 diagonal matrix stored as a NxN matrix, and NC is 
%   the differentiation term of the nonlinear part (and hence is linear).
%
% Remark: DISCRETIZE will fail to discretize SPINOP2 objects which are not of 
%         the right form. See HELP/SPINOP2.
%
% See also SPINOP2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Set-up:
 
% Get the domain DOM, the linear part LFUN, the nonlinear part NFUN, and the 
% number of variables NVARS from S:
dom = S.domain;
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
 
% Create a CHEBOPPREF object with TRIGSPEC discretization:
pref = cheboppref();
pref.discretization = @trigspec;

%% Discretize the linear part:

% Second-order Fourier differentiation matrix with TRIGSPEC (sparse diagonal 
% matrix):
D2 = trigspec.diffmat(N,2)*(2*pi/(dom(2) - dom(1)))^2;
if ( mod(N,2) == 0 )
    D2 = fftshift(D2);
else
    D2 = ifftshift(D2);
end

% Look for 'laplacian'/'lap', 'biharmonic'/'biharm', 'triharmonic'/'triharm',
% 'quadharmonic'/'quadharm' or 'quintharmomic'/'quintharm':
strL = func2str(funcL);
isLap = isempty(strfind(strL,'laplacian')) && isempty(strfind(strL,'lap'));
isLap = ~isLap;
isBih = isempty(strfind(strL,'biharmonic')) && isempty(strfind(strL,'biharm'));
isBih = ~isBih;
isTrih = isempty(strfind(strL,'triharmonic')) ...
    && isempty(strfind(strL,'triharm'));
isTrih = ~isTrih;
isQuadh = isempty(strfind(strL,'quadharmonic')) ...
    && isempty(strfind(strL,'quadharm'));
isQuadh = ~isQuadh;
isQuinth = isempty(strfind(strL,'quintharmonic')) ...
    && isempty(strfind(strL,'quintharm'));
isQuinth = ~isQuinth;

% NxN identity matrix for the Kronecker products:
I = eye(N);

% Construct the Laplacian operator -- needed for all the operators:
if ( isLap || isBih || isTrih || isQuadh || isQuinth )
    
    % Compute the N^2xN^2 Laplacian with KRON:
    lapmat = kron(I, D2) + kron(D2, I);
    
    % Create a NxN matrix with the diagonal of the N^2xN^2 Laplacian:
    lapmat = reshape(full(diag(lapmat)), N, N);
    
else
    lapmat = 0;
end

% The linear part has a B*biharmonic(u) term:
if ( isBih == 1 )
    
    % Pointwise multiplication since we only store the diagonal elements:
    bihmat = lapmat.^2;
    
else
    bihmat = 0;
end

% The linear part has a C*triharmonic(u) term:
if ( isTrih == 1 )
    
    % Pointwise multiplication since we only store the diagonal elements:
    trihmat = lapmat.^3;
    
else
    trihmat = 0;
end

% The linear part has a D*quadharmonic(u) term:
if ( isQuadh == 1 )
    
    % Pointwise multiplication since we only store the diagonal elements:
    quadhmat = lapmat.^4;
    
else
    quadhmat = 0;
end

% The linear part has a E*quintharmonic(u) term:
if ( isQuinth == 1 )
    
    % Pointwise multiplication since we only store the diagonal elements:
    quinthmat = lapmat.^5;
    
else
    quinthmat = 0;
end

% Convert to a string and initialize L:
strL = func2str(funcL);
L = [];

% Get the constants A in front of the Laplacians:
str = strrep(strL, 'laplacian', '');
str = strrep(str, 'lap', '');
str = strrep(str, 'biharmonic', '0*');
str = strrep(str, 'biharm', '0*');
str = strrep(str, 'triharmonic', '0*');
str = strrep(str, 'triharm', '0*');
str = strrep(str, 'quadharmonic', '0*');
str = strrep(str, 'quadharm', '0*');
str = strrep(str, 'quintharmonic', '0*');
str = strrep(str, 'quintharm', '0*');
func = eval(str);
inputs = cell(1, nVars);
for k = 1:nVars
   inputs{k} = 1; 
end
A = feval(func, inputs{:}); 

% Get the constants B in front of the biharmonic operators:
str = strrep(strL, 'laplacian', '0*');
str = strrep(str, 'lap', '0*');
str = strrep(str, 'biharmonic', '');
str = strrep(str, 'biharm', '');
str = strrep(str, 'triharmonic', '0*');
str = strrep(str, 'triharm', '0*');
str = strrep(str, 'quadharmonic', '0*');
str = strrep(str, 'quadharm', '0*');
str = strrep(str, 'quintharmonic', '0*');
str = strrep(str, 'quintharm', '0*');
func = eval(str);
B = feval(func, inputs{:}); 

% Get the constants C in front of the triharmonic operators:
str = strrep(strL, 'laplacian', '0*');
str = strrep(str, 'lap', '0*');
str = strrep(str, 'biharmonic', '0*');
str = strrep(str, 'biharm', '0*');
str = strrep(str, 'triharmonic', '');
str = strrep(str, 'triharm', '');
str = strrep(str, 'quadharmonic', '0*');
str = strrep(str, 'quadharm', '0*');
str = strrep(str, 'quintharmonic', '0*');
str = strrep(str, 'quintharm', '0*');
func = eval(str);
C = feval(func, inputs{:}); 

% Get the constants D in front of the quadharmonic operators:
str = strrep(strL, 'laplacian', '0*');
str = strrep(str, 'lap', '0*');
str = strrep(str, 'biharmonic', '0*');
str = strrep(str, 'biharm', '0*');
str = strrep(str, 'triharmonic', '0*');
str = strrep(str, 'triharm', '0*');
str = strrep(str, 'quadharmonic', '');
str = strrep(str, 'quadharm', '');
str = strrep(str, 'quintharmonic', '0*');
str = strrep(str, 'quintharm', '0*');
func = eval(str);
D = feval(func, inputs{:}); 

% Get the constants E in front of the quintharmonic operators:
str = strrep(strL, 'laplacian', '0*');
str = strrep(str, 'lap', '0*');
str = strrep(str, 'biharmonic', '0*');
str = strrep(str, 'biharm', '0*');
str = strrep(str, 'triharmonic', '0*');
str = strrep(str, 'triharm', '0*');
str = strrep(str, 'quadharmonic', '');
str = strrep(str, 'quadharm', '');
str = strrep(str, 'quintharmonic', '');
str = strrep(str, 'quintharm', '');
func = eval(str);
E = feval(func, inputs{:}); 

% Compute L:
for k = 1:nVars
    L = [L; A(k)*lapmat + B(k)*bihmat + C(k)*trihmat + D(k)*quadhmat + ...
        E(k)*quinthmat]; 
end

%% Disretize the differentiation term of the nonlinear part:

% We only support no differentiation, i.e., NC = 1.
% [TODO]: Support diff order > 1.
Nc = 1;

end