function [L, Nc] = discretize(S, N)
%DISCRETIZE   Discretize a SPINOP.
%   [L, NC] = DISCRETIZE(S, N) uses a Fourier spectral method in coefficient 
%   space to discretize the SPINOP S with N grid points. L is the linear part, a 
%   diagonal matrix stored as a vector, and NC is the diffenriation term of the 
%   nonlinear part (and hence is linear).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Set-up:

% Get the domain DOM, the linear part LFUN, the nonlinear part NFUN, and the 
% number of variables NVARS from S:
dom = S.domain;
funcL = S.linearPart;
funcNc = S.nonlinearPartCoeffs;
nVars = nargin(funcL);

% Create a CHEBOPPREF object with TRIGSPEC discretization:
pref = cheboppref();
pref.discretization = @trigspec;

%% Discretize the linear part:

% Add the dependent variable 'x': 
strL = func2str(funcL);
strL = strrep(strL, '@(', '@(x,');
funcL = str2func(strL);

% USE CHEBOP, LINOP and MATRIX with TRIGSPEC discretization:
matL = matrix(linop(chebop(funcL, dom)), N, pref);
L = [];
for l = 1:nVars
    idx = (l-1)*N + 1;
    temp = full(diag(matL(idx:idx+N-1,idx:idx+N-1)));
    if ( mod(N, 2) == 0 )
        temp = fftshift(temp);
    else
        temp = ifftshift(temp);
    end
    L = [L; temp]; %#ok<*AGROW>
end

%% Disretize the differentiation term of the nonlinear part:

% If NC if a FUNCTION HANDLE, discretize it with TRIGSPEC:
if ( isa(funcNc, 'function_handle') == 1 )
    
    % Add the dependent variable 'x':
    strNc = func2str(funcNc);
    strNc = strrep(strNc, '@(', '@(x,');
    funcNc = str2func(strNc);
    
    % USE CHEBOP, LINOP and MATRIX with TRIGSPEC discretization:
    matNc = matrix(linop(chebop(funcNc, dom)), N, pref);
    Nc = [];
    for l = 1:nVars
        idx = (l-1)*N + 1;
        temp = full(diag(matNc(idx:idx+N-1,idx:idx+N-1)));
        if ( mod(N, 2) == 0 )
            temp = fftshift(temp);
        else
            temp = ifftshift(temp);
        end
        Nc = [Nc; temp]; %#ok<*AGROW>
    end

% NC might be 1 and in that case, nothing to do:
elseif ( isa(funcNc, 'double') == 1 )
    Nc = funcNc;
end

end
