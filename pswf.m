function [P, V] = pswf(N, c, dom, flag)
%PSWF   Prolate spheroidal wave functions.
% P = PSWF(N, C) computes a CHEBFUN representing the Nth prolate spheroidal
% wavefunction (PSWF) with bandwidth C on the interval [-1,1]. C must be a
% scalar but N may be a vector of integers, in which case the output is
% an array-valued CHEBFUN with NUMEL(N) columns. 
%
% P = PSWF(N) is as above, but with assumption that C = max(N).
%
% P = PSWF(N, C, DOM) computes the PSWFs as above, but scaled to the interval 
% DOM, which must be a finite two-vector.
%
% [P, V] = PSWF(...) returns also a the matrix of Legendre coefficients for
% the computed PSWF. V = PSWF(N, C, DOM, 'coeffs') returns only the
% coefficient matrix and not the CHEBFUN P.
%
% See also PSWFPTS.

% Copyright 2020 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: The approach is to compute the (approximate) normalised
% Legendre coefficients of the PSWFs by solving an eigenvalue problem [1].
% The Legendre coefficients are then converted to Chebyshev via LEG2CHEB,
% and a Chebfun constructed.
%
% [1] H Xiao, V Rokhlin and N Yarvin, Prolate spheroidal wavefunctions,
% quadrature and interpolation, Inverse Problems, 17 (2001) 805–838.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults:
if ( nargin < 2 )
    c = max(N);
end
if ( nargin < 3 )
    dom = [-1 1];
end
if ( nargin < 4 )
    flag = 'chebfun';
end

% Set discretisation size. Heuristic estimates for initialisation.
M = max(ceil([2*sqrt(c)*N, 2*c, 20]));

% Increase discretisation size until the trailing Legendre coefficients are
% sufficiently small:
ishappy = 0;
tol = 1e-14;
count = 0;

while ( ~ishappy )

    % Construct the matrix (see Xiao et al):
    j = (0:M).';
    Asub = c^2*j.*(j-1)./((2*j-1).*sqrt((2*j-3).*(2*j+1)));
    Adia = j.*(j+1) + c^2*(2*j.*(j+1)-1)./((2*j+3).*(2*j-1));
    Asup = c^2*(j+2).*(j+1)./((2*j+3).*sqrt((2*j+5).*(2*j+1)));
    A = diag(Asub(3:end), -2) + diag(Adia, 0) + diag(Asup(1:end-2), 2);
    
    % Split in to even and odd parts for efficiency/accuracy. 
    Ae = A(1:2:end,1:2:end);
    Ao = A(2:2:end,2:2:end);
    
    % Compute (sorted) eigenvectors:
    [Ve, De] = eig(Ae);
    [~, idx] = sort(diag(De), 'ascend');
    Ve = Ve(:,idx);
    [Vo, Do] = eig(Ao);
    [~, idx] = sort(diag(Do), 'ascend');
    Vo = Vo(:,idx);
    
    % Reassemble full V:
    V = zeros(M+1,M+1);
    V(1:2:end,1:2:end) = Ve;
    V(2:2:end,2:2:end) = Vo;
    
    % Check discretisation size was large enough;
    ishappy = sum(abs(V(end-3:end,N)))/(2*numel(N)) < tol;
    if ( ~ishappy )
        M = 2*M;
    end
    
    % Failsafe:
    count = count + 1;
    if ( count > 10 )
        break
    end
    
end

% Extract required columns and unnormalise:
V = bsxfun(@times, V(:,N), sqrt((0:M)'+1/2) );

% Quit now if only coefficients are required:
if ( strcmpi(flag, 'coeffs') )
    P = V;
    return
end

% Convert to Chebyshev coeffs:
W = leg2cheb(V);

% Enforce even/oddness (which is lost in leg2cheb):
idx = abs(W(1,:)) < tol;
W(1:2:end,idx) = 0;
W(2:2:end,~idx) = 0;

% Create a Chebfun from the Chebyshev  coefficients:
P = chebfun(W, dom, 'coeffs');
P = simplify(P);    

end