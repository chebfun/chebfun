function [P, lam] = pswf(N, c, dom, output_type)
%PSWF   Prolate spheroidal wave functions.
% P = PSWF(N, C) returns a CHEBFUN P representing the Nth prolate spheroidal
% wave function with bandwidth C on the interval [-1,1], i.e., the Nth
% eigenfunction of the differential eigenvalue problem
%
%    [(1-x^2)*P(x)')' + (LAM - C^2*x^2)*P(x) = 0,
%
% with N = 0,1,2,....  C must be a positive scalar but N may be a vector of
% non-negative integers, in which case the output is an array-valued CHEBFUN
% with LENGTH(N) columns. P is scaled so that P'*P = 2/(2N+1), with the sign
% such that sign(P(0)) = (-1)^(N/2) if N is even and sign(P'(0)) =
% (-1)^((N-1)/2) if N is odd (see [3], eq. (30.4.1)).
%
% [P, LAM] = PSWF(N, C) also returns the eigenvalue(s).
%
% Examples:
%    f = pswf(2,pi); f(0.3)
%    plot(pswf(0:2:6, 100))
%
% [1] H. Xiao, V. Rokhlin and N. Yarvin, Prolate spheroidal wavefunctions,
% quadrature and interpolation, Inverse Problems, 17 (2001), 805-838.
%
% [2] https://reference.wolfram.com/language/ref/SpheroidalPS.html
% 
% [3] https://dlmf.nist.gov/30.4
%
% See also PSWFPTS, LEGPOLY.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: The approach is to compute the (approximate) normalised
% Legendre coefficients of the PSWFs by solving an eigenvalue problem [1].
% The Legendre coefficients are then converted to Chebyshev via LEG2CHEB,
% and a Chebfun constructed.
%
% There is functionality for scaling the domain, but this is currently
% undocumented, since we are not sure at present whether this choice
% is the right one.
%
% Copyright 2020 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults:
if ( nargin < 3 )
    dom = [-1 1];
    output_type = 'chebfun';
end
if ( nargin == 3 )
    if ( isnumeric(dom) )   
        output_type = 'chebfun';
    else
        output_type = dom;
        dom = [-1 1];
    end
end

% Parse inputs:
assert( nargin >= 2, 'PSWF requires at least two input arguments.')
assert( all(round(N)==N) && all(N>=0) && isvector(N), ...
    'Input N must be vector of non-negative integers.');
assert( (numel(c)==1) && (c>0) , ...
    'Input C must be a positive scalar.');
assert( numel(dom)==2 && all(isfinite(dom)) , ...
    'Domain must be a finite two-vector.');

% Set discretisation size. Heuristic estimates for initialisation.
M = max(ceil([2*sqrt(c)*N, 2*c, 20]));

% Increase discretisation size until the trailing Legendre coefficients are
% sufficiently small:
ishappy = 0;
tol = eps;
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
    [lame, idx] = sort(diag(De), 'ascend');
    Ve = Ve(:,idx);
    [Vo, Do] = eig(Ao);
    [lamo, idx] = sort(diag(Do), 'ascend');
    Vo = Vo(:,idx);
    
    % Reassemble full V and eigenvalues:
    V = zeros(M+1,M+1);
    V(1:2:end,1:2:end) = Ve;
    V(2:2:end,2:2:end) = Vo;
    lam = zeros(M+1,1);
    lam(1:2:end) = lame;
    lam(2:2:end) = lamo;
    
    % Check discretisation size was large enough;
    ishappy = sum(abs(V(end-3:end,N+1)))/(2*length(N)) < tol;
    if ( ~ishappy )
        M = 2*M;
    end
    
    % Failsafe:
    count = count + 1;
    if ( count > 10 )
        break
    end
    
end

% Extract required columns and unnormalise the Legendre coeffients:
V = bsxfun(@times, V(:,N+1), sqrt((0:M)'+1/2) );
lam = lam(N+1);

% Trim trailing coefficients that are below machine precision:
M = max(abs(V), [], 2);
idx = find(M > eps, 1, 'last');
V = V(1:idx,:);

% Scale as per Wolfram Alpha definition [1]:
V = (1./sqrt(N+0.5)).*V;

%%

% Enforce P_N(0) > 0 for N even and P'_N(0) > 0 for N odd.
m = 0:(idx-1)/2; idx = ~mod(N,2);
if ( any(idx) )
    L0 = (-1).^m./(beta(m,.5).*m); L0(1) = 1;
    V0 = L0*V(1:2:end,idx);
    currentSign = sign(V0);
    desiredSign = 1 - mod(N(idx),4);
    scl = currentSign.*desiredSign;
    V(:,idx) = scl.*V(:,idx);
end
if ( ~all(idx) )
    Lp0 = 2*(-1).^m./beta(m+1,.5); Lp0 = Lp0(1:length(V(2:2:end,1)));
    Vp0 = Lp0*V(2:2:end,~idx);
    currentSign = sign(Vp0);
    desiredSign = 1 - mod(N(~idx)-1,4);
    scl = currentSign.*desiredSign;
    V(:,~idx) = scl.*V(:,~idx);
end

% Quit now if only coefficients are required:
if ( strcmpi(output_type, 'coeffs') )
    P = V;
    return
end

% Convert Legendre coeffs to Chebyshev coeffs:
W = leg2cheb(V);

% Enforce even/oddness (which is lost in leg2cheb):
idx = logical(mod(N,2));
W(1:2:end,idx) = 0;
W(2:2:end,~idx) = 0;

% Create a Chebfun from the Chebyshev  coefficients:
P = chebfun(W, dom, 'coeffs');

% The coefficients are trimmed in V, so simplifying should not be necessary.
% P = simplify(P);    

end
