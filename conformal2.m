function [f, finv, rho, pol, polinv] = conformal2(C1, C2, varargin)
%% CONFORMAL2  Doubly-connected conformal map to annulus
%   [F, FINV] = CONFORMAL2(C1, C2) computes a conformal map F of the smooth
%   annular region bounded by the two complex periodic chebfuns C1 (outer
%   boundary) and C2 (inner boundary, enclosing the origin) to a circular
%   annulus and its inverse FINV.  Both maps are represented by function
%   handles evaluating rational functions.  The circular annulus has radii 1
%   and RHO < 1, the conformal modulus, which is determined as part of the
%   calculation.
%
%   CONFORMAL2(..., 'tol', tol) uses tolerance tol instead of the default 1e-6.
%
%   CONFORMAL2(..., 'plots') produces plots of the map and its inverse.
%
%   CONFORMAL2(..., 'numbers') prints various quantities.
%
%   [F, FINV, RHO, POL, POLINV] = CONFORMAL2(...) returns the conformal modulus
%   and the poles of F and FINV.
%
%   This experimental code is only good for very simple regions.  Easy to break.
%
%   Examples:
%
%   circle = chebfun('exp(1i*pi*z)','trig');
%   ellipse2 = real(circle) + .6i*imag(circle);
%   ellipse1 = (2+1i)*ellipse2 + .5;
%   conformal2(ellipse1, ellipse2, 'plots');        
%
%   z = chebfun('exp(1i*pi*z)','trig');
%   C1 = z*abs(1+.1*z^4); C2 = .5*z*abs(1+.2*z^3);
%   conformal2(C1, C2, 'plots');        
%
%   z = chebfun('exp(pi*1i*t)','trig');
%   C1 = z*(1+.1*randnfun(.5,'trig'));
%   C2 = z*(.5+.1*randnfun(1,'trig'));
%   [f, finv, rho] = conformal2(C1, C2, 'plots', 'numbers'); 
%
% See also CONFORMAL.

%%
%   Algorithm:
%     (1) Use least-squares to find constant rho < 1 and harmonic u s.t.
%         u(z) = -log(abs(z)) on C1 and  u(z) = -log(abs(z)) + log(rho) on C2
%     (2) Set f = z*exp(u(z)+iv(z));
%     (3) Use AAA to approximate f and its inverse by rational functions
%
%   Although in principle these algorithms should be imbedded in the
%   Chebfun constructor, for simplicity in this numerically challenging
%   area we have not done that.
%
% This code was written by L. N. Trefethen mainly in October 2019.  For
% information about the use of AAA approximation, see Gopal and Trefethen,
% Representation of conformal maps by rational functions, Numer. Math.
% 142 (2019), 359-382.  See also Trefethen, Numerical conformal mapping
% with rational functions, Comp. Meth. Funct. Thy. 20 (2020), 369-387.

t1 = tic;
[tol, plots, numbers, poly] = parseinputs(C1, C2, varargin{:});
err = Inf;

w1 = warning('off', 'MATLAB:rankDeficientMatrix');
logn = 4;
dom1 = domain(C1); dom1 = dom1([1 end]);    % endpoints of parameter range
dom2 = domain(C2); dom2 = dom2([1 end]);
errvec = [];
while err > tol
    n = round(2^logn);                      % degree of polynomial
    M = 8*n;                                % # of sample points on each curve
    MM = (1:M)';
    Z1 = C1(dom1(1) + MM*diff(dom1)/M);     % sample points, outer curve
    Z2 = C2(dom2(1) + MM*diff(dom2)/M);     % sample points, inner curve
    Z = [Z1; Z2];
    H = -log(abs(Z));                       % RHS for Dirichlet problem
    H(M+MM) = H(M+MM)+1;
    rvec = [zeros(M,1); ones(M,1)];
    [Hes, P] = VAorthog(Z,n);               % orthogonalize nonnegative powers
    [Hes2, P2] = VAorthog(Z.^(-1),n);       % orthogonalize negative powers
    A = [real(P) real(P2) -imag(P) -imag(P2) rvec];
    N = size(A,2);
    c = A\H;                                % soln of least-squares problem
    logrho = 1 - c(end);
    rho = exp(logrho);                      % conformal modulus
    err = norm(A*c-H, inf);                 % max error
    errvec = [errvec err];
    c(end) = [];
    c = reshape(c, [], 2)*[1; 1i];
    F = [P P2]*c;
    W = Z.*(exp(F));
    logn = logn + 1/2;
    if (err > tol) && (logn > 8)
        warning('CONFORMAL2 did not converge')
        break
    end
end
warning(w1.state, 'MATLAB:rankDeficientMatrix')
%semilogy(errvec,'.-','linewidth',1,'markersize',10), grid on, pause

w2 = warning('off', 'CHEBFUN:aaa:Froissart');
[f, pol] = aaa(W, Z, 'tol', tol);           % forward map
[finv, polinv] = aaa(Z, W, 'tol', tol);     % inverse map
warning(w2.state, 'CHEBFUN:aaa:Froissart')

%inC = inpolygon(real(pol), imag(pol), ...  % check for poles of f in region
%                real(Z), imag(Z));
%if max(inC) > 0
%    warning('CONFORMAL: pole in region')
%end
%if min(abs(polinv)) < 1                    % check for poles of finv in disk
%    warning('CONFORMAL: pole in disk')
%end
tcomp = toc(t1);

% Plot solution if requested
if plots
    t2 = tic;
    clf
    LW = 'linewidth'; PO = 'position'; MS = 'markersize';
    FW = 'fontweight'; NO = 'normal'; XT = 'xtick'; YT = 'ytick';
    circ = exp(2i*pi*(0:300)'/300);
  
    h1 = axes(PO, [.04 .38 .45 .53]);       % plot C1, C2 and poles of f
    plot([C1 C2], 'b', LW, 1)
    axis(1.2*norm(C1,inf)*[-1 1 -1 1])
    axis square, hold on, plot(pol, '.r', MS, 8)
    title([num2str(length(pol)) ' poles'],FW,NO)
  
    h2 = axes(PO, [.52 .38 .45 .53]);       % plot annulus and poles of finv
    plot([circ rho*circ], 'b', LW, 1), axis(1.5*[-1 1 -1 1])
    axis square, hold on, set(gca,XT,-1:1,YT,-1:1)
    plot(polinv, '.r', MS, 8)
    title([num2str(length(polinv)) ' poles'],FW,NO)
  
    ncirc = 8;                  % plot concentric circles and their images
    for r = rho + (1-rho)*(1:ncirc-1)/ncirc
        axes(h1), plot(finv(r*circ), '-k', LW, .5)
        axes(h2), plot(r*circ, '-k', LW, .5)
    end

    nrad = 16;                  % plot radii and their images
    ray = chebpts(101,[rho 1]);
    for k = 1:nrad
        axes(h1), plot(finv(ray*exp(2i*pi*k/nrad)), '-k', LW, .5)
        axes(h2), plot(ray*exp(2i*pi*k/nrad), '-k', LW, .5)
    end

    axes(h1), hold off
    axes(h2), hold off
    tplot = toc(t2);

end

% Print various quantities if requested
if numbers
    disp(' ')
    fprintf('                 computation time in seconds:%6.2f\n', tcomp)
    if plots
        fprintf('                    plotting time in seconds:%6.2f\n', tplot)
    end
    fprintf('umber of sample points Z on each boundary, M:  %d\n', M)
    fprintf('              numbers of poles of f and finv:  %d, %d\n', ...
        length(pol), length(polinv))
    fprintf('       boundary error norm(Z-finv(f(Z)),inf):  %6.1e\n',...
        norm(Z-finv(f(Z)), inf))
    fprintf('        inverse error norm(W-f(finv(W)),inf):  %6.1e\n',...
        norm(W-f(finv(W)), inf))
    fprintf('                 interior inverse error norm:  %6.1e\n',...
    norm(sqrt(rho)*W(1:M)-f(finv(sqrt(rho)*W(1:M))), inf))
    fprintf('     max error of least-squares problem, err:  %6.1e\n', err)
    fprintf('                        polynomial degree, n:  %d\n', n)
    fprintf('        number of real degrees of freedom, N:  %d\n', 2*n+1)
    disp(' ')
end

end   % end of conformal2

function [tol, plots, numbers, poly] = parseinputs(C1, C2, varargin)
tol = 1e-6; 
plots = 0;
numbers = 0;
poly = 0;
j = 2;
while j < nargin
    j = j+1;
    v = varargin{j-2};
    if strcmp(v, 'tol')
        j = j+1;
        tol = varargin{j-2};
    elseif strcmp(v, 'plots')
        plots = 1;
    elseif strcmp(v, 'numbers')
        numbers = 1;
    elseif strcmp(v, 'poly')
        poly = 1;
    end
end
%winding_number = sum(diff(C)/C)/(2i*pi);
%if abs(winding_number - 1) > tol
%    error('CONFORMAL:parseinputs','C must wind once counterclockwise around 0')
%end
end   % end of parseinputs
