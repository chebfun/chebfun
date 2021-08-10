function [f, finv, pol, polinv] = conformal(C, varargin)
%% CONFORMAL  Conformal map to unit disk
%   [F, FINV] = CONFORMAL(C, ctr) computes a conformal map F of the region
%   bounded by the complex periodic chebfun C to the unit disk and its inverse
%   FINV, with F(ctr) = 0 and F'(ctr) > 0.  Both maps are represented by
%   function handles evaluating rational functions.  If ctr is omitted it is
%   set to 0.
%
%   CONFORMAL(..., 'tol', tol) uses tolerance tol instead of the default 1e-5.
%
%   CONFORMAL(..., 'plots') produces plots of the map and its inverse.
%
%   CONFORMAL(..., 'numbers') prints various quantities.
%
%   CONFORMAL(..., 'poly') uses a less robust algorithm based on polynomials
%                          instead of the Kerzman-Stein integral equation.
%
%   [F, FINV, POL, POLINV] = CONFORMAL(...) returns the poles of F and FINV.
%
%   This experimental code is good for smooth simple regions, but easy to break.
%
%   Examples:
%
%   C = chebfun('exp(pi*1i*t)','trig')*(1+.15*randnfun(.2,'trig'));
%   [f, finv] = conformal(C, 'plots');                               % random
%
%   C = chebfun('2*cos(t)+1i*(sin(t)+2*cos(t).^3)',[0 2*pi],'trig'); % Ellacott
%   conformal(C, 'poly', 'plots');                                   % blade
%
%   C = chebfun('exp(pi*1i*t)*(1+.3*cos(6*pi*t))','trig');
%   conformal(C, 'plots', 'tol', 1e-8, 'numbers');                   % snowflake
%
%   s = chebfun('s'); C = join(s-.5i,1+.5i*s,.5i-s,-1-.5i*s);
%   conformal(C, 'poly', 'plots');                                   % rectangle
%
% References:
% 
% A. Gopal and L. N. Trefethen, Representation of conformal maps by
% rational functions, Numer. Math. 142 (2019), 359--382.
%
% L. N. Trefethen, Numerical conformal mapping with rational 
% functions, Comp. Meth. Funct. Th. 20 (2020), 369-387.
%
% See also CONFORMAL2.

%%
%   Default algorithm:
%     (1) Solve discretized Kerzman-Stein integral equation
%     (2) Use AAA to approximate f and its inverse by rational functions
%
%   Alternative algorithm if 'poly' is specified:
%     (1) Use least-squares to find harmonic u s.t. u(z) = -log(abs(z-ctr)) on C
%     (2) Set w = u+iv (v = harmonic conjugate of u)
%     (3) Set f = z*exp(w(z-ctr));
%     (4) Use AAA to approximate f and its inverse by rational functions
%
%   Although in principle these algorithms should be embedded in the
%   Chebfun constructor, for simplicity in this numerically challenging
%   area we have not done that.
%
% This code was written by L. N. Trefethen in September 2019.  The
% Kerzman-Stein part originates with Anne Greenbaum and Trevor Caldwell.

t1 = tic;
[ctr, tol, plots, numbers, poly] = parseinputs(C, varargin{:});
err = Inf;
scl = norm(C-ctr, inf);

if poly == 0               % DEFAULT ALGORITHM: KERZMAN-STEIN INTEGRAL EQUATION

    M = 300;
    while err > tol
        M = M + 300;
        [g, Z, W] = kerzstein((C-ctr)/scl, M, 0);
        Z = Z*scl + ctr;
        gc = trigcoeffs(g);
        err = norm(gc([1:10 end-9:end]));           % a crude error measure
        if (err > tol) && (M >= 1200)
            warning('CONFORMAL did not converge')
            break
        end
    end

else                       % ALTERNATIVE ALGORITHM: POLYNOMIAL EXPANSION

    w1 = warning('off', 'MATLAB:rankDeficientMatrix');
    logn = 4;
    dom = domain(C); 
    dom = dom([1 end]);
    while err > tol
        n = round(2^logn);                    % degree of polynomial
        M = 8*n;                              % number of sample points
        logn = logn + 1/2;
        Z = C(dom(1) + (1:M)'*diff(dom)/M);   % sample points
        Zscl = (Z-ctr)/scl;                   % rescale to improve conditioning
        G = -log(abs(Zscl));                  % RHS for Dirichlet problem
        n1 = 0:n;                             % exponents for real part
        n2 = 1:n;                             % of exponents for imag part
        Q = ones(size(Zscl));                 % Arnoldi.  Q has orthog cols of norm sqrt(M).
        H = zeros(n+1,n);  
        for k = 1:n
            v = Zscl.*Q(:,k);
            v = v - Q*(Q'*v)/M;               % or execute twice for better orthogonality!
            H(k+1,k) = norm(v)/sqrt(M);       % At end, Zscl.*Q(:,1:n)/M = Q*H
            Q = [Q v/H(k+1,k)];               % Q has orthog cols of norm sqrt(M)
        end
        A = [real(Q) imag(Q(:,2:end))];
        c = A\G;                              % soln of least-squares problem
        err = norm(A*c-G, inf);               % max error
        cc = c(n1+1)-1i*[0; c(n+1+n2)];       % coeffs for harmonic -> analytic
        W = Zscl.*exp(Q*cc);                  % images on boundary
        if (err > tol) && (logn >= 9.5)
            warning('CONFORMAL did not converge')
            break
        end
    end
    warning(w1.state, 'MATLAB:rankDeficientMatrix')

end

w2 = warning('off', 'CHEBFUN:aaa:Froissart');
[f0, pol] = aaa(W, Z, 'tol', tol);            % forward map
zz  = 1e-4*scl*[1 1i -1 -1i];                 % finite diff for simplicity
dwdz = sum(f0(ctr+zz)./zz);                   % derivative at ctr
rot = exp(-1i*angle(dwdz));                   % rotation for f'(ctr) > 0
f = @(z) rot*f0(z);                           % rotate mapping function
W = rot*W;                                    % rotate points on circle
[finv, polinv] = aaa(Z, W, 'tol', tol);       % inverse map
warning(w2.state, 'CHEBFUN:aaa:Froissart')

inC = inpolygon(real(pol), imag(pol), ...     % check for poles of f in region
	real(Z), imag(Z));
if max(inC) > 0
warning('CONFORMAL: pole in region')
end
if min(abs(polinv)) < 1                       % check for poles of finv in disk
    warning('CONFORMAL: pole in disk')
end
tcomp = toc(t1);

% Plot solution if requested
if plots
    t2 = tic;
    clf
    LW = 'linewidth'; PO = 'position'; MS = 'markersize';
    FW = 'fontweight'; NO = 'normal'; XT = 'xtick'; YT = 'ytick';
    circ = exp(2i*pi*(0:300)/300);
  
    h1 = axes(PO, [.04 .38 .45 .53]);         % plot C and poles of f
    plot(C, 'b', LW, 1)
    axis([real(ctr)+1.4*scl*[-1 1] imag(ctr)+1.4*scl*[-1 1]])
    axis square, hold on, plot(pol, '.r', MS, 8)
    title([num2str(length(pol)) ' poles'],FW,NO)
  
    h2 = axes(PO, [.52 .38 .45 .53]);         % plot disk and poles of finv
    plot(circ, 'b', LW, 1), axis(1.6*[-1 1 -1 1])
    axis square, hold on, set(gca,XT,-1:1,YT,-1:1)
    plot(polinv, '.r', MS, 8)
    title([num2str(length(polinv)) ' poles'],FW,NO)
  
    ncirc = 8;                  % plot concentric circles and their images
    for r = (1:ncirc-1)/ncirc
        axes(h1), plot(finv(r*circ), '-k', LW, .5)
        axes(h2), plot(r*circ, '-k', LW, .5)
    end

    nrad = 16;                  % plot radii and their images
    ray = chebpts(301);
    ray = ray(ray>=0);
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
    fprintf('                         computation time in seconds:%6.2f\n', tcomp)
    if plots
        fprintf('                            plotting time in seconds:%6.2f\n', tplot)
    end
    fprintf('            number of sample points Z on boundary, M:  %d\n', M)
    fprintf('                      numbers of poles of f and finv:  %d, %d\n', ...
        length(pol), length(polinv))
    fprintf('back-and-forth boundary error norm(Z-finv(f(Z)),inf):  %6.1e\n',...
        norm(Z-finv(f(Z)), inf))
    fprintf(' inverse back-and-forth error norm(W-f(finv(W)),inf):  %6.1e\n',...
        norm(W-f(finv(W)), inf))
    fprintf(' interior inverse error norm(.9*W-f(finv(.9*W)),inf):  %6.1e\n',...
        norm(.9*W-f(finv(.9*W)), inf))
    if poly
        fprintf('             max error of least-squares problem, err:  %6.1e\n', err)
        fprintf('                                polynomial degree, n:  %d\n', n)
        fprintf('                number of real degrees of freedom, N:  %d\n', 2*n+1)
        fprintf('          condition number of least-squares matrix A:  %6.1e\n', cond(A))
    else
        fprintf('                            rough error measure, err:  %6.1e\n', err)
    end
    disp(' ')
end

end   % of conformal

function [ctr, tol, plots, numbers, poly] = parseinputs(C, varargin)
ctr = 0;                        % defaults
tol = 1e-5; 
plots = 0;
numbers = 0;
poly = 0;
j = 1;
while j < nargin
    j = j+1;
    v = varargin{j-1};
    if isnumeric(v)
        ctr = v;
    elseif strcmp(v, 'tol')
        j = j+1;
        tol = varargin{j-1};
    elseif strcmp(v, 'plots')
        plots = 1;
    elseif strcmp(v, 'numbers')
        numbers = 1;
    elseif strcmp(v, 'poly')
        poly = 1;
    end
end
winding_number = sum(diff(C)/(C-ctr))/(2i*pi);
if abs(winding_number - 1) > tol
    error('CONFORMAL:parseinputs','C must wind once counterclockwise around ctr')
end
end   % end of parseinputs

function [g, Z, W] = kerzstein(C, M, ctr)
%
% Given a smooth periodic chebfun C defining the boundary of a region Omega 
% in the counterclockwise direction, this function computes the boundary
% correspondence function g for a conformal mapping from Omega to the unit
% disk, mapping ctr to 0.  It first finds M points (Z) on C that are equally
% spaced in arclength.  It then solves an integral equation to determine the
% images of these points (W) on the unit circle.  Finally, a periodic
% chebfun g is defined to map the arclengths associated with Z to their
% images on the unit circle.  The original version of this code cones
% from Anne Greenbaum and Trevor Caldwell.
% Reference:  Kerzman and Trummer, "Numerical conformal mapping via the Szego
% kernel," Journal of Computational and Applied Mathematics 14 (1986), 111-123.

% Express C as a function of arclength and compute its derivative

S = arcLength(C);
s = cumsum(abs(diff(C)));
dC_arg = angle(diff(C));
dC_arg = unwrap(dC_arg - dC_arg(0));
u = inv(s);
C = newDomain(C, minandmax(u));
D = C(u);
Dprime = diff(D);

%  Set up and solve the integral equation

ds = S/M;
svec = [0:M-1]'*ds;
Dvec = D(svec);
Z = Dvec;
Dprimevec = Dprime(svec);
gamdotvec = Dprimevec ./ abs(Dprimevec);   % unit tangents
d = 1/(2i*pi);
gvec = d*conj(gamdotvec./(ctr - Dvec));
A = eye(M);
for j = 1:M
    w = Dvec(j);
    for i = 1:M
        z = Dvec(i);
        if i ~= j
            A(i,j) = A(i,j) - d*(conj(gamdotvec(i)/(z-w)) + gamdotvec(j)/(w-z))*ds;
        end
    end
end
fvec = A\gvec;

%  Compute boundary correspondence function g

Rprimevec = fvec.^2;
Rvec = -1i*gamdotvec.*(Rprimevec./abs(Rprimevec));
W = Rvec;
g = chebfun(Rvec, [0,S], 'trig');

end   % of kerzstein
