function [x, w, v] = lagpts(n, int, meth, alf)
%LAGPTS  Laguerre points and Gauss-Laguerre quadrature weights.
%   LAGPTS(N) returns N Laguerre points X in (0,inf).
%
%   [X, W] = LAGPTS(N) returns also a row vector W of weights for Gauss-Laguerre
%   quadrature. [X, W, V] = LAGPTS(N) returns in addition a column vector V
%   of the barycentric weights corresponding to X.
%
%   LAGPTS(N, D) scales the nodes and weights for the semi-infinite domain D.
%   D can be either a domain object or a vector with two components.
%
%   [X, W] = LAGPTS(N, METHOD) allows the user to select which method to use.
%   METHOD = 'GW' will use the traditional Golub-Welsch eigenvalue method,
%   which is best for when N is small. 
%   METHOD = 'RH' will use Newton iterations on asymptotics of Laguerre polynomials. 
%   METHOD = 'RHW' is the same as 'RH' except it costs only O(sqrt(n)) 
%   operations. This is because it stops when the weights fall below 
%   realmin. The command [X, W] = lagpts(round( (n/17)^2), 'RHW') 
%   returns about n nodes and weights above realmin for large n.
%   METHOD = 'FAST' will use the Glaser-Liu-Rokhlin fast algorithm.
%   METHOD = 'REC' or 'RECW' uses Newton iterations on the recurrence relation.
%   METHOD = 'EXP' or 'EXPW' uses explicit expansions of nodes and weights.
%
%   [X, W] = LAGPTS(N, alpha) or LAGPTS(N, [0, inf], METHOD, alpha) will return
%   the nodes and weights for the generalised Laguerre polynomials with 
%   parameter alpha. Inaccurate results can be expected when alpha and n are both large.
%
%   By default, LAGPTS uses 'GW' with replacement of inaccurate small weights by
%   weights from the recurrence relation when N < 128, and else 'EXP'.
%
% See also CHEBPTS, LEGPTS, HERMPTS, and JACPTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTES AND REFERENCES:
% 
% Methods:
% 'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
% 'RH' by Peter Opsomer, June 2016 - algorithm adapted from [3], based on [5].
% 'REC' and 'EXP' by Peter Opsomer, August 2017 - algorithms adapted from [4].
% 
% References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature rules",
%       Math. Comp. 23:221-230, 1969.
%   [2] A. Glaser, X. Liu and V. Rokhlin, "A fast algorithm for the
%       calculation of the roots of special functions", SIAM Journal
%       on Scientific Computing", 29(4):1420-1438:, 2007.
%   [3] D. Huybrechs and P. Opsomer, Construction and implementation of asymptotic
%       expansions for Laguerre-type orthogonal polynomials", IMA J. Numer. Anal., 
%       published online on 29 July 2017, https://doi.org/10.1093/imanum/drx030
%   [4] D. Huybrechs and P. Opsomer, "Arbitrary-order asymptotic expansions of
%       generalized Gaussian quadrature rules", (in preparation).
%   [5] M. Vanlessen, "Strong asymptotics of Laguerre-Type orthogonal
%       polynomials and applications in Random Matrix Theory", Constr. Approx.,
%       25:125-175, 2007.
% 
% Historical note:
%   March 2009 - GW [1] algorithm.
%   April 2009 - GLR [2] added for N >= 129.
%   July 2016 - RH and RHW [3] replace GLR [2].
%   August 2017 - Added forward recurrence and explicit expansions [4].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults:
method = 'default';
interval = [0, inf];
alpha = 0;

if ( n < 0 )
    error('CHEBFUN:lagpts:n', 'First input should be a positive integer.');
end

% Return empty vector if n = 0.
if ( n == 0 )
    [x, w, v] = deal([]);
    return
end

% Check the inputs.
if ( nargin > 1 )
    if ( nargin == 4 )
        interval = int;
        method = meth;
        alpha = alf;
    elseif ( nargin == 3 )
        interval = int;
        method = meth;
    elseif ( nargin == 2 )
        if ( ischar(int) )
            method = int;
        elseif ( length(int) == 1 )
            alpha = int;
        else
            interval = int;
        end
    end
    if ( ~any(strcmpi(method, {'default', 'GW', 'fast', 'RH', 'RHW', 'REC', 'RECW', 'EXP', 'EXPW'})) )
        error('CHEBFUN:lagpts:inputs', 'Unrecognised input string %s.', method);
    elseif ( imag(alpha) ~= 0 ) || ( alpha < -1 )
        error('CHEBFUN:lagpts:inputs', ['alpha = ' num2str(alpha) ...
            ' is not allowed.']);
    end
    if ( numel(interval) > 2 )
        warning('CHEBFUN:lagpts:domain',...
            'Piecewise intervals not supported and will be ignored.');
        interval = interval([1, end]);
    end
end

if ( sum(isinf(interval)) ~= 1 )
    error('CHEBFUN:lagpts:inf', 'LAGPTS only supports semi-infinite domains.');
end

% Decide to use GW or RH.

if strcmpi(method,'GW')
    % GW, see [1]
    alph = 2*(1:n)-1 + alpha;  % 3-term recurrence coeffs
    beta = sqrt( (1:n-1).*(alpha + (1:n-1) ) );
    T = diag(beta,1) + diag(alph) + diag(beta,-1);  % Jacobi matrix
    [V, D] = eig(T);                      % eigenvalue decomposition
    [x, indx] = sort(diag(D));            % Laguerre points
    w = V(1,indx).^2;                     % Quadrature weights
    
elseif any(strcmpi(method, {'RH', 'RHW'}))
    % RH, see [3] and [5]
    [x, w] = newton(n, strcmpi(method, 'RHW'), alpha, 1);
    
elseif ( n < 128 ) && strcmpi(method,'default')
    % GW with forward recurrence for small weights
    [x,w] = lagpts(n, [0, inf], 'GW', alpha);
    ixs = find(w < 1e-10);
    for k = ixs
        [~, pder] = lagpnRecDer(n, alpha, x(k));
        w(k) = (n^2 +alpha*n)^(-1/2)/lagpnRecDer(n-1, alpha, x(k))/pder;
    end
    
elseif strcmpi(method,'fast')
    % Fast, see [2]
    [x, ders] = alg0_Lag(n);
    w = transpose(exp(-x)./(x.*ders.^2));
    
elseif any(strcmpi(method,{'rec', 'recW'}))
    % Forward recurrence
    [x, w] = newton(n, strcmpi(method, 'recW'), alpha, 0);
            
else
    % Explicit expansions, see [4]
    [x,w] = laguerreExp(n, strcmpi(method,'expW'), alpha);
    
end
w = (gamma(alpha+1)/sum(w))*w;    % Normalise so that sum(w) = Gamma(alpha+1)
v = (-1).^(0:n-1)'.*sqrt(w'.*x);  % Barycentric weights
v = v./max(abs(v));

% Nonstandard interval
if ( ~all(interval == [0, inf]) )
    a = interval(1);
    b = interval(2);
    if ( isinf(b) )
        x = x + a;
        w = w*exp(-a);
    else
        x = -x + b;
        w = w*exp(b);
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Newton Routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rh = 0 if forward recurrence, 1 if RH and 2 if RH with correction iso error
function [x, w] = newton(n, compRepr, alpha, rh)

if compRepr
    % Get a heuristic for the indices where the weights are about above realmin.
    mn = min(n,ceil(17*sqrt(n)));
else
    mn = n;
end
w = zeros(1, mn);

if ~rh % Use forward recurrence
    bes = besselroots(alpha, min(mn, 7) ).^2; 
    % [Tricomi 1947 pg. 296]
    bes = bes/(4*n + 2*alpha+2).*(1 + (bes + 2*(alpha^2 - 1) )/(4*n + 2*alpha+2)^2/3 );
    x = [ bes; zeros(mn - length(bes), 1)];
else
    itric = round(3.6*n^0.188);
    % Heuristics to switch between Bessel, extrapolation and Airy initial guesses.
    igatt = round(mn + 1.31*n^0.4 - n);
    
    bes = besselroots(alpha, itric).^2; % [Tricomi 1947 pg. 296]
    bes = bes/(4*n + 2*alpha+2).*(1 + (bes + 2*(alpha^2 - 1) )/(4*n + 2*alpha+2)^2/3 );
    
    ak = [-13.69148903521072; -12.828776752865757; -11.93601556323626;...
        -11.00852430373326; -10.04017434155809; -9.02265085340981; -7.944133587120853;...
        -6.786708090071759; -5.520559828095551; -4.08794944413097; -2.338107410459767];
    t = 3*pi/2*( (igatt:-1:12)'-0.25); % [DLMF (9.9.6)]
    ak = [-t.^(2/3).*(1 + 5/48./t.^2 - 5/36./t.^4 + 77125/82944./t.^6 ...
        -10856875/6967296./t.^8); ak(max(1,12-igatt):end)];
    nu = 4*n+2*alpha+2; % [Gatteshi 2002 (4.9)]
    air = (nu+ak*(4*nu)^(1/3)+ ak.^2*(nu/16)^(-1/3)/5 + (11/35-alpha^2-12/175*...
        ak.^3)/nu + (16/1575*ak+92/7875*ak.^4)*2^(2/3)*nu^(-5/3) -(15152/3031875*...
        ak.^5+1088/121275*ak.^2)*2^(1/3)*nu^(-7/3));
    
    x = [ bes; zeros(mn - itric -max(igatt,0), 1) ; air];
    fact = zeros(2,1);
    for k = 1:2
        a = alpha + k - 1;
        fact(k) = fact(k) + (1/3840*a^10 - 5/2304*a^9 + 11/2304*a^8 + 7/1920*a^7 ...
            - 229/11520*a^6 + 107/34560*a^5 + 2653/103680*a^4 - 989/155520*a^3 - ...
            3481/311040*a^2 + 139/103680*a + 9871/6531840)/(n - k + 1)^5;
        fact(k) = fact(k) + (1/384*a^8 - 1/96*a^7 + 1/576*a^6 + 43/1440*a^5 - 5/384*...
            a^4 - 23/864*a^3 + 163/25920*a^2 + 31/6480*a - 139/155520)/(n - k + 1)^4;
        fact(k) = fact(k) + (1/48*a^6 - 1/48*a^5 - 1/24*a^4 + 5/144*a^3 + ...
            1/36*a^2 - 1/144*a - 31/6480)/(n - k + 1)^3;
        fact(k) = fact(k) + (1/8*a^4 + 1/12*a^3 - 1/24*a^2 + 1/72)/(n - k + 1)^2;
        fact(k) = fact(k) + (1/2*a^2 + 1/2*a + 1/6)/(n - k + 1)^1 +1;
    end
    % We factored out some constants from the ratio or product of the asymptotics.
    factorx = sqrt(fact(1)/fact(2) )/2/(1 - 1/n)^(1+alpha/2);
    factorw =  -(1 - 1/(n + 1) )^(n + 1+ alpha/2)*(1 - 1/n)^(1 + alpha/2)*...
        exp(1 + 2*log(2) )*4^(1+alpha)*pi*n^alpha*sqrt(prod(fact))*...
        (1 + 1/n)^(alpha/2);
    
    % This is a heuristic for the number of terms in the expansions that follow.
    T = ceil(34/log(n) );
    if ( alpha^2/n > 1 )
        warning('CHEBFUN:lagpts:inputs',['A large alpha may lead to inaccurate ' ...
            'results because the weight is low and R(z) is not close to identity.']);
    end
end

noUnderflow = 1;
for k = 1:mn
    if ( x(k) == 0 ) % Use sextic extrapolation for the initial guesses.
        x(k) = 7*x(k-1)-21*x(k-2)+35*x(k-3)-35*x(k-4)+21*x(k-5)-7*x(k-6)+x(k-7);
    end
    step = x(k);
    l = 0; % Newton-Raphson iteration number
    ov = inf; % Previous/old value
    ox = x(k); % Old x
    % The constants multiplying eps or the maximal l can be increased to avoid
    % the error below, at the cost of possibly lower accuracy.
    while ( ( abs(step) > eps*40*x(k) ) && ( l < 5) )
        l = l + 1;
        if rh
            pe = polyAsyRH(n, x(k), alpha, T);
            % poly' = (p*exp(-Q/2) )' = exp(-Q/2)*(p' -p/2) with orthonormal p
            step = pe/(polyAsyRH(n-1, x(k), alpha+1, T)*factorx - pe/2);
        else
            [pe, pder] = lagpnRecDer(n, alpha, x(k));
            step = pe/pder;
        end
        if (abs(pe) >= abs(ov)*(1-50*eps) )
            % The function values do not decrease enough any more due to
            % roundoff errors, so set to the previous value and quit.
            x(k) = ox;
            break
        end
        ox = x(k);
        x(k) = x(k) -step;
        ov = pe;
    end
    if ( x(k) < 0 ) || ( x(k) > 4*n + 2*alpha + 2 ) || ( l == 5 ) || ...
            ( ( k ~= 1 ) && ( x(k - 1) >= x(k) ) )
        if (rh == 2) && (alpha == 0)
            [x, ders] = alg0_Lag(n);
            w = transpose(exp(-x)./(x.*ders.^2));
            return
        elseif rh == 2
           [x,w] = laguerreExp(n, 0, alpha); 
           return
        end
        error('CHEBFUN:lagpts:converge','Newton method may not have converged.');
    end
    if noUnderflow && rh
        w(k) = factorw/polyAsyRH(n-1, x(k), alpha + 1, T)/...
            polyAsyRH(n+1, x(k), alpha, T)/exp(x(k));
    elseif noUnderflow
        w(k) = (n^2 +alpha*n)^(-1/2)/lagpnRecDer(n-1, alpha, x(k))/pder;
    end
    if noUnderflow && ( w(k) == 0 ) && ( k > 1 ) && ( w(k - 1) > 0 )
        if compRepr
            w = w(1:k-1);
            x = x(1:k-1);
            return;
        end
        noUnderflow = 0;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% Routines for RH algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute the expansion of the orthonormal polynomial without e^(x/2) nor a
% constant factor based on some heuristics

function p = polyAsyRH(np, y, alpha, T)
% We could avoid these tests by splitting the loop k=1:mn into three parts
% with heuristics for the bounding indices. np + alpha is always n, so we
% do not risk dividing by the derivative of an expansion in another region.
if y < sqrt(np+alpha)
    % The fixed delta in the RHP would mean this bound has to be
    % proportional to n, but x(1:k) are O(1/n) so choose the bound in between
    % them to make more use of the (cheap) expansion in the bulk.
    p = asyBessel(np, y, alpha, T);
    return
elseif y > 3.7*(np+alpha)
    % Use the expansion in terms of the (expensive) Airy function, although
    % the corresponding weights will start underflowing for n >= 186.
    p = asyAiry(np, y, alpha, T);
    return
end
p = asyBulk(np, y, alpha, T);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the expansion of the orthonormal polynomial in the bulk without e^(x/2)
function p = asyBulk(np, y, alpha, T)
z = y/4/np;
mnxi = 2*np*( sqrt(z).*sqrt(1 - z) - acos(sqrt(z) ) ); % = -n*xin/i
d = z - 1;
if T == 1
    p = real( 2/z^(1/4 + alpha/2)/(1 - z)^(1/4)*cos(acos(2*z - 1)*...
        (1/2 + alpha/2) - mnxi - pi/4) );
    return;
end
% Getting the higher order terms is hard-coded for speed and code length, but
% can be made to get arbitrary orders for even more general weight functions.
R1 = 0.0;
R2 = 0.0;
if (alpha == 0)
    if ( T >= 7 )
        R1 = R1 + (+6.975493821675553e-05*z^(-1) +2.308068495174802e-05*z^(-2) )/np^6;
        R2 = R2 + (-8.294589753829749e-05*z^(-1) -0.0001387411504294034*z^(-2) )/np^6;
        R1 = R1 + (+3.465393092483285e-06*z^(-3) )/np^6;
        R2 = R2 + (-0.0001432362478226424*z^(-3) )/np^6;
        R1 = R1 + (-6.975493821675553e-05*d^(-1) +2.944141339480553e-05*d^(-2) )/np^6;
        R2 = R2 + (+1.401453805748335e-05*d^(-1) -4.189865305987762e-06*d^(-2) )/np^6;
        R1 = R1 + (-3.343882101991324e-05*d^(-3) +2.53485975588054e-05*d^(-4) )/np^6;
        R2 = R2 + (-1.993525079283811e-06*d^(-3) -9.103419242066082e-06*d^(-4) )/np^6;
        R1 = R1 + (+1.652740321445771e-05*d^(-5) -8.823225788038024e-05*d^(-6) )/np^6;
        R2 = R2 + (+0.0001048908434950015*d^(-5) -0.005941442455789839*d^(-6) )/np^6;
        R1 = R1 + (-0.001224011071160469*d^(-7) -0.002439398363695296*d^(-8) )/np^6;
        R2 = R2 + (-0.04919568821745263*d^(-7) -0.09585400629108581*d^(-8) )/np^6;
        R1 = R1 + (-0.001482771554403024*d^(-9) )/np^6;
        R2 = R2 + (-0.05337977595850887*d^(-9) )/np^6;
    end
    if ( T >= 6 )
        R1 = R1 + (-6.169113234727972e-05*z^(-1) -0.0001287894944349925*z^(-2) )/np^5;
        R2 = R2 + (-0.0002803040936650563*z^(-1) -0.000105773409207662*z^(-2) )/np^5;
        R1 = R1 + (-0.0001108925789594651*z^(-3) )/np^5;
        R2 = R2 + (+0.0001108925789594651*z^(-3) )/np^5;
        R1 = R1 + (+6.169113234728014e-05*d^(-1) -2.735243584977629e-05*d^(-2) )/np^5;
        R2 = R2 + (-9.749909808489811e-05*d^(-1) +7.023189392787169e-05*d^(-2) )/np^5;
        R1 = R1 + (+2.342076600182224e-05*d^(-3) -7.31386404171321e-05*d^(-4) )/np^5;
        R2 = R2 + (-6.224353508098193e-05*d^(-3) +0.0001481373301066077*d^(-4) )/np^5;
        R1 = R1 + (+0.0001987962105882525*d^(-5) -0.001539873713001479*d^(-6) )/np^5;
        R2 = R2 + (-0.004849931367753502*d^(-5) -0.03068677518409481*d^(-6) )/np^5;
        R1 = R1 + (-0.01282539667078748*d^(-7) -0.01377542605380878*d^(-8) )/np^5;
        R2 = R2 + (-0.04227630754444762*d^(-7) -0.01377542605380878*d^(-8) )/np^5;
    end
    if ( T >= 5 )
        R1 = R1 + (-0.0001123610837959949*z^(-1) -1.490116119384768e-05*z^(-2) )/np^4;
        R2 = R2 + (+0.0002556506498360341*z^(-1) +0.000452995300292969*z^(-2) )/np^4;
        R1 = R1 + (+0.0001123610837959949*d^(-1) -4.159894009185921e-05*d^(-2) )/np^4;
        R2 = R2 + (-3.220671979488066e-05*d^(-1) -1.006970189726202e-05*d^(-2) )/np^4;
        R1 = R1 + (+0.0001014470072930732*d^(-3) -0.0001985441019505633*d^(-4) )/np^4;
        R2 = R2 + (+0.000150605189947433*d^(-3) -0.004846322487411191*d^(-4) )/np^4;
        R1 = R1 + (-0.0008181122595390668*d^(-5) -0.000793068006696034*d^(-6) )/np^4;
        R2 = R2 + (-0.02270678924434961*d^(-5) -0.01903363216070482*d^(-6) )/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (+0.0002585517035590279*z^(-1) +0.0005722045898437501*z^(-2) )/np^3;
        R2 = R2 + (+0.0009774102105034725*z^(-1) -0.0005722045898437501*z^(-2) )/np^3;
        R1 = R1 + (-0.0002585517035590275*d^(-1) -1.465597270447586e-05*d^(-2) )/np^3;
        R2 = R2 + (+0.000218577443817516*d^(-1) +8.356541763117039e-05*d^(-2) )/np^3;
        R1 = R1 + (+0.0003698466736593355*d^(-3) -0.002262821903935187*d^(-4) )/np^3;
        R2 = R2 + (-0.006447629575376164*d^(-3) -0.02017682864342207*d^(-4) )/np^3;
        R1 = R1 + (-0.008014160909770449*d^(-5) )/np^3;
        R2 = R2 + (-0.008014160909770449*d^(-5) )/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (+0.0001627604166666672*z^(-1) )/np^2;
        R2 = R2 + (-0.004557291666666669*z^(-1) )/np^2;
        R1 = R1 + (-0.0001627604166666668*d^(-1) -0.0007052951388888891*d^(-2) )/np^2;
        R2 = R2 + (+0.001085069444444443*d^(-1) -0.01323784722222223*d^(-2) )/np^2;
        R1 = R1 + (-0.001898871527777778*d^(-3) )/np^2;
        R2 = R2 + (-0.02278645833333334*d^(-3) )/np^2;
    end
    R1 = R1 + (-0.015625*z^(-1) )/np^1;
    R2 = R2 + (+0.015625*z^(-1) )/np^1;
    R1 = R1 + (+0.015625*d^(-1) -0.02604166666666667*d^(-2) )/np^1 + 1;
    R2 = R2 + (-0.05729166666666667*d^(-1) -0.02604166666666667*d^(-2) )/np^1;
elseif (alpha == 1)
    if ( T >= 7 )
        R1 = R1 + (+4.421050825524612e-05*z^(-1) -9.845300465477299e-05*z^(-2) )/np^6;
        R2 = R2 + (-0.0004130488626772792*z^(-1) -0.0002886747975868209*z^(-2) )/np^6;
        R1 = R1 + (-0.0001284762402065099*z^(-3) )/np^6;
        R2 = R2 + (+0.0002936599776148798*z^(-3) )/np^6;
        R1 = R1 + (-4.42105082552552e-05*d^(-1) +3.097987510232154e-05*d^(-2) )/np^6;
        R2 = R2 + (-3.368568855349101e-05*d^(-1) +2.181357967052824e-05*d^(-2) )/np^6;
        R1 = R1 + (-2.025154726419684e-05*d^(-3) +1.804009902887319e-06*d^(-4) )/np^6;
        R2 = R2 + (-5.460950725967714e-05*d^(-3) +0.0004144046920969578*d^(-4) )/np^6;
        R1 = R1 + (+8.838002681793751e-05*d^(-5) -0.005903692556356525*d^(-6) )/np^6;
        R2 = R2 + (-0.02360578463834522*d^(-5) -0.2085459891532966*d^(-6) )/np^6;
        R1 = R1 + (-0.04927174711746236*d^(-7) -0.09714545248363038*d^(-8) )/np^6;
        R2 = R2 + (-0.4808177876586013*d^(-7) -0.4040791642450564*d^(-8) )/np^6;
        R1 = R1 + (-0.05486254751291189*d^(-9) )/np^6;
        R2 = R2 + (-0.1067595519170177*d^(-9) )/np^6;
    end
    if ( T >= 6 )
        R1 = R1 + (-0.0001335134292826253*z^(-1) -5.663062135378524e-05*z^(-2) )/np^5;
        R2 = R2 + (+1.53532918588624e-05*z^(-1) +0.0007243330279986067*z^(-2) )/np^5;
        R1 = R1 + (+0.0001355353742837906*z^(-3) )/np^5;
        R2 = R2 + (-0.0001355353742837906*z^(-3) )/np^5;
        R1 = R1 + (+0.0001335134292826214*d^(-1) -0.0001154726233560362*d^(-2) )/np^5;
        R2 = R2 + (-0.0001697125535676412*d^(-1) +7.79156868957573e-05*d^(-2) )/np^5;
        R1 = R1 + (+0.0001007355757637324*d^(-3) +1.598066688109512e-05*d^(-4) )/np^5;
        R2 = R2 + (+0.0003025175562278589*d^(-3) -0.01889981810881747*d^(-4) )/np^5;
        R1 = R1 + (-0.004761903800468913*d^(-5) -0.03071486300933196*d^(-6) )/np^5;
        R2 = R2 + (-0.1327588315026021*d^(-5) -0.2290429970088575*d^(-6) )/np^5;
        R1 = R1 + (-0.04132627816142633*d^(-7) -0.01377542605380878*d^(-8) )/np^5;
        R2 = R2 + (-0.1277789520163642*d^(-7) -0.01377542605380878*d^(-8) )/np^5;
    end
    if ( T >= 5 )
        R1 = R1 + (-4.915484675654808e-05*z^(-1) +0.0003713369369506837*z^(-2) )/np^4;
        R2 = R2 + (+0.001377543696650754*z^(-1) -0.0009346008300781254*z^(-2) )/np^4;
        R1 = R1 + (+4.915484675654708e-05*d^(-1) -6.56338876166932e-05*d^(-2) )/np^4;
        R2 = R2 + (+4.188788771141544e-05*d^(-1) +0.0004643114505971865*d^(-2) )/np^4;
        R1 = R1 + (+0.000152441503579723*d^(-3) -0.004777727892369407*d^(-4) )/np^4;
        R2 = R2 + (-0.01907248320402925*d^(-3) -0.100409803665224*d^(-4) )/np^4;
        R1 = R1 + (-0.0228570547614078*d^(-5) -0.01982670016740085*d^(-6) )/np^4;
        R2 = R2 + (-0.1208802603890376*d^(-5) -0.03806726432140963*d^(-6) )/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (+0.0005976359049479174*z^(-1) -0.0008010864257812502*z^(-2) )/np^3;
        R2 = R2 + (-0.003878275553385417*z^(-1) +0.0008010864257812502*z^(-2) )/np^3;
        R1 = R1 + (-0.0005976359049479167*d^(-1) +0.0008296636887538567*d^(-2) )/np^3;
        R2 = R2 + (+0.001602040985484176*d^(-1) -0.02611068913966053*d^(-2) )/np^3;
        R1 = R1 + (-0.007273111225646224*d^(-3) -0.01923398618344908*d^(-4) )/np^3;
        R2 = R2 + (-0.09346341733579284*d^(-3) -0.07109032148196376*d^(-4) )/np^3;
        R1 = R1 + (-0.008014160909770449*d^(-5) )/np^3;
        R2 = R2 + (-0.008014160909770449*d^(-5) )/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (-0.00048828125*z^(-1) )/np^2;
        R2 = R2 + (+0.007812500000000007*z^(-1) )/np^2;
        R1 = R1 + (+0.0004882812499999989*d^(-1) -0.01177300347222222*d^(-2) )/np^2;
        R2 = R2 + (-0.0529513888888889*d^(-1) -0.1154513888888889*d^(-2) )/np^2;
        R1 = R1 + (-0.02468532986111112*d^(-3) )/np^2;
        R2 = R2 + (-0.04557291666666667*d^(-3) )/np^2;
    end
    R1 = R1 + (+0.046875*z^(-1) )/np^1;
    R2 = R2 + (-0.046875*z^(-1) )/np^1;
    R1 = R1 + (-0.046875*d^(-1) -0.02604166666666667*d^(-2) )/np^1 + 1;
    R2 = R2 + (-0.2447916666666667*d^(-1) -0.02604166666666667*d^(-2) )/np^1;
else
    if ( T >= 5 )
        R1 = R1 + (5/49152*alpha^8 - 35/49152*alpha^7 + 67/49152*alpha^6 + ...
            11/36864*alpha^5 - 1529/589824*alpha^4 + 1891/2359296*alpha^3 + ...
            26827/26542080*alpha^2 - 109/524288*alpha - 190867/1698693120)*z^(-1)/np^4;
        R2 = R2 + (-1/3072*alpha^8 + 89/49152*alpha^7 - 421/294912*alpha^6 - ...
            1465/294912*alpha^5 + 2525/589824*alpha^4 + 10907/2359296*alpha^3 - ...
            424361/212336640*alpha^2 - 229/262144*alpha + 13571/53084160)*z^(-1)/np^4;
        R1 = R1 + (-1/32768*alpha^8 + 5/49152*alpha^7 + 5/18432*alpha^6 - 175/...
            196608*alpha^5 - 1253/2359296*alpha^4 + 1295/786432*alpha^3 + ...
            415/2359296*alpha^2 - 375/1048576*alpha - 125/8388608)*z^(-2)/np^4;
        R2 = R2 + (1/24576*alpha^8 - 5/49152*alpha^7 - 143/294912*alpha^6 + ...
            175/196608*alpha^5 + 2107/1179648*alpha^4 - 1295/786432*alpha^3 - ...
            10517/4718592*alpha^2 + 375/1048576*alpha + 475/1048576)*z^(-2)/np^4;
        R1 = R1 + (-5/49152*alpha^8 + 35/49152*alpha^7 - 67/49152*alpha^6 - ...
            11/36864*alpha^5 + 1529/589824*alpha^4 - 1891/2359296*alpha^3 - ...
            26827/26542080*alpha^2 + 109/524288*alpha + 190867/1698693120)*d^(-1)/np^4;
        R2 = R2 + (-1/3072*alpha^8 + 13/16384*alpha^7 + 293/294912*alpha^6 - ...
            3683/1474560*alpha^5 - 605/589824*alpha^4 + 14383/7077888*alpha^3 + ...
            30179/70778880*alpha^2 - 34231/106168320*alpha - 5129/159252480)*d^(-1)/np^4;
        R1 = R1 + (-1/32768*alpha^8 - 1/6144*alpha^7 + 145/147456*alpha^6 - ...
            1999/2949120*alpha^5 - 327/262144*alpha^4 + 2897/3538944*alpha^3 + ...
            3127/7077888*alpha^2 - 63391/424673280*alpha - 423983/10192158720)*d^(-2)/np^4;
        R2 = R2 + (-1/24576*alpha^8 - 17/24576*alpha^7 - 107/147456*alpha^6 + ...
            1751/983040*alpha^5 + 1535/1179648*alpha^4 - 853/884736*alpha^3 - ...
            7583/26542080*alpha^2 + 14131/141557760*alpha - 12829/1274019840)*d^(-2)/np^4;
        R1 = R1 + (-1/16384*alpha^7 - 35/49152*alpha^6 + 1127/1474560*alpha^5 + ...
            973/1179648*alpha^4 - 1051/7077888*alpha^3 - 2093/3538944*alpha^2 - ...
            533/21233664*alpha + 258491/2548039680)*d^(-3)/np^4;
        R2 = R2 + (-1/16384*alpha^7 - 49/32768*alpha^6 - 4249/737280*alpha^5 - ...
            2009/294912*alpha^4 - 9077/2359296*alpha^3 - 94099/70778880*alpha^2 + ...
            18133/212336640*alpha + 31979/212336640)*d^(-3)/np^4;
        R1 = R1 + (-1001/589824*alpha^5 - 5005/2359296*alpha^4 - 1001/589824*...
            alpha^3 + 5005/7077888*alpha^2 + 32461/141557760*alpha - ...
            674531/3397386240)*d^(-4)/np^4;
        R2 = R2 + (-1001/589824*alpha^5 - 5005/393216*alpha^4 - 101101/...
            3538944*alpha^3 - 695695/21233664*alpha^2 - 8411117/424673280*...
            alpha - 6174311/1274019840)*d^(-4)/np^4;
        R1 = R1 + (-85085/7077888*alpha^3 - 425425/42467328*alpha - 833833/1019215872)*d^(-5)/np^4;
        R2 = R2 + (-85085/7077888*alpha^3 - 1616615/42467328*alpha^2 - ...
            85085/1769472*alpha - 1446445/63700992)*d^(-5)/np^4;
        R1 = R1 + (-1616615/84934656*alpha - 1616615/2038431744)*d^(-6)/np^4;
        R2 = R2 + (-1616615/84934656*alpha - 1616615/84934656)*d^(-6)/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (1/1024*alpha^6 - 5/1536*alpha^5 + 31/24576*alpha^4 + 7/2048*...
            alpha^3 - 13/9216*alpha^2 - 1/1536*alpha + 305/1179648)*z^(-1)/np^3;
        R2 = R2 + (-1/384*alpha^6 + 19/3072*alpha^5 + 41/8192*alpha^4 - 11/1024*...
            alpha^3 - 737/147456*alpha^2 + 113/49152*alpha + 1153/1179648)*z^(-1)/np^3;
        R1 = R1 + (-1/6144*alpha^6 + 35/24576*alpha^4 - 259/98304*alpha^2 + ...
            75/131072)*z^(-2)/np^3;
        R2 = R2 + (1/6144*alpha^6 - 35/24576*alpha^4 + 259/98304*alpha^2 ...
            - 75/131072)*z^(-2)/np^3;
        R1 = R1 + (-1/1024*alpha^6 + 5/1536*alpha^5 - 31/24576*alpha^4 - ...
            7/2048*alpha^3 + 13/9216*alpha^2 + 1/1536*alpha - 305/1179648)*d^(-1)/np^3;
        R2 = R2 + (-1/384*alpha^6 - 1/1024*alpha^5 + 133/24576*alpha^4 + ...
            19/9216*alpha^3 - 287/147456*alpha^2 - 83/147456*alpha + 11603/53084160)*d^(-1)/np^3;
        R1 = R1 + (-1/6144*alpha^6 - 1/512*alpha^5 + 5/2048*alpha^4 + 23/18432*...
            alpha^3 - 151/294912*alpha^2 - 1/4608*alpha - 389/26542080)*d^(-2)/np^3;
        R2 = R2 + (-1/6144*alpha^6 - 13/3072*alpha^5 - 145/12288*alpha^4 - 83/9216*...
            alpha^3 - 51/32768*alpha^2 + 83/147456*alpha + 1109/13271040)*d^(-2)/np^3;
        R1 = R1 + (-35/8192*alpha^4 - 35/18432*alpha^3 - 217/147456*alpha^2 + ...
            19633/53084160)*d^(-3)/np^3;
        R2 = R2 + (-35/8192*alpha^4 - 35/1536*alpha^3 - 2611/73728*alpha^2 - ...
            3619/147456*alpha - 114089/17694720)*d^(-3)/np^3;
        R1 = R1 + (-5005/294912*alpha^2 - 1001/442368)*d^(-4)/np^3;
        R2 = R2 + (-5005/294912*alpha^2 - 5005/147456*alpha - 107107/5308416)*d^(-4)/np^3;
        R1 = R1 + (-85085/10616832)*d^(-5)/np^3;
        R2 = R2 + (-85085/10616832)*d^(-5)/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (1/128*alpha^4 - 1/128*alpha^3 - 1/384*alpha^2 + 1/512*alpha + 1/6144)*z^(-1)/np^2;
        R2 = R2 + (-1/64*alpha^4 + 1/128*alpha^3 + 17/768*alpha^2 - 1/512*alpha - 7/1536)*z^(-1)/np^2;
        R1 = R1 + (-1/128*alpha^4 + 1/128*alpha^3 + 1/384*alpha^2 - 1/512*alpha - 1/6144)*d^(-1)/np^2;
        R2 = R2 + (-1/64*alpha^4 - 11/384*alpha^3 - 3/256*alpha^2 + 1/512*alpha + 5/4608)*d^(-1)/np^2;
        R1 = R1 + (-5/384*alpha^3 + 1/512*alpha - 13/18432)*d^(-2)/np^2;
        R2 = R2 + (-5/384*alpha^3 - 35/768*alpha^2 - 67/1536*alpha - 61/4608)*d^(-2)/np^2;
        R1 = R1 + (-35/1536*alpha - 35/18432)*d^(-3)/np^2;
        R2 = R2 + (-35/1536*alpha - 35/1536)*d^(-3)/np^2;
    end
    R1 = R1 + (1/16*alpha^2 - 1/64)*z^(-1)/np^1;
    R2 = R2 + (-1/16*alpha^2 + 1/64)*z^(-1)/np^1;
    R1 = R1 + (-1/16*alpha^2 + 1/64)*d^(-1)/np^1;
    R2 = R2 + (-1/16*alpha^2 - 1/8*alpha - 11/192)*d^(-1)/np^1;
    R1 = R1 + (-5/192)*d^(-2)/np^1 + 1;
    R2 = R2 + (-5/192)*d^(-2)/np^1;
end
p = real( 2/z^(1/4 + alpha/2)*(cos(acos(2*z-1)*(1/2+alpha/2) - mnxi-pi/4)*R1 ...
    -cos(acos(2*z-1)*(-1/2+alpha/2)-mnxi-pi/4)*R2)/(1 - z)^(1/4) );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the expansion of the orthonormal polynomial near zero without e^(x/2)
function p = asyBessel(np, y, alpha, T)
z = y/4/np;
npb = 2*np*(pi/2 + sqrt(z).*sqrt(1 - z) - acos(sqrt(z) ) ); % = 2i n sqrtphitn

if ( T == 1 )
    p = real( sqrt(2*pi)*(-1)^np*sqrt(npb)/z^(1/4)/(1 - z)^(1/4)*...
        z^(-alpha/2)*(sin( (alpha + 1)/2*acos(2*z - 1) - pi*alpha/2)*...
        besselj(alpha,npb) + cos( (alpha + 1)/2*acos(2*z - 1) - ...
        pi*alpha/2)*(besselj(alpha-1,npb) - alpha/(npb)*...
        besselj(alpha, npb) ) ) );
    return
end
% Use the series expansion of R because it is faster and we use pl only very
% close to zero to have less calls to besselj.
R1 = 0.0;
R2 = 0.0;
if ( alpha == 0 )
    if ( T >= 7 )
        R1 = R1 + (-0.01196102063075393*z^2 -0.00277571228121701*z^1 )/np^6;
        R2 = R2 + (+0.1904949571852236*z^2 +0.01937961523964843*z^1 )/np^6;
        R1 = R1 + (-0.0003486406879182943 )/np^6;
        R2 = R2 + (+0.0004183688255019682 )/np^6;
    end
    if ( T >= 6 )
        R1 = R1 + (-0.6637971470977652*z^3 -0.1730340465711127*z^2 )/np^5;
        R2 = R2 + (+0.3913221799067465*z^3 +0.1312375428294408*z^2 )/np^5;
        R1 = R1 + (-0.03190107388066272*z^1 -0.003136156886880289 )/np^5;
        R2 = R2 + (+0.03048894498589488*z^1 +0.003920196108600354 )/np^5;
    end
    if ( T >= 5 )
        R1 = R1 + (+0.001501174033247706*z^4 +0.006473841918350666*z^3 )/np^4;
        R2 = R2 + (-0.9524842232974917*z^4 -0.3550440568543063*z^3 )/np^4;
        R1 = R1 + (+0.005174830777647138*z^2 +0.002518294998040369*z^1 )/np^4;
        R2 = R2 + (-0.1020198851574237*z^2 -0.01827557013031548*z^1 )/np^4;
        R1 = R1 + (+0.0006884162808641975 )/np^4;
        R2 = R2 + (-0.0009178883744855927 )/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (+0.8486152325952079*z^5 +0.4599905388515398*z^4 )/np^3;
        R2 = R2 + (-0.02362629629145646*z^5 -0.07442783509439153*z^4 )/np^3;
        R1 = R1 + (+0.2227888285411434*z^3 +0.0915386169900059*z^2 )/np^3;
        R2 = R2 + (-0.07472004547236026*z^3 -0.05208755878894767*z^2 )/np^3;
        R1 = R1 + (+0.02883322310405644*z^1 +0.005362654320987655 )/np^3;
        R2 = R2 + (-0.02605544532627866*z^1 -0.008043981481481482 )/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (+0.03274424783742833*z^6 +0.02138984808385162*z^5 )/np^2;
        R2 = R2 + (+0.5345928880817182*z^6 +0.3890960516718894*z^5 )/np^2;
        R1 = R1 + (+0.01205177881103808*z^4 +0.004772192827748389*z^3 )/np^2;
        R2 = R2 + (+0.2664715497122905*z^4 +0.1667621987066432*z^3 )/np^2;
        R1 = R1 + (-0.0003747795414462069*z^2 -0.003240740740740739*z^1 )/np^2;
        R2 = R2 + (+0.09005731922398588*z^2 +0.03657407407407406*z^1 )/np^2;
        R1 = R1 + (-0.003472222222222222 )/np^2;
        R2 = R2 + (+0.006944444444444446 )/np^2;
    end
    R1 = R1 + (-0.2113083089153498*z^7 -0.1842728591546934*z^6 )/np^1;
    R2 = R2 + (-0.1378214216323702*z^7 -0.1106475197021934*z^6 )/np^1;
    R1 = R1 + (-0.1569653787325745*z^5 -0.129241088129977*z^4 )/np^1;
    R2 = R2 + (-0.08313723470337227*z^5 -0.05509753620864732*z^4 )/np^1;
    R1 = R1 + (-0.1008289241622575*z^3 -0.07116402116402117*z^2 )/np^1;
    R2 = R2 + (-0.02615520282186949*z^3 +0.0044973544973545*z^2 )/np^1;
    R1 = R1 + (-0.03888888888888889*z^1 -3.469446951953614e-18 )/np^1 + 1;
    R2 = R2 + (+0.0388888888888889*z^1 +0.08333333333333333 )/np^1;
elseif ( alpha == 1 )
    if ( T >= 7 )
        R1 = R1 + (+0.1888634820675584*z^2 +0.01806148723662363*z^1 )/np^6;
        R2 = R2 + (-0.2473524833501428*z^2 -0.02552670106483584*z^1 )/np^6;
        R1 = R1 + (+6.972813758367975e-05 )/np^6;
        R2 = R2 + (+0.0003920196108599757 )/np^6;
    end
    if ( T >= 6 )
        R1 = R1 + (+0.2500575428553918*z^3 +0.07916956633223091*z^2 )/np^5;
        R2 = R2 + (+0.381735599523504*z^3 +0.03570599642395061*z^2 )/np^5;
        R1 = R1 + (+0.01519177406366366*z^1 +0.0007840392217200875 )/np^5;
        R2 = R2 + (-0.006645317399079066*z^1 -0.0001147360468107504 )/np^5;
    end
    if ( T >= 5 )
        R1 = R1 + (-0.9842194497752994*z^4 -0.3638886994207212*z^3 )/np^4;
        R2 = R2 + (+0.4445716141609543*z^4 +0.2691722827926532*z^3 )/np^4;
        R1 = R1 + (-0.1025708529296493*z^2 -0.01716940402704293*z^1 )/np^4;
        R2 = R2 + (+0.1099354240152853*z^2 +0.02321226484420927*z^1 )/np^4;
        R1 = R1 + (-0.0002294720936213962 )/np^4;
        R2 = R2 + (-0.001340663580246942 )/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (+0.07190356645934227*z^5 -0.008994495608252214*z^4 )/np^3;
        R2 = R2 + (-1.04475270413498*z^5 -0.5663204149612885*z^4 )/np^3;
        R1 = R1 + (-0.03246736545347658*z^3 -0.0269951499118166*z^2 )/np^3;
        R2 = R2 + (-0.2509383517716854*z^3 -0.07455357142857155*z^2 )/np^3;
        R1 = R1 + (-0.01297949735449737*z^1 -0.002681327160493828 )/np^3;
        R2 = R2 + (-0.004497354497354576*z^1 +0.001736111111111079 )/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (+0.5818842080253015*z^6 +0.422774248874778*z^5 )/np^2;
        R2 = R2 + (+0.4633747833589102*z^6 +0.2662342736628452*z^5 )/np^2;
        R1 = R1 + (+0.2885276040831597*z^4 +0.1792151675485009*z^3 )/np^2;
        R2 = R2 + (+0.1157407407407409*z^4 +0.01269841269841278*z^3 )/np^2;
        R1 = R1 + (+0.09497354497354499*z^2 +0.03611111111111111*z^1 )/np^2;
        R2 = R2 + (-0.04087301587301589*z^2 -0.03888888888888888*z^1 )/np^2;
        R1 = R1 + (+0.003472222222222224 )/np^2;
        R2 = R2 + (+0.04166666666666669 )/np^2;
    end
    R1 = R1 + (-0.147667867682724*z^7 -0.1203555135830268*z^6 )/np^1;
    R2 = R2 + (+0.04286875484833783*z^7 +0.06814331582585553*z^6 )/np^1;
    R1 = R1 + (-0.09264242400750337*z^5 -0.06428731762065096*z^4 )/np^1;
    R2 = R2 + (+0.0927247285342524*z^5 +0.1158922558922559*z^4 )/np^1;
    R1 = R1 + (-0.03481481481481481*z^3 -0.003174603174603177*z^2 )/np^1;
    R2 = R2 + (+0.1358730158730159*z^3 +0.1476190476190476*z^2 )/np^1;
    R1 = R1 + (+0.03333333333333333*z^1 +0.08333333333333333 )/np^1 + 1;
    R2 = R2 + (+0.1333333333333333*z^1 +0.5000000000000001 )/np^1;
else
    if ( T >= 5 )
        R1 = R1 + (-1/70761600*alpha^17 + 1/9580032*alpha^16 + 8209/6793113600*alpha^15 - 12557/1277337600*alpha^14 - 4529803/149448499200*alpha^13 + 32533/106444800*alpha^12 + 322620017/1207084032000*alpha^11 - 14482007/3359232000*alpha^10 - 16413757/36578304000*alpha^9 + 2685618697/86220288000*alpha^8 -...
            1951839107/1207084032000*alpha^7 - 78927123169/747242496000*alpha^6 - 21286749097/196151155200*alpha^5 + 328901198447/1810626048000*alpha^4 + 8020862411/8915961600*alpha^3 - 7637940403/33949238400*alpha^2 - 1621342262117/980755776000*alpha + 2650113187/1765360396800)*z^4/np^4;
        R2 = R2 + (1/70761600*alpha^17 + 1/19160064*alpha^16 - 15797/10674892800*alpha^15 - 4901/1277337600*alpha^14 + 8666507/149448499200*alpha^13 + 565729/5748019200*alpha^12 - 1337582357/1207084032000*alpha^11 - 1630981/1469664000*alpha^10 + 3726654553/329204736000*alpha^9 +...
            474130753/86220288000*alpha^8 - 81309310193/1207084032000*alpha^7 + 70322427437/3138418483200*alpha^6 + 882127042519/3923023104000*alpha^5 - 2530098562949/5884534656000*alpha^4 - 462027031849/1176906931200*alpha^3 + 14073935210333/8826801984000*alpha^2 + 1266579893857/2942267328000*alpha -...
            8407389631931/8826801984000)*z^4/np^4;
        R1 = R1 + (1/2534400*alpha^15 - 53/21772800*alpha^14 - 5671/273715200*alpha^13 + 8191/58060800*alpha^12 + 329687/1149603840*alpha^11 - 145877/58060800*alpha^10 - 224891/174182400*alpha^9 + 179041/9331200*alpha^8 + 1620097/522547200*alpha^7 - 208582621/3448811520*alpha^6 -...
            223498109/2874009600*alpha^5 + 240670031/2463436800*alpha^4 + 1343705443/2874009600*alpha^3 - 14957779603/129330432000*alpha^2 - 503372099/718502400*alpha + 209316193/32332608000)*z^3/np^4;
        R2 = R2 + (-1/2534400*alpha^15 - 1/870912*alpha^14 + 50257/1916006400*alpha^13 + 2819/58060800*alpha^12 - 3597151/5748019200*alpha^11 - 16583/24883200*alpha^10 + 1212661/174182400*alpha^9 + 1830587/522547200*alpha^8 - 299023/6967296*alpha^7 + 272776067/17244057600*alpha^6 +...
            26321231/179625600*alpha^5 - 1018402829/4311014400*alpha^4 - 488108053/1724405760*alpha^3 + 92447911573/129330432000*alpha^2 + 646668943/2155507200*alpha - 11479500313/32332608000)*z^3/np^4;
        R1 = R1 + (-1/138240*alpha^13 + 1/27648*alpha^12 + 59/276480*alpha^11 - 6833/5806080*alpha^10 - 128221/87091200*alpha^9 + 184301/17418240*alpha^8 + 19547/4354560*alpha^7 - 2028469/65318400*alpha^6 - 23797/483840*alpha^5 + 437701/9331200*alpha^4 + 368551/1741824*alpha^3 - 10051211/195955200*alpha^2 ...
            - 5386001/21772800*alpha + 202807/39191040)*z^2/np^4;
        R2 = R2 + (1/138240*alpha^13 + 29/1935360*alpha^12 - 79/276480*alpha^11 - 131/414720*alpha^10 + 334501/87091200*alpha^9 + 193/103680*alpha^8 - 11381/435456*alpha^7 + 549599/52254720*alpha^6 + 3977899/43545600*alpha^5 - 14527739/130636800*alpha^4 - 2290789/13063680*alpha^3 +...
            101496679/391910400*alpha^2 + 2072993/13063680*alpha - 19991327/195955200)*z^2/np^4;
        R1 = R1 + (5/64512*alpha^11 - 7/23040*alpha^10 - 1087/967680*alpha^9 + 241/51840*alpha^8 + 1861/483840*alpha^7 - 229841/17418240*alpha^6 - 151051/5806080*alpha^5 + 173/9072*alpha^4 + 29773/387072*alpha^3 - 4726201/261273600*alpha^2 - 19021/290304*alpha + 164491/65318400)*z^1/np^4;
        R2 = R2 + (-5/64512*alpha^11 - 1/11520*alpha^10 + 1633/967680*alpha^9 + 143/207360*alpha^8 - 7153/483840*alpha^7 + 112871/17418240*alpha^6 + 305863/5806080*alpha^5 - 118261/2903040*alpha^4 - 1547069/17418240*alpha^3 + 16784701/261273600*alpha^2 + 262201/4354560*alpha - 170533/9331200)*z^1/np^4;
        R1 = R1 + (-1/2560*alpha^9 + 5/4608*alpha^8 + 1/512*alpha^7 - 13/3840*alpha^6 - 1351/138240*alpha^5 + 769/138240*alpha^4 + 505/27648*alpha^3 - 5201/1244160*alpha^2 - 139/13824*alpha + 571/829440)/np^4;
        R2 = R2 + (1/2560*alpha^9 + 1/4608*alpha^8 - 11/1536*alpha^7 + 13/3840*alpha^6 + 707/27648*alpha^5 - 1399/138240*alpha^4 - 101/3072*alpha^3 + 8417/1244160*alpha^2 + 139/10368*alpha - 571/622080)/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (1/389188800*alpha^17 - 29/8172964800*alpha^16 - 89/233513280*alpha^15 + 464407/653837184000*alpha^14 + 91103/4670265600*alpha^13 - 1016209/27433728000*alpha^12 - 1740359/3772137600*alpha^11 + 21624227/27433728000*alpha^10 + 9561221/1714608000*alpha^9 - 258861697/33530112000*alpha^8 -...
            112625083/3143448000*alpha^7 + 1434417007/43110144000*alpha^6 + 1837963583/17513496000*alpha^5 + 74826019453/980755776000*alpha^4 - 2257518019/14859936000*alpha^3 - 48709416163/53495769600*alpha^2 + 2669768693/24518894400*alpha + 624213218227/735566832000)*z^5/np^3;
        R2 = R2 + (-1/389188800*alpha^17 - 61/2043241200*alpha^16 + 37/166795200*alpha^15 + 24173/8072064000*alpha^14 - 59231/9340531200*alpha^13 - 33351151/301771008000*alpha^12 + 5837/86220288*alpha^11 + 52760503/27433728000*alpha^10 - 49691/326592000*alpha^9 - 190468991/11176704000*alpha^8 -...
            4253279/3592512000*alpha^7 + 7916648611/100590336000*alpha^6 - 16738865651/490377888000*alpha^5 - 20185524839/140107968000*alpha^4 + 48175779259/98075577600*alpha^3 - 122433856319/2942267328000*alpha^2 - 1153342121/851350500*alpha - 3475743983/147113366400)*z^5/np^3;
        R1 = R1 + (-1/9979200*alpha^15 + 13/119750400*alpha^14 + 233/23950080*alpha^13 - 14051/958003200*alpha^12 - 109799/359251200*alpha^11 + 117469/261273600*alpha^10 + 136891/32659200*alpha^9 - 63661/12441600*alpha^8 - 60631/2177280*alpha^7 + 86313929/3919104000*alpha^6 + 6990713/89812800*alpha^5 +...
            1018395439/15324309000*alpha^4 - 6588979/59875200*alpha^3 - 140683954849/245188944000*alpha^2 + 6929479/89812800*alpha + 112784594471/245188944000)*z^4/np^3;
        R2 = R2 + (1/9979200*alpha^15 + 17/17107200*alpha^14 - 629/119750400*alpha^13 - 11833/191600640*alpha^12 + 11993/143700480*alpha^11 + 350363/261273600*alpha^10 - 1009/2419200*alpha^9 - 1128817/87091200*alpha^8 + 527/21772800*alpha^7 + 33932833/559872000*alpha^6 - 22218013/718502400*alpha^5 -...
            115334458321/980755776000*alpha^4 + 62394391/179625600*alpha^3 + 28057787633/490377888000*alpha^2 - 1703363/2138400*alpha - 18248882291/245188944000)*z^4/np^3;
        R1 = R1 + (1/362880*alpha^13 - 1/453600*alpha^12 - 71/435456*alpha^11 + 1361/7257600*alpha^10 + 625/217728*alpha^9 - 21569/7257600*alpha^8 - 1337/64800*alpha^7 + 291101/21772800*alpha^6 + 9329/170100*alpha^5 + 3224237/59875200*alpha^4 - 272927/3628800*alpha^3 - 6823639/20528640*alpha^2 +...
            138793/2721600*alpha + 40018577/179625600)*z^3/np^3;
        R2 = R2 + (-1/362880*alpha^13 - 41/1814400*alpha^12 + 5/62208*alpha^11 + 2897/3628800*alpha^10 - 13891/21772800*alpha^9 - 33493/3628800*alpha^8 + 4913/3628800*alpha^7 + 98863/2177280*alpha^6 - 24343/870912*alpha^5 - 22858097/239500800*alpha^4 + 99343/435456*alpha^3 + 10059743/102643200*alpha^2 -...
            94511/226800*alpha - 13421633/179625600)*z^3/np^3;
        R1 = R1 + (-1/20160*alpha^11 + 13/483840*alpha^10 + 391/241920*alpha^9 - 23/17920*alpha^8 - 3401/241920*alpha^7 + 481/69120*alpha^6 + 865/24192*alpha^5 + 95807/2419200*alpha^4 - 2143/45360*alpha^3 - 370619/2177280*alpha^2 + 787/25920*alpha + 498263/5443200)*z^2/np^3;
        R2 = R2 + (1/20160*alpha^11 + 31/96768*alpha^10 - 5/6912*alpha^9 - 2767/483840*alpha^8 + 43/16128*alpha^7 + 6733/207360*alpha^6 - 437/17280*alpha^5 - 553631/7257600*alpha^4 + 6949/51840*alpha^3 + 21173/217728*alpha^2 - 29/160*alpha - 283523/5443200)*z^2/np^3;
        R1 = R1 + (1/1920*alpha^9 - 1/5760*alpha^8 - 23/2880*alpha^7 + 257/103680*alpha^6 + 233/11520*alpha^5 + 36343/1451520*alpha^4 - 893/34560*alpha^3 - 20681/290304*alpha^2 + 131/8640*alpha + 10463/362880)*z^1/np^3;
        R2 = R2 + (-1/1920*alpha^9 - 7/2880*alpha^8 + 1/288*alpha^7 + 2173/103680*alpha^6 - 29/1280*alpha^5 - 84391/1451520*alpha^4 + 455/6912*alpha^3 + 103489/1451520*alpha^2 - 9/160*alpha - 1891/72576)*z^1/np^3;
        R1 = R1 + (-1/384*alpha^7 + 1/128*alpha^5 + 67/5760*alpha^4 - 1/96*alpha^3 - 17/864*alpha^2 + 1/192*alpha + 139/25920)/np^3;
        R2 = R2 + (1/384*alpha^7 + 1/96*alpha^6 - 7/384*alpha^5 - 217/5760*alpha^4 + 7/288*alpha^3 + 125/3456*alpha^2 - 1/128*alpha - 139/17280)/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (-1/5108103000*alpha^17 - 1/681080400*alpha^16 + 113/2918916000*alpha^15 + 257/1362160800*alpha^14 - 1/364000*alpha^13 - 1133/122472000*alpha^12 + 430603/4715172000*alpha^11 + 190201/857304000*alpha^10 - 654749/428652000*alpha^9 - 11059/3969000*alpha^8 + 40110733/3143448000*alpha^7 +...
            512587/28066500*alpha^6 - 6115957717/122594472000*alpha^5 - 128380519/2554051500*alpha^4 - 604989479/12259447200*alpha^3 + 5824114243/91945854000*alpha^2 + 2071465391/3405402000*alpha + 3010697831/91945854000)*z^6/np^2;
        R2 = R2 + (1/5108103000*alpha^17 + 1/227026800*alpha^16 + 1/416988000*alpha^15 - 383/817296480*alpha^14 - 461/265356000*alpha^13 + 3011/168399000*alpha^12 + 102953/1178793000*alpha^11 - 9277/30618000*alpha^10 - 373693/214326000*alpha^9 + 108959/47628000*alpha^8 + 16662799/1047816000*alpha^7 -...
            2972201/449064000*alpha^6 - 8029748039/122594472000*alpha^5 + 951900073/40864824000*alpha^4 + 337197743/15324309000*alpha^3 - 60473961767/183891708000*alpha^2 + 1174327333/4378374000*alpha + 49153599637/91945854000)*z^6/np^2;
        R1 = R1 + (1/97297200*alpha^15 + 1/14968800*alpha^14 - 89/64864800*alpha^13 - 499/89812800*alpha^12 + 5573/89812800*alpha^11 + 1381/8164800*alpha^10 - 70177/57153600*alpha^9 - 1301/544320*alpha^8 + 1063813/95256000*alpha^7 + 1123/68040*alpha^6 - 10035559/224532000*alpha^5 -...
            116387263/2554051500*alpha^4 - 139623829/3143448000*alpha^3 + 7118462903/122594472000*alpha^2 + 118854433/261954000*alpha + 655569283/30648618000)*z^5/np^2;
        R2 = R2 + (-1/97297200*alpha^15 - 1/4989600*alpha^14 - 1/4324320*alpha^13 + 149/11226600*alpha^12 + 169/3592512*alpha^11 - 1199/4082400*alpha^10 - 75409/57153600*alpha^9 + 1409/544320*alpha^8 + 1325917/95256000*alpha^7 - 22781/2721600*alpha^6 - 845351/14033250*alpha^5 +...
            356801723/13621608000*alpha^4 + 86755349/3143448000*alpha^3 - 34602825671/122594472000*alpha^2 + 125152271/785862000*alpha + 11925256253/30648618000)*z^5/np^2;
        R1 = R1 + (-1/2494800*alpha^13 - 1/453600*alpha^12 + 29/855360*alpha^11 + 43/388800*alpha^10 - 2441/2721600*alpha^9 - 3491/1814400*alpha^8 + 2843/302400*alpha^7 + 1409/97200*alpha^6 - 107323/2721600*alpha^5 - 268879/6652800*alpha^4 - 26591/680400*alpha^3 + 586921/11226600*alpha^2 +...
            73009/226800*alpha + 270601/22453200)*z^4/np^2;
        R2 = R2 + (1/2494800*alpha^13 + 1/151200*alpha^12 + 61/5987520*alpha^11 - 97/388800*alpha^10 - 317/388800*alpha^9 + 19/6720*alpha^8 + 10373/907200*alpha^7 - 28913/2721600*alpha^6 - 148493/2721600*alpha^5 + 904117/29937600*alpha^4 + 45917/1360800*alpha^3 - 330854/1403325*alpha^2 +...
            49711/680400*alpha + 5983139/22453200)*z^4/np^2;
        R1 = R1 + (1/90720*alpha^11 + 1/20160*alpha^10 - 7/12960*alpha^9 - 11/8064*alpha^8 + 89/12096*alpha^7 + 157/12960*alpha^6 - 11/324*alpha^5 - 62869/1814400*alpha^4 - 30127/907200*alpha^3 + 124667/2721600*alpha^2 + 5363/25200*alpha + 3247/680400)*z^3/np^2;
        R2 = R2 + (-1/90720*alpha^11 - 1/6720*alpha^10 - 23/90720*alpha^9 + 17/6048*alpha^8 + 163/20160*alpha^7 - 347/25920*alpha^6 - 619/12960*alpha^5 + 32693/907200*alpha^4 + 1361/33600*alpha^3 - 73703/388800*alpha^2 + 2183/226800*alpha + 22693/136080)*z^3/np^2;
        R1 = R1 + (-1/5040*alpha^9 - 1/1440*alpha^8 + 11/2240*alpha^7 + 79/8640*alpha^6 - 241/8640*alpha^5 - 113/4032*alpha^4 - 23/864*alpha^3 + 1739/45360*alpha^2 + 91/720*alpha - 17/45360)*z^2/np^2;
        R2 = R2 + (1/5040*alpha^9 + 1/480*alpha^8 + 23/6720*alpha^7 - 139/8640*alpha^6 - 329/8640*alpha^5 + 2717/60480*alpha^4 + 103/2160*alpha^3 - 1867/12960*alpha^2 - 67/2160*alpha + 817/9072)*z^2/np^2;
        R1 = R1 + (1/480*alpha^7 + 1/192*alpha^6 - 59/2880*alpha^5 - 29/1440*alpha^4 - 11/576*alpha^3 + 253/8640*alpha^2 + 1/16*alpha - 7/2160)*z^1/np^2;
        R2 = R2 + (-1/480*alpha^7 - 1/64*alpha^6 - 61/2880*alpha^5 + 7/120*alpha^4 + 31/576*alpha^3 - 173/1728*alpha^2 - 7/144*alpha + 79/2160)*z^1/np^2;
        R1 = R1 + (-1/96*alpha^5 - 1/96*alpha^4 - 1/96*alpha^3 + 5/288*alpha^2 + 1/48*alpha - 1/288)/np^2;
        R2 = R2 + (1/96*alpha^5 + 7/96*alpha^4 + 5/96*alpha^3 - 17/288*alpha^2 - 1/24*alpha + 1/144)/np^2;
    end
    R1 = R1 + (1/10216206000*alpha^16 - 1/56756700*alpha^14 + 103/84199500*alpha^12 - 179/4286520*alpha^10 + 6617/8930250*alpha^8 - 202/30375*alpha^6 + 10369976/383107725*alpha^4 + 332327847221/7815397590000*alpha^2 - 37533146551/177622672500)*z^7/np^1;
    R2 = R2 + (-1/10216206000*alpha^16 - 1/638512875*alpha^15 + 1/170270100*alpha^14 + 1/5212350*alpha^13 + 17/84199500*alpha^12 - 61/7016625*alpha^11 - 2371/107163000*alpha^10 + 1957/10716300*alpha^9 + 2609/4465125*alpha^8 - 302/165375*alpha^7 - 454/70875*alpha^6 + 56108/7016625*alpha^5 +...
        11109664/383107725*alpha^4 - 3287552/273648375*alpha^3 + 276425325499/7815397590000*alpha^2 + 69963511/547296750*alpha - 269282301619/1953849397500)*z^7/np^1;
    R1 = R1 + (-1/170270100*alpha^14 + 1/1403325*alpha^12 - 163/5103000*alpha^10 + 263/396900*alpha^8 - 1388/212625*alpha^6 + 1124/40095*alpha^4 + 91465217/2189187000*alpha^2 - 1411927117/7662154500)*z^6/np^1;
    R2 = R2 + (1/170270100*alpha^14 + 1/12162150*alpha^13 - 1/5613300*alpha^12 - 1/155925*alpha^11 - 47/5103000*alpha^10 + 89/510300*alpha^9 + 349/793800*alpha^8 - 599/297675*alpha^7 - 181/30375*alpha^6 + 2032/212625*alpha^5 + 41906/1403325*alpha^4 - 3008/200475*alpha^3 + 9755509/294698250*alpha^2 +...
        164284139/1277025750*alpha - 77072581/696559500)*z^6/np^1;
    R1 = R1 + (1/3742200*alpha^12 - 1/48600*alpha^10 + 5/9072*alpha^8 - 59/9450*alpha^6 + 1231/42525*alpha^4 + 419761691/10216206000*alpha^2 - 400897661/2554051500)*z^5/np^1;
    R2 = R2 + (-1/3742200*alpha^12 - 1/311850*alpha^11 + 1/340200*alpha^10 + 1/6804*alpha^9 + 11/45360*alpha^8 - 61/28350*alpha^7 - 29/5670*alpha^6 + 164/14175*alpha^5 + 1289/42525*alpha^4 - 824/42525*alpha^3 + 309830809/10216206000*alpha^2 + 24307/187110*alpha - 212336779/2554051500)*z^5/np^1;
    R1 = R1 + (-1/113400*alpha^10 + 1/2520*alpha^8 - 23/4050*alpha^6 + 4/135*alpha^4 + 15199/374220*alpha^2 - 241823/1871100)*z^4/np^1;
    R2 = R2 + (1/113400*alpha^10 + 1/11340*alpha^9 - 2/945*alpha^7 - 29/8100*alpha^6 + 19/1350*alpha^5 + 4/135*alpha^4 - 74/2835*alpha^3 + 100889/3742200*alpha^2 + 3743/28350*alpha - 103093/1871100)*z^4/np^1;
    R1 = R1 + (1/5040*alpha^8 - 1/216*alpha^6 + 4/135*alpha^4 + 9257/226800*alpha^2 - 5717/56700)*z^3/np^1;
    R2 = R2 + (-1/5040*alpha^8 - 1/630*alpha^7 - 1/1080*alpha^6 + 1/60*alpha^5 + 7/270*alpha^4 - 1/27*alpha^3 + 5263/226800*alpha^2 + 257/1890*alpha - 1483/56700)*z^3/np^1;
    R1 = R1 + (-1/360*alpha^6 + 1/36*alpha^4 + 65/1512*alpha^2 - 269/3780)*z^2/np^1;
    R2 = R2 + (1/360*alpha^6 + 1/60*alpha^5 + 1/72*alpha^4 - 1/18*alpha^3 + 79/3780*alpha^2 + 13/90*alpha + 17/3780)*z^2/np^1;
    R1 = R1 + (1/48*alpha^4 + 37/720*alpha^2 - 7/180)*z^1/np^1;
    R2 = R2 + (-1/48*alpha^4 - 1/12*alpha^3 + 23/720*alpha^2 + 1/6*alpha + 7/180)*z^1/np^1;
    R1 = R1 + (1/12*alpha^2)/np^1 + 1;
    R2 = R2 + (1/6*alpha^2 + 1/4*alpha + 1/12)/np^1;
end

p = real( sqrt(2*pi)*(-1)^np*sqrt(npb)/z^(1/4)/ ...
    (1 - z)^(1/4)*z^(-alpha/2)*( (sin( (alpha + 1)/2*acos(2*z - 1) - ...
    pi*alpha/2)*R1 -sin( (alpha - 1)/2*acos(2*z - 1) - pi*alpha/2)*...
    R2)*besselj(alpha, npb) + (cos( (alpha + 1)/2*acos(2*z - ...
    1)- pi*alpha/2)*R1 - cos( (alpha - 1)/2*acos(2*z - 1) - pi*alpha/2)*...
    R2)*(besselj(alpha-1, npb) - alpha/npb*besselj(alpha, npb) ) ) );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the expansion of the orthonormal polynomial near 4n without e^(x/2)
function p = asyAiry(np, y, alpha, T)
z = y/4/np;
fn = (np*3i*( sqrt(z).*sqrt(1 - z) - acos(sqrt(z) ) ))^(2/3);
d = z - 1;
if T == 1
    p = real( 4*sqrt(pi)/z^(1/4)/d^(1/4)*z^(-alpha/2)* ...
        (cos( (alpha + 1)/2*acos(2*z - 1) )*fn^(1/4)*airy(0,fn) + ...
        -1i*sin( (alpha + 1)/2*acos(2*z - 1) )*fn^(-1/4)*airy(1,fn) ) );
    return
end
R1 = 0.0;
R2 = 0.0;
if ( alpha == 0 )
    if ( T >= 7 )
        R1 = R1 + (+1.309628097160176e-05*d^2 -6.505325803676646e-06*d^1 )/np^6;
        R2 = R2 + (-0.001326131973043531*d^2 +0.0007635409517350945*d^1 )/np^6;
        R1 = R1 + (+2.100224847639528e-06 )/np^6;
        R2 = R2 + (-0.000344405038138915 )/np^6;
    end
    if ( T >= 6 )
        R1 = R1 + (+0.001509301781329904*d^3 -0.0009676950522129104*d^2 )/np^5;
        R2 = R2 + (-0.0001802687708021319*d^3 -0.0001238977582304063*d^2 )/np^5;
        R1 = R1 + (+0.00053567814019176*d^1 -0.0002131630878045043 )/np^5;
        R2 = R2 + (+0.0003181332313542115*d^1 -0.0004025677940905445 )/np^5;
    end
    if ( T >= 5 )
        R1 = R1 + (-2.796022042244435e-05*d^4 +2.011295619003682e-05*d^3 )/np^4;
        R2 = R2 + (+0.002475382540803074*d^4 -0.002022143985881445*d^3 )/np^4;
        R1 = R1 + (-1.295611489560738e-05*d^2 +6.757898749300978e-06*d^1 )/np^4;
        R2 = R2 + (+0.001569093104321315*d^2 -0.001116765488720164*d^1 )/np^4;
        R1 = R1 + (-2.026867991649663e-06 )/np^4;
        R2 = R2 + (+0.000666664945504233 )/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (-0.003253020840537458*d^5 +0.002696059409107524*d^4 )/np^3;
        R2 = R2 + (+0.001993718633874486*d^5 -0.001436399596567729*d^4 )/np^3;
        R1 = R1 + (-0.002141460647635705*d^3 +0.001590456473343028*d^2 )/np^3;
        R2 = R2 + (+0.0008821041836667568*d^3 -0.0003328986988216683*d^2 )/np^3;
        R1 = R1 + (-0.001045363705958944*d^1 +0.0005110818194151536 )/np^3;
        R2 = R2 + (-0.000206871642466881*d^1 +0.0007267403892403877 )/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (+7.18695260864394e-05*d^6 -6.692939666830364e-05*d^5 )/np^2;
        R2 = R2 + (-0.004480794599931599*d^6 +0.004460948839100688*d^5 )/np^2;
        R1 = R1 + (+6.109774830863335e-05*d^4 -5.407654074320793e-05*d^3 )/np^2;
        R2 = R2 + (-0.0044332342004771*d^4 +0.004392654646940364*d^3 )/np^2;
        R1 = R1 + (+4.5413316841889e-05*d^2 -3.439153439153486e-05*d^1 )/np^2;
        R2 = R2 + (-0.004329315657887089*d^2 +0.004221019721019722*d^1 )/np^2;
        R1 = R1 + (+1.984126984127046e-05 )/np^2;
        R2 = R2 + (-0.004007936507936511 )/np^2;
    end
    R1 = R1 + (+0.01272750016968636*d^7 -0.01255318649195314*d^6 )/np^1;
    R2 = R2 + (-0.01249169963638477*d^7 +0.01227111943654102*d^6 )/np^1;
    R1 = R1 + (+0.01234301148871137*d^5 -0.01208266968103703*d^4 )/np^1;
    R2 = R2 + (-0.0119973041110328*d^5 +0.01164522908686174*d^4 )/np^1;
    R1 = R1 + (+0.01174829614829615*d^3 -0.01129622758194187*d^2 )/np^1;
    R2 = R2 + (-0.01117002997002997*d^3 +0.01048155019583591*d^2 )/np^1;
    R1 = R1 + (+0.01063492063492063*d^1 -0.009523809523809525 )/np^1 + 1;
    R2 = R2 + (-0.009365079365079364*d^1 +0.007142857142857144 )/np^1;
elseif ( alpha == 1 )
    if ( T >= 7 )
        R1 = R1 + (-0.001116872178853068*d^2 +0.0006145274231811699*d^1 )/np^6;
        R2 = R2 + (+0.0003914779895445831*d^2 +0.0001817422654701506*d^1 )/np^6;
        R1 = R1 + (-0.0002421391621335997 )/np^6;
        R2 = R2 + (-0.0004609315000727462 )/np^6;
    end
    if ( T >= 6 )
        R1 = R1 + (-0.001217631725961165*d^3 +0.0007064505509263123*d^2 )/np^5;
        R2 = R2 + (-0.001284213690985701*d^3 +0.00112752157728142*d^2 )/np^5;
        R1 = R1 + (-0.0003330216038310101*d^1 +9.757071356383969e-05 )/np^5;
        R2 = R2 + (-0.0008332130005899818*d^1 +0.0004030585196629511 )/np^5;
    end
    if ( T >= 5 )
        R1 = R1 + (+0.001882808665107142*d^4 -0.001505607402838352*d^3 )/np^4;
        R2 = R2 + (-0.003079494756738815*d^4 +0.0021480231321534*d^3 )/np^4;
        R1 = R1 + (+0.001127942655094661*d^2 -0.0007500605012473177*d^1 )/np^4;
        R2 = R2 + (-0.00121916734519941*d^2 +0.0002978082751158517*d^1 )/np^4;
        R1 = R1 + (+0.0003729230795475838 )/np^4;
        R2 = R2 + (+0.0006008030043149 )/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (+0.004744560105086139*d^5 -0.003932667625282309*d^4 )/np^3;
        R2 = R2 + (-0.00143856736376185*d^5 +0.0006338239758711377*d^4 )/np^3;
        R1 = R1 + (+0.003122487606361571*d^3 -0.002315841844741004*d^2 )/np^3;
        R2 = R2 + (+0.0001615072603597489*d^3 -0.0009363363590450462*d^2 )/np^3;
        R1 = R1 + (+0.001516718193622956*d^1 -0.000734979989146656 )/np^3;
        R2 = R2 + (+0.001661552988457739*d^1 -0.002244249269249276 )/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (-0.0006821631748549443*d^6 +0.0006790103374176761*d^5 )/np^2;
        R2 = R2 + (+0.007247642006208204*d^6 -0.007153198600144216*d^5 )/np^2;
        R1 = R1 + (-0.0006715318212236979*d^4 +0.0006562001666763589*d^3 )/np^2;
        R2 = R2 + (+0.007017323836035319*d^4 -0.006807624756196191*d^3 )/np^2;
        R1 = R1 + (-0.0006256188256188247*d^2 +0.0005622895622895622*d^1 )/np^2;
        R2 = R2 + (+0.00644986864986863*d^2 -0.005737854737854754*d^1 )/np^2;
        R1 = R1 + (-0.0004166666666666676 )/np^2;
        R2 = R2 + (+0.003888888888888907 )/np^2;
    end
    R1 = R1 + (-0.04460324252123175*d^7 +0.04443907573247043*d^6 )/np^1;
    R2 = R2 + (+0.04348512118855518*d^7 -0.04306174737502545*d^6 )/np^1;
    R1 = R1 + (-0.04423440188249792*d^5 +0.04396981497716192*d^4 )/np^1;
    R2 = R2 + (+0.04248039440813511*d^5 -0.04163157791565959*d^4 )/np^1;
    R1 = R1 + (-0.04361026909598338*d^3 +0.04308472479901051*d^2 )/np^1;
    R2 = R2 + (+0.04027808382094096*d^3 -0.03780416408987839*d^2 )/np^1;
    R1 = R1 + (-0.04222222222222222*d^1 +0.04047619047619047 )/np^1 + 1;
    R2 = R2 + (+0.0320634920634921*d^1 +0.1571428571428572 )/np^1;
else
    if ( T >= 5 )
        R1 = R1 + (- 1013879/11823903759738470400*alpha^21 + 45064039/7882602506492313600*alpha^20 - 688951/4863802451558400*alpha^19 + 157459079/103718454032793600*alpha^18 - 3648427763/1067689967984640000*alpha^17 - 2550629833/36976275947520000*alpha^16 + ...
            8527600735679/16015349519769600000*alpha^15 + 55893547/10266249692160000*alpha^14 - 341936294746721/27102899187302400000*alpha^13 + 565487075439283/18068599458201600000*alpha^12 + 59655162402274583/587229482391552000000*alpha^11 - 2008247669052787/4448708199936000000*alpha^10 -...
            18993666326329/116729952768000000*alpha^9 + 31237554475352881/13346124599808000000*alpha^8 - 1078157996924747183/850815443237760000000*alpha^7 - 6134986008321159443/2268841181967360000000*alpha^6 - 103232383455348087349/395156505859315200000000*alpha^5 - 3502648344145848473/8781255685762560000000*alpha^4 +...
            180559269933337639/33329664799200000000*alpha^3 + 16610189119023418049/8672845121740800000000*alpha^2 - 12207497503256717345429/9572483411615124000000000*alpha - 19666701992772351781487/47126072180259072000000000)*d^4/np^4;
        R2 = R2 + (- 1013879/11823903759738470400*alpha^21 + 177334499/39413012532461568000*alpha^20 - 1808203/24892428967870464*alpha^19 + 81757531/518592270163968000*alpha^18 + 110848813709/18150729455738880000*alpha^17 - 3505768471/81347807084544000*alpha^16 -...
            2067463006801/16015349519769600000*alpha^15 + 4831136071/2566562423040000*alpha^14 - 330226749484013/352337689434931200000*alpha^13 - 612287562559613/18068599458201600000*alpha^12 + 35929915687839863/587229482391552000000*alpha^11 + 1232444348557757/4448708199936000000*alpha^10 -...
            4050296304548881/5719767685632000000*alpha^9 - 91760498774033/83413278748800000*alpha^8 + 3849766087526241787/850815443237760000000*alpha^7 - 5891116878044670803/2268841181967360000000*alpha^6 - 2517421413160487142949/395156505859315200000000*alpha^5 + 15500210153740590841/1254465097966080000000*alpha^4 +...
            535448090870817469/1176061029343200000000*alpha^3 - 876031664564464102199/78055606095667200000000*alpha^2 + 143449595956451079161/870225764692284000000000*alpha + 98306009664670335668113/47126072180259072000000000)*d^4/np^4;
        R1 = R1 + (-37/9654373048320*alpha^19 + 1/4311014400*alpha^18 - 22163/4234374144000*alpha^17 + 26401/498161664000*alpha^16 - 1514633/8717829120000*alpha^15 - 2369/2286144000*alpha^14 + 188962531/19615115520000*alpha^13 - 101881/8382528000*alpha^12 - 3191929607/30177100800000*alpha^11 +...
            1195568051/3621252096000*alpha^10 + 21979238557/78460462080000*alpha^9 - 224264056069/112086374400000*alpha^8 + 13216977279629/13730580864000000*alpha^7 + 176875487699/78460462080000*alpha^6 + 51777130063/260513253000000*alpha^5 + 2416907662549/10003708915200000*alpha^4 -...
            101457798071933/22627436832000000*alpha^3 - 27341574857/111152321280000*alpha^2 + 1450813856961629/1385930505960000000*alpha + 546571755859/27175107960000000)*d^3/np^4;
        R2 = R2 + (-37/9654373048320*alpha^19 + 239/1270312243200*alpha^18 - 1759/604910592000*alpha^17 + 569/57480192000*alpha^16 + 3530881/26153487360000*alpha^15 - 553031/523069747200*alpha^14 - 16827929/19615115520000*alpha^13 + 153704773/6035420160000*alpha^12 - 977315567/30177100800000*alpha^11 -...
            21967814521/90531302400000*alpha^10 + 122432162953/235381386240000*alpha^9 + 220166487539/196151155200000*alpha^8 - 55659560808911/13730580864000000*alpha^7 + 405710676571/217945728000000*alpha^6 + 55836859575139/8336424096000000*alpha^5 - 52372063912321/5001854457600000*alpha^4 -...
            943973862545603/475176173472000000*alpha^3 + 3900102001763159/369581468256000000*alpha^2 + 246033374218259/1385930505960000000*alpha - 254777367043331/125993682360000000)*d^3/np^4;
        R1 = R1 + (-19/169374965760*alpha^17 + 179/28466380800*alpha^16 - 541/4151347200*alpha^15 + 14221/11496038400*alpha^14 - 684823/145297152000*alpha^13 - 935293/160944537600*alpha^12 + 33475339/335301120000*alpha^11 - 1450837/7315660800*alpha^10 - 24868537/60963840000*alpha^9 +...
            4497268903/2682408960000*alpha^8 - 3960305773/7264857600000*alpha^7 - 27136359533/13076743680000*alpha^6 - 1210955477/18162144000000*alpha^5 + 42792983/65383718400000*alpha^4 + 67569898733/20583763200000*alpha^3 + 4899641437/32691859200000*alpha^2 - 23594745318577/30798455688000000*alpha -...
            21001491079/1620971352000000)*d^2/np^4;
        R2 = R2 + (-19/169374965760*alpha^17 + 263/49816166400*alpha^16 - 1189/14944849920*alpha^15 + 26069/74724249600*alpha^14 + 531437/326918592000*alpha^13 - 3264001/201180672000*alpha^12 + 1203443/111767040000*alpha^11 + 370159/1905120000*alpha^10 - 66717731/182891520000*alpha^9 -...
            707915987/670602240000*alpha^8 + 120338278349/32691859200000*alpha^7 - 8763140129/5943974400000*alpha^6 - 565047950849/81729648000000*alpha^5 + 14221072519/1634592960000*alpha^4 + 1788611634253/555761606400000*alpha^3 - 7700210671649/926269344000000*alpha^2 - 216475631431/473822395200000*alpha +...
            144976933331359/92395367064000000)*d^2/np^4;
        R1 = R1 + (-289/149448499200*alpha^15 + 13/127733760*alpha^14 - 20783/10674892800*alpha^13 + 715/41803776*alpha^12 - 967949/14370048000*alpha^11 + 226711/5225472000*alpha^10 + 1280737/2612736000*alpha^9 - 436663/348364800*alpha^8 + 248247547/2011806720000*alpha^7 + 4880461/2737152000*alpha^6 -...
            87915739/4670265600000*alpha^5 - 352344007/1868106240000*alpha^4 - 6273086233/2971987200000*alpha^3 - 8775229/133436160000*alpha^2 + 114086343257/231567336000000*alpha + 276160343/40864824000000)*d^1/np^4;
        R2 = R2 + (-289/149448499200*alpha^15 + 149/1660538880*alpha^14 - 101021/74724249600*alpha^13 + 79841/11496038400*alpha^12 + 37673/7185024000*alpha^11 - 103223/746496000*alpha^10 + 80941/373248000*alpha^9 + 240781/248832000*alpha^8 - 6650525873/2011806720000*alpha^7 + 105546299/95800320000*alpha^6 ...
            + 1186883177/166795200000*alpha^5 - 3232076083/467026560000*alpha^4 - 144767064173/32691859200000*alpha^3 + 198377506939/32691859200000*alpha^2 + 24026725061/33081048000000*alpha - 775819227479/694702008000000)*d^1/np^4;
        R1 = R1 + (-23/1660538880*alpha^13 + 5/6967296*alpha^12 - 8473/638668800*alpha^11 + 2143/19353600*alpha^10 - 25037/58060800*alpha^9 + 37643/58060800*alpha^8 + 72713/290304000*alpha^7 - 43651/34836480*alpha^6 + 6671/228096000*alpha^5 + 23093/85536000*alpha^4 + 9359089/9434880000*alpha^3 +...
            13679/4790016000*alpha^2 - 210525163/908107200000*alpha - 22087361/10897286400000)/np^4;
        R2 = R2 + (-23/1660538880*alpha^13 + 29/42577920*alpha^12 - 7213/638668800*alpha^11 + 12581/174182400*alpha^10 - 14911/174182400*alpha^9 - 49019/58060800*alpha^8 + 852013/290304000*alpha^7 - 217939/290304000*alpha^6 - 111019/15206400*alpha^5 + 4086851/798336000*alpha^4 +...
            27653519/4942080000*alpha^3 - 18286033/4790016000*alpha^2 - 121127887/123832800000*alpha + 201801079/302702400000)/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (137/292416106752000*alpha^20 - 1/54486432000*alpha^19 + 41269/200074178304000*alpha^18 + 551/1389404016000*alpha^17 - 410869/19615115520000*alpha^16 + 971/12770257500*alpha^15 + 96338279/147113366400000*alpha^14 - 5741261/1362160800000*alpha^13 - 423142409/62240270400000*alpha^12 +...
            3076720159/36306824400000*alpha^11 - 2188394969/107270163000000*alpha^10 - 328873402133/429080652000000*alpha^9 + 1894875073367/2574483912000000*alpha^8 + 55167609839/17878360500000*alpha^7 - 3293297140077977/875324530080000000*alpha^6 - 33208583791723/18235927710000000*alpha^5 -...
            54475015980190937/15245235565560000000*alpha^4 - 2088551585106851/7622617782780000000*alpha^3 + 360008253055089407/25408725942600000000*alpha^2 + 1674345981481/11144178045000000*alpha - 216246444678354844343/66475579247327250000000)*d^5/np^3;
        R2 = R2 + (137/292416106752000*alpha^20 - 10589/950352346944000*alpha^19 - 107/9527341824000*alpha^18 + 823/490377888000*alpha^17 - 93629/19615115520000*alpha^16 - 26567/272432160000*alpha^15 + 8284187/21016195200000*alpha^14 + 29148803/10508097600000*alpha^13 - 21923549/1728896400000*alpha^12 -...
            5882440051/145227297600000*alpha^11 + 1278649279/6501222000000*alpha^10 + 13873987321/47675628000000*alpha^9 - 3837595190047/2574483912000000*alpha^8 - 11567507071/12872419560000*alpha^7 + 1983516634032821/291774843360000000*alpha^6 - 9849992297539/1870351560000000*alpha^5 -...
            14337677195357843/3811308891390000000*alpha^4 + 2926860258278599/282319177140000000*alpha^3 - 21018445670230627/2823191771400000000*alpha^2 - 25277539613069/11634032025000000*alpha + 397600803128989077941/199426737741981750000000)*d^5/np^3;
        R1 = R1 + (47/1818856166400*alpha^18 - 1/1089728640*alpha^17 + 19099/1961511552000*alpha^16 - 19/35026992000*alpha^15 - 5621123/9807557760000*alpha^14 + 5682269/2451889440000*alpha^13 + 21133337/2263282560000*alpha^12 - 12509719/188606880000*alpha^11 - 156750961/6601240800000*alpha^10 +...
            655310503/943034400000*alpha^9 - 10629710969/21454032600000*alpha^8 - 64279006669/21454032600000*alpha^7 + 219466363153/64362097800000*alpha^6 + 15643485311/8581613040000*alpha^5 + 157729235747/60786425700000*alpha^4 + 341371466027/1458874216800000*alpha^3 - 11863501792245991/1016349037704000000*alpha^2 -...
            17939/126360000*alpha + 540816589349851/200595204810000000)*d^4/np^3;
        R2 = R2 + (47/1818856166400*alpha^18 - 71/123502579200*alpha^17 + 257/980755776000*alpha^16 + 13943/245188944000*alpha^15 - 1807283/9807557760000*alpha^14 - 925027/445798080000*alpha^13 + 438493/51438240000*alpha^12 + 6548873/188606880000*alpha^11 - 4144229479/26404963200000*alpha^10 -...
            508264153/1886068800000*alpha^9 + 224468429441/171632260800000*alpha^8 + 9318491023/10727016300000*alpha^7 - 6775508879891/1029793564800000*alpha^6 + 904295467573/171632260800000*alpha^5 + 524074996199/108064756800000*alpha^4 - 15068497679749/1458874216800000*alpha^3 +...
            5000521151686019/1016349037704000000*alpha^2 + 91000970429207/42347876571000000*alpha - 5474562553987601/3811308891390000000)*d^4/np^3;
        R1 = R1 + (1/958003200*alpha^16 - 1/29937600*alpha^15 + 853/2594592000*alpha^14 - 3503/7783776000*alpha^13 - 5783/574801920*alpha^12 + 60817/1397088000*alpha^11 + 429649/6531840000*alpha^10 - 3417521/5715360000*alpha^9 + 1060721/4656960000*alpha^8 + 30056657/10478160000*alpha^7 -...
            97555360009/32691859200000*alpha^6 - 13586567/7484400000*alpha^5 - 2736174359/1634592960000*alpha^4 - 761324623/4086482400000*alpha^3 + 425829704527/46313467200000*alpha^2 + 2134801/16216200000*alpha - 16488420215951/7699613922000000)*d^3/np^3;
        R2 = R2 + (1/958003200*alpha^16 - 17/778377600*alpha^15 + 79/2223936000*alpha^14 + 20609/15567552000*alpha^13 - 4481/958003200*alpha^12 - 38813/1397088000*alpha^11 + 755969/6531840000*alpha^10 + 345937/1428840000*alpha^9 - 1107929/997920000*alpha^8 - 8689063/10478160000*alpha^7 +...
            68986687457/10897286400000*alpha^6 - 14357592647/2724321600000*alpha^5 - 690902021/116756640000*alpha^4 + 41976041801/4086482400000*alpha^3 - 111444353899/46313467200000*alpha^2 - 12252648061/5789183400000*alpha + 452790776881/513307594800000)*d^3/np^3;
        R1 = R1 + (179/6227020800*alpha^14 - 1/1209600*alpha^13 + 887/119750400*alpha^12 - 631/39916800*alpha^11 - 37531/381024000*alpha^10 + 8327/18144000*alpha^9 + 54611/762048000*alpha^8 - 56563/21168000*alpha^7 + 205608943/83825280000*alpha^6 + 3118943/1746360000*alpha^5 +...
            228067699/272432160000*alpha^4 + 11665859/90810720000*alpha^3 - 6123651893/908107200000*alpha^2 - 29489/249480000*alpha + 10230493571/6432426000000)*d^2/np^3;
        R2 = R2 + (179/6227020800*alpha^14 - 1789/3113510400*alpha^13 + 703/479001600*alpha^12 + 211/10886400*alpha^11 - 12461/169344000*alpha^10 - 5599/27216000*alpha^9 + 1363321/1524096000*alpha^8 + 737/952560*alpha^7 - 253419281/41912640000*alpha^6 + 220829699/41912640000*alpha^5 +...
            1264927661/181621440000*alpha^4 - 2773666603/272432160000*alpha^3 - 27020101/302702400000*alpha^2 + 3730921/1801800000*alpha - 6424038737/19297278000000)*d^2/np^3;
        R1 = R1 + (13/27371520*alpha^12 - 1/80640*alpha^11 + 8681/87091200*alpha^10 - 2671/10886400*alpha^9 - 11503/29030400*alpha^8 + 529/226800*alpha^7 - 21883/12441600*alpha^6 - 30839/18144000*alpha^5 - 20563/149688000*alpha^4 - 34523/598752000*alpha^3 + 5359207/1235520000*alpha^2 + 451/4536000*alpha -...
            237325577/227026800000)*d^1/np^3;
        R2 = R2 + (13/27371520*alpha^12 - 151/15966720*alpha^11 + 2801/87091200*alpha^10 + 3361/21772800*alpha^9 - 19051/29030400*alpha^8 - 4981/7257600*alpha^7 + 498563/87091200*alpha^6 - 38179/7257600*alpha^5 - 9576269/1197504000*alpha^4 + 5999563/598752000*alpha^3 + 22045619/8648640000*alpha^2 -...
            6495887/3243240000*alpha - 46965407/227026800000)*d^1/np^3;
        R1 = R1 + (1/290304*alpha^10 - 1/11520*alpha^9 + 313/483840*alpha^8 - 131/80640*alpha^7 + 1949/2419200*alpha^6 + 191/134400*alpha^5 - 1067/3628800*alpha^4 - 29/1209600*alpha^3 - 29929/14784000*alpha^2 - 1/14400*alpha + 1657561/3243240000)/np^3;
        R2 = R2 + (1/290304*alpha^10 - 1/12960*alpha^9 + 187/483840*alpha^8 + 25/48384*alpha^7 - 13031/2419200*alpha^6 + 1583/302400*alpha^5 + 2053/226800*alpha^4 - 35129/3628800*alpha^3 - 654361/133056000*alpha^2 + 63127/33264000*alpha + 1571329/2162160000)/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (-1/808122744000*alpha^19 + 1/45972927000*alpha^18 + 83/868377510000*alpha^17 - 271/81729648000*alpha^16 + 1013/357567210000*alpha^15 + 421769/2145403260000*alpha^14 - 74489/153243090000*alpha^13 - 409201/70727580000*alpha^12 + 411925007/22691765250000*alpha^11 +...
            246248521/2750517000000*alpha^10 - 5376311051/17878360500000*alpha^9 - 2386255199/3352192593750*alpha^8 + 182033236859/78217827187500*alpha^7 + 53058329293/20113155562500*alpha^6 - 433281906517/56987274093750*alpha^5 + 25206947887853/6078642570000000*alpha^4 + 5317008498553/17437360941000000*alpha^3 -...
            5868217447/4169165000000*alpha^2 - 154015299939251/617573199993750000*alpha + 133154079621709/1852719599981250000)*d^6/np^2;
        R2 = R2 + (-1/808122744000*alpha^19 + 1/1042053012000*alpha^18 + 37/133596540000*alpha^17 - 47/306486180000*alpha^16 - 1648/67043851875*alpha^15 + 10319/1072701630000*alpha^14 + 28459/25540515000*alpha^13 - 661/2182950000*alpha^12 - 418609067/15127843500000*alpha^11 +...
            20925739/4125775500000*alpha^10 + 359262577/957769312500*alpha^9 - 2363908969/53635081500000*alpha^8 - 611929740043/234653481562500*alpha^7 + 5921049677/33521925937500*alpha^6 + 383863454012/47489395078125*alpha^5 - 144527997487591/9117963855000000*alpha^4 - 3265529483516039/4446527039955000000*alpha^3 +...
            516521870325379/23526598095000000*alpha^2 + 245634249243059/617573199993750000*alpha - 2767218659594471/617573199993750000)*d^6/np^2;
        R1 = R1 + (-19/231567336000*alpha^17 + 1/778377600*alpha^16 + 23/6810804000*alpha^15 - 29/224532000*alpha^14 + 47/272432160*alpha^13 + 953/196465500*alpha^12 - 957653/78586200000*alpha^11 - 610493/7144200000*alpha^10 + 137743259/550103400000*alpha^9 + 8171983/11226600000*alpha^8 -...
            3206699123/1489863375000*alpha^7 - 355798627/127702575000*alpha^6 + 2371222339/319256437500*alpha^5 - 279087041369/71513442000000*alpha^4 - 207669760061/607864257000000*alpha^3 + 5271198011/3972969000000*alpha^2 + 5099280384049/17644948571250000*alpha - 2302077509/34395611250000)*d^5/np^2;
        R2 = R2 + (-19/231567336000*alpha^17 + 1/13621608000*alpha^16 + 263/20432412000*alpha^15 - 79/10216206000*alpha^14 - 15767/20432412000*alpha^13 + 2413/7858620000*alpha^12 + 1770277/78586200000*alpha^11 - 10307/1786050000*alpha^10 - 185430361/550103400000*alpha^9 + 28964977/550103400000*alpha^8 +...
            2783759599/1117397531250*alpha^7 - 134429569/638512875000*alpha^6 - 2560071383/319256437500*alpha^5 + 567439200521/35756721000000*alpha^4 + 176850964159/202621419000000*alpha^3 - 295639943713/13508094600000*alpha^2 - 53515639897/115326461250000*alpha + 4630188990877/1037938151250000)*d^5/np^2;
        R1 = R1 + (-17/4086482400*alpha^15 + 1/17463600*alpha^14 + 1/14594580*alpha^13 - 29/8164800*alpha^12 + 92093/15717240000*alpha^11 + 111263/1428840000*alpha^10 - 135391/714420000*alpha^9 - 702449/952560000*alpha^8 + 265276031/137525850000*alpha^7 + 11624803/3929310000*alpha^6 -...
            229827869/31925643750*alpha^5 + 1051670021/291891600000*alpha^4 + 2736245641/7151344200000*alpha^3 - 139318793/113513400000*alpha^2 - 1123775113/3377023650000*alpha + 12136973/198648450000)*d^4/np^2;
        R2 = R2 + (-17/4086482400*alpha^15 + 1/227026800*alpha^14 + 8/18243225*alpha^13 - 2/7016625*alpha^12 - 263437/15717240000*alpha^11 + 1171/178605000*alpha^10 + 208139/714420000*alpha^9 - 3869/59535000*alpha^8 - 80752051/34381462500*alpha^7 + 5142581/19646550000*alpha^6 +...
            92325203/11609325000*alpha^5 - 16238539043/1021620600000*alpha^4 - 7603702313/7151344200000*alpha^3 + 2472711463/113513400000*alpha^2 + 1834530377/3377023650000*alpha - 1361012431/307002150000)*d^4/np^2;
        R1 = R1 + (-1/6486480*alpha^13 + 1/544320*alpha^12 + 1/5544000*alpha^11 - 97/1512000*alpha^10 + 1829/15876000*alpha^9 + 947/1296000*alpha^8 - 18497/11340000*alpha^7 - 21491/6804000*alpha^6 + 53689/7796250*alpha^5 - 11260331/3492720000*alpha^4 - 746231/1746360000*alpha^3 + 640867/582120000*alpha^2 ...
            + 1198759/3153150000*alpha - 219229/4054050000)*d^3/np^2;
        R2 = R2 + (-1/6486480*alpha^13 + 1/4989600*alpha^12 + 173/16632000*alpha^11 - 7/972000*alpha^10 - 11173/47628000*alpha^9 + 1327/15876000*alpha^8 + 8191/3780000*alpha^7 - 3929/11340000*alpha^6 - 10939/1386000*alpha^5 + 83477897/5239080000*alpha^4 + 1408177/1047816000*alpha^3 -...
            490811131/22702680000*alpha^2 - 6044693/9459450000*alpha + 13850699/3153150000)*d^3/np^2;
        R1 = R1 + (-13/3326400*alpha^11 + 1/25200*alpha^10 - 1/37800*alpha^9 - 59/86400*alpha^8 + 89/73500*alpha^7 + 341/100800*alpha^6 - 4859/756000*alpha^5 + 2747/1008000*alpha^4 + 54961/116424000*alpha^3 - 1411/1512000*alpha^2 - 809969/1891890000*alpha + 2203/48510000)*d^2/np^2;
        R2 = R2 + (-13/3326400*alpha^11 + 1/151200*alpha^10 + 73/453600*alpha^9 - 17/151200*alpha^8 - 10147/5292000*alpha^7 + 383/756000*alpha^6 + 743/94500*alpha^5 - 3457/216000*alpha^4 - 629603/349272000*alpha^3 + 354743/16632000*alpha^2 + 475907/630630000*alpha - 744599/171990000)*d^2/np^2;
        R1 = R1 + (-11/181440*alpha^9 + 1/1920*alpha^8 - 173/302400*alpha^7 - 31/8640*alpha^6 + 49/8640*alpha^5 - 1231/604800*alpha^4 - 13/25920*alpha^3 + 211/302400*alpha^2 + 1931/4158000*alpha - 13/378000)*d^1/np^2;
        R2 = R2 + (-11/181440*alpha^9 + 1/6720*alpha^8 + 457/302400*alpha^7 - 19/21600*alpha^6 - 343/43200*alpha^5 + 2447/151200*alpha^4 + 349/129600*alpha^3 - 6269/302400*alpha^2 - 3679/4158000*alpha + 17551/4158000)*d^1/np^2;
        R1 = R1 + (-1/2240*alpha^7 + 1/288*alpha^6 - 1/240*alpha^5 + 1/960*alpha^4 + 1/2240*alpha^3 - 1/2880*alpha^2 - 11/25200*alpha + 1/50400)/np^2;
        R2 = R2 + (-1/2240*alpha^7 + 1/480*alpha^6 + 1/120*alpha^5 - 49/2880*alpha^4 - 103/20160*alpha^3 + 11/576*alpha^2 + 1/1050*alpha - 101/25200)/np^2;
    end
    R1 = R1 + (1/781539759000*alpha^18 - 43/153243090000*alpha^16 + 6613/268175407500*alpha^14 - 49487/44204737500*alpha^12 + 29136593/1031443875000*alpha^10 - 1331657837/3352192593750*alpha^8 + 74404812682/25141444453125*alpha^6 - 14883503848/1439072578125*alpha^4 -...
        10205946194812597/205857733331250000*alpha^2 + 178640068357141/14035754545312500)*d^7/np^1;
    R2 = R2 + (1/781539759000*alpha^18 + 1/43418875500*alpha^17 - 1/11787930000*alpha^16 - 1/290233125*alpha^15 - 1367/268175407500*alpha^14 + 1271/6385128750*alpha^13 + 28183/44204737500*alpha^12 - 6929/1227909375*alpha^11 - 23873287/1031443875000*alpha^10 + 2833531/34381462500*alpha^9 +...
        654133979/1676096296875*alpha^8 - 331517866/558698765625*alpha^7 - 80301183938/25141444453125*alpha^6 + 555843532/310388203125*alpha^5 + 68957024/5922109375*alpha^4 - 6724463104/5276599453125*alpha^3 + 9847410960789173/205857733331250000*alpha^2 - 2903284233728/4288702777734375*alpha -...
        1928634729450727/154393299998437500)*d^7/np^1;
    R1 = R1 + (1/10216206000*alpha^16 - 19/1277025750*alpha^14 + 863/982327500*alpha^12 - 5737/223256250*alpha^10 + 13538081/34381462500*alpha^8 - 245565074/79814109375*alpha^6 + 186118004/16930265625*alpha^4 + 35146182161/721586250000*alpha^2 - 129208525865897/10292886666562500)*d^6/np^1;
    R2 = R2 + (1/10216206000*alpha^16 + 1/638512875*alpha^15 - 2/638512875*alpha^14 - 2/13030875*alpha^13 - 73/245581875*alpha^12 + 124/22325625*alpha^11 + 2029/111628125*alpha^10 - 2062/22325625*alpha^9 - 25521383/68762925000*alpha^8 + 6122537/8595365625*alpha^7 + 264991141/79814109375*alpha^6 -...
        177095444/79814109375*alpha^5 - 2341022746/186232921875*alpha^4 + 296762176/186232921875*alpha^3 - 1310256809111/28141863750000*alpha^2 + 486442496/586288828125*alpha + 42101747210723/3430962222187500)*d^6/np^1;
    R1 = R1 + (1/170270100*alpha^14 - 1/1701000*alpha^12 + 653/29767500*alpha^10 - 45539/119070000*alpha^8 + 2618237/818606250*alpha^6 - 62690297/5320940625*alpha^4 - 15770340667/331080750000*alpha^2 + 1929751931/156343687500)*d^5/np^1;
    R2 = R2 + (1/170270100*alpha^14 + 1/12162150*alpha^13 - 1/18711000*alpha^12 - 23/4677750*alpha^11 - 109/9922500*alpha^10 + 101/992250*alpha^9 + 39451/119070000*alpha^8 - 12989/14883750*alpha^7 - 2807183/818606250*alpha^6 + 1165336/409303125*alpha^5 + 24376991/1773646875*alpha^4 -...
        1224632/591215625*alpha^3 + 1651448867/36786750000*alpha^2 - 7241408/6897515625*alpha - 16881324883/1407093187500)*d^5/np^1;
    R1 = R1 + (1/3742200*alpha^12 - 1/60750*alpha^10 + 2833/7938000*alpha^8 - 3523/1063125*alpha^6 + 1044779/81860625*alpha^4 + 5683373/122850000*alpha^2 - 200016967/16554037500)*d^4/np^1;
    R2 = R2 + (1/3742200*alpha^12 + 1/311850*alpha^11 + 1/850500*alpha^10 - 1/9450*alpha^9 - 1997/7938000*alpha^8 + 361/330750*alpha^7 + 14753/4252500*alpha^6 - 2717/708750*alpha^5 - 1251361/81860625*alpha^4 + 25874/9095625*alpha^3 - 402946889/9459450000*alpha^2 + 1184/853125*alpha +...
        192775559/16554037500)*d^4/np^1;
    R1 = R1 + (1/113400*alpha^10 - 23/75600*alpha^8 + 641/189000*alpha^6 - 142/10125*alpha^4 - 718433/16170000*alpha^2 + 264601/22522500)*d^3/np^1;
    R2 = R2 + (1/113400*alpha^10 + 1/11340*alpha^9 + 1/10800*alpha^8 - 13/9450*alpha^7 - 619/189000*alpha^6 + 523/94500*alpha^5 + 2467/141750*alpha^4 - 43/10125*alpha^3 + 1900291/48510000*alpha^2 - 1976/1010625*alpha - 27953/2502500)*d^3/np^1;
    R1 = R1 + (1/5040*alpha^8 - 1/300*alpha^6 + 11/700*alpha^4 + 5267/126000*alpha^2 - 9133/808500)*d^2/np^1;
    R2 = R2 + (1/5040*alpha^8 + 1/630*alpha^7 + 1/450*alpha^6 - 2/225*alpha^5 - 257/12600*alpha^4 + 23/3150*alpha^3 - 467/14000*alpha^2 + 8/2625*alpha + 25423/2425500)*d^2/np^1;
    R1 = R1 + (1/360*alpha^6 - 13/720*alpha^4 - 947/25200*alpha^2 + 67/6300)*d^1/np^1;
    R2 = R2 + (1/360*alpha^6 + 1/60*alpha^5 + 17/720*alpha^4 - 1/60*alpha^3 + 523/25200*alpha^2 - 1/175*alpha - 59/6300)*d^1/np^1;
    R1 = R1 + (1/48*alpha^4 + 7/240*alpha^2 - 1/105)/np^1 + 1;
    R2 = R2 + (1/48*alpha^4 + 1/12*alpha^3 + 7/240*alpha^2 + 1/60*alpha + 1/140)/np^1;
end

p = real( 4*sqrt(pi)/z^(1/4)/d^(1/4)*z^(-alpha/2)* ...
    ( (R1*cos( (alpha + 1)/2*acos(2*z - 1) ) -cos( (alpha - 1)/2* ...
    acos(2*z - 1) )*R2)*fn^(1/4)*airy(0,fn) + ...
    1i*(-sin( (alpha + 1)/2*acos(2*z - 1) )*R1 +sin( (alpha - 1)/2*...
    acos(2*z - 1) )*R2)*fn^(-1/4)*airy(1,fn) ) );

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% Routines for FAST algorithm %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, ders] = alg0_Lag(n)
ders = zeros(n, 1);
xs = 1/(2*n+1);
n1 = 20;
n1 = min(n1, n);
x = zeros(n, 1);
for k = 1:n1
    [xs, ders(k)] = alg3_Lag(n, xs);
    x(k) = xs;
    xs = 1.1*xs;
end
[x, ders] = alg1_Lag(x, ders, n, n1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [roots, ders] = alg1_Lag(roots, ders, n, n1)

% number of terms in Taylor expansion
m = 30;

% initialise
hh1 = ones(m+1, 1); 
zz = zeros(m, 1); 
u = zeros(1, m+1); 
up = zeros(1, m+1);

x = roots(n1);
for j = n1:(n - 1)
    
    % initial approx
    h = rk2_Lag(pi/2, -pi/2, x, n) - x;
    
    % scaling:
    M = 1/h; 
    M2 = M^2; 
    M3 = M^3; 
    M4 = M^4;
    
    % recurrence relation for Laguerre polynomials
    r = x*(n + .5 - .25*x);  
    p = x^2;
    u(1:2) = [0 ; ders(j)/M];
    u(3) = -.5*u(2)/(M*x) - (n + .5 - .25*x)*u(1)/(x*M^2);
    u(4) = -u(3)/(M*x) + ( -(1+r)*u(2)/6/M^2 - (n+.5-.5*x)*u(1)/M^3 ) / p;
    up(1:3) = [u(2) ; 2*u(3)*M ; 3*u(4)*M];
    
    for k = 2:(m - 2)
        u(k+3) = ( -x*(2*k+1)*(k+1)*u(k+2)/M - (k*k+r)*u(k+1)/M2 - ...
                   (n+.5-.5*x)*u(k)/M3 + .25*u(k-1)/M4 ) / (p*(k+2)*(k+1));
        up(k+2) = (k+2)*u(k+3)*M;
    end
    up(m+1) = 0;
    
    % Flip for more accuracy in inner product calculation.
    u = u(m+1:-1:1);  
    up = up(m+1:-1:1);
    
    % Newton iteration
    hh = hh1; 
    hh(end) = M;    
    step = inf;  
    l = 0;
    if ( M == 1 )
        Mhzz = (M*h) + zz;
        hh = [M ; cumprod(Mhzz)];
        hh = hh(end:-1:1);
    end
    while ( (abs(step) > eps) && (l < 10) )
        l = l + 1;
        step = (u*hh)/(up*hh);
        h = h - step;
        Mhzz = (M*h) + zz;
        % Powers of h (This is the fastest way!)
        hh = [M ; cumprod(Mhzz)];     
        % Flip for more accuracy in inner product
        hh = hh(end:-1:1);          
    end
    
    % Update
    x = x + h;
    roots(j+1) = x;
    ders(j+1) = up*hh;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1, d1] = alg3_Lag(n, xs)
[u, up] = eval_Lag(xs, n);
theta = atan(sqrt(xs/(n + .5 - .25*xs))*up/u);
x1 = rk2_Lag(theta, -pi/2, xs, n);

% Newton iteration
step = inf;  
l = 0;
while ( (abs(step) > eps || abs(u) > eps) && (l < 200) )
    l = l + 1;
    [u, up] = eval_Lag(x1, n);
    step = u/up;
    x1 = x1 - step;
end

[~, d1] = eval_Lag(x1, n);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluate Laguerre polynomial via recurrence
function [L, Lp] = eval_Lag(x, n)
Lm2 = 0; 
Lm1 = exp(-x/2); 
Lpm2 = 0; 
Lpm1 = 0;
for k = 0:n-1
    L = ( (2*k+1-x).*Lm1 - k*Lm2 ) / (k + 1);
    Lp = ( (2*k+1-x).*Lpm1 - Lm1 - k*Lpm2 ) / (k + 1);
    Lm2 = Lm1; 
    Lm1 = L;
    Lpm2 = Lpm1; 
    Lpm1 = Lp;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Runge-Kutta for Laguerre Equation
function x = rk2_Lag(t, tn, x, n)
m = 10; 
h = (tn - t)/m;
for j = 1:m
    f1 = (n + .5 - .25*x);
    k1 = -h/( sqrt(f1/x) + .25*(1/x-.25/f1)*sin(2*t) );
    t = t + h;  
    x = x + k1;   
    f1 = (n + .5 - .25*x);
    k2 = -h/( sqrt(f1/x) + .25*(1/x-.25/f1)*sin(2*t) );
    x = x + .5*(k2 - k1);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%% Routine for recurrence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orthonormal associated Laguerre polynomial with positive leading coefficient and derivative
function [pn, pnd] = lagpnRecDer(n,alpha,x)
    pnprev = 0*alpha;
    pn= 1/sqrt(gamma(alpha+1) );
    pndprev = 0*alpha;
    pnd = 0*alpha;
    for k=1:n
        pnold = pn;
        pn = (x -2*k -alpha+1)/sqrt(k*(alpha+k))*pnold-sqrt((k-1+alpha)*(k-1)/k/(k+alpha))*pnprev;
        pnprev = pnold;
        pndold = pnd;
        pnd = (pnold+(x-2*k-alpha+1)*pndold)/sqrt(k*(alpha+k)) -sqrt((k-1+alpha)*(k-1)/k/(alpha+k))*pndprev;
        pndprev = pndold;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%% Routine for explicit expansion %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w] = laguerreExp(n, compRepr, alpha)

if compRepr
    % Get a heuristic for the indices where the weights are about above realmin.
    mn = min(n,ceil(17*sqrt(n)));
else
    mn = n;
end
w = zeros(1, mn);

ibes = max(round(sqrt(n)), 7);
% Heuristics to switch between Bessel, extrapolation and Airy initial guesses.
iair = floor(0.9*n);

% This is a heuristic for the number of terms in the expansions that follow.
T = ceil(50/log(n) );
if ( alpha^2/n > 1 )
    warning('CHEBFUN:lagpts:inputs',['A large alpha may lead to inaccurate ' ...
        'results because the weight is low and R(z) is not close to identity.']);
end
d = 1/(4*n+2*alpha+2);

pt = (4*n -4*transpose((ibes+1):min(mn, iair-1)) +3)*d;
t = pi^2/16*(pt -1).^2;
% This t is not a very good initial approximation of the inverse function 
% of f(y) = (4*n -4*k +3)*d*pi +2*sqrt(y)*sqrt(1-y) -acos(2*y-1);
for it = 1:6
    t = t - (pt*pi +2*sqrt(t-t.^2) -acos(2*t-1) ).*sqrt(t./(1-t))/2;
end
% Another option is to use fzero or initial approximations with barycentric interpolation.

jak = besselroots(alpha, ibes);
bes = jak*0;
wbes = 0*bes;

bulk = t*0;
wbulk = t*0;

ak = [-13.69148903521072; -12.828776752865757; -11.93601556323626;  -11.00852430373326;...
    -10.04017434155809; -9.02265085340981; -7.944133587120853;    -6.786708090071759; ...
    -5.520559828095551; -4.08794944413097; -2.338107410459767];
tair = 3*pi/2*( transpose((mn-iair+1):-1:12)-0.25);
% [DLMF (9.9.6)]
ak = [-tair.^(2/3).*(1 + 5/48./tair.^2 - 5/36./tair.^4 + 77125/82944./tair.^6 ...
    -10856875/6967296./tair.^8); ak(max(1,12-mn+iair-1):11) ];
air = 0*ak;

if (T >= 9)
    % These higher order terms in the left and bulk region are derived in [4]
    bes = bes + (10644*jak.^8 + 60*(887*alpha^2 - 2879)*jak.^6 + (125671*alpha^4 -729422*alpha^2 + 1456807)*jak.^4 ...
        + 3*(63299*alpha^6 - 507801*alpha^4 + 1678761*alpha^2 - 2201939)*jak.^2 ...
        + 2*(107959*alpha^8 - 1146220*alpha^6 + 5095482*alpha^4 -10087180*alpha^2 + 6029959) )*d^8/42525;
    wbes = wbes + (215918*alpha^8 + 53220*jak.^8 + 240*(887*alpha^2 - 2879)*jak.^6 -2292440*alpha^6 ...
        + 3*(125671*alpha^4 - 729422*alpha^2 + 1456807)*jak.^4 + 10190964*alpha^4 ...
        + 6*(63299*alpha^6 - 507801*alpha^4 + 1678761*alpha^2 -2201939)*jak.^2 - 20174360*alpha^2 + 12059918)/42525*d^8;
    
    bulk = bulk + d^7/10886400.*(1-t).^3./t.^3.*(43222750000*(1-t).^(-14) - 241928673000*(1-t).^(-13) + 566519158800*(1-t).^(-12)...
        -714465642135*(1-t).^(-11) + 518401904799*(1-t).^(-10) + 672*(12000*alpha^4 - 24000*alpha^2 +64957561)*(1-t).^(-8) ...
        - 212307298152*(1-t).^(-9) + 24883200*alpha^8 - 192*(103425*alpha^4 -206850*alpha^2 + 15948182)*(1-t).^(-7) ...
        + 3360*(4521*alpha^4 - 9042*alpha^2 - 7823)*(1-t).^(-6) -232243200*alpha^6 - 1792*(3375*alpha^6 - 13905*alpha^4 ...
        + 17685*alpha^2 - 1598)*(1-t).^(-5)+ 16128*(450*alpha^6 - 2155*alpha^4 + 2960*alpha^2 - 641)*(1-t).^(-4) ...
        + 812851200*alpha^4 -768*(70875*alpha^8 - 631260*alpha^6 + 2163630*alpha^4 - 2716980*alpha^2 +555239)*(1-t).^(-3) ...
        + 768*(143325*alpha^8 - 1324260*alpha^6 + 4613070*alpha^4 -5826660*alpha^2 + 1193053)*(1-t).^(-2) - 1028505600*alpha^2 ...
        - 5806080*(15*alpha^8 -140*alpha^6 + 490*alpha^4 - 620*alpha^2 + 127)*(1-t).^(-1) + 210677760);
    %Computing this term for wbulk takes too much memory
end

if (T >= 7)
    bes = bes + (657*jak.^6 +36*jak.^4*(73*alpha^2-181) +2*jak.^2*(2459*alpha^4 -10750*alpha^2 +14051) ...
        + 4*(1493*alpha^6 -9303*alpha^4 +19887*alpha^2 - 12077) )*d^6/2835;
	wbes = wbes + (1493*alpha^6 + 657*jak.^6 + 27*(73*alpha^2 - 181)*jak.^4 - 9303*alpha^4 ...
        + (2459*alpha^4 -10750*alpha^2 + 14051)*jak.^2 + 19887*alpha^2 - 12077)*4/2835*d^6;

	bulk = bulk -d^5/181440.*(1-t).^2./t.^2.*(10797500*(1-t).^(-10) - 43122800*(1-t).^(-9) + 66424575*(1-t).^(-8) ...
        -48469876*(1-t).^(-7) + 193536*alpha^6 + 16131880*(1-t).^(-6) + 80*(315*alpha^4 - 630*alpha^2 -221)*(1-t).^(-4) ...
        - 1727136*(1-t).^(-5) - 967680*alpha^4 - 320*(63*alpha^4 - 126*alpha^2 +43)*(1-t).^(-3) ...
        + 384*(945*alpha^6 - 4620*alpha^4 + 6405*alpha^2 - 1346)*(1-t).^(-2) +1354752*alpha^2 ...
        - 23040*(21*alpha^6 - 105*alpha^4 + 147*alpha^2 - 31)*(1-t).^(-1) -285696);
    wbulk = wbulk - (1-t).^3./90720./t.^3.*d^6.*(43190000*(1-t).^(-12) -204917300*(1-t).^(-11) + 393326325*(1-t).^(-10) ...
        - 386872990*(1-t).^(-9) + 201908326*(1-t).^(-8) +80*(315*alpha^4 - 630*alpha^2 + 53752)*(1-t).^(-6) ...
        - 50986344*(1-t).^(-7) - 320*(189*alpha^4 -378*alpha^2 - 89)*(1-t).^(-5) + 480*(63*alpha^4 - 126*alpha^2 + 43)*(1-t).^(-4) ...
        -384*(315*alpha^6 - 1470*alpha^4 + 1995*alpha^2 - 416)*(1-t).^(-3) + 2304*(21*alpha^6 -105*alpha^4 + 147*alpha^2 - 31)*(1-t).^(-2) );
end
if (T >= 5)
    bes = bes + (11*jak.^4 +3*jak.^2.*(11*alpha^2-19) +46*alpha^4 -140*alpha^2 +94)*d^4/45;
    wbes = wbes + (46*alpha^4 + 33*jak.^4 +6*jak.^2*(11*alpha^2 -19) -140*alpha^2 +94)/45*d^4;
    
    air = air -(15152/3031875*ak.^5+1088/121275*ak.^2)*2^(1/3)*d^(7/3);
    %[Gatteshi 2002 (4.9)], gives an O(n^{-4}) relative error
    
    bulk = bulk +d^3.*(1-t)./t./720.*(1600*(1-t).^(-6) - 3815*(1-t).^(-5) + 480*alpha^4 +2814*(1-t).^(-4) - 576*(1-t).^(-3) ...
        - 960*alpha^2 - 48*(15*alpha^4 - 30*alpha^2 + 7)*(1-t).^(-1) -16*(1-t).^(-2) + 224);
	wbulk = wbulk +(1-t).^2/720./t.^2.*d^4.*(8000*(1-t).^(-8) - 24860*(1-t).^(-7) + 27517*(1-t).^(-6) - 12408*(1-t).^(-5) ...
        + 1712*(1-t).^(-4) +16*(15*alpha^4 - 30*alpha^2 + 7)*(1-t).^(-2) + 32*(1-t).^(-3));
end
if (T >= 3)
    bes = bes + (jak.^2 + 2*alpha^2 - 2)*d^2/3;
    wbes = wbes + (alpha^2 + jak.^2 -1)*2/3*d^2;
    air = air +  ak.^2*(d*16)^(1/3)/5 + (11/35-alpha^2-12/175*ak.^3)*d + ...
        (16/1575*ak+92/7875*ak.^4)*2^(2/3)*d^(5/3);
        
    bulk = bulk -d/12*(12*alpha^2 + 5*(1-t).^(-2) - 4*(1-t).^(-1) - 4);
    wbulk = wbulk  + d^2/6*(2*t + 3)./(t-1).^3;
end
bes = jak.^2*d.*(1 + bes );
air = 1/d +ak*(d/4)^(-1/3) + air;
bulk = bulk + t/d;

w = [ 4*d*bes.^alpha.*exp(-bes)./(besselj(alpha-1, jak)).^2.*(1+ wbes); ...
    bulk.^alpha.*exp(-bulk)*2*pi.*sqrt(t./(1-t)).*(1 +wbulk); ...
    4^(1/3)*air.^(alpha+1/3).*exp(-air)./(airy(1,ak)).^2]';
% For the Airy region, only O(n^{-4}) relative error for x and O(n^{-2/3}) for w as the latter are 
% extremely small or even underflow, and the formula for the expansion would be thrice as long.
x = [ bes; bulk ; air];

if ( min(x) < 0.0 ) || ( max(x) > 4*n + 2*alpha + 2 ) ||  ( min(diff(x)) <= 0.0 ) || (min(w) < 0.0)
    error('Wrong node or weight.')
end

if compRepr
    k = find(w, 1, 'last');
    x = x(1:k)';
    w = w(1:k);
end

end
