function h = conv(f, g, pt)
%CONV   Convolution of BNDFUN objects.
%   H = CONV(F, G) produces the convolution of BNDFUN objects F and G:
%                     - 
%                    /
%           H(x) =   |    F(t) G(x-t) dt,  x in [a + c, b + d]
%                    /
%                   -
%   where domain(F) is [a, b] and domain(G) is [c, d]. The integral is taken
%   over all t for which the integrand is defined: max(a, x - d) <= t <= min(b,
%   x - c).  The breakpoints of H are all pairwise sums of the breakpoints of F
%   and G.
%
%   Note that CONV only supports piecewise-smooth functions on bounded domains.
%
% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
%
% Nick Hale and Alex Townsend, 2014
% Kuan Xu, 2017

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%
% For further details, see Hale and Townsend, "An algorithm for the convolution
% of Legendre series", SIAM Journal on Scientific Computing, Vol. 36, No. 3,
% pages A1207-A1220, 2014.
% 
% For convolution based on all classic orthogonal polynomials, including 
% Chebyshev polynomials, Gegenbauer polynomials, Jacobi polynomials, and 
% Laguerre polynomials, see
%
% [1]. Xu and Loureiro, "Spectral approximation of convolution operator" (2017).
% [2]. Loureiro and Xu, "convolution of Chebyshev polynomials" (2017).
% 
% In the following, it is assumed that the length of the domain of g is greater
% than the length of the domain of f. If this is not the case, then simply
% compute h = conv(g, f), which is equivalent. We assume f and g are polynomials
% of degree M and N, resepectively. If f and g are piecewise-defined, one can
% use the bilinearity of convolution and convolve each each the FUNs
% individually.
%
% The general convolution domain (for smooth functions) is as follows:
%           .___________________________
%          /|                  |      /
%        /  |                  |    /
%      /    |                  |  /
%    /______|__________________|/
%  a+c     b+c                a+d     b+d 
%
% The triangular pieces at end are dealt with using a convolution theorem for
% Legendre polynomials, which leads to a convenient recurrence relation. The
% cost of this is O(m*n) and results in a polynomial of degree m+n. See
% legConv() for details.
%
% Similarly one can show the interior rectangle results in a polynomial of
% degree n. This can be computed by patching with R = floor[(d-c) / (b-a)]
% parallelograms (although triangle Z is not used). Each parallelogram requires
% restricting g to a suitable subdomain, which costs O(n^2) operations. Total
% complexity ratio*(m*n + n*n)
%            ___________________________
%          /       /       /:     /   /
%        /       /       /  : Z /   /          <-- R patches
%      /       /       /    : /   /           
%    /_______/_______/______/__ /
%  a+c     b+c             fl  a+d     b+d
%
% The final piece is computed via further parallelogram subdivision starting
% from the right (B), and a smaller subdivision in both f and g (C).
% Contributions from D and E are discarded, as they have already been counted
% above. Complexity O(m^2 + n^2)
%   ________________
%   : /E|C/:      / 
%   :   /  :    /
%   : /D|B :  /
%   /___|__:/
%  a+c fl a+d     b+d
%
% Rather than make a BNDFUN corresponding to the each of the patches, we
% instead evaluate directly on a corresonding Chebyshev grid of appropriate
% size, which turns out to be far more efficient. 
%
% Total complexity: O( R*(m*n + n*n) + m*m )
%
% [TODO]: It's possible this should be pushed further to the chebtech level.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    h = bndfun();
    return
end

% Extract the domain:
domF = f.domain;
a = domF(1);
b = domF(2);

domG = g.domain;
c = domG(1);
d = domG(2);

% Ensure that g is the signal (i.e., on the larger domain) and f is the filter:
if ( (b - a) > (d - c) )
    h = conv(g, f, pt);
    return
end

% Method options:
if ( strcmpi(pt, 'cheb') )
    pt = 1;
elseif ( strcmpi(pt, 'leg') )
    pt = 2;
else
    error('CHEBFUN:BNDFUN:conv:validFlag', 'Unrecognizable flag.')
end
    
% Useful things:
N = length(g);                               % Length of g
numPatches = floor((d - c) / (b - a));       % Number of patches required
x = chebpts(N, [b+c, a+d], 1);               % Chebyshev grid for interior piece
y = 0*x;                                     % Initialise values in interior
map = @(x, a, b) (x-a)/(b-a) - (b-x)/(b-a);  % Map from [a, b] --> [-1, 1]

fc = get(f, 'coeffs');
if ( pt == 2 )
    fc = cheb2leg(fc);  % Legendre coefficients of f
end

if ( numPatches > 100 )
    error('CHEBFUN:BNDFUN:conv:tooManyPatches', ...
        ['Max number of patches (100) exceeded. ', ...
         'Small interval? Check breakpoints.']);
end

% Restrict g:
doms = c + (b - a)*(0:numPatches);
g_restricted = restrict(g, doms);
if ( ~iscell(g_restricted) )
    % If doms happened to be domain(g), restrict would return a cell.
    g_restricted = {g_restricted};
end

% Loop over the patches:
for k = 1:numPatches      
                          %         _____
    dk = doms([k, k+1]);  %       /|    /
    dk_left  = a + dk(1); %     /  |  /
    dk_mid   = a + dk(2); %   /____|/
    dk_right = b + dk(2); %  dkl  dkm   dkr
    gk = g_restricted{k};                          % g on this subdomain
    gk = simplify(gk);                             % Simplify for efficiency
    
    gkc = get(gk, 'coeffs');
    
    if ( pt == 1 )
        hLc = chebConv(fc, gkc, 1);      % Convolution on this domain
        hRc = chebConv(fc, gkc, 0);
    elseif ( pt == 2 )
        gkc = cheb2leg(gkc);  % Legendre coefficients
        [hLc, hRc] = legConv(fc, gkc);      % Convolution on this domain
    end
    
    
    % The left triangle for the kth patch:
    ind = (dk_left <= x) & (x < dk_mid); % Locate the grid values in [dkl, dkr]:
    if ( k == 1 ) % First piece:
        if ( pt == 2 )
            hLc = leg2cheb(hLc);           % Cheb. coeffs of left tri.
        end
        data.domain = [dk_left, dk_mid];
        h_left = bndfun({[], hLc}, data);        % Make BNDFUN from coeffs
    else          % Subsequent left pieces
        z = map(x(ind), dk_left, dk_mid);          % Map grid points to [-1, 1]
        if ( pt == 1 )
            tmp = clenshawCheb(z, hLc);          % Evaluate via recurrence
        elseif ( pt == 2 )
            tmp = clenshawLeg(z, hLc);          % Evaluate via recurrence
        end
        y(ind) = y(ind) + tmp;                     % Append
    end
    
    % The right triangle for the kth patch:
    if ( k < numPatches )                          % Not needed for final patch!
        % Locate the grid values in [dkl, dkr]:
        ind = (dk_mid <= x) & (x < dk_right);
        z = map(x(ind), dk_mid, dk_right);
        if ( pt == 1 )
            tmp = clenshawCheb(z, hRc);          
        elseif ( pt == 2 )
            tmp = clenshawLeg(z, hRc);
        end
        y(ind) = y(ind) + tmp;
    end
end

if ( abs((b-a)-(d-c)) < 10*eps(norm([a b c d], inf)) )
    % If there's only one patch, then we already have all the information reqd.
    
    if ( pt == 2)
        hRc = leg2cheb(hRc);       % Cheb coeffs of right tri.
    end
    data.domain = d + [a b];
    h_right = bndfun({[], hRc}, data);   % Make BNDFUN from coeffs
    h_mid = bndfun();
    
else  
    % Final right right triangle:
    %  ________________
    %  : /E|C/:      / 
    %  :   /  : A  /
    %  : /D|B :  /
    %  /___|__:/
    %     fl a+d     b+d

    finishLocation = a + c + numPatches*(b - a);    % Where patches got to. (fl) 
    gk = restrict(g, d-[(b-a) 0]);                  % g on appropriate domain   
    gk = simplify(gk);                              % Simplify for efficiency
    
    gkc = get(gk, 'coeffs');
    if ( pt == 1 )
        hLc = chebConv(fc, gkc, 1);       % Conv on A and B
        hRc = chebConv(fc, gkc, 0);       % Conv on A and B
    elseif ( pt == 2 )
        gkc = cheb2leg(gkc);           % Legendre coeffs
        [hLc, hRc] = legConv(fc, gkc);       % Conv on A and B
        hRc = leg2cheb(hRc);                % Cheb coeffs on A
    end
    
    data.domain = [d+a, d+b];
    h_right = bndfun({[], hRc}, data);            % Make BNDFUN from coeffs
    
    % Remainder piece: (between fl and a+d)
    remainderWidth = d + a - finishLocation; % b+d-fl-(b-a)
    if ( remainderWidth > 0 )
        ind = finishLocation <= x;           % Discard D and E
        
        % B: (Coeffs were computed above)
        z = map(x(ind), d - b + 2*a, d + a); % Map grid to [-1, 1]
        if ( pt == 1 )
            tmp = clenshawCheb(z, hLc);    % Evaluate via recurrence
        elseif ( pt == 2 )
            tmp = clenshawLeg(z, hLc);    % Evaluate via recurrence
        end
        y(ind) = tmp;                        % Store
        
        % C: 
        domfk = b + [-remainderWidth, 0];               % Domain of fk
        domfk(1) = max(domfk(1), f.domain(1));          % Ensure domfk is a
        domfk(end) = min(domfk(end), f.domain(end));    %  valid subdomain
        fk = restrict(f, domfk);                        % Restrict f
        fk = simplify(fk);                              % Simplify f
        
        fkc = get(fk, 'coeffs');
        if ( pt == 2 )
            fkc = cheb2leg(fkc);   % Legendre coeffs
        end
        
        domgk = [finishLocation, d + a] - b;            % Domain of gk
        domgk(1) = max(domgk(1), g.domain(1));          % Ensure domgk is a
        domgk(end) = min(domgk(end), g.domain(end));    %  valid subdomain
        gk = restrict(g, domgk);                        % Restrict g
        gk = simplify(gk);                              % Simplify g
        
        gkc = get(gk, 'coeffs');
        if ( pt == 1 )
            hRc = chebConv(fkc, gkc, 0);    % Conv
        elseif ( pt == 2 )
            gkc = cheb2leg(gkc);   % Legendre coeffs
            [ignored, hRc] = legConv(fkc, gkc);    % Conv
        end
        z = map(x(ind), finishLocation, d + a);         % Map to [-1, 1]
        
        if ( pt == 1 )
            tmp = clenshawCheb(z, hRc);               % Eval via recurrence
        elseif ( pt == 2 )
            tmp = clenshawLeg(z, hRc);               % Eval via recurrence
        end
        y(ind) = y(ind) + tmp*remainderWidth/(b - a);   % Scale and append
    end
    % Convert values to coeffs (we don't want to construct a chebtech1)
    y = chebtech1.vals2coeffs(y);
    % Construct BNDFUN of the interior (rectangle) using coefficients:
    data.domain = [b+c, a+d];
    h_mid = bndfun({[], y}, data);
    
end

% h_mid can be empty so return the three or two pieces as a cell array:
if ( isempty(h_mid) )
    h = {h_left*(b-a)/2, h_right*(b-a)/2};
else    
    h = {h_left*(b-a)/2, h_mid*(b-a)/2, h_right*(b-a)/2};
end

% Ensure endpoints of BNDFUNs for adjacent subintervals match exactly.
%
% NB:  This is a bit inefficient if consecutive endpoints get adjusted, as some
% BNDFUNs may have their maps changed twice, but the map-change operation is
% fast enough that this shouldn't matter.
for (n = 2:1:numel(h))
    end_left = h{n-1}.domain(end);
    end_right = h{n}.domain(1);
    if ( end_left ~= end_right )
        hs_left = norm(h{n-1}.domain, Inf);  % hscale for left piece.
        hs_right = norm(h{n}.domain, Inf);   % hscale for right piece.

        % If there's a mismatch, it should be small because the BNDFUNs should
        % come out in order corresponding to consecutive adjacent subintervals.
        if ( abs(end_left - end_right) < 2*eps*max(hs_left, hs_right) )
            new_end = (end_left + end_right)/2;
            h{n-1} = changeMap(h{n-1}, [h{n-1}.domain(1) new_end]);
            h{n} = changeMap(h{n}, [new_end h{n}.domain(end)]);
        else
            % We should only get here if a programmer error made elsewhere
            % causes the BNDFUNs to get out of order.
            error('CHEBFUN:BNDFUN:conv:nonConsecutive', ...
                'Pieces do not belong to consecutive subintervals.');
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gammaL, gammaR] = legConv(alpha, beta)
% Convolution using Legendre expansions and the analoguous convolution theorem.
% See Hale and Townsend, "An algorithm for the convolution of Legendre series",
% SIAM Journal on Scientific Computing, Vol. 36, No. 3, pages A1207-A1220, 2014.

% Better computational efficiency is achieved when g has the lower degree:
if ( length(beta) > length(alpha) )
    tmp = alpha;
    alpha = beta;
    beta = tmp;
end

% Maximum degree of result:
MN = length(alpha) + length(beta);

% Pad to make length n + 1.
alpha = [ alpha ; zeros(MN - length(alpha), 1) ];

% S represents multiplication by 1/z in spherical Bessel space:
e = [[1 ; 1./(2*(1:(MN-1)).'+1)], [1 ; zeros(MN-1, 1)], -1./(2*(0:MN-1).'+1)];
S = spdiags(e, -1:1, MN, MN);

gammaL = rec(S, alpha, beta, -1); % Chebyshev coeffs for the left piece
S(1,1) = -1;                      % Update S
gammaR = rec(S, -alpha, beta, 1); % Chebyshev coeffs for the right piece

%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX FREE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function gamma = rec(S, alpha, beta, sgn)
        % Compute the Legendre coefficients of the convolution on L/R piece.
        % See Theorem 4.1 of paper.
        
        % Initialise scl:
        N = length(beta);
        scl = 1./(2*(1:N).'-1);
        scl(2:2:end) = -scl(2:2:end);
        
        % First column of B:
        vNew = S*alpha;
        v = vNew;
        gamma = beta(1)*vNew;
        beta_scl = scl.*beta;
        beta_scl(1) = 0;
        gamma(1) = gamma(1) + vNew(1:N).'*beta_scl;
        
        % The scalar case is trivial:
        if ( length(beta) == 1 )
            return
        end
        
        % Second column of B:
        vNew = S*v + sgn*v;
        vOld = v;
        v = vNew;
        vNew(1) = 0;

        gamma = gamma + beta(2)*vNew;
        beta_scl = -beta_scl*((2 - 0.5)/(2 - 1.5));
        beta_scl(2) = 0;
        gamma(2) = gamma(2) + vNew(1:N).'*beta_scl;
        
        % Loop over remaining columns using recurrence:
        for n = 3:N
            vNew = (2*n-3) * (S * v) + vOld; % Recurrence
            vNew(1:n-1) = 0;                 % Zero terms 
            gamma = gamma + vNew*beta(n);    % Append to g
            
            % Recurrence is unstable for j < k. Correct for upper-tri part:
            beta_scl = -beta_scl*((n-.5)/(n-1.5));
            beta_scl(n) = 0;
            gamma(n) = gamma(n) + vNew(1:N).'*beta_scl;
            
            vOld = v;
            v = vNew;
            
        end
        
        ag = abs(gamma);
        mg = max(ag);
        loc = find(ag > eps*mg, 1, 'last');
        gamma = gamma(1:loc);

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = clenshawLeg(x, alpha) 
% Evaluate a Legendre expansion with coefficient alpha at x. 

n = length(alpha); 
b_old = 0; 
b_cur = 0; 
for k = (n-1):-1:1
  b_new = alpha(k+1) + (2*k + 1)/(k + 1)*x.*b_cur - (k + 1)/(k + 2)*b_old;
  b_old = b_cur; 
  b_cur = b_new; 
end
val = alpha(1) + x.*b_cur - .5*b_old; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ch = chebConv(cf, cg, lr)
% Convolution using Chebyshev expansions based a unified framework for all 
% classic orthogonal polynomials.
%
% See [1]. Xu and Loureiro, "Spectral approximation of convolution operator" (2017),
% and [2]. Loureiro and Xu, "Convolution of Chebyshev polynomials" (2017).

% Lengths:
M = length(cf)-1;
N = length(cg)-1;

% Better computational efficiency is achieved when g has the lower degree:
if ( (N < M && N > 1) || (M < N && M < 2) )
    tmp = cf;
    cf = cg;
    cg = tmp;
    tmp = M;
    M = N;
    N = tmp;
end

ch = zeros(M+N+2, 1);
v = ones(1, M+1);
v(1:2:end) = -1;
cr = zeros(2*M+5, 1);

%% First three columns:

c0 = int(cf, lr);
ch(1:M+2) = cg(1)*c0;
if ( N == 0 )
    return
end
cr(1) = c0(2);

c1 = timesx(c0) + (-1)^lr*[c0; 0] - int(timesx(cf), lr);
ch(1:M+3) = ch(1:M+3) + cg(2)*c1(1:M+3);
if ( N == 1 )
    return
end

c0 = int(cf, 2);
c2 = (-1)^(lr+1)*[c0; 0; 0] + 4*int(c1, 2);
ch(3:M+4) = ch(3:M+4) + cg(3)*c2(3:M+4);
if ( N == 2 )
    ch(1:2) = ch(1:2) + cg(3)*c2(1:2);
    return
end

c1 = [c1(2:M+3); 0; 0];
c2 = [c2(3:M+4); 0; 0];
c3 = zeros(M+4, 1);

cr(2:3) = c1(1:2);
cr(4:5) = c2(1:2);

%% Recursion below the diagonal

for n = 2:M-1  % n is the true order index
    m = n + 2; % m is the matrix index
    c3(1:M+2) = (1+2/(n-1))*c1(3:M+4) + ...
        (n+1)*(c2(1:M+2)-c2(3:M+4))./(((n+1):(n+M+2)).');
    c3(1:M-m+3) = c3(1:M-m+3) + (2*(-1)^(lr*n)/(n-1))*c0(m:M+2);
    
    ch(m:m+M+1) = ch(m:m+M+1) + cg(m)*c3(1:M+2);
    
    cr(2*n+2:2*n+3) = c3(1:2);
    c1 = c2;
    c2 = c3;
end

r = zeros(2, N+2);

for n = M:M+1 % n is the true order index
    
    % Recursion:
    m = n + 2; % m is the matrix index
    c3(1:M+2) = (1+2/(n-1))*c1(3:M+4) + ...
        (n+1)*(c2(1:M+2)-c2(3:M+4))./(((n+1):(n+M+2)).');
    c3(1:M-m+3) = c3(1:M-m+3) + (2*(-1)^(lr*n)/(n-1))*c0(m:M+2);
    
    if ( n < N )
        ch(m:m+M+1) = ch(m:m+M+1) + cg(m)*c3(1:M+2);
    end
    
    cr(2*n+2:2*n+3) = c3(1:2);
    c1 = c2;
    c2 = c3;
    
    % Reflection:
    r(M+2-n, 3:M+3) = v.*(m:m+M).*c3(2:M+2).'/(n+1);
    ch(m) = ch(m) + sum(r(M+2-n, 3:N-n+1).'.*cg(m+1:N+1));
    
end

cr(2*(M+2)+1) = [];

r(1, 1:2) = cr(2*(M+1)+1:2*(M+2));
r(2, 1:2) = cr(2*M+1:2*(M+1));

for n = M+2:N-M-2  % n is the true order index
    
    % Recursion:
    m = n + 2; % m is the matrix index
    c3(1:M+2) = (1+2/(n-1))*c1(3:M+4) + ...
        (n+1)*(c2(1:M+2)-c2(3:M+4))./(((n+1):(n+M+2)).');
    ch(m:m+M+1) = ch(m:m+M+1) + cg(m)*c3(1:M+2);
    c1 = c2;
    c2 = c3;
    
    % Reflection:
    %     m = n + 1; % m is the matrix index
    ch(m) = ch(m) + (v.*(m:m+M).*c3(2:M+2).'/(n+1))*cg(m+1:m+M+1);
end

for n = max(M+2, N-M-1):N-1  % n is the true order index
    
    % Recursion:
    m = n + 2; % m is the matrix index
    c3(1:M+2) = (1+2/(n-1))*c1(3:M+4) + ...
        (n+1)*(c2(1:M+2)-c2(3:M+4))./(((n+1):(n+M+2)).');
    ch(m:m+M+1) = ch(m:m+M+1) + cg(m)*c3(1:M+2);
    c1 = c2;
    c2 = c3;
    
    % Reflection:
    %     m = n + 1; % m is the matrix index
    ch(m) = ch(m) + (v(1:N-n-1).*(m:N).*c3(2:N-n).'/(n+1))*cg(m+1:N+1);
end

%% Recursion for the top rows

r1 = r(1, :);
r2 = r(2, :);
r3 = zeros(1, N+2);

for m = M+1:-1:2
    n = m:m+N-1;
    r3(1:2) = cr(2*(m-1)-1:2*(m-1));
    r3(3:N+2) = m*r2(1:N)./(1-n) + m*r2(3:N+2)./(n+1) + r1(1:N)...
        -2*m*(-1).^(lr*n).*c0(m+1)./(n.^2-1);
    
    ch(m) = ch(m) + r3(3:N-m+3)*cg(m+1:N+1);
    r1 = r2;
    r2 = r3;
end

n = 2:N;
ch(1) = ch(1) + (r2(n)./(1-n) + r2(4:N+2)./(n+1) + r1(n)...
    -2*(-1).^(lr*n).*c0(2)./(n.^2-1))*cg(3:N+1)/2;


%% integration:

    function b = int(c, mode)
        
        k = length(c);
        c = [c; 0; 0];          % Pad with zeros
        b = zeros(k+1, 1);   % Initialize vector b = {b_r}
        
        % Compute b_(2) ... b_(n+1):
        b(3:k+1) = (c(2:k)-c(4:k+2))./(2*(2:k).');
        b(2) = c(1)-c(3)/2;        % Compute b_1
        
        if ( mode == 1 )
            w = ones(1, k);
            w(2:2:end) = -1;
            b(1) = w*b(2:end);             % Compute b_0 (satisfies f(-1) = 0)
        elseif ( mode == 0 )
            b(1) = sum(b(2:end));
            b(2:end) = -b(2:end);
        elseif( mode == 2 )
            b(1) = -sum(b(2:end));
        end
        
    end

%% product with x:

    function b = timesx(c)
        
        j = length(c);
        c = [c; 0; 0];       % Pad with zeros
        b = zeros(j+1, 1);   % Initialize vector b = {b_r}
        
        % Compute b_(2) ... b_(j+1):
        b(3:j+1) = (c(2:j)+c(4:j+2))/2;
        b(2) = c(1)+c(3)/2;  % Compute b_1
        b(1) = c(2)/2;
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = clenshawCheb(x, alpha)
% Evaluate a Chebyshev expansion with coefficient alpha at x.

bk1 = 0*x;
bk2 = bk1;
x = 2*x;
n = size(alpha,1)-1;
for k = (n+1):-2:3
    bk2 = alpha(k) + x.*bk1 - bk2;
    bk1 = alpha(k-1) + x.*bk2 - bk1;
end
if ( mod(n, 2) )
    [bk1, bk2] = deal(alpha(2) + x.*bk1 - bk2, bk1);
end
val = alpha(1) + .5*x.*bk1 - bk2;

end
