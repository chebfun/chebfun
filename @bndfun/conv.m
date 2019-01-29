function h = conv(f, g)
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%
% For further details, see Hale and Townsend, "An algorithm for the convolution
% of Legendre series", SIAM Journal on Scientific Computing, Vol. 36, No. 3,
% pages A1207-A1220, 2014.
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
% EASYCONV() for details.
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
% However, in the case when R is too big we revert to a Clenshaw-Curtis
% quadrature-based approach to compute the convolution in the inner rectangle.
% This approach has a complexity O( (m + n)^3 ). We naively compare this with
% the big-O term above to determine which approach to use.

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
    h = conv(g, f);
    return
end
    
% Useful things:
M = length(f);                               % Length of F
N = length(g);                               % Length of g
numPatches = floor((d - c) / (b - a));       % Number of patches required
x = chebpts(N, [b+c, a+d], 1);               % Chebyshev grid for interior piece
y = 0*x;                                     % Initialise values in interior
map = @(x, a, b) (x-a)/(b-a) - (b-x)/(b-a);  % Map from [a, b] --> [-1, 1]

% If there are too many patches then the HT approach is too slow. In such a case
% we resort to the standard quadrature-based approach (but still use the HT
% approach for the two triangular domains at the ends).
coeffsConvCost = numPatches*((M+N)*N);
quadConvCost = (M+N)^3;
if ( numPatches > 1 && coeffsConvCost > quadConvCost )
    h_left = conv(f, restrict(g, c+[0, b-a]));     % Left triangle
    h_right = conv(f, restrict(g, d-[(b-a) 0]));   % Right tirangle
    % Middle:
    [t, w] = legpts((M+N+5)/2, [a,b]);             % Legendre grid
    [tt, xx] = meshgrid(t,x);                      % Cheb/Leg grid to evaluate g
    ft = feval(f,t);                               % Evaluate f
    gxmt = feval(g, xx-tt);                        %    and g (expensive)
    y = gxmt*(w'.*ft);                             % Compute integral
    y = chebtech1.vals2coeffs(y);                  % Convert values to coeffs 
    % Trim small coefficients:
    ay = abs(y); my = max(ay); loc = max(find(ay > 10*eps*my, 1, 'last'),1);
    if ( isempty(loc) ), loc = 1; end              % Deal with case when y = 0
    y = y(1:loc);
    data.domain = [b+c, a+d];
    h_mid = bndfun({[],y}, data);                  % Make h_mid from coeffs
    h = {h_left{1}, h_mid, h_right{end}};          % Combine three pieces
    h = fixMaps(h);                                % Ensure domain ends match
    return
end

% Trim small coefficients in f:
f_cheb = get(f, 'coeffs'); af = abs(f_cheb); mf = max(af); 
loc = find(af > 10*eps*mf, 1, 'last'); if ( isempty(loc) ), loc = 1; end
f_cheb = f_cheb(1:loc);
f_leg = cheb2leg(f_cheb);                   % Legendre coefficients of f

% Restrict g:
doms = c + (b - a)*(0:numPatches);
g_restricted = restrict(g, doms);
doms = c + (b - a)*(0:numPatches);           % Subdomains
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
    gk_leg = cheb2leg(get(gk, 'coeffs'));          % Legendre coefficients
    [hLegL, hLegR] = easyConv(f_leg, gk_leg);      % Convolution on this domain
        
    % The left triangle for the kth patch:
    ind = (dk_left <= x) & (x < dk_mid); % Locate the grid values in [dkl, dkr]:
    if ( k == 1 ) % First piece:
        hLegL = leg2cheb(hLegL);                   % Cheb. coeffs of left tri.
        data.domain = [dk_left, dk_mid];
        h_left = bndfun({[], hLegL}, data);        % Make BNDFUN from coeffs
    else          % Subsequent left pieces
        z = map(x(ind), dk_left, dk_mid);          % Map grid points to [-1, 1]
        tmp = clenshawLegendre(z, hLegL);          % Evaluate via recurrence
        y(ind) = y(ind) + tmp;                     % Append
    end
    
    % The right triangle for the kth patch:
    if ( k < numPatches )                          % Not needed for final patch!
        % Locate the grid values in [dkl, dkr]:
        ind = (dk_mid <= x) & (x < dk_right);
        z = map(x(ind), dk_mid, dk_right);
        tmp = clenshawLegendre(z, hLegR);
        y(ind) = y(ind) + tmp;
    end
    
end

if ( abs((b-a)-(d-c)) < 10*eps(norm([a b c d], inf)) )
    % If there's only one patch, then we already have all the information reqd.
    hLegR = leg2cheb(hLegR);                        % Cheb coeffs of right tri.
    data.domain = d + [a b];
    h_right = bndfun({[], hLegR}, data);            % Make BNDFUN from coeffs
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
    gk_leg = cheb2leg(get(gk, 'coeffs'));           % Legendre coeffs
    [hLegL, hLegR] = easyConv(f_leg, gk_leg);       % Conv on A and B
    hLegR = leg2cheb(hLegR);                        % Cheb coeffs on A
    data.domain = [d+a, d+b];
    h_right = bndfun({[], hLegR}, data);            % Make BNDFUN from coeffs

    % Remainder piece: (between fl and a+d)
    remainderWidth = d + a - finishLocation; % b+d-fl-(b-a)
    if ( remainderWidth > 0 )
        ind = finishLocation <= x;           % Discard D and E

        % B: (Coeffs were computed above)
        z = map(x(ind), d - b + 2*a, d + a); % Map grid to [-1, 1]
        tmp = clenshawLegendre(z, hLegL);    % Evaluate via recurrence
        y(ind) = tmp;                        % Store

        % C: 
        domfk = b + [-remainderWidth, 0];               % Domain of fk
        domfk(1) = max(domfk(1), f.domain(1));          % Ensure domfk is a
        domfk(end) = min(domfk(end), f.domain(end));    %  valid subdomain
        fk = restrict(f, domfk);                        % Restrict f
        fk = simplify(fk);                              % Simplify f
        fk_leg = cheb2leg(get(fk, 'coeffs'));           % Legendre coeffs
        domgk = [finishLocation, d + a] - b;            % Domain of gk
        domgk(1) = max(domgk(1), g.domain(1));          % Ensure domgk is a
        domgk(end) = min(domgk(end), g.domain(end));    %  valid subdomain
        gk = restrict(g, domgk);                        % Restrict g
        gk = simplify(gk);                              % Simplify g
        gk_leg = cheb2leg(get(gk, 'coeffs'));           % Legendre coeffs
        [~, hLegR] = easyConv(fk_leg, gk_leg);          % Conv 
        z = map(x(ind), finishLocation, d + a);         % Map to [-1, 1]
        tmp = clenshawLegendre(z, hLegR);               % Eval via recurrence
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
h = fixMaps(h);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = fixMaps(h)
% Ensure endpoints of BNDFUNs for adjacent subintervals match exactly. NB: This
% is a bit inefficient if consecutive endpoints get adjusted, as some BNDFUNs
% may have their maps changed twice, but the map-change operation is fast enough
% that this shouldn't matter.
for n = 2:1:numel(h)
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

function [gammaL, gammaR] = easyConv(alpha, beta)
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

function val = clenshawLegendre(x, alpha) 
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
