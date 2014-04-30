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
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
%
% Nick Hale and Alex Townsend, 2014

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Devoloper note:
%
% For further details, see Hale and Townsend, "The convolution of compactly
% supported functions", (To appear in SISC)
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
%   : ??/  :    /
%   : /D|B :  /
%   /___|__:/
%  a+c fl a+d     b+d
%
% Rather than make a BNDFUN corresponding to the each of the patches, we
% instead evaluate directly on a corresonding Chebyshev grid of appropriate
% size, which turns out to be far more efficient. 
%
% Total complexity: O( R*(m*n + n*n) + m*m )
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
N = length(g);                               % Length of g
numPatches = floor((d - c) / (b - a));       % Number of patches required
x = chebpts(N, [b+c, a+d], 1);               % Chebyshev grid for interior piece
y = 0*x;                                     % Initialise values in interior
map = @(x, a, b) (x-a)/(b-a) - (b-x)/(b-a);  % Map from [a, b] --> [-1, 1]
f_leg = chebtech.cheb2leg(get(f, 'coeffs')); % Legendre coefficients of f

% Restrict g:
doms = c + (b - a)*(0:numPatches);
g_restricted = restrict(g, doms);
if ( ~iscell(g_restricted) )
    % If doms happened to be domain(g), restrict would return a cell.
    g_restricted = {g_restricted};
end

% Loop over the patches:
for k = 1:numPatches      
    
    dk = doms([k, k+1]);  %       /|????/
    dk_left  = a + dk(1); %     /  |  /
    dk_mid   = a + dk(2); %   /____|/
    dk_right = b + dk(2); %  dkl  dkm   dkr
    gk = g_restricted{k};                          % g on this subdomain
    gk = simplify(gk);                             % Simplify for efficiency
    gk_leg = chebtech.cheb2leg(get(gk, 'coeffs')); % Its Legendre coefficients
    [hLegL, hLegR] = easyConv(f_leg, gk_leg);      % Convolution on this domain
    
    % The left triangle for the kth patch:
    ind = (dk_left <= x) & (x < dk_mid); % Locate the grid values in [dkl, dkr]:
    if ( k == 1 ) % First piece:
        hLegL = chebtech.leg2cheb(flipud(hLegL));  % Cheb. coeffs of left tri.
        h_left = bndfun({[], hLegL}, [dk_left, dk_mid]); % Make BNDFUN from coeffs
    else          % Subsequent left pieces
        z = map(x(ind), dk_left, dk_mid);          % Map grid points to [-1, 1]
        tmp = clenshawLegendre(z, hLegL);          % Evaluate via recurrence
        y(ind) = y(ind) + tmp;                     % Append
    end
    
    % The right triangle for the kth patch:
    if ( k < numPatches )                         % Not needed for final patch!
        % Locate the grid values in [dkl, dkr]:
        ind = (dk_mid <= x) & (x < dk_right);
        z = map(x(ind), dk_mid, dk_right);
        tmp = clenshawLegendre(z, hLegR);
        y(ind) = y(ind) + tmp;
    end
end

if ( numPatches == 1 )
    % If there's only one patch, then we already have all the information reqd.
    hLegR = chebtech.leg2cheb(flipud(hLegR));       % Cheb coeffs of right tri.
    h_right = bndfun({[], hLegR}, d + [a, b]); % Make BNDFUN from coeffs
    h_mid = bndfun();
    
else  
    % Final right right triangle:
    %  ________________
    %  : /E|C/:      / 
    %  : ??/  : A  /
    %  : /D|B :  /
    %  /___|__:/
    %     fl a+d     b+d

    finishLocation = a + c + numPatches*(b - a);    % Where patches got to. (fl) 
    gk = restrict(g, d-[(b-a) 0]);                  % g on appropriate domain   
    gk = simplify(gk);                              % Simplify for efficiency
    gk_leg = chebtech.cheb2leg(get(gk, 'coeffs'));  % Legendre coeffs
    [hLegL, hLegR] = easyConv(f_leg, gk_leg);       % Conv on A and B
    hLegR = chebtech.leg2cheb(flipud(hLegR));       % Cheb coeffs on A
    h_right = bndfun({[], hLegR}, [d+a, d+b]);      % Make BNDFUN from coeffs
    
    % Remainder piece: (between fl and a+d)
    remainderWidth = d + a - finishLocation; % b+d-fl-(b-a)
    if ( remainderWidth > 0 )
        ind = finishLocation <= x;           % Discard D and E
        
        % B: (Coeffs were computed above)
        z = map(x(ind), d - b + 2*a, d + a); % Map grid to [-1, 1]
        tmp = clenshawLegendre(z, hLegL);    % Evaluate via recurrence
        y(ind) = tmp;                        % Store
        
        % C: 
        fk = restrict(f, b + [-remainderWidth, 0]);     % Restrict f
        fk = simplify(fk);                              % Simplify f
        fk_leg = chebtech.cheb2leg(get(fk, 'coeffs'));  % Legendre coeffs
        gk = restrict(g, [finishLocation, d + a] - b);  % Restrict g
        gk = simplify(gk);                              % Simplify g
        gk_leg = chebtech.cheb2leg(get(gk, 'coeffs'));  % Legendre coeffs
        [ignored, hLegR] = easyConv(fk_leg, gk_leg);    % Conv 
        z = map(x(ind), finishLocation, d + a);         % Map to [-1, 1]
        tmp = clenshawLegendre(z, hLegR);               % Eval via recurrence
        y(ind) = y(ind) + tmp*remainderWidth/(b - a);   % Scale and append
    end
    
    % Convert values to coeffs (we don't want to construct a chebtech1)
    y = chebtech1.vals2coeffs(y);
    % Construct BNDFUN of the interior (rectangle) using coefficients:
    h_mid = bndfun({[], y}, [b+c, a+d]);
    
end

% h_mid can be empty so return the three or two pieces as a cell array:
if ( isempty(h_mid) )
    h = {h_left*(b-a)/2, h_right*(b-a)/2};
else    
    h = {h_left*(b-a)/2, h_mid*(b-a)/2, h_right*(b-a)/2};
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gammaL, gammaR] = easyConv(alpha, beta)
% Convolution using Legendre expansions and the analoguous convolution theorem.
% See Hale and Townsend, "The convolution of compactly supported functions", (In
% prep)

% TODO: Document the simple case a little more fully? (Wait until paper done).

% Better computational efficiency is achieved when g has the lower degree:
if ( length(beta) > length(alpha) )
    tmp = alpha;
    alpha = beta;
    beta = tmp;
end

% Flip, as per convention:
alpha = flipud(alpha);
beta = flipud(beta);

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
        % TODO: Document further once paper is complete.
        
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