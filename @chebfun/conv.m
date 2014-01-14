function h = conv(f, g)
%CONV   Convolution of CHEBFUN objects.
% H = CONV(F, G) produces the convolution of CHEBFUN objects F and G:
%                   - 
%                  /
%         H(x) =   |    F(t) G(x-t) dt,  x in [a + c, b + d]
%                  /
%                 -
% where domain(F) is [a, b] and domain(G) is [c, d]. The integral is taken over
% all t for which the integrand is defined: max(a, x - d) <= t <= min(b, x - c).
% The breakpoints of H are all pairwise sums of the breakpoints of F and G.
%
% Note that CONV only supports piecewise-smooth functions on bounded domains.
%
% Example:
%   f = chebfun(1/2); g = f;
%   subplot(2, 2, 1), plot(f)
%   for j = 2:4, g = conv(f, g); subplot(2, 2, j), plot(g), end
%   figure, for j = 1:4, subplot(2,2,j), plot(g), g = diff(g); end

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
%
% Nick Hale and Alex Townsend, 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Devoloper note:
%
% For further details, see Hale and Townsend, "The convolution of compactly
% supported functions", (In prep)
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
% from the right (B), and a smaller subdivisin in both f and g (C).
% Contributions from D and E are discarded, as they have already been counted
% above. Complexity O(m^2 + n^2)
%   ________________
%   : /E|C/:      / 
%   : ??/  :    /
%   : /D|B :  /
%   /___|__:/
%  a+c fl a+d     b+d
%
% Rather than make a CHEBFUN corresponding to the each of the patches, we
% instead evaluate directly on a corresonding Chebyshev grid of appropriate
% size, which turns out to be far more efficient. 
%
% Total complexity: O( R*(m*n + n*n) + m*m )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Force size(B) = [n, min(m, n)] ifor the interior convolutions?
% TODO: Refactor so that most of this lives at the FUN and/or CHEBTECH level?

% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    h = chebfun();
    return
end

% No support for quasimatrices:
if ( min(size(f)) > 1 || min(size(g)) > 1 )
    error('CHEBFUN:conv:quasi', 'No support for array-valued CHEBFUN objects.');
end

% Check transpose state:
if ( xor(f(1).isTransposed, g(1).isTransposed) )
    error('CHEBFUN:conv:transposed', 'CHEBFUN dimensions do not agree.');
end

% Extract the domain:
[a, b] = domain(f);
[c, d] = domain(g);

% Ensure that g is the signal (i.e., on the larger domain) and f is the filter:
if ( b-a > d-c )
    h = conv(g, f);
    return
end
    
% Deal with piecewise CHEBFUN objects.
if ( numel(f.funs) > 1 || numel(g.funs) > 1 )
    h = 0;
    % Loopover each of the interactions:
    for j = 1:numel(f.funs)
        for k = 1:numel(g.funs)
            % TODO: Tidy this!
            h = myplus(h, conv(chebfun(f.funs(j)), chebfun(g.funs(k))));
        end
    end
    return
else
    f = f.funs{1};
    g = g.funs{1};
end
% Note, for simplicity we work with the FUNs, rather than the CHEBFUNs.

% TODO: Replace this with a call to CONV_OLD()?
if ( ~isa(f.onefun, 'chebtech') || ~isa(g.onefun, 'chebtech') )
    error('CHEBFUN:conv:notsmooth', ...
        'CONV() only supports smooth CHEBFUN objects on bounded domains.')
end

% Useful things..
m = length(f); n = length(g);                % Lengths of f anf g
numPatches = floor((d - c) / (b - a));       % Number of patches required
x = chebpts(n, [b + c, a + d], 1);           % Chebyshev grid for interior piece
y = 0*x;                                     % Initialise values in interior
map = @(x, a, b) (x-a)/(b-a) - (b-x)/(b-a);  % Map from [a, b] --> [-1, 1]
f_leg = chebtech.cheb2leg(get(f, 'coeffs'));  % Legendre coefficients of f


% Restrict g:
doms = c + (b-a)*(0:numPatches);
g_restricted = restrict(g, doms);
if ( ~iscell(g_restricted) )
    g_restricted = {g_restricted};
end

% Loop over the patches:
for k = 1:numPatches      
    
    dk = doms([k, k+1]);  %       /|????/
    dk_left  = a + dk(1); %     /  |  /
    dk_mid   = a + dk(2); %   /____|/
    dk_right = b + dk(2); %  dkl  dkm   dkr
    gk = g_restricted{k};                          % g on this subdomain
    gk_leg = chebtech.cheb2leg(get(gk, 'coeffs')); % Its Legendre coefficients
    [hLegL, hLegR] = easyConv(f_leg, gk_leg);      % Convolution on this domain
    
    % The left triangle for the kth patch:
    idx = dk_left <= x & x < dk_mid; % Locate the grid values in [dkl, dkr]:
    if ( k == 1 ) % First piece:
        hLegL = chebtech.leg2cheb(flipud(hLegL));  % Chebyshev coeffs of left tri.
        h_left = chebfun(hLegL, [dk_left, dk_mid], 'coeffs'); % Make CHEBFUN
    else          % Subsequent left pieces
        z = map(x(idx), dk_left, dk_mid);          % Map grid points to [-1, 1]
        tmp = clenshawLegendre(z, hLegL);          % Evaluate via recurrence
        y(idx) = y(idx) + tmp;                     % Append
    end
    
    % The right triangle for the kth patch:
    if ( k < numPatches ) % Not needed for the final patch!
        % Locate the grid values in [dkl, dkr]:
        idx = dk_mid <= x & x < dk_right;
        z = map(x(idx), dk_mid, dk_right);
        tmp = clenshawLegendre(z, hLegR);
        y(idx) = y(idx) + tmp;
    end
end

if ( numPatches == 1 )
    % If there's only one patch, then we already have all the information reqd.
    hLegR = chebtech.leg2cheb(flipud(hLegR));       % Cheb coeffs of right tri.
    h_right = chebfun(hLegR, d + [a, b], 'coeffs'); % Make CHEBFUN
    h_mid = chebfun();
    
else  
    % Final right right triangle:
    %  ________________
    %  : /E|C/:      / 
    %  : ??/  : A  /
    %  : /D|B :  /
    %  /___|__:/
    %     fl a+d     b+d

    finishLocation = a + c + numPatches*(b - a);   % Where patches got to. (fl) 
    gk = restrict(g, d-[(b-a) 0]);                 % g on appropriate domain        
    gk_leg = chebtech.cheb2leg(get(gk, 'coeffs')); % Legendre coeffs
    [hLegL, hLegR] = easyConv(f_leg, gk_leg);      % Conv on A and B
    hLegR = chebtech.leg2cheb(flipud(hLegR));      % Cheb coeffs on A
    h_right = chebfun(hLegR, [d + a, d + b], 'coeffs'); % Make CHEBFUN
    
    % Remainder piece: (between fl and a+d)
    remainderWidth = d + a - finishLocation; % b+d-fl-(b-a)
    if ( remainderWidth > 0 )
        idx = finishLocation <= x;           % Discard D and E
        
        % B: (Coeffs were computed above)
        z = map(x(idx), d - b + 2*a, d + a); % Map grid to [-1, 1]
        tmp = clenshawLegendre(z, hLegL);    % Evaluate via recurrence
        y(idx) = tmp;                        % Store
        
        % C: 
        fk = restrict(f, b + [-remainderWidth, 0]);     % Restrict f
        fk_leg = chebtech.cheb2leg(get(fk, 'coeffs'));  % Legendre coeffs
        gk = restrict(g, [finishLocation, d + a] - b);  % Restrict g
        gk_leg = cheb2leg(get(gk, 'coeffs'));           % Legendre coeffs
        [ignored, hLegR] = easyConv(fk_leg, gk_leg);    % Conv 
        z = map(x(idx), finishLocation, d + a);         % Map to [-1, 1]
        tmp = clenshawLegendre(z, hLegR);               % Eval via recurrence
        y(idx) = y(idx) + tmp*remainderWidth/(b - a);   % Scale and append
    end
    
    % Convert values to coeffs (we don't want to construct a chebteh1)
    y = chebtech1.vals2coeffs(y);
    % Construct CHEBFUN of the interior (rectangle).
    h_mid = chebfun(y, [b+c, a+d], 'coeffs');
    
end

% Join the three pieces:
h = join(h_left, h_mid, h_right)*(b-a)/2;

% % TODO: Remove this when the time comes:
% f = chebfun({f});
% g = chebfun({g});
% tic
% exact = conv(f, g);
% toc
% err = norm(exact - h)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
N = length(alpha) + length(beta);

% Pad to make length n + 1.
alpha = [ alpha ; zeros(N - length(alpha), 1) ];

% M represents multiplication by 1/z in spherical Bessel space:
e = [[1 ; 1./(2*(1:(N-1)).'+1)], [1 ; zeros(N-1, 1)], -1./(2*(0:N-1).'+1)];
M = spdiags(e, -1:1, N, N);

gammaL = rec(M, alpha, beta, -1); % Chebyshev coeffs for the left piece
M(1,1) = -1;                      % Update M
gammaR = rec(M, -alpha, beta, 1); % Chebyshev coeffs for the right piece

%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX FREE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function gamma = rec(M, alpha, beta, sgn)
        
        % Initialise scl:
        nb = length(beta);
        scl = 1./(2*(1:nb).'-1);
        scl(2:2:end) = -scl(2:2:end);
        
        % First column of B:
        vNew = M*alpha; v = vNew;
        gamma = beta(1)*vNew;
        beta_scl = scl.*beta; beta_scl(1) = 0;
        gamma(1) = gamma(1) + vNew(1:nb).'*beta_scl;
        
        % The scalar case is trivial:
        if ( length(beta) == 1 ), return, end
        
        % Second column of B:
        vNew = M*v + sgn*v; vOld = v; v = vNew; vNew(1) = 0;
        gamma = gamma + beta(2)*vNew;
        beta_scl = -beta_scl*((2-.5)/(2-1.5)); beta_scl(2) = 0;
        gamma(2) = gamma(2) + vNew(1:nb).'*beta_scl;
        
        % Loop over remaining columns using recurrence:
        for k = 3:nb
            vNew = (2*k-3) * (M * v) + vOld; % Recurrence
            vNew(1:k-1) = 0;                 % Zero terms 
            gamma = gamma + vNew*beta(k);    % Append to g
            
            % Recurrence is unstable for j < k. Correct for upper-tri part:
            beta_scl = -beta_scl*((k-.5)/(k-1.5)); beta_scl(k) = 0;
            gamma(k) = gamma(k) + vNew(1:nb).'*beta_scl;
            
            vOld = v;
            v = vNew;
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%% FORMING B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO: Delete me when we're sure I'm not needed anymore.
    function gamma = rec2(M, alpha, beta, sgn)
        % Initialise B and scl:
        nb = length(beta);
        B = zeros(numel(alpha), nb);
        scl = 1./(2*(1:nb).'-1); scl(2:2:end) = -scl(2:2:end);
        
        % First 2 columns of B:
        vNew = M*alpha; v = vNew;
        B(:,1) = vNew;
        idx2 = 1:nb;
        B(1,:) = vNew(idx2).*scl(idx2);
        
        % The scalar case is trivial:
        if ( length(beta) == 1 ), return, end
        
        vNew = M*v + sgn*v; vOld = v; v = vNew;  B(:,2) = vNew;
        scl = -scl*(2-.5)/(2-1.5);
        idx2 = 3:nb;
        B(2,idx2) = vNew(idx2).*scl(idx2);
        
        % Loop over remaining columns using recurrence:
        for k = 3:nb
            vNew = (2*k-3) * ( M * v ) + vOld;
            idx1 = (k-1):numel(alpha);
            B(idx1, k) = vNew(idx1);
            vNew(1:k-2) = 0;
            % Recurrence is unstable for j < k. Correct for upper-tri part:
            scl = -scl*(k-.5)/(k-1.5);
            idx2 = (k+1):nb;
            B(k,idx2) = vNew(idx2).*scl(idx2);
            
            vOld = v;
            v = vNew;
        end
        
        gamma = B * beta;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = clenshawLegendre(x, alpha) 
% Evaluate a Legendre expansion with coefficient alpha at x. 
n = length(alpha); 
b_old = 0; 
b_cur = 0; 
for k = (n-1):-1:1
  b_new = alpha(k+1) + (2*k+1)/(k+1)*x.*b_cur - (k+1)/(k+2)*b_old;
  b_old = b_cur; 
  b_cur = b_new; 
end
val = alpha(1) + x.*b_cur - .5*b_old; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = myplus(f, g)
% Modified PLUS() which pads with zeros to fulfil domain requirements.
if ( isnumeric(f) )
    h = f + g;
else
    [a, b] = domain(f);
    [c, d] = domain(g);
    dom = union([a, b], [c, d]);
    h = chebfun(0, dom);
    h = defineInterval(h, [a, b], f);        % h{a, b} = f;
    hTmp = restrict(h, [c, d]);
    h = defineInterval(h, [c, d], hTmp + g); % h{c, d} = h{c, d} + g;
end
end
