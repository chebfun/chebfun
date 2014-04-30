function h = conv(f, g)
%CONV   Convolution of DELTAFUN objects.
%   H = CONV(F, G) produces the convolution of DELTAFUN objects F and G:
%                     - 
%                    /
%           H(x) =   |    F(t) G(x-t) dt,  x in [a + c, b + d]
%                    /
%                   -
%   Example:
%     f = chebfun(1/2); g = f;
%     subplot(2, 2, 1), plot(f)
%     for j = 2:4, g = conv(f, g); subplot(2, 2, j), plot(g), end
%     figure, for j = 1:4, subplot(2,2,j), plot(g), g = diff(g); end

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
%
% Nick Hale and Alex Townsend, 2014

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    h = deltafun();
    return
end

if ( ~isa(f, 'deltafun') )
    % Then g must be a deltafun
    h = conv(g, f);
    return;
end

hFun = conv(f.funPart, g.funPart);
            
% Extract the domain:
[a, b] = domain(f.funPart);
[c, d] = domain(g.funPart);

if ( any(isinf([a b c d])) )
    error('CHEBFUN:conv:bounded', ...
        'CONV only supports CHEBFUN objects on bounded domains.');
end
    
% Deal with piecewise CHEBFUN objects.
if ( (numel(f.funs) > 1) || (numel(g.funs) > 1) )
    h = 0;
    % Loop over each of the interactions:
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
        h_left = chebfun(hLegL, [dk_left, dk_mid], 'coeffs'); % Make CHEBFUN
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

    finishLocation = a + c + numPatches*(b - a);    % Where patches got to. (fl) 
    gk = restrict(g, d-[(b-a) 0]);                  % g on appropriate domain   
    gk = simplify(gk);                              % Simplify for efficiency
    gk_leg = chebtech.cheb2leg(get(gk, 'coeffs'));  % Legendre coeffs
    [hLegL, hLegR] = easyConv(f_leg, gk_leg);       % Conv on A and B
    hLegR = chebtech.leg2cheb(flipud(hLegR));       % Cheb coeffs on A
    h_right = chebfun(hLegR, [d+a, d+b], 'coeffs'); % Make CHEBFUN
    
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
    % Construct CHEBFUN of the interior (rectangle).
    h_mid = chebfun(y, [b+c, a+d], 'coeffs');
    
end

% Join the three pieces:
h = join(h_left, h_mid, h_right)*(b-a)/2;

if ( transState )
    h = h.';
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = oldConv(f, g)

% Find all breakpoints in the convolution:
[A, B] = meshgrid(f.domain, g.domain);
dom = unique(A(:) + B(:)).';

% Coalesce breaks that are close due to roundoff:
dom(diff(dom) < 10*eps*max(abs(dom([1,end])))) = [];
dom(isnan(dom)) = [];

% Combine vertical and horizontal scales:
hs = max(hscale(f), hscale(g));
vs = 2*max([vscale(f), vscale(g)]);

% Avoid resampling for speed up:
p = chebpref();
p.enableBreakpointDetection = false;
p.enableSingularityDetection = false;
p.techPrefs.extrapolate = true;
p.techPrefs.resampling = false;
p.techPrefs.sampletest = false;

% Construct FUNS:
funs = cell(1, length(dom)-1);
for k = 1:length(dom)-1  
    newFun = bndfun(@(x) convIntegral(x, f, g), dom(k:k+1), vs, hs, p);
    vs = max(get(newFun, 'vscale'), vs); 
    funs{k} = newFun;
end

% Construct CHEBFUN:
h = chebfun(funs);
h = simplify(h);
h.isTransposed = f.isTransposed;

end

function out = convIntegral(x, f, g)
%CONVINTEGRAL   Evaluate convolution integral.
%   Y = CONVINTEGRAL(X, F, G) evaluates the convolution of the CHEBFUNs F and G
%   at the points X.

a = f.domain(1);
b = f.domain(end);
c = g.domain(1);
d = g.domain(end);

out = 0*x;
for k = 1:length(x)
    A = max(a, x(k) - d); 
    B = min(b, x(k) - c);
    if ( A < B )
        ends = union(x(k) - g.domain, f.domain);
        dom = [A, ends((A < ends) & (ends < B)), B];
        for j = 1:length(dom)-1
            out(k) = out(k) + integral(@(t) feval(f, t).*feval(g, x(k) - t), ...
                dom(j), dom(j+1), 'AbsTol', 1e-15, 'RelTol', 1e-15);
        end
    end
end

end

