function h = conv(f, g)
%CONV   Convolution of CHEBFUN objects.
%   H = CONV(F, G) produces the convolution of CHEBFUN objects F and G:
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
%   Example:
%     f = chebfun(1/2); g = f;
%     subplot(2, 2, 1), plot(f)
%     for j = 2:4, g = conv(f, g); subplot(2, 2, j), plot(g), end
%     figure, for j = 1:4, subplot(2,2,j), plot(g), g = diff(g); end

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
% Rather than make a CHEBFUN corresponding to the each of the patches, we
% instead evaluate directly on a corresonding Chebyshev grid of appropriate
% size, which turns out to be far more efficient. 
%
% Total complexity: O( R*(m*n + n*n) + m*m )
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Force size(B) = [n, min(m, n)] for the interior convolutions?
% TODO: Refactor so that most of this lives at the FUN and/or CHEBTECH level?
% TODO: Support for delta functions?

% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    h = chebfun();
    return
end

% No support for quasimatrices:
if ( numColumns(f) > 1 || numColumns(g) > 1 )
    error('CHEBFUN:conv:quasi', 'No support for array-valued CHEBFUN objects.');
end

% Check transpose state:
if ( xor(f(1).isTransposed, g(1).isTransposed) )
    error('CHEBFUN:conv:transposed', 'CHEBFUN dimensions do not agree.');
end
transState = f(1).isTransposed;

% Extract the domain:
[a, b] = domain(f);
[c, d] = domain(g);

if ( any(isinf([a b c d])) )
    error('CHEBFUN:conv:bounded', ...
        'CONV only supports CHEBFUN objects on bounded domains.');
end

if ( issing(f) || issing(g) )
    % Call the old (and slow) version of CONV if we are not based on CHETECHS.
    h = oldConv(f, g);
    return
end

% Ensure that g is the signal (i.e., on the larger domain) and f is the filter:
if ( (b - a) > (d - c) )
    h = conv(g, f);
    return
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

end

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

