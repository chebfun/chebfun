function h = conv(f, g, flag)
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
%   If F and G are simple, in the sense that their FUNS are CHEBTECH objects, a
%   fast algorithm due to Hale and Townsend is used [1]. Otherwise, the integral
%   is computed by brute force. CONV(F, G, 'old') forces the brute force
%   approach, even when the fast algorithm may be used.
%
%   Note that CONV only supports piecewise-smooth functions on bounded domains.
%
%   Example:
%     f = chebfun(1/2); g = f;
%     subplot(2, 2, 1), plot(f)
%     for j = 2:4, g = conv(f, g); subplot(2, 2, j), plot(g), end
%     figure, for j = 1:4, subplot(2,2,j), plot(g), g = diff(g); end

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
%
% [1] N. Hale and A. Townsend, "An algorithm for the convolution of Legendre
% series", (To appear in SISC)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% No support for unbounded domains:s
if ( any(isinf([a b c d])) )
    error('CHEBFUN:conv:bounded', ...
        'CONV only supports CHEBFUN objects on bounded domains.');
end

if ( issing(f) || issing(g) || nargin == 3 )
    % Call the old (and slow) version of CONV if we are not based on CHEBTECHS.
    h = oldConv(f, g);
    return
end

% Ensure that g is the signal (i.e., on the larger domain) and f is the filter:
if ( (b - a) > (d - c) )
    h = conv(g, f);
    return
end

% Initialize the output:
h = chebfun(0, [a + c, b + d]);
% Deal with piecewise CHEBFUN objects by looping over each of the interactions:
for j = 1:numel(f.funs)
    for k = 1:numel(g.funs)
        % Compute the contribution of jth fun of f with kth fun of g:
        hjk = conv(f.funs{j}, g.funs{k});  
        % Add this contribution:
        for i = 1:numel(hjk)
            h = myplus(h, chebfun(hjk(i)));
        end
    end
end

if ( transState )
    h = h.';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = myplus(f, g)
% Modified PLUS() which pads with zeros to fulfil domain requirements.

% Tidy the domains:
[f, g] = tweakDomain(f, g);

% f is always on the largest possible domain. g is on a subdomain of f:
[c, d] = domain(g);    
fTmp = restrict(f, [c, d]); % f{c, d}
f = defineInterval(f, [c, d], fTmp + g); % f{c, d} = f{c, d} + g;

% Make sure that point values are not added twice:
dom = domain(f);
intDom = dom(2:end-1);
f.pointValues(2:end-1) = 1/2*(feval(f, intDom.', 'left') + feval(f, intDom.', 'right'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = oldConv(f, g)
% The old convolution algorithm based on quadrature. See section 2 of [1] for a
% description of the algorithm. 

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
p = chebfunpref();
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
