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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
%
% Nick Hale and Alex Townsend, 2014

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Force size(B) = [n, min(m, n)] for the interior convolutions?

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

% Initialize output:
h = 0;
% Deal with piecewise CHEBFUN objects.
% Loop over each of the interactions:
for j = 1:numel(f.funs)
    for k = 1:numel(g.funs)
        % Compute the contribution of jth fun of with kth fun of g:
        hjk = chebfun(conv(f.funs{j}, g.funs{k}));        
        % Add this contribution:
        h = myplus(h, hjk);
    end
end

if ( transState )
    h = h.';
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
