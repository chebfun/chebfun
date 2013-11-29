function h = conv(f, g)
%CONV   Convolution of CHEBFUN objects.
% H = CONV(F, G) produces the convolution of CHEBFUN objects F and G:
%                   -
%                  /
%         H(x) =   |    F(t) G(x-t) dt,
%                  /
%                 -
% defined for x in [a+c, b+d], where F.domain is [a,b] and G.domain is [c, d].
% The integral is taken over all t for which the integrand is defined: max(a,
% x-d) <= t <= min(b, x-c), and the breakpoints of H are all pairwise sums of
% the breakpoints of F and G.
%
% Example:
%   f = chebfun(1/2); g = f;
%   subplot(2, 2, 1), plot(f)
%   for j = 2:4, g = conv(f, g); subplot(2, 2, j), plot(g), end
%   figure, for j = 1:4, subplot(2, 2, j), plot(g), g = diff(g); end

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: f and g may be defined on different domains!
h = chebfun();
if ( isempty(f) || isempty(g) )
    return
end

if ( xor(f.isTransposed, g.isTransposed) )
    error('CHEBFUN:conv:transposed', 'CHEBFUN dimensions do not agree.');
end

if ( min(size(f)) > 1 || min(size(g)) > 1 ) % TODO
    error('CHEBFUN:conv:quasi', 'No support for array-valued CHEBFUN objects.');
end

% Remove deltas from f:
fImps = f.impulses(:,:,2:end);
if ( size(f.impulses, 3) > 1 )
    f.impulses = f.impulses(:,:,1); 
end
gImps = g.impulses(:,:,2:end);
% Remove deltas from g:
if ( size(g.impulses, 3) > 1 )
    g.impulses = g.impulses(:,:,1); 
end

% Find all breakpoints in the convolution:
[A, B] = meshgrid(f.domain, g.domain);
dom = unique( A(:) + B(:) ).';

% Coalesce breaks that are close due to roundoff:
dom( diff(dom) < 10*eps*max(abs(dom([1,end]))) ) = [];
dom(isnan(dom)) = [];

if ( any(isinf(dom)) )
    error('CHEBFUN:conv:bounded', ...
        'CONV only supports CHEBFUN objects on bounded intervals.');
end

% Combine vertical and horizontal scales:
hs = max(hscale(f), hscale(g));
vs = 2*max([vscale(f), vscale(g)]);

% Avoid resampling for speed up:
p = chebpref;
p.enableBreakpointDetection = false;
p.enableSingularityDetection = false;
p.techPrefs.extrapolate = true;
p.techPrefs.resampling = false;
p.techPrefs.sampletest = false;

% For convenience:
a = f.domain(1); 
b = f.domain(end); 
c = g.domain(1); 
d = g.domain(end);

% Construct FUNS:
funs = cell(1, length(dom)-1);
for k = 1:length(dom)-1  
    newFun = bndfun(@(x) convIntegral(x, a, b, c, d, f, g), dom(k:k+1), vs, hs, p);
    vs = max(get(newFun, 'vscale'), vs); 
    funs{k} = newFun;
end

% Construct CHEBFUN:
h = chebfun(funs);
h = simplify(h);
h.isTransposed = f.isTransposed;

% TODO: Delta functions.

end

function out = convIntegral(x, a, b, c, d, f, g)
out = 0*x;
for k = 1:length(x)
    A = max(a, x(k) - d); 
    B = min(b, x(k) - c);
    if ( A < B )
        ends = union(x(k) - g.domain, f.domain);
        dom = [A, ends(A < ends & ends < B),  B];
        for j = 1:length(dom)-1
            out(k) = out(k) + integral(@(t) feval(f, t).*feval(g, x(k) - t), ...
                dom(j), dom(j+1), 'AbsTol', 1e-15, 'RelTol', 1e-15);
        end
    end
end
end