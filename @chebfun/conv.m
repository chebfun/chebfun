function h = conv(f, g, varargin)
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
%   H = CONV(F, G, 'same') will truncate the domain of H so that it is the same
%   as F. This is useful when F and G represent rapidly decaying functions on
%   large but finite intervals which are used to approximate infinity.
%
%   If F and G are simple, in the sense that their FUNS are CHEBTECH objects, a
%   fast algorithm due to Hale and Townsend is used [1]. Otherwise, the integral
%   is computed by brute force. CONV(F, G, 'old') forces the brute force
%   approach, even when the fast algorithm may be used.
%
%   Note that CONV only supports piecewise-smooth functions on bounded domains.
%
% Examples:
%     cheb.x; f = 0.8 - abs(x-0.2);
%     phi = @(t) chebfun(@(x) exp(-x^2/(4*t))/sqrt(4*pi*t));
%     fsmooth = conv(f,phi(1e-2),'same');
%     plot(f,'b',fsmooth','r')
%
%     f = chebfun(1/2); g = f;
%     subplot(2, 2, 1), plot(f)
%     for j = 2:4, g = conv(f, g); subplot(2, 2, j), plot(g), end
%     figure, for j = 1:4, subplot(2,2,j), plot(g), g = diff(g); end
%
% Reference:
%   [1] N. Hale and A. Townsend, "An algorithm for the convolution of Legendre
%   series", SIAM Journal on Scientific Computing, Vol. 36, No. 3,
%   pp. A1207-A1220, 2014.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check transpose state:
if ( xor(f(1).isTransposed, g(1).isTransposed) )
    error('CHEBFUN:CHEBFUN:conv:transposed', ...
        'CHEBFUN dimensions do not agree.');
end
transState = f(1).isTransposed;

% Support for quasimatrices:
nf = numColumns(f);
ng = numColumns(g);
if ( nf > 1 || ng > 1 )
    if ( nf == ng )
        f = mat2cell(f);
        g = mat2cell(g);
        h = f;
        for k = 1:nf
            h{k} = conv(f{k}, g{k}, varargin{:});
        end
    elseif ( nf == 1 )
        g = mat2cell(g);
        h = g;
        for k = 1:ng
            h{k} = conv(f, g{k}, varargin{:});
        end
    elseif ( ng == 1 )
        f = mat2cell(f);
        h = f;
        for k = 1:nf
            h{k} = conv(f{k}, g, varargin{:});
        end
    else
        error('CHEBFUN:CHEBFUN:conv:dimagree', 'CHEBFUN dimensions must agree.');
    end
    if ( ~transState )
        h = horzcat(h{:});
    else
        h = vertcat(h{:});
    end
    return
end

% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    h = chebfun();
    return
end

% Parse inputs:
oldMethod = false;
same = false;
for k = 1:numel(varargin)
    vk = varargin{k};
    if ( strcmpi(vk, 'old') )
        oldMethod = true;
    elseif ( strcmpi(vk, 'same') )
        same = true;
    elseif ( strcmpi(vk, 'full') )
        % Do nothing.
    elseif ( strcmpi(vk, 'valid') )
        % TODO: Supoprt 'valid'. Presumably where domains of f and g overlap?
        error('CHEBFUN:CHEBFUN:conv:validFlag', '''valid'' is not yet supprted.');
    else
        error('CHEBFUN:CHEBFUN:conv:badInput', 'Unknown input option %s.', vk);
    end
end

% Return a warning if F and G have too many pieces (the computation is probably
% going to be very slow):
if ( ( numel(f.funs) + numel(g.funs) ) > 50 ) 
    % Give a warning and proceed. 
   warning('CHEBFUN:CHEBFUN:conv:piecewise',...
       ['Convolving CHEBFUNs with many pieces can be very slow.\n', ...
        'Try calling MERGE() on the inputs before calling CONV().']);
end

% Extract the domain:
[a, b] = domain(f);
[c, d] = domain(g);

% No support for unbounded domains:
if ( any(isinf([a b c d])) )
    error('CHEBFUN:CHEBFUN:conv:bounded', ...
        'CONV only supports CHEBFUN objects on bounded domains.');
end

if ( oldMethod || issing(f) || issing(g) )
    % Call the old (and slow) version of CONV if we are not based on CHEBTECHS.
    h = oldConv(f, g);
    
else

    % Ensure g is the signal (i.e., on the larger domain) and f is the filter:
    if ( (b - a) > (d - c) )
        h = conv(g, f);
        return
    end

    % Initialize the output:
    h = chebfun(0, [a + c, b + d]);
    % Deal with piecewise CHEBFUN objects by looping over each interaction:
    for j = 1:numel(f.funs)
        for k = 1:numel(g.funs)
            % Compute the contribution of jth fun of f with kth fun of g:
            hjk = conv(f.funs{j}, g.funs{k});  
            % Add this contribution:            
            h = myplus(h, chebfun(hjk));
        end
    end
    
    % Make sure that point values are not added twice:
    dom = domain(h);
    intDom = dom(2:end-1);
    h.pointValues(2:end-1) = 1/2*(feval(h, intDom.', 'left') + ...
        feval(h, intDom.', 'right'));
    
end

% Truncate:
if ( same )
    h = restrict(h, [a, b]);
end

% Transpose:
if ( transState )
    h = h.';
end

% Simplify the result.
h = simplify(h);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = myplus(f, g)
% Modified PLUS() which pads with zeros to fulfil domain requirements.
%  Note f is always on the largest possible domain. g is on a subdomain of f

% Tidy the domains:
[f, g] = tweakDomain(f, g);
[c, d] = domain(g);    

fTmp = restrict(f, [c, d]); % f{c, d}
hTmp = fTmp + g;            % h{c, d} = f{c, d} + g;
h = defineInterval(f, [c, d], hTmp); 

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

% Set preferences:
%
% TODO:  CHEBFUN is not supposed to set the refinementFunction preference
% because it doesn't belong to the list of "abstract" preferences required of
% all techs.  Do we really need to alter it here?
p = chebfunpref();
p.splitting = false;
p.blowup = false;
p.techPrefs.extrapolate = true;
p.techPrefs.refinementFunction = 'nested';
p.techPrefs.sampleTest = false;

% Construct FUNS:
funs = cell(1, length(dom)-1);
for k = 1:length(dom)-1  
    data.domain = dom(k:k+1);
    data.vscale = vs;
    data.hscale = hs;
    newFun = bndfun(@(x) convIntegral(x, f, g), data, p);
    vs = max(get(newFun, 'vscale'), vs); 
    funs{k} = newFun;
end

% Construct CHEBFUN:
h = chebfun(funs);
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
            % INTEGRAL is not available in versions of MATLAB prior to R2012a,
            % so if we're running on an older version, fall back to QUADGK.
            integrand = @(t) feval(f, t).*feval(g, x(k) - t);
            if ( verLessThan('matlab', '7.14') )
                out(k) = out(k) + quadgk(integrand, dom(j), dom(j+1), ...
                    'AbsTol', 1e-15, 'RelTol', 100*eps);
            else
                out(k) = out(k) + integral(integrand, dom(j), dom(j+1), ...
                    'AbsTol', 1e-15, 'RelTol', 1e-15);
            end
        end
    end
end

end
