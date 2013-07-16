function [g gsing Cm] = cumsum(g)
% CUMSUM	Indefinite integral
%
% For functions with exponents, things are more complicated. We switch to
% a Jacobi polynomial representation with the correct weights. We can then
% integrate all the terms for r > 0 exactly.
%
% In these cases [F1 F2] = cumsum(G) will return two funs, the first F1 will
% be the a smooth part (or a smooth part with exponents), whilst the 2nd F2
% will contain the terms which are harder to represent (using a sing-type
% map). (This is very experimental!)
%
% There is limited support for functions whose indefinite integral is also
% unbounded (i.e. G.exps <=1). In particular, the form of the blow up must
% be an integer, and G may only have exponents at one end of it's interval.
%
% Functions with both exponents and nonlinear maps can only be dealt with
% by switching to and from a linear map, and are therefore often very slow.

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

ends = g.map.par(1:2);
Cm = 0;

% linear map (simplest case)
if strcmp(g.map.name,'linear')
    
    if isempty(g), return, end
    
    if ~any(g.exps)
        g.vals = g.vals*g.map.der(0); % From change of variables to [-1,1]
        g.coeffs = g.coeffs*g.map.der(0);
        g = cumsum_unit_interval(g);
        g = g - g.vals(1);
        gsing = fun(0,g.map.par(1:2));
    elseif any(g.exps<=-1)
        if nargout > 1
            [g gsing Cm] = unbdnd(g);
        else
            g = unbdnd(g);
        end
    else
        if nargout > 1
            [g gsing Cm] = jacsum(g);
        else
            g = jacsum(g);
        end
    end
    
else
    
    map = g.map;
    if any(g.exps)
        warning('FUN:cumsum:dblexp',['Cumsum does not fully support functions ', ...
            'with both maps and exponents. Switching to a linear map (which may be VERY slow!)']);
        pref = chebfunpref;
        pref.splitting = false;
        pref.resampling = false;
        pref.blowup = false;
        % make the map linear
        exps = g.exps; g.exps = [0 0];
        g = fun(@(x) feval(g,x),linear(g.map.par(1:2)),pref);
        g.exps = exps;
        % do cumsum in linear case
        g = cumsum(g);
        % change the map back
        exps = g.exps; g.exps = [0 0];
        g = fun(@(x) feval(g,x),map,pref);
        g.exps = exps;
    else
        g.map = linear([-1 1]);
        g = cumsum_unit_interval(g.*fun(map.der,g.map));
        g.map = map;
    end
    
end

end

function g = cumsum_unit_interval(g)

n = g.n;
c = [0;0;chebpoly(g)];                        % obtain Cheb coeffs {c_r}
cout = zeros(n-1,1);                          % initialize vector {C_r}
cout(1:n-1) = (c(3:end-1)-c(1:end-3))./...    % compute C_(n+1) ... C_2
    (2*(n:-1:2)');
cout(n,1) = c(end) - c(end-2)/2;              % compute C_1
v = ones(1,n); v(end-1:-2:1) = -1;
cout(n+1,1) = v*cout;                         % compute C_0

% Trim small coeffs, as suggested in #128
tol = chebfunpref('eps')/norm(cout,inf);
idx = find(abs(cout)>tol,1);
cout(1:idx-1) = [];
if ~isempty(idx)
    n = n-idx+1;
end

% Recover vals and return the new fun
g.vals = chebpolyval(cout);
g.coeffs = cout;
g.scl.v = max(g.scl.v, norm(g.vals,inf));
g.n = n+1;

end

function [f G const] = jacsum(f)
% for testing - delete this eventually
% h = f; h.exps = [0 0];

% Get the exponents
ends = f.map.par(1:2);
exps = f.exps;
a = exps(2); b = exps(1);

% Compute Jacobi coefficients of F
j = jacpoly(f,a,b).';

if abs(j(end)) < chebfunpref('eps'), j(end) = 0; end

% Integrate the nonconstant terms exactly to get new coefficients
k = (length(j)-1:-1:1).';
jhat = -.5*j(1:end-1)./k;

% Convert back to Chebyshev series
c = jac2cheb2(a+1,b+1,jhat);

% Construct fun
f.vals = chebpolyval(c);
f.coeffs = c;
f.n = length(f.vals);
f.exps = f.exps + 1;
f = f*diff(ends)/2;
f.scl.v = max(f.scl.v, norm(f.vals,inf));

% Deal with the constant part
if j(end) == 0
    G = 0;
    const = 0;
elseif exps(2)
    const = j(end)*2^(a+b+1)*beta(b+1,a+1)*(diff(ends)/2);
    
    % Choose the right sing map
    mappar = [b a];
    mappar(mappar<=0) = mappar(mappar<=0)+1;
    mappar(mappar>1) = mappar(mappar>1)-floor(mappar(mappar>1)) ;
    map = maps(fun,{'sing',mappar},ends);
    
    pref = chebfunpref;
    if all(mappar), pref.exps = [mappar(1) 0]; end
    G = fun(@(x) const*betainc((x-ends(1))/diff(ends),b+1,a+1),map,pref,f.scl);
else
    G = fun(j(end)/(1+exps(1)),f.map.par(1:2));
    const = (2/diff(ends)).^exps(1);
    G = const*setexps(G,[exps(1)+1 0]);
    
end

% Add together smooth and singular terms
if nargout == 1 || ~exps(2)
    f = f + G;
end

f = replace_roots(f);

end

function cheb = jac2cheb2(a,b,jac)
N = length(jac)-1;

if ~N, cheb = jac; return, end

% Chebyshev-Gauss-Lobatto nodes
x = chebpts(N+1);

apb = a + b;

% Jacobi Vandermonde Matrix
P = zeros(N+1,N+1);
P(:,1) = 1;
P(:,2) = 0.5*(2*(a+1)+(apb+2)*(x-1));
for k = 2:N
    k2 = 2*k;
    k2apb = k2+apb;
    q1 =  k2*(k + apb)*(k2apb - 2);
    q2 = (k2apb - 1)*(a*a - b*b);
    q3 = (k2apb - 2)*(k2apb - 1)*k2apb;
    q4 =  2*(k + a - 1)*(k + b - 1)*k2apb;
    P(:,k+1) = ( (q2+q3*x).*P(:,k) - q4*P(:,k-1) ) / q1;
end

f = fun;
f.vals = P*flipud(jac(:)); f.n = length(f.vals);
cheb = chebpoly(f);

end



function [u M Cm] = unbdnd(f)
% If only one output is asked for, u is the whole function.
% For two outputs, u is the smooth part, and M contains the log singularity

% Get the exponents. At least one of these is less than -1.
exps = f.exps;
% Get the ends
ends = f.map.par(1:2);

% If the numerator is constant, things are simple.
if f.n == 1
    M = fun(0, f.map);
    if exps(1) && ~exps(2)
        if exps(1) == -1
            % Construct a representation of log
            f = setexps(f,[0 0]);
            u = makelog(f.vals(1),2,ends,f.scl.v);
            if nargout == 2, M = u; u = fun(0,ends); end
        else
            u = setexps(f,[f.exps(1)+1 0])/(f.exps(1)+1); 
            if round(exps(1))==exps(1), u = u - feval(u,ends(2)); end
        end
    elseif exps(2) && ~exps(1)
        if exps(2) == -1
            % Construct a representation of log
            f = setexps(f,[0 0]);
            u = makelog(f.vals(1),2,ends,f.scl.v);
            % Flip it
            u = flipud(u);
            if strcmp(u.map.name,'sing')
                u.map = maps(fun,{'sing',u.map.par([4 3])},u.map.par(1:2));
            end
            if nargout == 2, M = u; u = fun(0,ends); end
        else
            u = -setexps(f,[0 f.exps(2)+1])/(f.exps(2)+1);  
            if round(exps(2))==exps(2), u = u - feval(u,ends(1)); end
        end
    else
        error('FUN:cumsum:both1',['cumsum does not yet support functions whose ', ...
            'definite integral diverges and has exponents at both boundaries.'])
    end
    return
end

% We can only do linear maps at the moment
if ~strcmpi(f.map.name,'linear')
    error('FUN:cumsum:exps','cumsum does not yet support exponents <= 1 with arbitrary maps.');
end

% We may need to flip so that singularity is at the left
flip = false;
if exps(2)~=0
    if ~exps(1)  % Flip so singularity is on the left
        f = flipud(f); exps = f.exps;  flip = true;
    else
        error('FUN:cumsum:both2',['cumsum does not yet support functions whose ', ...
            'definite integral diverges and has exponents at both boundaries.']);
    end
end

% Remove exps from f
f = setexps(f,[0 0]);
% The order of the pole
a = -exps(1);
ra = max(round(a),1);               % ra = [a]

x = fun('x',ends);
xf = (x-ends(1)).*f;
N = length(xf)-1;
oldN = N;
if N < ra+2
    % Put some padding in
    N = ra+2;
    xf = prolong(xf,N+1);
end
aa = flipud(chebpoly(xf));


% The is the recurrence to get the coefficients for u'
dd = zeros(N,1);
dd(N) = 2*aa(N+1)./(1-a./N);
dd(N-1) = 2*(aa(N)-dd(N))./(1-a./(N-1));
for k = N-2:-1:ra+1
    dd(k) = 2*(aa(k+1)-dd(k+1)-dd(k+2)*.5*(1+a./k))./(1-a./k);
end

%%%%%
% % Subtract out a suitable term
% Cm = 2^(ra-1)*(aa(ra+1)-dd(ra+1)-dd(ra+2)*.5*(1+a./ra));
% % Adjustment for arbitrary intervals
de = diff(ends);
K = (4/de)^(ra-1)/de;
Cm = K*(2*aa(ra+1)-2*dd(ra+1)-dd(ra+2)*(1+a./ra));

% Find the new coefficients
M1 = fun(@(x)(x-ends(1)).^ra,ends);
aa(1:ra+1) = aa(1:ra+1) - Cm*flipud(chebpoly(M1));

% % Some testing
% df0 = feval(diff(f,ra-1),ends(1));% f^([a]-1)(-1). This is just for testing.
% % f = f - Cm*fun(@(x) (x-ends(1)).^(ra-1),ends);
% % xf = (x-ends(1)).*f;
% test1 = aa(ra+1)-dd(ra+1)-dd(ra+2)*.5*(1+a./ra);
% test2 = abs(df0-Cm);
% fprintf(' This should be zero:              \t %12.12g \n',test1)
% fprintf(' This should go to zero as a-->[a]:\t %12.12g \n',test2)

for k = ra-1:-1:1
    dd(k) = 2*(aa(k+1)-dd(k+1)-dd(k+2)*.5*(1+a./k))./(1-a./k);
end

% Convert coefficients for u' to those of u
kk = (1:N)';
dd = .5*dd; dd1 = dd./kk; dd2 = -dd(3:end)./kk(1:end-2);
cc = [0 ; dd1 + [dd2 ; 0 ; 0]];

% Choose first coefficient so that u(ends(1)) = (x-ends(1))*f(ends(1)) = 0;
cc(1) = sum(cc(2:2:end))-sum(cc(3:2:end));

% Reove some padding we put in
if N > oldN+2, cc = cc(1:oldN+2); end
% Construct the chebfun of the solution
u = fun(chebpolyval(flipud(cc)),ends);

% Plot for testing
% plot(xf - Cm*M1,'-b'); hold on
% plot((x-ends(1)).*diff(u)-a*u,'--r'); hold off

% Adding in the log term
tol = chebfunpref('eps');
if abs(ra-a) > tol*u.scl.v && 0     % No log term
    if nargout == 1
        u = u+(Cm/(ra-a))*M1;
        u = setexps(u,exps);
    else
        u = setexps(u,exps);
        M = Cm*(setexps(M1,exps)-1)/(ra-a);
    end
else                                % Log term
    u = setexps(u,exps);
    if abs(Cm) < 1e3*tol*u.scl.v
        % Contribution is small, so ignore it.
        M = fun(0,ends);
    else
        % Construct a representation of log
        M = makelog(Cm,ra,ends,u.scl.v);
        if nargout == 1, u = u+M;  end
    end
end

if strcmpi(u.map.name,'linear')
%     u = extract_roots(u);
    val = get(u,'lval');
    if isinf(val) || isnan(val), val = get(u,'rval'); end
    u = u - val; 
end

% Flip back so singularity is on the right
if flip
    u = -flipud(u);
    if strcmp(u.map.name,'sing')
        u.map = maps(fun,{'sing',u.map.par([4 3])},u.map.par(1:2));
    end
    if nargout == 2 && ~isempty(M)
        M = -flipud(M);
        if strcmp(M.map.name,'sing')
            M.map = maps(fun,{'sing',M.map.par([4 3])},M.map.par(1:2));
        end
    end
end

end

function M = makelog(Cm,ra,ends,scl)
% Constuct a representation of log(x-ends(1)) on interval ends
if ra == 1, ra = 2; end
if ra == 2, map = maps(fun,{'sing',[.125 1]},ends);
else        map = maps(fun,{'sing',[.25 1]},ends); end
pref = chebfunpref; pref.extrapolate = 1; pref.scl = scl;
M = fun(@(x) Cm*(x-ends(1)).^(ra-1).*log(x-ends(1)),map,pref);
M = setexps(M,[1-ra,0]);


% % plots for testing
% MM = chebfun(M,ends);
% xx = linspace(ends(1),ends(2),1e5);
% close all,
% plot(MM); hold on
% plot(xx,Cm*log(xx-ends(1)),'--r'); hold off
% figure
end