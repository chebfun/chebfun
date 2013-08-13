function f = cumsum(f)
%CUMSUM	Indefinite integral of a SINGFUN.
% CUMSUM(F) is the indefinite integral of the SINGFUN F.
%
% For functions with exponents, things are more complicated. We switch to
% a Jacobi polynomial representation with the correct weights. 
%
% WARNING: The current version of cumsum@singfun can cover limited cases. The
% case that the antiderivative has logarithm (including atanh) and the case that
% the integrand has singularities at both ends are not supported. 

% [TODO]: Improvement on the algorithm to handle the case with singularities at 
% both the end points. Improvement on the whole singfun class to handle 
% logarithmic singularity, otherwise the singfun class is not closed. Here, by 
% 'closed' we mean that the result of an operator on an arbitrary singfun (unary
% operator) or any two arbitrary singfuns (binary operator) is representable by 
% a singfun object. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% If f has no singularity at any endpoint, then just integrate its smooth part
% and return.

if all( ~f.exponents ) % no exponents, only integrate the smooth part.
    f = f.smoothPart.cumsum();
    return
elseif ( any( f.exponents < 0 ) && any( f.exponents == 0 ) ) % one singularity
    f = singintegral(f);
else % error message thrown for other cases
    error('SINGFUN:cumsum:nosupport',['cumsum does not support the case ', ...
        'with the current exponents.'])
end

end

%%
% singintegral was originally named unbdnd in Chebfun v4. It tryies to find a
% singfun representation for the indefinite integral of a singfun. It now can
% only handle integrand with singularity at one end of the domain, but not at 
% both. The algorithm was due to Nick Hale and Sheehan Olver. Some minor
% modification are made during Chebfun v5 refactorization. See the workig note
% by Nick Hale and Sheehan Olver for the detail about the algorithm.

function u = singintegral(f)

% When the singularity is at the right end, we flip f so that the singularity is
% at the left end of the domain.

flip = false;
if ( f.exponents(2) ~= 0 )
    f = flipud(f);
    flip = true;
end

% Get the exponents.
exps = f.exponents;

% Get the smooth part of f
s = f.smoothPart;

% Get the order of the singularity and negate it.
a = -exps(1);

% Compute the rounded integer of a, that is ra = [a].
ra = max(round(a),1);      

% Compute (x+1)*s. 
xs = smoothfun(@(x) x+1).*s;

% If the length of xs is less than ra+2, we pad the length of xs by prolonging
% it. This will save us from branch out for different cases when computing the 
% coefficients c_k of the smooth part of the u. The system for c_k is indicated
% by (*) below and c_k is solve by recursive substitution.

N = length(xs)-1;
oldN = N;
if N < ra+2
    N = ra+2;
    xs = prolong(xs, N+1);
end

% We flip up and down to have the coefficients of xs ordered with ascending 
% indice. 
aa = flipud(xs.coeffs);

% The recurrence to solve for the coefficients for u', i.e. c_k. (*)
c = zeros(N,1);
c(N) = 2*aa(N+1)./(1-a./N);
c(N-1) = 2*(aa(N)-c(N))./(1-a./(N-1));
for k = N-2:-1:ra+1
    c(k) = 2*(aa(k+1)-c(k+1)-c(k+2)*.5*(1+a./k))./(1-a./k);
end

% Compute Cm
Cm = (2^(ra-1))*(2*aa(ra+1)-2*c(ra+1)-c(ra+2)*(1+a./ra));

% Compute the smoothfun representation for (x+1)^[a]
xa = smoothfun(@(x)(x-ends(1)).^ra);

% Intermediate results for temporary use.
aa(1:ra+1) = aa(1:ra+1) - Cm*flipud(xa.coeffs);

% Some testing:
% df0 = feval(diff(f,ra-1),ends(1));% f^([a]-1)(-1). This is just for testing.
% % f = f - Cm*fun(@(x) (x-ends(1)).^(ra-1),ends);
% % xf = (x-ends(1)).*f;
% test1 = aa(ra+1)-dd(ra+1)-dd(ra+2)*.5*(1+a./ra);
% test2 = abs(df0-Cm);
% fprintf(' This should be zero:              \t %12.12g \n',test1)
% fprintf(' This should go to zero as a-->[a]:\t %12.12g \n',test2)

% Compute the rest of the coefficients c_k.
for k = ra-1:-1:1
    c(k) = 2*(aa(k+1)-c(k+1)-c(k+2)*.5*(1+a./k))./(1-a./k);
end

% Convert coefficients for u' to those of u
kk = (1:N)';
c = .5*c; 
dd1 = c./kk; 
dd2 = -c(3:end)./kk(1:end-2);
cc = [0 ; dd1 + [dd2 ; 0 ; 0]];

% Choose first coefficient so that u(-1) = (x+1)*f(-1) = 0.
cc(1) = sum(cc(2:2:end))-sum(cc(3:2:end));

% Remove the padding we put in.
if ( N > (oldN + 2) )
    cc = cc(1:oldN+2);
end

% Construct u as a singfun object.
u = singfun;


% Plot for testing
% plot(xf - Cm*M1,'-b'); hold on
% plot((x-ends(1)).*diff(u)-a*u,'--r'); hold off

% Construct the singfun object of the solution.
tol = singfun.pref.singfun.eps;
if abs(ra-a) > tol*u.smoothPart.vscale     % No log term
    u.smoothPart = f.smoothPart.make({[],cc}) + (Cm/(ra-a))*xa;
    u.exponents = exps;
else                                % Log term
    if abs(Cm) < tol*u.smoothPart.vscale
        % Contribution is small, so ignore it.
        u.smoothPart = f.smoothPart.make({[],cc});
    else
        % Construct a representation of log
        error('SINGFUN:cumsum:nolog',['cumsum does not support the case ', ...
            'in which the indefinite integral has a logarithm term.'])
    end
end

% Flip back so singularity is on the right for the case with singularity at the
% right end of the domain.

if flip
    u = flipud(u);
end

end

%% 
% jacsum is inherited from Chebfun v4. It computes the indefinite integral of a 
% singular function f(x) = s(x)*(1+x)^a*(1-x)^b, where -1 < a, b < 0. Since its
% main functionalities can be covered by singintegral, it serves now only as a
% historical legacy for reference. 

% function [f G const] = jacsum(f)
% % for testing - delete this eventually
% % h = f; h.exps = [0 0];
% 
% % Get the exponents
% ends = f.map.par(1:2);
% exps = f.exps;
% a = exps(2); b = exps(1);
% 
% % Compute Jacobi coefficients of F
% j = jacpoly(f,a,b).';
% 
% if abs(j(end)) < chebfunpref('eps'), j(end) = 0; end
% 
% % Integrate the nonconstant terms exactly to get new coefficients
% k = (length(j)-1:-1:1).';
% jhat = -.5*j(1:end-1)./k;
% 
% % Convert back to Chebyshev series
% c = jac2cheb2(a+1,b+1,jhat);
% 
% % Construct fun
% f.vals = chebpolyval(c);
% f.coeffs = c;
% f.n = length(f.vals);
% f.exps = f.exps + 1;
% f = f*diff(ends)/2;
% f.scl.v = max(f.scl.v, norm(f.vals,inf));
% 
% % Deal with the constant part
% if j(end) == 0
%     G = 0;
%     const = 0;
% elseif exps(2)
%     const = j(end)*2^(a+b+1)*beta(b+1,a+1)*(diff(ends)/2);
%     
%     % Choose the right sing map
%     mappar = [b a];
%     mappar(mappar<=0) = mappar(mappar<=0)+1;
%     mappar(mappar>1) = mappar(mappar>1)-floor(mappar(mappar>1)) ;
%     map = maps(fun,{'sing',mappar},ends);
%     
%     pref = chebfunpref;
%     if all(mappar), pref.exps = [mappar(1) 0]; end
%     G = fun(@(x) const*betainc((x-ends(1))/diff(ends),b+1,a+1),map,pref,f.scl);
% else
%     G = fun(j(end)/(1+exps(1)),f.map.par(1:2));
%     const = (2/diff(ends)).^exps(1);
%     G = const*setexps(G,[exps(1)+1 0]);
%     
% end
% 
% % Add together smooth and singular terms
% if nargout == 1 || ~exps(2)
%     f = f + G;
% end
% 
% f = replace_roots(f);
% 
% end

%% 
% jac2cheb2 is inherited from Chebfun v4. It is called by jacsum and converts 
% Jacobi coefficients to Chebyshev coefficients. Since jacsum is fully covered
% by singintegral, jac2cheb2 is to be removed eventually.

% function cheb = jac2cheb2(a,b,jac)
% N = length(jac)-1;
% 
% if ~N, cheb = jac; return, end
% 
% % Chebyshev-Gauss-Lobatto nodes
% x = chebpts(N+1);
% 
% apb = a + b;
% 
% % Jacobi Vandermonde Matrix
% P = zeros(N+1,N+1);
% P(:,1) = 1;
% P(:,2) = 0.5*(2*(a+1)+(apb+2)*(x-1));
% for k = 2:N
%     k2 = 2*k;
%     k2apb = k2+apb;
%     q1 =  k2*(k + apb)*(k2apb - 2);
%     q2 = (k2apb - 1)*(a*a - b*b);
%     q3 = (k2apb - 2)*(k2apb - 1)*k2apb;
%     q4 =  2*(k + a - 1)*(k + b - 1)*k2apb;
%     P(:,k+1) = ( (q2+q3*x).*P(:,k) - q4*P(:,k-1) ) / q1;
% end
% 
% f = fun;
% f.vals = P*flipud(jac(:)); f.n = length(f.vals);
% cheb = chebpoly(f);
% 
% end

%% 
% makelog is inherited from Chebfun v4. It builds a Chebfun v4 fun 
% representation for logarithm functions via a sing map. Since sing map is no
% longer used in Chebfun v5, this function only serves for references.

% function M = makelog(Cm,ra,ends,scl)
% % Constuct a representation of log(x-ends(1)) on interval ends
% if ra == 1, ra = 2; end
% if ra == 2, map = maps(fun,{'sing',[.125 1]},ends);
% else        map = maps(fun,{'sing',[.25 1]},ends); end
% pref = chebfunpref; pref.extrapolate = 1; pref.scl = scl;
% M = fun(@(x) Cm*(x-ends(1)).^(ra-1).*log(x-ends(1)),map,pref);
% M = setexps(M,[1-ra,0]);
% 
% 
% % % plots for testing
% % MM = chebfun(M,ends);
% % xx = linspace(ends(1),ends(2),1e5);
% % close all,
% % plot(MM); hold on
% % plot(xx,Cm*log(xx-ends(1)),'--r'); hold off
% % figure
% end
