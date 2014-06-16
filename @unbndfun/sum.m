function out = sum(g)
%SUM   Definite integral of a UNBNDFUN on singly- or doubly-infinite domain.
%    SUM(G) is the definite integral of G, whose domain can be either singly-
%    infinite or doubly-infinite.
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the domain.
dom = g.domain;

% Cancel vanishing boundary values with negative exponents:
g.onefun = cancelExponents(g.onefun);

% Get the function values at the end of the domain. Note that the end point of
% the domain can be infinite.
vends = [get(g, 'lval'), get(g, 'rval')];

% Get the epslevel and vscale of the function g.
tol = get(g, 'epslevel')*get(g, 'vscale');

% A dirty checklist:

% Check 1: Check if the function values are vanishing at infinity/ies.
unbounded = [];
if ( isinf(dom(1)) && abs(vends(1)) > 1e4*tol )
    unbounded(1) = sign(vends(1))*Inf;
end
if ( isinf(dom(2)) && abs(vends(2)) > 1e4*tol )
    unbounded(2) = sign(vends(2))*Inf;
end
if ( ~isempty(unbounded) )
    out = sum(unbounded);
    return
end

% Check 2: Check the speed of decay at infinity/ties. The integrand is
% integrable only when it decays faster than 1/x towards infinity/ties.
if ( any(~isdecay(g.onefun) & isinf(dom)) )
    warning('CHEBFUN:UNBNDFUN:sum:slowDecay', ...
        ['Result may not be accurate ' ...
         'as the function decays slowly at infinity.'])
end

% If we reach here, the function decays sufficiently fast. Construct a ONEFUN
% for the integrand and integrate it.

tech = get(g.onefun, 'tech');
if ( isa(tech(), 'chebtech') )
    techPrefs.exactLength = length(g);
    % TODO: Using an exact-length construction here is a hack. It only works if
    % g.onefun is a CHEBTECH, and, even then, the idea that it "works" is merely
    % heuristic; there's no a priori reason that it should. We really would like
    % to do an adaptive construction, but the fact that the nonlinear map
    % compresses very wide intervals near +/-Inf into very small ones near +/-1,
    % means that the filtered function produced by unbndfunIntegrand() will
    % appear to exhibit sharp transitions to zero when sampled on a Chebyshev
    % grid of the usual sizes, causing the constructor to fail. Until we can
    % solve this problem there doesn't seem to be much else we can do here.
else
    techPrefs = [];
end
integrand = tech(@(x) unbndfunIntegrand(x, g), [], techPrefs);
out = sum(integrand);

end

function y = unbndfunIntegrand(x, g)
%UNBNDFUNINTEGRAND   Compute the integrand in the integral of an UNBNDFUN.
%   Y = UNBNDFUNINTEGRAND(X, G) evaluates at X the mapped integrand with the
%   chain rule applied when one goes to compute the integral of G, i.e., it
%   computes G(X).*G.MAPPING.DER(X).  Values of G less than a tolerance are
%   set to zero.  This filtering is applied to avoid problems with round-off
%   error in G near infinite endpoints being amplified by the derivative of the
%   map, which is large there.
%
%   This function is meant to be used to create a ONEFUN for the integrand,
%   which can then be integrated with the ONEFUN's implementation of SUM.

tol = 10*get(g, 'epslevel')*get(g, 'vscale');
y = feval(g, g.mapping.for(x));
y(abs(y) < tol) = 0;
y = y.*g.mapping.der(x);

end
