function out = sum(g, pref)
% SUM	Definite integral of a UNBNDFUN on singly- or doubly-infinite domain.
%    SUM(G) is the definite integral of G, whose domain can be either singly-
%    infinite or doubly-infinite.
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the domain.
dom = g.domain;

% Get the fun class preference if no preference is passed.
if ( nargin < 2 )
    pref = chebfunpref();
end

% Get the function values at the end of the domain. Note that the end point of
% the domain can be infinite.
vends = [get(g, 'lval'), get(g, 'rval')];

% Get the epslevel and vscale of the function g.
tol = 1e1*get(g, 'epslevel')*get(g, 'vscale');

if ( ~issing(g) )
    
    % Cases for ONEFUN has no singularities. The function may be decaying or
    % approach a constant value at infinities.
    
    % A dirty checklist:

    % Check 1: Check if the function values are vanishing at infinity/ies.
    unbounded = [];
    if ( isinf(dom(1)) && abs(vends(1)) > tol )
        unbounded(1) = sign(vends(1))*Inf;
    end
    if ( isinf(dom(2)) && abs(vends(2)) > tol )
        unbounded(2) = sign(vends(2))*Inf;
    end
    if ( ~isempty(unbounded) )
        out = sum(unbounded);
        return
    end
    
    % Check 2: Check the speed of decay at infinity/ties. The integrand is
    % integrable only when it decays faster than 1/x towards infinity/ties.
    
    % Call ISDECAY of ONEFUN to check if G decays fast enough to be integrable:
    decay = isdecay(g.onefun);
    
    maskInf = isinf(dom);
    if ( any(~decay & maskInf) )
        warning('CHEBFUN:UNBNDFUN:sum:slowdecay', ...
            ['Result may not be accurate ' ...
            'as the function decays slowly at infinity.'])
    end
    
    % If we reach here, the function decays sufficiently fast.
    
    % Construct the ONEFUN presentation of the derivative of the forward map.
    pref.singPrefs.exponents = g.mapping.forDerExps;
    forDer = onefun.constructor(@(x) g.mapping.forDer(x), [], [], pref);
    
    % Peel off the boundary roots to cancel the negative exponents of forDer:
    h = extractBoundaryRoots(g.onefun, -g.mapping.forDerExps.');
    
    % Form the new integrand:
    integrand = h.*forDer.smoothPart;
    
    % Call the sum at onefun level.
    out = sum(integrand);
    
elseif ( issing(g) ) % Cases for ONEFUN has singularities at the end points.
    
    % Construct the ONEFUN presentation of the derivative of the forward map.
    pref.singPrefs.exponents = g.mapping.forDerExps;
    forDer = onefun.constructor(@(x) g.mapping.forDer(x), [], [], pref);
    
    % Form the new integrand:
    integrand = g.onefun.*forDer;
    
    % Simplify the exponents:
    if ( isa(integrand, 'singfun') )
        integrand = cancelExponents(integrand);
    end
    
    % Call the sum at ONEFUN level.
    out = sum(integrand);
    
else
    error('CHEBFUN:UNBNDFUN:IrrecognizableInput', ...
        'The input can not be recognized.');
end

end