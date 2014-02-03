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
    pref = chebpref();
end

% Get the function values at the end of the domain. Note that the end point of
% the domain can be infinite.
vends = [get(g, 'lval'), get(g, 'rval')];

% Get the epslevel and vscale of the function g.
tol = 1e1*get(g, 'epslevel')*get(g, 'vscale');

%%
% Cases for onefun has no singularities.

if ( ~issing(g) || ( issing(g) && ( ~any( get(g, 'exponents') ) ) ) )
    
    % A dirty checklist to see if the integrand is integrable or not. This
    % checklist may potentially be removed, since SINGFUN is supposed to be able
    % to handle functions with positive singularities.
    
    % TODO: how is this possible is ~issing(g)?
    % Check 1: Check if the function values are vanishing at infinity/ties.
    unbounded = [];
    if ( isinf(dom(1)) && abs(vends(1)) > tol )
        unbounded(1) = sign(vends(1))*inf;
    end
    if ( isinf(dom(2)) && abs(vends(2)) > tol )
        unbounded(2) = sign(vends(2))*inf;
    end
    if ( ~isempty(unbounded) )
        out = sum(unbounded);
        return
    end
    
    % Check 2: Check for the decay speed at infinity/ties.
    % The integrand is integrable only when it decays faster than 1/x towards
    % infinity/ties.
    
    % [NOTE]: Temporarily this is commented out since the following snippet
    % which is copied from Chebfun v4 is not doing what it is supposed to do!
    
    %     % Get the points at which the function is sampled.
    %     y = g.onefun.points();
    %     unbounded = [];
    %     if ( isinf(dom(2)) )
    %         gval = g.onefun.values;
    %         gval = gval./(1-y);   % peel off a factor of 1/x
    % %       gval = extrapolate(gval);
    %         if ( abs(gval(end)) > 1e3*tol && ...
    %                 diff(gval((end-1:end))/diff(y(end-1:end))) > -g.onefun.vscale/g.onefun.hscale )
    %             unbounded(1) = inf*sign(g.onefun.values(end-1));
    %         elseif ( abs(gval(end)) > tol )
    %             warning('FUN:sum:slowdecay', ...
    %                 'Result is likely inaccurate. (Slow decay).')
    %         end
    %     end
    %     if ( isinf(dom(1)) )
    %         gval = g.onefun;
    %         gval.values = gval.values./(1+y);
    % %         gval = extrapolate(gval);
    %         if ( abs(gval(1)) > 1e3*tol && ...
    %                 diff(gval(1:2)./diff(y(1:2))) < g.onefun.vscale/g.onefun.hscale )
    %             unbounded(2) = inf*sign(g.onefun.values(2));
    %         elseif ( abs(gval(1)) > tol )
    %             warning('FUN:sum:slowdecay', ...
    %                 'Result is likely inaccurate. (Slow decay).')
    %         end
    %     end
    %     if ( ~isempty(unbounded) )
    %         out = sum(unbounded);
    %         return
    %     end
    
    % If we reach here, the function decays sufficiently fast that we can integrate
    % by the onefun sum.
    
    % Construct the ONEFUN presentation of the derivative of the map.
    pref.singPrefs.exponents = g.mapping.forDerExps;
    forDer = onefun.constructor(@(x) g.mapping.forDer(x), [], [], pref);
    
    % Peel off the boundary roots to cancel off the negative exponents of forDer:
    h = extractBoundaryRoots(g.onefun, -g.mapping.forDerExps.');
    
    % Form the new integrand:
    integrand = h.*forDer.smoothPart;
    
    % Call the sum at onefun level.
    out = sum(integrand);
    
    return
    
    %%
    % Cases for onefun has singularities at the end points.
elseif ( issing(g) )
    
    % Now assume that g.onefun is a SINGFUN: g.onefun = (1+x)^a*(1-x)^b*s(x),
    % where s(x) is the smooth part and a and b are the powers of the zeros or
    % orders of the poles/fractional poles. The UNBNDFUN g is integrable only if
    % both of the following conditions are met: 
    %  1. For a finite boundary, the corresponding exponent of g.onefun is
    %     larger than -1.
    %  2. For an infinite boundary, the corresponding exponent of g.onefun is
    %     larger than 1.
    
    % TODO: The what of the what!?
    % Get the exponents of the singfun.
    exps = g.exponents;
    
    % Locate the finite boundary and the infinite boundary/ies.
    infMask = isinf(dom);
    fntMask = isfinite(dom);
    
    if ( ( isempty(fntMask) && all( exps > 1) ) || ...
            ( ~isempty(fntMask) && exps(fntMask) > -1 && exps(infMask) > 1 ) )
        
        % This if branch covers the following cases, in which g is integrable:
        % 1. singly-infinite domain [d inf] with exponents [a b] where a > -1,
        %    b > 1.
        % 2. singly-infinite domain [inf d] with exponents [a b] where a > 1,
        %    b > -1.
        % 3. doubly-infinite domain [-inf inf] with exponents [a b] where a > 1,
        %    b > 1.
        
        % Construct the onefun presentation of the derivative of the map.
        pref.singPrefs.exponents = g.mapping.forDerExps;
        forDer = onefun.constructor(@(x) g.mapping.forDer(x), [], [], pref);
        
        % Peel off the boundary roots to cancel off the negative exponents of forDer:
        h = extractBoundaryRoots(g.onefun, -g.mapping.forDerExps.');
        
        % Form the new integrand:
        integrand = h.*forDer.smoothPart;
        
        % Call the onefun sum.
        out = sum(integrand);
        
    elseif ( isempty(fntMask) && ~all( exps <= 1 ) )
        
        % This elseif condition covers the following case for which the
        % integral is infinite due to the non-integrablity at one of the
        % boundaries at infinity:
        
        % 1. doubly-infinite domain [-inf inf] with exponents [a b] where a > 1,
        %    b <= 1 or a <= 1, b > 1.
        
        % Find the boundary at which the integrand is not integrable.
        blowMask = exps <= 1;
        
        % Set the infinite the correct sign.
        out = sign(vends(blowMask))*Inf;
        
    elseif ( ~isempty(fntMask) && ( exps(fntMask) <= -1 ) && ( exps(infMask) > 1 ) )
        
        % This elseif condition covers the following cases for which the
        % integral is infinite due to non-integrability at the finite boundary:
        % 1. singly-infinite domain [d inf] with exponents [a b] where a <= -1,
        %    b > 1.
        % 2. singly-infinite domain [inf d] with exponents [a b] where a > 1,
        %    b <= -1.
        
        % Set the infinite the correct sign.
        out = sign(vends(fntMask))*Inf;
        
    elseif ( ~isempty(fntMask) && ( exps(fntMask) > -1 ) && ( exps(infMask) <= 1 ) )
        
        % This elseif condition covers the following cases for which the
        % integral is infinite due to non-integrability at the infinite boundary:
        % 1. singly-infinite domain [d inf] with exponents [a b] where a > -1,
        %    b <= 1.
        % 2. singly-infinite domain [inf d] with exponents [a b] where a <= 1,
        %    b > -1.
        
        % Set the infinite the correct sign.
        out = sign(vends(infMask))*Inf;
        
    elseif ( sign(vends(1)) == sign(vends(2)) )
        % This elseif condition covers the following cases for which the
        % integral is infinite:
        % 1. singly-infinite domain [d inf] with exponents [a b] where a <= -1,
        %    b <= 1 and sign(s(-1)) = sign(s(1)).
        % 2. singly-infinite domain [inf d] with exponents [a b] where a <= 1,
        %    b <= -1 and sign(s(-1)) = sign(s(1)).
        % 3. doubly-infinite domain [-inf inf] with exponents [a b] where a <= 1,
        %    b <= 1 and sign(s(-1)) = sign(s(1)).
        
        out = sign(vends(1))*Inf;
        
    else
        % This else condition covers the following cases for which the
        % integral doesn't exist due to the opposite signs of the function values
        % at the boundaries where the integrand is not integrable:
        % 1. singly-infinite domain [d inf] with exponents [a b] where a <= -1,
        %    b <= 1 and sign(s(-1)) ~= sign(s(1)).
        % 2. singly-infinite domain [inf d] with exponents [a b] where a <= 1,
        %    b <= -1 and sign(s(-1)) ~= sign(s(1)).
        % 3. doubly-infinite domain [-inf inf] with exponents [a b] where a <= 1,
        %    b <= 1 and sign(s(-1)) ~= sign(s(1)).
        
        out = NaN;
        
    end
    
    return
    
else
    error('CHEBFUN:UNBNDFUN:IrrecognizableInput',...
        'The input can not be recognized.');
end

end