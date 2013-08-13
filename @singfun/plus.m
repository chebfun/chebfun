function s = plus(f,g)
%PLUS Add SINGFUNS F and G.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% if one of the arguments is a double
if ( isa(f,'double') )
    aDouble = f;
    f = zeroSingFun();
    f.smoothPart = chebtech.constructor(aDouble, [], [], []);
elseif isa(g,'double')
    aDouble = g;
    g = zeroSingFun();
    g.smoothPart = chebtech.constructor(aDouble, [], [], []);
end

fExps = f.exponents;
gExps = g.exponents;
tol = singfun.pref.singfun.eps;

if ( all(abs(fExps-gExps) < tol ) )
    % Case 1: Exponents exactly alike. Just add the smooth parts.
    s = f;
    s.smoothPart = f.smoothPart + g.smoothPart;
    if ( iszero(s.smoothPart) )
       s = singfun.zeroSingFun();     
    end
elseif ( all(abs(round(fExps-gExps) - (fExps-gExps) ) < tol) )
    % Case 2: Both exponents differ by integers. Factor out the common
    % singular parts to leave the sum of smooth quotients.
    
    % At each endpoint, the smaller exponent will be factored out
    % of both summands.
    
    % Start off each compensating quotient as 1.
    factorF = @(x) 1;
    factorG = @(x) 1;
    for side = 1:2
        % The smaller of the two exponents is the exponent of the sum.
        [e,k] = sort([ fExps(side),gExps(side)] );
        newExps(side) = e(1);
        
        % The quotient factor is the difference in the exponents.
        if ( side == 1 )
            newFactor = @(x) (1+x).^diff(e);
        else
            newFactor = @(x) (1-x).^diff(e);
        end
        
        % Who had the smaller exponent? The other one gets the factor.
        if ( k(1) == 1 )
            factorG = @(x) factorG(x).*newFactor(x);
        else
            factorF = @(x) factorF(x).*newFactor(x);
        end
    end
    
    % Construct the function handle for the new smooth fun
    sF = f.smoothPart;
    sG = g.smoothPart;
    smoothOp = @(x) feval(sF, x).*factorF(x) + feval(sG, x).*factorG(x);
    
    % Construct the new smooth fun
    s = singfun.zeroSingFun();
    smoothPrefs = chebtech.pref('tech', 'cheb1', 'extrapolate', false);
    vscale = [];
    hscale = [];
    s.smoothPart = chebtech.constructor(smoothOp, vscale, hscale, smoothPrefs);
    
    % Assing new exponents
    s.exponents = newExps;
    
else
    % Case 3: Nontrivial exponent difference.       
    % Form a new function handle for the sum from F and G.    
    
    % Retrieve function handle of F
    s1 = f.smoothPart;
    a1 = f.exponents(1);
    b1 = f.exponents(2);    
    op1 = @(x) feval(s1, x).*(1+x).^a1.*(1-x).^b1;
    
    % Retrieve function handle of G
    s2 = g.smoothPart;
    a2 = g.exponents(1);
    b2 = g.exponents(2);    
    op2 = @(x) feval(s2, x).*(1+x).^a2.*(1-x).^b2;
    
    % Define a function handle for the sum
    op = @(x) op1(x) + op2(x);
    
    % construct a new SINGFUN for sum
    s = singfun( op, [], {'sing', 'sing'}, singfun.pref );
end

end

%% LEGACY CODE
% This is the solution used for case 3 in Chebfun v4. It relies on
% coordinate mappings. It's included here so that it becomes part of the
% version history before deletion. 
% 
%     pref = chebfunpref;
%     if pref.splitting
%         pref.splitdegree = 8*pref.splitdegree;
%     end
%     pref.resampling = false;
%     pref.blowup = 1;
%     
%     scl.h = max(g1.scl.h,g2.scl.h);
%     scl.v = max(g1.scl.v,g2.scl.v);
%     
%     % Choose the correct singmap
%     dexps = fExps - gExps;
%     newexps = [0 0];
%     pows = [0 0 ];       % will be the powers in the sing map
%     lr = 0;
%     % left
%     if round(dexps(1)) ~= dexps(1) % then trouble at the left
%         lr = -1;        % flag for sing.m
%         expsl = sort([fExps(1) gExps(1)]);
%         if expsl(1) < 0 % ==> blow up, so use exponents in new representation
%             newexps(1) = expsl(1);
%             pows(1) = expsl(2)-expsl(1);
%         else            % ==> no blow up, so use largest power
%             pows(1) = expsl(2);
%         end
%     else
%         newexps(1) = min(fExps(1),gExps(1));
%         % We could probably do something cleverer like extracting out the blowup
%         % by hand, but for now we just let the constructor do it.
%     end
%     % right (as above)
%     if round(dexps(2)) ~= dexps(2)
%         lr = lr + 1;
%         expsr = sort([fExps(2) gExps(2)]);
%         if expsr(1) < 0
%             newexps(2) = expsr(1);
%             pows(2) = expsr(2)-expsr(1);
%         else
%             pows(2) = expsr(2);
%         end
%     else
%         newexps(2) = min(fExps(2),gExps(2));
%     end
%     pows = pows-floor(pows); % Should be < 1;
%     
%     % The new map
%     map = maps(fun,{'sing',pows},ends);
%     % The new exponents
%     pref.exps = [newexps(1) newexps(2)];
%     pref.sampletest = 0;
%     
%     % Call the fun constructor
%     g1 = fun(@(x) feval(g1,x)-feval(g2,x),map,pref,scl);
%     if ~g1.ish
%         warning('FUN:minus:failtoconverge','Operation may have failed to converge');
%     end
%     
%     if any(g1.exps < 0)
%         g1 = checkzero(g1);
%         g1 = extract_roots(g1);
%     end
%     
% end
% 
% 
% function g1 = checkzero(g1)
% % With exps, if the relative deifference is O(eps) we set it to zero.
% % Same goes for unbounded domains.
% if all(abs(g1.vals) < 100*g1.scl.v*chebfunpref('eps'))
%     g1.vals = 0;
%     g1.coeffs = 0;
%     g1.n = 1;
%     g1.exps = [0 0];
%     g1.scl.v = 0;
% end
% 
% end



