function s = plus(f, g)
%+   Addition of SINGFUN objects with SINGFUNs and SMOOTHFUNs.
%   F + G adds F and G, where F and G may be SINGFUN objects or scalars.
%
% See also MINUS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% If one of the arguments is empty:
if ( isempty(f) || isempty(g) )
    % Create an empty SINGFUN and return:
    s = singfun();
    return
end

% Check if inputs are other than SINGFUNS, SMOOTHFUNS or doubles:
if ( (~isa(f, 'singfun') && ~isa(f, 'smoothfun') && ~isa(f, 'double')) || ...
     (~isa(g, 'singfun') && ~isa(g, 'smoothfun') && ~isa(g, 'double')) )
    error('CHEBFUN:SINGFUN:plus:badInput' , ...
        'Input can only be a SINGFUN, a SMOOTHFUN or a double')
end
% One of the arguments i.e. f or g is necessarily a SINGFUN object. Otherwise, 
% this overloaded plus would not have been called.

% If one of the arguments is a double, upgrade it to a SINGFUN:
if ( isa(f, 'double') )
    
    if ( f == 0 )
        % If F is zero, return G.
        s = g;
        return
    end
    
    % Make a SMOOTHFUN of the double f:
    f = g.smoothPart.make(f);
    % Convert f to a SINGFUN:
    f = singfun.smoothFun2SingFun(f);
elseif ( isa(g, 'double') )
    
    if ( g == 0 )
        % If G is zero, return F.
        s = f;
        return
    end
    
    % Make a SMOOTHFUN of the double g:
    g = f.smoothPart.make(g);
    % Convert g to a SINGFUN:
    g = singfun.smoothFun2SingFun(g);    
end

% If one of the arguments is a SMOOTHFUN, upgrade it to a SINGFUN:
if ( isa(f, 'smoothfun') )
    f = singfun.smoothFun2SingFun(f);    
elseif ( isa(g, 'smoothfun') )
    g = singfun.smoothFun2SingFun(g);    
end

fExps = f.exponents;
gExps = g.exponents;
tolExps = chebfunpref().blowupPrefs.exponentTol;
tolSmth = 1e2*eps;

%%
if ( all(abs(fExps - gExps) < tolExps) )
    % Case 1: Exponents exactly same. Just add the smooth parts.
    s = f;
    s.smoothPart = f.smoothPart + g.smoothPart;
    if ( normest(s.smoothPart) < tolSmth )
       s = singfun.zeroSingFun();     
    end
    
elseif ( all(abs(round(fExps - gExps) - (fExps - gExps)) < tolExps) )
    % Case 2: Both exponents differ by integers. Factor out the more singular
    % exponent to leave the sum of smooth quotients.
    
    % At each endpoint, the algebraically smaller exponent will be factored out
    % of both summands. For example, if we are adding (1+x).^-4 and (1+x).^-1,
    % we first factor out (1+x).^-4:
    %   (1+x).^-4 + (1+x).^-1 = (1+x).^-4 .* ( 1 + (1+x).^3 )
    
    % Start off each compensating quotient as 1.
    factorF = @(x) 1;
    factorG = @(x) 1;
    newExps = zeros(1, 2);    
    for side = 1:2
        % The algebraically smaller of the two exponents is the exponent of 
        % the sum.
        [e, k] = sort([fExps(side), gExps(side)]);
        newExps(side) = e(1);
        
        % The quotient factor is the difference in the exponents.
        if ( side == 1 )
            newFactor = @(x) (1 + x).^diff(e);
        else
            newFactor = @(x) (1 - x).^diff(e);
        end
        
        % Who had the algebraically smaller exponent? The other one gets the 
        % factor.
        if ( k(1) == 1 )
            factorG = @(x) factorG(x).*newFactor(x);
        else
            factorF = @(x) factorF(x).*newFactor(x);
        end
    end
    
    % Construct the function handle for the new smooth fun:
    sF = f.smoothPart;
    sG = g.smoothPart;
    smoothOp = @(x) feval(sF, x).*factorF(x) + feval(sG, x).*factorG(x);
    
    % Construct the new smooth fun:
    s = singfun.zeroSingFun();
    s.smoothPart = singfun.constructSmoothPart(smoothOp, [], []);
    
    % Assign new exponents:
    s.exponents = newExps;
    
else
    % Case 3: Nontrivial difference in the exponents of F and G. Form a new
    % function handle for the sum from F and G.
    
    warning('CHEBFUN:SINGFUN:plus:exponentDiff', ...
        ['Non-integer difference in the exponents of the two SINGFUN ' ...
        'objects: The result may not be accurate.']);
    
    % Define a function handle for the sum:
    op = @(x) feval(f, x) + feval(g, x);
    
    % The new scales for the sum:
    vScale = get(f, 'vscale') + get(g, 'vscale');
    hScale = get(f, 'hscale');
    
    % Take the smallest exponents to be those for the summation:
    exps = [get(f, 'exponents'); get(g, 'exponents')];
    exps = min(exps);
    
    % Construct a new SINGFUN for the sum:
    data.exponents = exps;
    data.vscale = vScale;
    data.hscale = hScale;
    s = singfun(op, data, []);
end

%%
% Check if after addition s has become smooth:
if ( issmooth(s) )
    s = s.smoothPart;
end

end
