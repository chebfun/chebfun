function g = restrict(f, s)
%RESTRICT Restrict an UNBNDFUN to a subinterval.
%   RESCTRICT(F, S) returns a FUN object that is restricted to the subinterval
%   [S(1), S(2)] of the domain of F.
%
%   If length(S) > 2, i.e., S = [S1, S2, ..., Sn], then RESCTRICT(F, S) returns
%   an array of FUN objects, where the cells contain F restricted to each of 
%   the subintervals defined by S. If there is only one FUN to be returned,
%   that is, if length(S) == 2, then the FUN object g is returned. This 
%   facilitates the use of the result by other functions, e.g. plot etc.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    g = f;
    return
end

% Check if the argument s is actually a subinterval:
if ( s(1) < f.domain(1) || s(end) > f.domain(2) || any(diff(s) <= 0) )
    error('CHEBFUN:UNBNDFUN:restrict:badInterval', 'Not a valid interval.')
elseif ( numel(s) == 2 && all(s == f.domain) )
    % Nothing to do here!
    g = f;
    return
end

% Number of subdomains:
numSubDom = numel(s) - 1;

% Preallocate the output cell:
g = cell(1, numSubDom);

% Any exponents?
exps = [];

if ( issing(f) )
    
    % Grab the exponents:
    exps = get(f, 'exponents');
        
    % If there is a non-trivial exponent:
    if ( any(exps) )
        
        % Interior exponents:
        interiorExps = zeros(1, numSubDom-1);  
        
        % Insert the interior exponents:
        exps = [exps(1) interiorExps exps(2)];
        
        % Negate the exponents for infinite endpoint, since the exponents stored
        % in SINGFUN are the negated values of those supplied to the UNBNDFUN
        % constructor:
        
        if ( s(1) == -Inf )
            exps(1) = -exps(1);
        elseif ( f.domain(1) == -Inf )
            exps(1) = 0;
        end
        
        if ( s(end) == Inf )
            exps(end) = -exps(end);
        elseif ( f.domain(end) == Inf )
            exps(end) = 0;
        end

    end
    
end

% Grab the vscale:
vscale = get(f, 'vscale');

% Loop over each subdomain:
for k = 1:numSubDom
    
    % Initialize preferences:
    pref = chebfunpref();
    
    % Pass the information about exponents through preference:
    if ( ~isempty(exps) && any(exps(k:k+1)) )
        data.exponents = exps(k:k+1);
        
        if ( (k == 1) || (k == numSubDom) )
            pref.techPrefs.extrapolate = 1;
        end
    else
        data.exponents = [];
    end
    
    data.domain = s(k:k+1);
    data.vscale = vscale;
    g{k} = classicfun.constructor(@(x) feval(f, x), data, pref);
end

% When there is only one cell, return the UNBNDFUN instead:
if ( numSubDom == 1 )
    g = g{1};
end

end
