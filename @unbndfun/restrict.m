function g = restrict(f, s)
%RESTRICT Restrict an UNBNDFUN to a subinterval.
%   RESCTRICT(F, S) returns an FUN object that is restricted to the subinterval
%   [S(1), S(2)] of F.domain.
%
%   If length(S) > 2, i.e., S = [S1, S2, S3, ...], then RESCTRICT(F, S) returns
%   an array of FUN objects, where the cells contain F restricted to each of 
%   the subintervals defined by S. If there is only one FUN to be returned,
%   that is, length(S) == 2, then the FUN object g is returned. This 
%   facilitates the use of the result by other functions, e.g. plot etc.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    g = f;
    return
end

% Check if subint is actually a subinterval:
if ( s(1) < f.domain(1) || s(end) > f.domain(2) || any(diff(s) <= 0) )
    error('CHEBFUN:UNBNDFUN:restrict:badinterval', 'Not a valid interval.')
elseif ( numel(s) == 2 && all(s == f.domain) )
    % Nothing to do here!
    g = f;
    return
end

% Number of subdomains:
numSubDom = numel(s) - 1;

% Preallocate the output cell:
g = cell(1, numSubDom);

% Loop over each subdomain:
for k = 1:numSubDom
    pref = chebpref();
    if ( issing(f) )
        ind = isinf(s(k:k+1));
        exps = get(f, 'exponents');
        
        % Negate the exponents for infinite endpoint, since the exponents stored
        % in SINGFUN are the negated value of those supplied to the BNDFUN 
        % constructor:
        exps(ind) = -exps(ind);
        
        if ( k == 1 ) && ( logical(exps(1)) )
            pref.singPrefs.exponents = [exps(1) 0];
            pref.techPrefs.extrapolate = 1;
        end
        
        if ( k == numSubDom ) && ( logical(exps(2)) )
            pref.singPrefs.exponents = [0 exps(2)];
            pref.techPrefs.extrapolate = 1;
        end
    end
    
    g{k} = fun.constructor(@(x) feval(f, x), s(k:k+1), [], [], pref);
end

end