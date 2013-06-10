function g = restrict(f, s)
%RESTRICT Restrict a BNDFUN to a subinterval.
%   RESCTRICT(F, S) returns a BNDFUN that is restricted to the subinterval
%   [S(1), S(2)] of F.domain.
%
%   If length(S) > 2, i.e., S = [S1, S2, S3, ...], then RESCTRICT(F, S) returns
%   a cell-array of BNDFUN objects, where the cells contain F restricted to each
%   of the subintervals defined by S.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    g = f;
    return
end

% Check if subint is actually a subinterval:
if ( s(1) < f.domain(1) || s(end) > f.domain(2) || any(diff(s) <= 0) )
    error('BNDFUN:restrict:badinterval', 'Not a valid interval.')
elseif ( numel(s) == 2 && all(s == f.domain) )
    % Nothing to do here!
    g = f;
    return
end

% Compute the breaks in [-1,1] space and restrict the onefun:
t = f.mapping.inv(s);

% Restrict the ONEFUN field of f.
restrictedOnefuns = restrict(f.onefun, t);

% In case S = [S1, S2], i.e., we are only restricting to a one subinterval,
% restrictedOnefuns will be a CHEBTECH. In case S = [S1, S2, ..., SN], i.e., we
% are restricting to multiple subinterval, restrictedOnefuns will be a
% cell-array of CHEBTECH objects. We need to treat the cases differently, as we
% also want to either return a BNDFUN or cell-array of BNDFUN objects depending
% on how many subintervals we are restricting to.
if ( length(s) == 2 )
    % Only restricting to one subinterval -- return a BNDFUN.
    
    % g = bndfun(restrictedOnefuns, s); [TODO]: This line is planned to be removed.
    g = bndfun();
    g.onefun = restrictedOnefuns;
    g.domain = s;
    g.mapping = bndfun.createMap(s);
else
    % Restricting to multiple subintervals -- return a cell-array of BNDFUN
    % objects.
    
    % Create a cell to be returned.
    g = cell(1, numel(s) - 1);
    
    % Loop over each of the new subintervals, make a bndfun with new mapping,
    % and store in the cell returned:
    for k = 1:(numel(s) - 1)
        g_tmp = bndfun();
        g_tmp.onefun = restrictedOnefuns{k};
        g_tmp.domain = s(k:k+1);
        g_tmp.mapping = bndfun.createMap(s(k:k+1));
        g{k} = g_tmp;
        % g{k} = bndfun(restrictedOnefuns{k}, s(k:k+1)); [TODO]: planned to be removed.
    end
end

end

