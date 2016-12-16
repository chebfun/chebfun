function g = restrict(f, s)
%RESTRICT   Restrict a BNDFUN to a subinterval.
%   RESTRICT(F, S) returns a BNDFUN that is restricted to the subinterval
%   [S(1), S(2)] of F.domain.
%
%   If LENGTH(S) > 2, i.e., S = [S1, S2, S3, ...], then RESTRICT(F, S) returns
%   a cell-array of BNDFUN objects, where the cells contain F restricted to each
%   of the subintervals defined by S.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    g = f;
    return
end

% Check if s is actually a subinterval:
hs = diff(f.domain)*eps; % Some wiggle room.
if ( s(1) < f.domain(1) )
    if ( abs(s(1) - f.domain(1)) < hs )
        s(1) = f.domain(1);
    else
        error('CHEBFUN:BNDFUN:restrict:badInterval', 'Not a valid interval.')
    end
end
if ( s(end) > f.domain(2) )
    if ( abs(s(end) - f.domain(2)) < hs )
        s(end) = f.domain(2);
    else
        error('CHEBFUN:BNDFUN:restrict:badInterval', 'Not a valid interval.')
    end
end
if ( any(diff(s) <= 0) )
    error('CHEBFUN:BNDFUN:restrict:badInterval', 'Not a valid interval.')
end

if ( (numel(s) == 2) && all( s == f.domain) )    
    % Nothing to do here!
    g = f;
    return
end

% Compute the breaks in [-1,1] space and restrict the onefun:
t = f.mapping.Inv(s);

% Restrict the ONEFUN field of f.
restrictedOnefuns = restrict(f.onefun, t);

% In case S = [S1, S2], i.e., we are only restricting to one subinterval,
% restrictedOnefuns will be a ONEFUN. In case S = [S1, S2, ..., SN], i.e., we
% are restricting to multiple subinterval, restrictedOnefuns will be a
% cell-array of ONEFUN objects. We need to treat the cases differently, as we
% also want to either return a BNDFUN or cell-array of BNDFUN objects depending
% on how many subintervals we are restricting to.
if ( length(s) == 2 )
    % Only restricting to one subinterval -- return a BNDFUN.
    
    % Create an empty BNDFUN, and assign fields directly. This is faster than
    % using the BNDFUN constructor.
    g = bndfun();
    g.onefun = restrictedOnefuns;
    g.domain = s;
    g.mapping = bndfun.createMap(s);
else
    % Restricting to multiple subintervals -- return a cell-array of BNDFUN
    % objects.
    
    % Create a cell to be returned.
    g = cell(1, numel(s) - 1);
    
    % Create an empty BNDFUN:
    emptyBndfun = bndfun();
    
    % Loop over each of the new subintervals, make a bndfun with new mapping,
    % and store in the cell returned:
    for k = 1:(numel(s) - 1)
        % Assign fields directly to an empty temporary BNDFUN. This is faster
        % than using the BNDFUN constructor.
        gTemp = emptyBndfun;
        gTemp.onefun = restrictedOnefuns{k};
        gTemp.domain = s(k:k+1);
        gTemp.mapping = bndfun.createMap(s(k:k+1));
        g{k} = gTemp;
    end
end

end

