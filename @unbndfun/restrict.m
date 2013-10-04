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

% Compute the breaks in [-1,1] space and restrict the onefun:
t = f.mapping.inv(s);

% Make sure -Inf and Inf are mapped to -1 and 1 respectively:
mask = isinf(x); 
z(mask) = sign(x(mask));

% Restrict the ONEFUN of f.
restrictedOnefuns = restrict(f.onefun, t);

if ( length(s) == 2 )
    % Only restricting to one subinterval - return a BNDFUN or an UNBNDFUN.
    
    % Create an empty BNDFUN or UNBNDFUN, and assign fields directly. This is
    % faster than calling the FUN constructor.
    if any( isinf(s) )
        g = unbndfun();
    else
        g = bndfun();
    end
    g.onefun = restrictedOnefuns;
    g.domain = s;
    g.mapping = g.createMap(s);
else
    % Restricting to multiple subintervals - return a cell-array of FUN objects.
    
    % Create the cell to be returned.
    g = cell(1, numel(s) - 1);
    
    % Create an empty FUN:
    emptyBndfun = bndfun();
    emptyUnbndfun = unbndfun();
    
    % Loop over each of the new subintervals, make a fun with new mapping,
    % and store in the cell returned:
    for k = 1:(numel(s) - 1)
        % Assign fields directly to an empty temporary BNDFUN or UNBNDFUN. This
        % is faster than calling the FUN constructor.
        if any( isinf(s(k:k+1)) )
            gTemp = emptyUnbndfun;
        else
            gTemp = emptyBndfun;
        end
        gTemp.onefun = restrictedOnefuns{k};
        gTemp.domain = s(k:k+1);
        gTemp.mapping = gTemp.createMap(s(k:k+1));
        g{k} = gTemp;
    end

end