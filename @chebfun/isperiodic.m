function out = isperiodic(f)
%ISEMPTY   Test for periodic CHEBFUN.
%   ISPERIODIC(F) returns logical true if F is based on a TRIGTECH object 
%   and false otherwise.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check the empty case:
if ( isempty(f) )
    out = true;
    return
end

if ( numel(f) == 1 )
    % Array-valued CHEBFUN case:
    out = columnIsPeriodic(f);
else
    % QUASIMATRIX case:
    f = mat2cell(f);
    % Loop over the columns:
    for k = 1:numel(f)
        out = columnIsequal(f{k});
        if ( ~out ), break, end
    end
end

end

function out = columnIsPeriodic(f)

% Periodic chebfuns should have just a sing fun.
% Check the number of FUNs and the periodicity:
if ( isperiodic(f.funs{1}) && numel(f.funs) < 2 )
    out = true;
else
    out = false;
end

end
