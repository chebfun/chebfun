function g = squeeze(f)
%SQUEEZE   Squeeze a CHEBFUN3 to two variables or to one variable, if 
%   possible.
%   G = squeeze(F) returns a CHEBFUN3 if F depends on all of the three 
%   variables x, y and z. If F depends only any two of the variables then
%   a CHEBFUN2 is returned. If F depends only on one of the variables x, y 
%   or z, then a column CHEBFUN is returned.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )  % Check for an empty CHEBFUN3.
    g = chebfun();    % Return an empty CHEBFUN because we are squeezing.
    return
end

% Get the low rank representation for f.
[core, cols, rows, tubes] = tucker(f);
dom = f.domain;
[r1, r2, r3] = rank(f);

% Simplify rows and cols:
cols = simplify(cols, [], 'globaltol');
rows = simplify(rows, [], 'globaltol');
tubes = simplify(tubes, [], 'globaltol');

if ( r1 == 1 && r2 == 1 && r3 == 1)
% If f is of rank 1, then it MAY be a function of just ONE variable:
    if ( length(cols) == 1 &&  length(rows) == 1 )
    % If cols and rows are constant then we have a function of z only.
        newCore = squeeze(chebfun3.txm(chebfun3.txm(core, cols.coeffs, 1), ...
             rows.coeffs, 1));
     % Note that cols.coeffs = cols.coeffs(1) is just a fixed constant scalar.
        newdomain = dom(5:6);
        g = chebfun(newCore*tubes, newdomain); 
    elseif ( length(cols) == 1 &&  length(tubes) == 1 )
        % If cols and tubes are constant then we have a function of y only.
        newCore = squeeze(chebfun3.txm(chebfun3.txm(core, cols.coeffs, 1), ...
             tubes.coeffs, 3));
        newdomain = dom(3:4);
        g = chebfun(newCore*rows, newdomain); 

    elseif ( length(rows) == 1 &&  length(tubes) == 1 )
        % If rows and tubes are constant then we have a function of x only.
        newCore = squeeze(chebfun3.txm(chebfun3.txm(core, rows.coeffs, 2), ...
             tubes.coeffs, 3));
        newdomain = dom(1:2);
        g = chebfun(newCore*cols, newdomain);
    end

elseif ( r1 == 1 && length(cols) == 1)
% If f has the x-rank 1, and cols are constant then we have a function of y
% and z.
    newCore = squeeze(chebfun3.txm(core, cols.coeffs, 1));
    newdomain = dom(3:6);
    g = chebfun2(rows*newCore*tubes.', newdomain);

elseif ( r2 == 1 && length(rows) == 1 )
% If f has the y-rank 1 and rows are constant then we have a function of x 
% and z.
    newCore = squeeze(chebfun3.txm(core, rows.coeffs, 2));
    newdomain = [dom(1:2) dom(5:6)];
    g = chebfun2(cols*newCore*tubes.', newdomain);

elseif ( r3 == 1 && length(tubes) == 1 )
% If f has z-rank 1 and tubes are constant then we have a function of x and
% y.
    newCore = squeeze(chebfun3.txm(core, tubes.coeffs, 3));
    newdomain = dom(1:4);
    g = chebfun2(cols*newCore*rows.', newdomain); 
else
    g = f;
end

end