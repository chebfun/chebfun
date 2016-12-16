function I = sum3(f)
%SUM3   Triple integral of a CHEBFUN3 over its domain.
%   I = SUM3(F) is the definite integral of F over its domain. The output
%   is a scalar.
%
% See also CHEBFUN3/SUM, CHEBFUN3/SUM2, CHEBFUN3/CUMSUM, CHEBFUN3/CUMSUM2 
% and CHEBFUN3/CUMSUM3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) ) 
    I = [];
    return
end

I = chebfun3.txm( ...
        chebfun3.txm( ...
            chebfun3.txm(f.core, sum(f.cols), 1), ... % integrals wrt x
        sum(f.rows), 2), ...                          % integrals wrt y
    sum(f.tubes), 3);                                 % integrals wrt z

end