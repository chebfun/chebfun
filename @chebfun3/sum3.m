function I = sum3(f)
%SUM3   Triple integral of a CHEBFUN3 over its domain.
%
%   I = sum3(F) is the definite integral of f over its domain. The output
%   is a scalar.
%
%   See also chebfun3/sum, chebfun3/sum2, chebfun3/cumsum, 
%   chebfun3/cumsum2 and chebfun3/cumsum3.

% Empty check: 
if ( isempty(f) ) 
    I = []; 
    return; 
end

I = chebfun3.txm(chebfun3.txm(chebfun3.txm(f.core, sum(f.cols), 1), ...
    sum(f.rows), 2), sum(f.tubes), 3);

end
