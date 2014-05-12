function f = fracCumSum(funs, fracM)
%FRACUMSUM   Indefinite integral of a BNDFUN for a fractional order.
%  FRACUMSUM(FUNS, FRACM) returns the fractional indefinite integral of order 
%  FRACM defined by the piecewise smooth CLASSICFUNs stored in the cell FUNS and 
%  its domain is same as the CLASSICFUN given in the last entry of FUNS. Note 
%  that FRACM needs to be non-integer (fractional).
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Grab the ONEFUNs and the MAPPINGs of FUNS:
numFuns = numel(funs);
oneFuns = cell(1, numFuns);
maps = cell(1, numFuns);

for k = 1:numFuns
    oneFuns{k} = funs{k}.onefun;
    maps{k} = funs{k}.mapping;
end

% Copy the last entry of FUNS to the result F, since F takes the same domain and 
% mapping:
f = funs{end};

% Pass the data of ONEFUNs and MAPPINGs to FRACCUMSUM@SINGFUN and rescale. The
% result is the ONEFUN of F:
f.onefun = (diff(funs{end}.domain)/2)^fracM*...
    singfun.fracCumSum(oneFuns, maps, fracM);

end