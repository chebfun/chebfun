function f = toFunctionIn(disc, coeffs)
%TOFUNCTION   Convert discrete values of a TRIGSPEC to a CHEBFUN.
%   F = TOFUNCTIONIN(DISC, COEFFS) converts the coeffs returned by TRIGSPEC to a
%   CHEBFUN. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Remove the FILPUD once Moshin's branch has been merged.

dom = disc.domain; % Domain.
c = mat2cell(flipud(coeffs), disc.dimension); 
funs = cell(numel(c), 1);  
for k = 1:numel(c)
    % Construct TRIGTECH object from each piece.
    ct = trigtech({[], c{k}});
    % Assign each piece to a subinterval with a BNDFUN.
    funs{k} = bndfun(ct, struct('domain', dom(k: k + 1)));
end
% Conver the FUNS cell-array to a CHEBFUN.
f = chebfun(funs);
f = real(f);

end
