function V = vscale( F )
%VSCALE   Vertical scale of a BALLFUN.
% 
% VSCL = VSCALE(F) returns the vertial scale of a BALLFUN as determined
% by evaluating on a coarse tensor-product grid. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( F )
    V = [];
    return
end

% Convert to CFF values
F = ballfun.coeffs2vals(coeffs3(F,65,65,65));

% Return the maximum of the absolute values of F
V = max(abs(F(:)));
end