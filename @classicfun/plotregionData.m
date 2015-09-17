function data = plotregionData(u, userEps, numpts)
%PLOTREGIONDATA   Useful data for plotting regions of analyticity.
%   PLOTREGIONDATA(U) returns a struct containing data that can be used for
%   plotting the estimate region of analyticity of U.
%
%   PLOTREGIONDATA(U, EPS) allows a user-specified EPS.
%
%   PLOTREGIONDATA(U, EPS, NUMPTS) uses N points to represent the boundary.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    numpts = 500;
end

if ( nargin < 2 )
    userEps = eps;
end

% Fetch the data from the ONEFUN of U:
data = plotregionData(u.onefun, userEps, numpts);

end
