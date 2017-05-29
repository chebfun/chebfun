function r = roots( f )
%ROOTS      Zero contours of a SPHEREFUN
%   R = ROOTS(F), returns the zero contours of F as a quasimatrix of 
%   array-valued chebfuns. Each column of R is one zero contour. This 
%   command only finds contours when there is a change of sign and it may 
%   group intersecting contours in a non-optimal way. Contours are computed
%   to, roughly, four digits of precision. In particular, this command 
%   cannot reliably compute isolated real roots of F. 
%
% See also CHEBFUN2/ROOTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If the spherefun is empty, then there are no roots.
if ( isempty( f ) )
    r = []; 
    return
end

% Convert f to a chebfun2 on a larger domain so that we can properly
% pick up roots that go accross -pi to pi in lambda.  This will allow us 
% to account for the inherent periodicity of the sphere.
rts = roots@separableApprox( chebfun2(f,[-pi-2*pi/100 pi 0 pi]) );

% If no roots were found, then there is nothing more to do
if ( isempty( rts ) )
    r = [];
    return
end

% Now make into a collection of array-valued chebfuns ready for plotting on
% the sphere.  Make sure we have enough points to sample the contours.
x = chebpts(max(length(rts),17) + 1);

vals = feval(rts, x);
r = cell(size(vals,2), 1);

% Go through each component and make it an array-valued chebfun: 
cnt = 1;
for k = 1:size(vals, 2)
    comp = feval(rts(:, k), x); 
    lam = real(comp);
    th = imag(comp);
    % We only want to return roots that are in the domain of [-pi pi 0 pi].
    idlam = lam <= pi & lam >= -pi;
    idth = th <= pi & th >= 0;
    id = idlam & idth;
    if any(id)
        % Convert to a curve on the sphere.
        AX = cos(lam).*sin(th);
        AY = sin(lam).*sin(th);
        AZ = cos(th);
        % If there is a contour right on the poles then don't include it
        % as a contour.
        if norm(abs(AZ)-1,inf) > 1e4*eps
            r{cnt} = chebfun([AX, AY, AZ]);
            cnt = cnt + 1;
        end
    end
end
% Remove the empty elements of r.
r = r(1:cnt-1);

end