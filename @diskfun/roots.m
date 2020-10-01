function r = roots(varargin)
%ROOTS   Zero contours of a DISKFUN.
%   R = ROOTS(F) returns the zero contours of F as a quasimatrix of 
%   array-valued chebfuns. Each column of R is one zero contour. This 
%   command only finds contours when there is a change of sign and it may 
%   group intersecting contours in a non-optimal way. Contours are computed
%   to, roughly, four digits of precision. In particular, this command 
%   cannot reliably compute isolated real roots of F or zero curves lying 
%   close to the boundary of the domain. 
%
%   R = ROOTS(F, G) returns the isolated points of F and G.
%
%
%   In the special case when F is of length 1 then the zero contours are found
%   to full precision.
%  
% See also DISKFUN2V/ROOTS, CHEBFUN2/ROOTS, CHEBFUN2V/ROOTS

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty(varargin{ 1 }) )
    r = [];
    return
end

% convert f to a polar chebfun2
f = varargin{1};
f = cart2pol(f);

%check for a second diskfun
if ( nargin > 1 )
    if isa(varargin{2}, 'diskfun')
        g = varargin{2};
        dom = [-pi pi 0 1];
        % we require f and g to be chebfuns. Cart2pol isn't enough because
        % it keeps the rows represented in the Fourier basis.
        f = chebfun2(@(t,r) feval(f, t, r), dom);
        g = chebfun2(@(t,r) feval(g, t, r, 'polar'), dom); 
        r = roots@separableApprox(f,g, varargin{3:end});
        
        %convert to cartesian
        if ~isempty(r)
            [x,y] = pol2cart(r(:,1),r(:,2));
            r = [x,y];
        end
        return
    end
end

rts = roots@separableApprox(f, varargin{2:end});


%convert from polar coords to complex-valued chebfuns: 
x = chebpts(max(length(rts),101));

r = feval(rts,x);
rvals = imag(r).*exp(1i*real(r)); %c = re^(i*x)

r = chebfun(rvals);

% The code below creates a real-valued parametrization of the curves. 
%
% Go through each component and make it an array-valued chebfun: 
% for k = 1:size(vals, 2)
%     comp = feval(rts(:, k), x); 
%     AX = imag(comp) .* cos(real(comp)); 
%     AY = imag(comp) .* sin(real(comp));
%     r{k} = chebfun([AX, AY ]);
% end

end