function f = randnfun(dt, dom)
%RANDNFUN   Random smooth function
%   F = RANDNFUN(DT) returns a smooth chebfun on [-1,1] with space
%   scale DT and standard normal distribution N(0,1) at each point.
%   It can be regarded as a sample path of a Gaussian process.
%   F is obtained by calling RANDNFUNTRIG on an interval of about
%   twice the length and restricting the result to [-1,1].
%
%   F = RANDNFUN(DT, DOM) returns a random chebfun on DOM = [A, B].
%
% Example:
%   f = randnfun(0.1);
%   std(f)
%   plot(f)
%   plotcoeffs(f, '.'), xlim([0 200])
%
% See also RANDNFUNTRIG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if nargin == 1
    dom = [-1 1];
end

% Call RANDNFUNTRIG on interval of approximately double length.
% and then restrict the result to the prescribed interval.

m = 2*round(diff(dom)/dt)+1;
dom2 = dom(1) + [0 m*dt];
f = randnfuntrig(dt, dom2);
f = f{dom(1),dom(2)};
