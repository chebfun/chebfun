function [f,fun] = gallery(N)
%GALLERY   Gallery of 1-dimensional functions.
%   GALLERY(N) returns an interesting 1D function as a CHEBFUN. N must
%   be between 1 and 7.
%
%   [F,FA] = GALLERY(N) also returns an anonymous function used to define
%   the CHEBFUN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% From ATAP, Chapter 2
F{1} = @(x) sin(6*x) + sign(sin(x+exp(2*x)));
domain{1} = [-1 1];

% From ATAP, Chapter 3
F{2} = @(x) sin(6*x) + sin(60*exp(x));
domain{2} = [-1 1];

% From ATAP, Chapter 3
F{3} = @(x) 1./(1+1000*(x+.5).^2) + 1./sqrt(1+1000*(x-.5).^2);
domain{3} = [-1 1];

% From ATAP, Chapter 5
F{4} = @(x) tanh(20*sin(12*x)) + .02*exp(3*x).*sin(300*x);
domain{4} = [-1 1];

% From ATAP, Chapter 9
F{5} = @(x) cos(x) - cos(3*x)/3 + cos(5*x)/5 - cos(7*x)/7 + cos(9*x)/9;
domain{5} = [-6 6];

% From ATAP, Chapter 13 (Runge function)
F{6} = @(x) 1./(1+25*x.^2);
domain{6} = [-1 1];

% From ATAP, Chapter 18
F{7} = @(x) exp(x).*sech(4*sin(40*x)).^exp(x);
domain{7} = [-1 1];

% From ATAP, Chapter 22
F{8} = @(x) cos(17*x)./(1+sin(100*x).^2);
domain{8} = [-1 1];

% Spike function from Example [quad/SpikeIntegral]
F{9} = @(x) sech(10*(x-0.2)).^2 + sech(100*(x-0.4)).^4 + ...
           sech(1000*(x-0.6)).^6 + sech(1000*(x-0.8)).^8;
domain{9} = [0 1];

% Needle on a corrugated surface [opt/Needle]
F{10} = @(x) 0.1*x.^2 + 0.1*sin(6*x) + 0.03*sin(12*x);
domain{10} = [-4 4];

% A complicated function from [opt/ExtremeExtrema]
F{11} = @(x) cos(x).*sin(exp(x));
domain{11} = [0 6];



f = chebfun2(F{N}, domain{N});
fun = F{N};

end