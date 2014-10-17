function [F,FA] = gallery(nn)
%GALLERY   Gallery of 1-dimensional functions.
%   GALLERY(N) returns interesting 1D functions as a CHEBFUN quasimatrix.
%   N must be a vector with integer entries taking values from 1 to 12.
%   All gallery functions have domain [-1, 1].
%
%   [F,FA] = GALLERY(N) also returns the anonymous functions used to define
%   the CHEBFUN. If N is an integer, FA is an anonymous function; if N is a
%   vector, FA is a cell array of anonymous functions.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% The current number of gallery functions:
totalFuns = 12;

% These cell arrays collect the functions, domains, and prefs.
funcs = cell(1, totalFuns);
prefs = cell(1, totalFuns);

% Parse the input:
indx_goodvalues = ismember(nn, 1:totalFuns);
if ( any(indx_goodvalues == 0) )
    % The user passed bad values (e.g. non-integers or too big a number),
    % so remove them.
    nn = nn(indx_goodvalues);

    % If the user passed at least one valid integer, just issue a warning
    % and move on. If the user did not issue any valid integers, give an
    % error.
    if ( length(nn) )
        warning('CHEBFUN:CHEB:gallery:input', 'Ignoring bad input values.')
    else
        error(['CHEBFUN:CHEB:gallery:input', 'Input value(s) are not valid ' ...
        'integers.'])
    end
end


%%
% Here are the functions for the gallery.

% From ATAP, Chapter 2
funcs{1} = @(x) sin(6*x) + sign(sin(x+exp(2*x)));
prefs{1} = {'splitting', 'on'};

% From ATAP, Chapter 3
funcs{2} = @(x) sin(6*x) + sin(60*exp(x));
prefs{2} = {};

% From ATAP, Chapter 3
funcs{3} = @(x) 1./(1+1000*(x+.5).^2) + 1./sqrt(1+1000*(x-.5).^2);
prefs{3} = {};

% From ATAP, Chapter 5
funcs{4} = @(x) tanh(20*sin(12*x)) + .02*exp(3*x).*sin(300*x);
prefs{4} = {};

% From ATAP, Chapter 9 (but made periodic since we can do that now)
funcs{5} = @(x) cos(2*pi*x) - cos(6*pi*x)/3 + cos(10*pi*x)/5 ...
              - cos(14*pi*x)/7 + cos(18*pi*x)/9;
prefs{5} = {'trig'};

% From ATAP, Chapter 13 (Runge function)
funcs{6} = @(x) 1./(1+25*x.^2);
prefs{6} = {};

% From ATAP, Chapter 18
funcs{7} = @(x) exp(x).*sech(4*sin(40*x)).^exp(x);
prefs{7} = {};

% From ATAP, Chapter 22
funcs{8} = @(x) cos(17*x)./(1+sin(100*x).^2);
prefs{8} = {};

% Spike function from Example [quad/SpikeIntegral]
funcs{9} = @(x) sech(5*(x+0.6)).^2 + sech(50*(x+0.2)).^4 + ...
           sech(500*(x-0.2)).^6 + sech(500*(x-0.6)).^8;
prefs{9} = {};

% Needle on a corrugated surface [opt/Needle]
funcs{10} = @(x) 0.1*(x*4).^2 + 0.1*sin(6*(x*4)) + 0.03*sin(12*(x*4));
prefs{10} = {};

% A complicated function from [opt/ExtremeExtrema]
funcs{11} = @(x) cos(3*x+3).*sin(exp(3*x+3));
prefs{11} = {};

% Create the function needed by funcs{12}:
if ( any(nn == 12) )
    g = chebfun(@(t) sign(sin(100*t./(2-t))), 'splitting', 'on');
    f = cumsum(g);
    funcs{12} = @(x) f(x);
    prefs{12} = {10000+1};
end

% NOTE: If you add a new function here, be sure to change the help text and the
% variable totalFuns above to reflect the largest acceptable N!


F = [];     % A quasimatrix of the chebfuns
FA = {};    % A cell array of the anonymous functions

% Only construct the chebfuns that the user wants.
for n = nn(:)'
    f = chebfun(funcs{n}, [-1 1], prefs{n}{:});
    F = [F, f];
    FA = {FA{:}, funcs{n}};
end

% If the user only wants one function, don't output the anon function as a
% cell array.
if ( length(nn) == 1 ),
    FA = FA{1};
end


end