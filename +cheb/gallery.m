function [F,FA] = gallery(nn)
%GALLERY   Gallery of 1-dimensional functions.
%   GALLERY(N) returns interesting 1D functions as a CHEBFUN quasimatrix.
%   N must be a vector with integer entries taking values from 1 to 11.
%   All gallery functions have domain [-1, 1].
%
%   [F,FA] = GALLERY(N) also returns the anonymous functions used to define
%   the CHEBFUN. If N is an integer, FA is a function handle; if N is a
%   vector, FA is a cell array of function handles.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% These cell arrays collect the functions, domains, and prefs.
funcs = cell(1);
prefs = cell(1);

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

% From ATAP, Chapter 13 (Runge function)
funcs{5} = @(x) 1./(1+25*x.^2);
prefs{5} = {};

% From ATAP, Chapter 18
funcs{6} = @(x) exp(x).*sech(4*sin(40*x)).^exp(x);
prefs{6} = {};

% From ATAP, Chapter 22
funcs{7} = @(x) cos(17*x)./(1+sin(100*x).^2);
prefs{7} = {};

% Spike function from Example [quad/SpikeIntegral]
funcs{8} = @(x) sech(5*(x+0.6)).^2 + sech(50*(x+0.2)).^4 + ...
           sech(500*(x-0.2)).^6 + sech(500*(x-0.6)).^8;
prefs{8} = {};

% Needle on a corrugated surface [opt/Needle]
funcs{9} = @(x) 0.1*(x*4).^2 + 0.1*sin(6*(x*4)) + 0.03*sin(12*(x*4));
prefs{9} = {};

% A complicated function from [opt/ExtremeExtrema]
funcs{10} = @(x) cos(3*x+3).*sin(exp(3*x+3));
prefs{10} = {};

% Sawtooth polynomial, from ATAP, Appendix.
if ( any(nn == 11) )
    g = chebfun(@(t) sign(sin(100*t./(2-t))), 10000);
    f = cumsum(g);
    funcs{11} = @(x) f(x);
    prefs{11} = {10000+1};
end

% DEVELOPER NOTE: If you add a new function here, be sure to change the help
% text to reflect the total number of functions available!


% Parse the input:
indx_goodvalues = ismember(nn, 1:length(funcs));
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
        error('CHEBFUN:CHEB:gallery:input', ...
            'Input value(s) are not valid integers.')
    end
end

% Assemble the output:
F  = [];    % A quasimatrix of the chebfuns.
FA = {};    % A cell array of the function handles.

% Construct only the chebfuns that the user wants.
for n = nn(:)'
    f = chebfun(funcs{n}, [-1 1], prefs{n}{:});
    F = [F, f];
    FA = {FA{:}, funcs{n}};
end

% If the user only wants one function, don't output
% the anonymous function as a cell array.
if ( length(nn) == 1 )
    FA = FA{1};
end

end
