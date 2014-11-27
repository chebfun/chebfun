function [F,FA] = gallerytrig(nn)
%GALLERY   Gallery of periodic 1-dimensional functions.
%   GALLERY(N) returns interesting periodic 1D functions as a CHEBFUN
%   quasimatrix. N must be a vector with integer entries taking values from
%   1 to 4. All gallery functions have domain [-1, 1].
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

% From ATAP, Chapter 9 (made periodic)
funcs{1} = @(x) cos(2*pi*x) - cos(6*pi*x)/3 + cos(10*pi*x)/5 ...
              - cos(14*pi*x)/7 + cos(18*pi*x)/9;
prefs{1} = {}; % 'trig' is included in the construction below

% From Hadrien's thesis (p. 18)
funcs{2} = @(x) cos(8*sin(pi*x+1/7));
prefs{2} = {};

% From Hadrien's thesis (p. 19)
funcs{3} = @(x) sin(6*pi*(x+1)) + sin(20*exp(sin(pi*(x+1))));
prefs{3} = {};

% From Hadrien's thesis (p. 27)
if ( any(nn == 4) )
    L = chebop(@(x,u) diff(u) + (1+sin(cos(10*pi*(x+1)))).*u, [-1 1], 'periodic');
    f = L \ chebfun(@(x) exp(sin(pi*(x+1))));
    funcs{4} = @(x) f(x);
    prefs{4} = {};
end

% % From Hadrien's thesis (p. 18)
% funcs{} = @(x) cos(8*sin(2*pi*x+1/7));
% prefs{} = {};

% % From Hadrien's thesis (p. 18)
% funcs{} = @(x) cos(8*sin(2*pi*x+1/7));
% prefs{} = {};

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
    f = chebfun(funcs{n}, [-1 1], 'trig', prefs{n}{:});
    F = [F, f];
    FA = {FA{:}, funcs{n}};
end

% If the user only wants one function, don't output
% the anonymous function as a cell array.
if ( length(nn) == 1 )
    FA = FA{1};
end

end
