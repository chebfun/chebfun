function [F,FA] = gallery2(nn)
%GALLERY2   Gallery of 2-dimensional functions.
%   GALLERY2(N) returns interesting 2D functions as a CHEBFUN2 cell array.
%   N must be a vector with integer entries taking values from 1 to 7.
%   All gallery functions have domain [-1, 1] x [-1, 1].
%
%   [F,FA] = GALLERY2(N) also returns the anonymous functions used to define
%   the CHEBFUN2 objects. If N is an integer, FA is a function handle; if N is
%   a vector, FA is a cell array of function handles.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% These cell arrays collect the functions, domains, and prefs.
funcs = cell(1);
prefs = cell(1);

%%
% Here are the functions for the gallery.

% PEAKS function
funcs{1} = @(x,y) 3*(1-3*x).^2.*exp(-(9*x.^2) - (3*y+1).^2) ...
         - 10*(3*x/5 - 27*x.^3 - (3*y).^5).*exp(-9*(x.^2+y.^2)) ...
         - 1/3*exp(-(3*x+1).^2 - (3*y).^2);
prefs{1} = {};

% Complicated function from the 100-digit challenge
funcs{2} = @(x,y) exp(sin(50*x)) + sin(60*exp(y)) + sin(70*sin(x)) +...
         sin(sin(80*y)) - sin(10*(x+y)) + (x.^2+y.^2)./4;
prefs{2} = {};

funcs{3} = @(x,y) 1./(1+100*(x.^2-y.^2).^2);
prefs{3} = {};

funcs{4} = @(x,y) 1./(1+1e3*((x.^2-.25).^2.*(y.^2-.25).^2));
prefs{4} = {};

funcs{5} = @(x,y) cos(10*(x.^2+y)).*sin(10*(x+y.^2));
prefs{5} = {};

funcs{6} = @(x,y) real(airy(5*(x+y.^2)).*airy(-5*(x.^2+y.^2)));
prefs{6} = {};

funcs{7} = @(x,y) tanh(10*x).*tanh(10*y)./tanh(10).^2+cos(5*x);
prefs{7} = {};

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
F  = {};    % A cell array of the chebfun2 objects.
FA = {};    % A cell array of the function handles.

% Construct only the chebfuns that the user wants.
for n = nn(:)'
    f = chebfun2(funcs{n}, [-1 1 -1 1], prefs{n}{:});
    F = {F{:}, f};
    FA = {FA{:}, funcs{n}};
end

% If the user only wants one function, don't output cell arrays.
if ( length(nn) == 1 )
    F = F{1};
    FA = FA{1};
end

end
