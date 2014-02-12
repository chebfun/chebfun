function ATAPformats(flag)
%ATAPFORMATS   Set default formats for Trefethen's book 
%              "Approximation Theory and Approximation Practice"
%
% ATAPFORMATS sets certain default properties and is called at the beginning of
% each section of the book so that the output from "publish" has a pleasing
% appearance. It is included in this directory as part of the Chebfun
% distribution so that readers of the book need only download Chebfun to be able
% to execute each section.
%
% ATAPFORMATS('factory') will reset these defaults to the MATLAB default state.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 0 )
    % Change to the ATAP settings.

    evalin('caller', 'clear all')
    close all 
    set(0, 'defaultfigureposition', [380 320 540 200]);
    set(0, 'defaultaxeslinewidth',  0.9);
    set(0, 'defaultaxesfontsize',   8);
    set(0, 'defaultlinelinewidth',  1.1);
    set(0, 'defaultpatchlinewidth', 1.1);
    set(0, 'defaultlinemarkersize', 15); 
    format compact
    format long
    chebpref.setDefaults('factory');
    
elseif ( any(strcmpi(flag, {'reset', 'factory'})) )
    % Revert to MATLAB factory values.

    close all
    set(0, 'defaultfigureposition', get(0, 'factoryfigureposition'));
    set(0, 'defaultaxeslinewidth',  get(0, 'factoryaxeslinewidth'));
    set(0, 'defaultaxesfontsize',   get(0, 'factoryaxesfontsize'));
    set(0, 'defaultlinelinewidth',  get(0, 'factorylinelinewidth'));
    set(0, 'defaultpatchlinewidth', get(0, 'factorypatchlinewidth'));
    set(0, 'defaultlinemarkersize', get(0, 'factorylinemarkersize')); 
    chebpref.setDefaults('factory');
    
end

end    
