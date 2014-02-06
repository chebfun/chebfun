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
% ATAPFORMATS('reset') will reset these defaults to the MATLAB default state.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

persistent defaultfigureposition
if ( isempty(defaultfigureposition) )
    defaultfigureposition = get(0, 'defaultfigureposition');
    mlock();
end

if ( nargin == 0 )
    % Change to the ATAP settings.

    evalin('caller', 'clear all')
    close all 
    set(0, 'defaultfigureposition', [380 320 540 200])
    set(0, 'defaultaxeslinewidth', 0.9)
    set(0, 'defaultaxesfontsize', 8)
    set(0, 'defaultlinelinewidth', 1.1)
    set(0, 'defaultpatchlinewidth', 1.1)
    set(0, 'defaultlinemarkersize', 15); 
    format compact
    format long
    chebpref.setDefaults('factory');
    
elseif ( any(strcmpi(flag, {'reset', 'factory'})) )
    % Revert to MATLAB defaults

    close all
    % The default options (harcoded). Matlab 2014a Mac OSX.
    set(0, 'defaultfigureposition', defaultfigureposition)
    set(0, 'defaultaxeslinewidth', 0.5)
    set(0, 'defaultaxesfontsize', 10)
    set(0, 'defaultlinelinewidth', 0.5)
    set(0, 'defaultpatchlinewidth', 0.5)
    set(0, 'defaultlinemarkersize', 6); 
    
end

end    
