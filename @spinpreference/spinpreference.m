classdef spinpreference
%SPINPREFERENCE   Abstract class for managing preferences for SPIN, SPIN2 and 
%SPIN3.
%   SPINPREFERENCE is an abstract class for mananig preferences when solving 
%   a time-dependent PDE defined by a SPINOPERATOR. SPINPREF (for SPIN in 1D), 
%   SPINPREF2 (for SPIN2 in 2D) and SPINPREF3 (for SPIN3 in 3D) are full 
%   implementations.
%
% See also SPINPREF, SPINPREF2, SPINPREF3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        dataToPlot = 'real';    % Which data to plot when complex values
        dealias = 'off';        % To use dealiasing with 2/3-rule
        dt                      % Time-step
        dtmin = 1e-10;          % Minimum time-step for apative time-grid
        dtmax                   % Maximum time-step for apative time-grid
        errTol                  % Desired accuracy on the solution
        iterPlot                % Plot every ITERPLOT iterations with 'movie'
        M                       % Number of points for complex means
        N                       % Number of points for space-grid
        Nmin                    % Min. number of points for apative space-grid
        Nmax                    % Max. number of points for apative space-grid
        plot  = 'movie';        % Plot options
        scheme = 'etdrk4';      % Time-stepping scheme
    end
   
end