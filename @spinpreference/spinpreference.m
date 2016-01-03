classdef spinpreference
%SPINPREFERENCE   Abstract class for managing preferences.
%
% See also SPINPREF, SPINPREF2, SPINPREF3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        dataToPlot = 'real';    % Which data to plot when complex values
        dealias = 0             % To use dealiasing with 2/3-rule
        dt                      % Timestep
        dtmin = 1e-10;          % Minimum timestep for apative time-grid
        dtmax                   % Maximum timestep for apative time-grid
        errTol                  % Desired accuracy on the solution
        iterPlot                % Plot every ITERPLOT iterations with 'movie'
        M = 64;                 % Number of points for complex means
        N                       % Number of points for space-grid
        Nmin                    % Min. number of points for apative space-grid
        Nmax                    % Max. number of points for apative space-grid
                                % discretization if using an adaptive space-grid
        plotting  = 'movie';    % Plotting options
        scheme = 'etdrk4';      % Time-stepping scheme
    end
   
end