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
        dealias = 0             % To use dealiasing with 2/3-rule
        dt                      % Timestep
        dtmin = 1e-10;          % Minimum timestep if using an apative time-grid
        errTol                  % Desired accuracy on the solution
        M = 64;                 % Number of points for complex means
        N                       % Number of points for spatial discretization if
                                % using a fixed space-gird
        Nmax                    % Max. number of points for space discretization 
                                % discretization if using an adaptive space-grid
        plotting  = 'movie';    % Plotting options
        scheme = 'etdrk4';      % Time-stepping scheme
    end
   
end