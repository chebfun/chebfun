classdef spinpreference
%SPINPREFERENCE   Abstract class for managing preferences for SPIN, SPIN2 and 
%SPIN3.
%   SPINPREFERENCE is an abstract class for managing preferences when solving 
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
        dataToPlot           % Which data to plot when complex values (STRING)
        dealias              % To use dealiasing with 2/3-rule (STRING)
        dt                   % Time-step (1x1 DOUBLE)
        dtmin                % Min. time-step for apative time-grid (1x1 DOUBLE)
        dtmax                % Max. time-step for apative time-grid (1x1 DOUBLE)
        errTol               % Desired accuracy on the solution (1x1 DOUBLE)
        iterPlot             % Plot every ITERPLOT iterations (1x1 INT)
        M                    % Number of points for complex means (1x1 INT)
        N                    % Number of points for space-grid (1x1 INT)
        Nmin                 % Min. # of pts for apative space-grid (1x1 INT)
        Nmax                 % Max. # of pts for apative space-grid (1x1 INT)
        Nplot                % Number of points for plotting grid (1x1 INT)
        plot                 % Plot options (STRING)
        scheme               % Time-stepping scheme (STRING)
    end
   
end
