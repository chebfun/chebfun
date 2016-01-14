classdef spinpref2 < spinpreference
%SPINPRE2   Class for managing SPIN2/SPINOP2 preferences.
%
% Available preferences ([] = defaults):
%
%   dataToPlot                * What data to plot when the solution is complex-
%     ['real']                  valued.
%      'imag'
%      'abs'
%
%   dealias                   * If 1, use the 2/3-rule to zero high wavenumbers.
%     [0]                       No dealiasing by default.
% 
%   dt                        * Timestep for time discretization. To switch to
%     [1]                       adaptive timestepping, set dt=[].
%
%   dtmax                     * Maximum timestep when using an apative grid in
%     [1]                       time.
%
%   dtmin                     * Minimum timestep when using an apative grid in
%     [1e-10]                   time.
%
%   errTol                    * Desired accuracy on the solution.
%     [1e-3]
%
%   M                         * Number of points for complex means to evaluate
%     [32]                      the phi-functions.
%
%   N                         * Number points in each direction for spatial 
%     [64]                      discretization. To switch to adaptive grid, set
%                               N=[].
%
%   Nmin                      * Minimum number of points in each direction when 
%     [64]                      using an adaptive grid in space.
%
%   Nmax                      * Maximum number of points in each direction when   
%     [512]                     using an adaptive grid in space.
%                                         
%   plotting                  * Plotting options: 'movie' for plotting a 
%     ['movie']                 movie of the solution, [] for no plotting.
%
%   scheme                    * Timestepping scheme. HELP/SPINPSCHEME for more
%     [@etdrk4]                 schemes.
%
%   view                      * Viewpoint specification when using 'movie'.
%     [0 90]   
%
% See also SPINPREF, SPINPREF3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        view = [0 90];        % Viewpoint of the plot
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function pref = spinpref2(varargin) 
            if ( nargin == 0 )
                pref.dt = 1;
                pref.dtmax = 1;
                pref.errTol = 1e-3;
                pref.iterPlot = 1;
                pref.M = 32;
                pref.N = 64;
                pref.Nmin = 64;
                pref.Nmax = 512;
            else
                pref = spinpref2();
                for k = 1:nargin/2
                    pref.(varargin{2*(k-1)+1}) = varargin{2*k};
                end
            end
        end
    end
    
end