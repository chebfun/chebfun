classdef spinpref2 < spinpreference
%SPINPRE2   Class for managing preferences when solving a 2D PDE with SPIN2.
%
% Available preferences ([] = defaults):
%
%   dataToPlot                * What data to plot when the solution is complex-
%     ['real']                  valued.
%      'imag'
%      'abs'
%
%   dealias                   * If it is 'on', use the 2/3-rule to zero high 
%     ['off']                   wavenumbers.
%      'on'
% 
%   dt                        * Time-step for time discretization. To switch to
%     [1]                       adaptive time-stepping, set dt=[].
%
%   dtmax                     * Maximum time-step when using an apative grid in
%     [5]                       time.
%
%   dtmin                     * Minimum time-step when using an apative grid in
%     [1e-10]                   time.
%
%   errTol                    * Desired accuracy on the solution.
%     [1e-4]
%
%   iterPlot                  * Plot the solution every ITERPLOT iterations of
%     [1]                       the time-stepping loop if 'plot' is 'movie'.
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
%   plot                      * Plot options: 'movie' to plot a movie of the
%     ['movie']                 the solution, 'off' otherwise. 
%      'off'
%
%   scheme                    * Time-stepping scheme. HELP/SPINPSCHEME for the
%     ['etdrk4']                list of available schemes.
%
%   view                      * Viewpoint specification when 'plot' is 'movie'.
%     [0 90]   
%
% Construction:
%
%   PREF = SPINPREF2() creates a SPINPREF2 object with the default values.
%
%   PREF = SPINPREF2(PROP, VALUE) creates a SPINPREF2 object with the property
%   PROP set to VALUE.
%
% See also SPIN2.

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
                pref.dtmax = 5;
                pref.errTol = 1e-4;
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