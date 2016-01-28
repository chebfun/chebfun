classdef spinpref3 < spinpreference
%SPINPREF3   Class for managing preferences when solving a 3D PDE with SPIN3.
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
%     [1e-2]
%
%   iterPlot                  * Plot the solution every ITERPLOT iterations of
%     [1]                       the time-stepping loop if 'plot' is 'movie'.
%
%   M                         * Number of points for complex means to evaluate
%     [32]                      the phi-functions.
%
%   N                         * Number points in each direction for spatial 
%     [32]                      discretization. To switch to adaptive grid, set
%                               N=[].
%
%   Nmin                      * Minimum number of points in each direction when 
%     [32]                      using an adaptive grid in space.
%
%   Nmax                      * Maximum number of points in each direction when   
%     [128]                     using an adaptive grid in space.
%                                         
%   plot                      * Plot options: 'movie' to plot a movie of the
%     ['movie']                 solution, 'off' otherwise.
%      'off'
%
%   scheme                    * Time-stepping scheme. HELP/SPINPSCHEME for the
%     ['etdrk4']                list of available schemes.
%
%   slices                    * Slices of the volumetric slice plot when 'plot'
%     []                        is 'movie'. Default is empty, i.e., 
%                               automatically chosen by the code.                                                                       
%                
% Construction:
%
%   PREF = SPINPREF3() creates a SPINPREF3 object with the default values.
%
%   PREF = SPINPREF3(PROP, VALUE) creates a SPINPRE32 object with the property
%   PROP set to VALUE.
%
% See also SPIN3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        slices = [];        % Slices of the volumetric slice plot
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function pref = spinpref3(varargin) 
            if ( nargin == 0 )
                pref.dt = 1;
                pref.dtmax = 5;
                pref.errTol = 1e-2;
                pref.iterPlot = 1;
                pref.M = 32;
                pref.N = 32;
                pref.Nmin = 32;
                pref.Nmax = 128;
            else
                pref = spinpref3();
                for k = 1:nargin/2
                    pref.(varargin{2*(k-1)+1}) = varargin{2*k};
                end
            end
        end
    end
    
end