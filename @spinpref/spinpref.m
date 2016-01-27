classdef spinpref < spinpreference
%SPINPREF   Class for managing preferences when solving a 1D PDE with SPIN.
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
%   dt                        * Time-step for time discretization. Default is 
%     []                        empty, i.e., adaptively chosen by the code to 
%                               achieve errTol. 
%
%   dtmax                     * Maximum time-step when using an apative grid in
%     [1e-1]                    time.
%
%   dtmin                     * Minimum time-step when using an apative grid in
%     [1e-10]                   time.
%
%   errTol                    * Desired accuracy on the solution.
%     [1e-6]
%
%   iterPlot                  * Plot the solution every ITERPLOT iterations of
%     [20]                      the time-stepping loop if 'plot' is 'movie'.
%
%   M                         * Number of points for complex means to evaluate
%     [64]                      the phi-functions.
%
%   N                         * Number points for spatial discretization. 
%     []                        Default is empty, i.e., adaptively chosen by the 
%                               code to achieve errTol.
%
%   Nmin                      * Minimum number of points when using an adaptive
%     [256]                     grid in space.
%
%   Nmax                      * Maximum number of points when using an adaptive
%     [4096]                    grid in space.
%                                         
%   plot                      * Plot options: 'movie' to plot a movie of the 
%     ['movie']                 solution, 'waterfall' to use the WATERFALL
%      'waterfall'              command, 'off' otherwise.
%      'off'
%
%   scheme                    * Time-stepping scheme. HELP/SPINPSCHEME for the
%     ['etdrk4']                list of available schemes.
%
%   Ylim                      * Limit of the y-axis when 'plot' is 'movie'.
%     []                        Default is empty, i.e., automatically chosen by 
%                               the code. 
%              
% Construction:
%   PREF = SPINPREF() creates a SPINPREF object with the default values.
%
%   PREF = SPINPREF(PROP, VALUE) creates a SPINPREF object with the property
%   PROP set to VALUE.
%
% See also SPIN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        Ylim         % Limit of the y-axis of the plot
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function pref = spinpref(varargin) 
            if ( nargin == 0 )
                pref.dtmax = 1e-1;
                pref.errTol = 1e-6;
                pref.iterPlot = 20;
                pref.M = 64;
                pref.Nmin = 256;
                pref.Nmax = 4096;
            else
                pref = spinpref();
                for k = 1:nargin/2
                    pref.(varargin{2*(k-1)+1}) = varargin{2*k};
                end
            end
        end
    end
    
end