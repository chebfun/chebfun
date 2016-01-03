classdef spinpref < spinpreference
%SPINPREF   Class for managing SPIN/SPINOP preferences.
%
% Available preferences ([] = defaults):
%
%   dealias                   * If 1, use the 2/3-rule to zero high wavenumbers.
%     [0]                       No dealiasing by default.
% 
%   dt                        * Timestep for time discretization. Default is 
%     []                        empty, i.e., adaptively chosen by the code to 
%                               achieve errTol. 
%
%   dtmax                     * Maximum timestep when using an apative grid in
%     [1e-1]                   time.
%
%   dtmin                     * Minimum timestep when using an apative grid in
%     [1e-10]                   time.
%
%   errTol                    * Desired accuracy on the solution.
%     [1e-6]
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
%   plotting                  * Plotting options: 'movie' for plotting a 
%     ['movie']                 movie of the solution, 'waterfall' to use the 
%      'waterfall'              CHEBFUN WATERFALL command. [] for no plotting.
%
%   scheme                    * Timestepping scheme.
%     [@etdrk4]  
%      @exprk5s8
%      @krogstad
%      @eglm433
%      @pecec433
%
%   Ylim                      * Limit of the y-axis when 'plotting' is 'movie'
%     []                        Default is empty, i.e., automatically chosen by 
%                               the code. 
%                                
% See also SPINPREF2, SPINPREF3.

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