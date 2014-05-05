function val = get(cg,propName)
% GET   Get chebgui properties.
%
%    'type' - 'bvp','pde','eig'
%    'domain' - spatial domain of BVP/PDE
%    'timedomain' - time domain of PDE
%    'de' - the differential operator or RHS F in u_t = F(x,t,u)
%    'lbc' - left boundary conditions
%    'rbc' - right boundary conditions
%    'bc' - general boundary conditions
%    'tol' - tolerance
%    'init' - intial condition/guess for nonlinear BVPs/PDEs
%    'sigma' - desired eigenvalues: 'LM','SM','LA','SA','LR','SR','LI','SI'
%    'options' - a structure containing the below
%      'numeigs' - number of desired eigenvalues
%      'damping' - damping in newton iteration [true/false]
%      'plotting' - plotting in nonlinear solves/PDEs [true/false]
%      'grid' - display a grid on these plots [true/false]
%      'pdeholdplot' - 
%      'fixn' - fixed spatial discretisation for PDEs (experimental)
%      'fixyaxislower' - fix y axis on plots (lower)
%      'fixyaxisupper' - fix y axis on plots (upper)


% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

switch lower(propName)
    case 'type'
        val = cg.type;
    case 'domain'
        val = cg.domain;
    case 'timedomain'
        val = cg.timedomain;
    case 'de'
        val = cg.DE;
    case 'lbc'
        val = cg.LBC;
    case 'rbc'
        val = cg.RBC;
    case 'bc'
        val = cg.BC;
    case 'tol'
        val = cg.tol;
    case 'init'
        val = cg.init;
    case 'sigma'
        val = cg.sigma;
    case 'options'
        val = cg.options;
    case 'damping'
        val = cg.options.damping;
    case 'plotting'
        val = cg.options.plotting;
    case 'grid'
        val = cg.options.grid;
    case 'pdeholdplot'
        val = cg.options.pdeholdplot;
    case 'fixn'
        val = cg.options.fixN;
    case 'fixyaxislower'
        val = cg.options.fixYaxisLower;
    case 'fixyaxisupper'
        val = cg.options.fixYaxisUpper;
    case 'numeigs'
        val = cg.options.numeigs;
    otherwise
        error('CHEBGUI:get:propname',...
      [propName,' is not a valid chebgui property.'])
end