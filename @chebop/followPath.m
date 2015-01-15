function [uquasi, lamvec, mvec] = followPath(H, A, g, BCstruct, u0, lam0, measure, direction, varargin)
% H  -- chebop
% bcFun -- anonymous function
% u0 -- initial solution on path
% lam0 -- initial value of lambda
% measure -- anonymous function which value we plot
% direction -- Go to the left or the right?
% varargin -- options

% Set default values
plotOn = 0; % Option for plotting
slmax = 4; slmin = 0.0001; % Maximum/min steplength
maxCounter = 7;

if ( nargin < 8 )
    direction = 1;
end

% Parse varargin
while ~isempty(varargin)  % Recurse
    if ~ischar(varargin{1}) && ~isnumeric(varargin{2})
        error('followpath:inputArgument','Incorrect input arguments');
    end
    val = varargin{2};
    switch varargin{1}
        case 'plotOn'
            plotOn = val;
        case 'sl0'
            sl0 = val;
        case 'maxCounter'
            maxCounter = val;
        case 'slmax'
            slmax = val;
    end
    
    % Throw away arguments and move on
    varargin(1:2) = [];
end

sl0 = slmax; % Initial steplength

% Store all the solutions to be returned
uquasi = u0;


% Constraint for tangent
J = @(u,lam) sum(u).^2+lam.^2;

% Iterate along path.
% Begin by finding a tangent direction, then set steplength, then compute
% Newton correction, and repeat.
counter = 1;
uold = u0; lamold = lam0;
if isa(lam0,'chebconst')
    lamvec = lam0.vals;
else
    lamvec = lam0;
end
mvec = measure(u0); sl = sl0;

% Before starting path following, plot initial information
% Obtain a nice string to set on ylabel of bifurcation diagram
mstring = func2str(measure);
mstring = mstring(min(strfind(mstring,')'))+1:end); % Throw away the @(u) part
% mstring = 
if plotOn
    subplot(1,2,1);
    plot(u0), title(['Solution for \lambda =' num2str(lamvec)]), xlabel('x'),ylabel('u(x)')
    subplot(1,2,2)
    plot(lamvec,mvec,'-*')
    title('Bifurcation diagram'), xlabel('\lambda'), ylabel(mstring)
    drawnow, shg
end

% Set up initial tangent and tau we want to be orthogonal to.
told = chebfun(0,domain(u0)); tauold = 1;

retract = 0; % retract == 1 if Newton told us to go back along the tangent.
fprintf('No. path iter    Newton iter   Steplength    Measure    Num. sol.\n')
fprintf('----------------------------------------------------------------\n')
numSols = 1;
while counter <= maxCounter
    % Find a tangent direction, but only if we were told by Newton not to
    % retract
    if ~retract
        [t, tau] = tangentBVP(H,A,g,BCstruct,uold, lamold, told,tauold);
        if counter == 1
            t = direction*t;
            tau = direction*tau;
        end
        % Move in the direction of the tangent
        uinit = uold+sl*t;
        laminit = lamold+sl*tau;
    end

    % Find a Newton correction
    [u, lam, iter, retract] = newtonBVP(H,A,g,BCstruct,uinit,laminit,t,tau);
    
    if retract % Newton told us we were trying to take too long tangent steps
        disp('retracted')
        % Move in the direction of the current tangent, but only with
        % quarter of the steplength
        sl = sl/4;
        uinit = uold+sl*t;
        laminit = lamold+sl*tau;       
        if sl < slmin
            disp('FAILED: sl < slmin')
            return
        end
        continue
    end
    
    if (measure(u) < -150)
        return
    end
    % Store values for plotting
    mvec = [mvec; measure(u)]; lamvec = [lamvec;lam];
    
    fprintf('%7i \t   %2i \t\t %6.4f       %6.4f \t   %i \n',counter,iter,sl, measure(u), numSols)    
    if (mvec(end)*mvec(end-1) < 0 )
        numSols = numSols + 1;
    end


    if plotOn
        subplot(1,2,1);
        plot(u),title(['Solution for \lambda =' num2str(lam)]), xlabel('x'),ylabel('u(x)')
        subplot(1,2,2)
        plot(lamvec,mvec,'-*'), 
        title('Bifurcation diagram'), xlabel('\lambda'), ylabel(mstring)
        drawnow, shg
    end
    
    counter = counter + 1;
    if iter >= 5
        sl = max(sl/2,slmin); % Half steplength
    else
        sl = min(sl*2,slmax); % Try to increase steplength
    end
    
    % If successful, update old values
    told = t; tauold = tau; uold = u; lamold = lam;
    
    % Update quasimatrix to be returned
    uquasi = [uquasi, u];
    
end

