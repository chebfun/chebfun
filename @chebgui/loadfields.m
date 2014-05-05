function initSuccess = loadfields(guifile,handles)

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Fill the String fields of the handles
set(handles.input_domain,'String',guifile.domain);
set(handles.input_DE,'String',guifile.DE);
set(handles.input_LBC,'String',guifile.LBC);
set(handles.input_RBC,'String',guifile.RBC);
set(handles.input_BC,'String',guifile.BC);
set(handles.input_GUESS,'String',guifile.init);

if strcmpi(guifile.type,'pde')
    set(handles.input_timedomain,'String',guifile.timedomain);
elseif strcmpi(guifile.type,'eig')
    sigma = guifile.sigma;
    switch sigma
        case '',   set(handles.popupmenu_sigma,'Value',1);
        case 'lm', set(handles.popupmenu_sigma,'Value',2);
        case 'sm', set(handles.popupmenu_sigma,'Value',3);
        case 'lr', set(handles.popupmenu_sigma,'Value',4);
        case 'sr', set(handles.popupmenu_sigma,'Value',5);
        case 'li', set(handles.popupmenu_sigma,'Value',6);
        case 'si', set(handles.popupmenu_sigma,'Value',7);
    end
    numeigs = guifile.options.numeigs;
    if isempty(numeigs), numeigs = 6; end
end


% Store the tolerance in the UserData of the tolerance menu object
set(handles.menu_tolerance,'UserData',guifile.tol);

% Change the checking of menu options
if strcmpi(guifile.type,'pde')
    set(handles.input_timedomain,'String',guifile.timedomain);
    
    if ~strcmp(guifile.options.plotting,'off')
        set(handles.menu_pdeplottingon,'Checked','On');
        set(handles.menu_pdeplottingoff,'Checked','Off');
    else
        set(handles.menu_pdeplottingon,'Checked','Off');
        set(handles.menu_pdeplottingoff,'Checked','On');
    end
    
    if guifile.options.pdeholdplot
        set(handles.menu_pdeholdploton,'Checked','On');
        set(handles.menu_pdeholdplotoff,'Checked','Off');
    else
        set(handles.menu_pdeholdploton,'Checked','Off');
        set(handles.menu_pdeholdplotoff,'Checked','On');
    end  
    
    if ~isempty(guifile.options.fixYaxisLower)
        set(handles.menu_pdefixon,'Checked','On');
        set(handles.menu_pdefixoff,'Checked','Off');
    else
        set(handles.menu_pdefixon,'Checked','Off');
        set(handles.menu_pdefixoff,'Checked','On');
    end

    if ~isempty(guifile.options.fixN)
        set(handles.menu_fixNon,'Checked','On');
        set(handles.menu_fixNoff,'Checked','Off');
    else
        set(handles.menu_fixNon,'Checked','Off');
        set(handles.menu_fixNoff,'Checked','On');
    end
    
elseif strcmpi(guifile.type,'bvp')
    if strcmp(guifile.options.damping,'1')
        set(handles.menu_odedampednewtonon,'Checked','On');
        set(handles.menu_odedampednewtonoff,'Checked','Off');
    else
        set(handles.menu_odedampednewtonon,'Checked','Off');
        set(handles.menu_odedampednewtonoff,'Checked','On');
    end
    
    if ~strcmp(guifile.options.plotting,'off')
        set(handles.menu_odeplottingon,'Checked','On');
        set(handles.menu_odeplottingoff,'Checked','Off');
        set(handles.menu_odeplottingpause,'UserData',guifile.options.plotting);
    else
        set(handles.menu_odeplottingon,'Checked','Off');
        set(handles.menu_odeplottingoff,'Checked','On');
    end
end

% Make sure that we enable the BCs fields again when we load a new demo
set(handles.input_LBC,'Enable','on');
set(handles.input_RBC,'Enable','on');

% Try to plot the initial guess/condition if one exist in the chebgui
% object. If an error is returned, we keep calm and carry on.
if ~isempty(guifile.init)
    try
        initString = guifile.init;
        % Obtain the name of the independend variable from the init field.
        % Need to do concatination if it's a cellstring
        if iscellstr(initString)
            allString = [];
            for initCounter = 1:length(initString)
                % Throw away everything left of = in init 
                equalSign = find(initString{initCounter}=='=');
                allString = [allString,',',initString{initCounter}(equalSign+1:end)];
            end
        else % Else wrap in a cell for later use
            % Throw away everything left of = in init
            equalSign = find(initString=='=');
            if isempty(equalSign), equalSign = 0; end
            allString = initString(equalSign+1:end);
            initString = cellstr(initString);
        end
        % Now obtain the name of the variables
        indVar = symvar(allString);
        
        % Create a domain and a temporary independent variable
        dom = [str2num(guifile.domain)];
        xTemp = chebfun('x',dom);
        % Only support one independent variable for initial
        % guesses/condition.
        if length(indVar) > 1
            return
        elseif length(indVar) == 0 % Only constants passed
            % Do nothing
        else
            % Assign the independent variable its correct name
            eval(sprintf('%s = xTemp;',indVar{1}))
        end
        
        initChebfun = chebfun;
        % Put the initString in a cell
        if ~iscell(initString), initString = cellstr(initString); end
        for initCounter = 1:length(initString)
            equalSign = find(initString{initCounter}=='=');
            if isempty(equalSign), equalSign = 0; end
            currInit = vectorize(initString{initCounter}(equalSign+1:end));
            if isempty(equalSign), equalSign = 0; end
            initChebfun = [initChebfun, chebfun(currInit,dom)];
        end
        axes(handles.fig_sol);
        plot(initChebfun,'LineWidth',2)
        
        if ~isempty(guifile.options.fixYaxisLower) % Fix y-axis
            ylim([str2num(guifile.options.fixYaxisLower),...
                str2num(guifile.options.fixYaxisUpper)]);
        end
        if guifile.options.grid, grid on, end
        
        initSuccess = 1;
    catch ME
        initSuccess = 0;
    end
else
    initSuccess = 0;
end

if strcmpi(guifile.LBC,'periodic')
        set(handles.input_RBC,'String','periodic');
        handles.guifile.RBC = '';
        set(handles.input_RBC,'Enable','off');
elseif strcmpi(guifile.RBC,'periodic')
        set(handles.input_LBC,'String','periodic');
        handles.guifile.LBC = 'periodic';
        handles.guifile.RBC = '';
        set(handles.input_RBC,'Enable','off');
end