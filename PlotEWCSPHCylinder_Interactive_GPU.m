function PlotEWCSPHCylinder_Flow_GPU
% Continuous polar-surface viewer with velocity arrows.
% Arrow directions use u,v from each CSV frame.
% Arrow lengths are automatically scaled; they do not show trajectories.
% FaceColor=interp smooths colours between fixed grid nodes.
% The surface covers fluid rings only; the cylinder interior is empty.
% Display interpolation does not change or stabilize the SPH solution.
%
% Expected C++ timestep format:
%   EWCSPH_hdp_<hcoef>_t_<time>.csv
%
% CSV particle columns:
%   ID,x,y,,rho,drhodt,pressure,,u,v,,dudt,dvdt,type
%
% The viewer:
%   - selects one h/dp case at a time
%   - caches all frames
%   - optionally uses the GPU for the cache
%   - clips colour limits using percentiles, so a few extreme values
%     do not dominate the colour scale
%   - provides playback, frame navigation, parameter selection,
%     particle data tips, and particle-history plots

%% =========================================================
% USER SETTINGS
% ==========================================================

simulationFolder = ...
    "EWCSPHCylinder_dr0_0.010000_Nr_20_Ntheta_31_Rr0_20.000000";

% Choose ONE h/dp case written by the C++ code:
% available in the current C++ code: 1.2, 1.5, 1.8, 2.0
hCoefficient = 2;

% Initial displayed variable:
% "rho", "drhodt", "pressure", "u", "v",
% "velocity", "dudt", or "dvdt"
plotVariable = "velocity";

% Leave empty for percentile-based automatic limits.
manualColorLimits = [];

% Global limits over all fluid particles and all frames.
% Example: [5 95] ignores the lowest/highest 5% when choosing colours.
colorPercentiles = [5 95];

playbackFPS = 60;

% Velocity arrows over the continuous field.
showVelocityArrows = true;
arrowStride = 3;              % Increase to show fewer arrows
arrowScale = 1.5;             % MATLAB automatic arrow-length scaling
arrowColor = [0 0 0];
if arrowStride < 1 || arrowStride ~= round(arrowStride)
    error("arrowStride must be a positive integer.");
end

% GPU is used only to cache/display the C++ output.
% No SPH quantity is recomputed on the GPU.
useGPU = true;
maximumGPUMemoryFraction = 0.60;


%% =========================================================
% FIND C++ TIMESTEP FILES
% ==========================================================

allFiles = dir(fullfile( ...
    simulationFolder, ...
    "EWCSPH_hdp_*_t_*.csv"));

if isempty(allFiles)
    error("No EWCSPH timestep CSV files found in:\n%s",simulationFolder);
end


%% Parse h/dp and time from every filename

allH = nan(numel(allFiles),1);
allTime = nan(numel(allFiles),1);

numberPattern = ...
    '[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?';

for k = 1:numel(allFiles)

    token = regexp( ...
        allFiles(k).name, ...
        ['_hdp_(' numberPattern ')_t_(' numberPattern ')\.csv$'], ...
        'tokens', ...
        'once');

    if isempty(token)
        error("Could not read h/dp and time from %s",allFiles(k).name);
    end

    allH(k) = str2double(token{1});
    allTime(k) = str2double(token{2});

end


%% Keep only the requested h/dp

hTolerance = 1e-10*max(1,abs(hCoefficient));

selected = abs(allH-hCoefficient) <= hTolerance;

if ~any(selected)

    availableH = unique(allH);

    error( ...
        "No files found for h/dp = %.6g.\nAvailable h/dp values: %s", ...
        hCoefficient, ...
        strjoin(string(availableH.'),", "));

end

files = allFiles(selected);
time = allTime(selected);

[time,order] = sort(time);
files = files(order);

numberOfFrames = numel(files);

fprintf("\nSelected h/dp = %.6g\n",hCoefficient);
fprintf("Found %d timestep files.\n",numberOfFrames);

firstFilename = fullfile(files(1).folder,files(1).name);


%% =========================================================
% EXACT C++ CSV FORMAT
% ==========================================================
%
% Line 1 : h/dp
% Line 2 : summary header
% Line 3 : summary values
% Line 4 : blank
% Line 5 : time
% Line 6 : particle header
% Line 7+: particle data

opts = delimitedTextImportOptions("NumVariables",14);

opts.DataLines = [7 Inf];
opts.Delimiter = ",";

opts.VariableNames = ...
    ["ID", ...
     "x", ...
     "y", ...
     "empty1", ...
     "rho", ...
     "drhodt", ...
     "pressure", ...
     "empty2", ...
     "u", ...
     "v", ...
     "empty3", ...
     "dudt", ...
     "dvdt", ...
     "Type"];

opts.VariableTypes = ...
    ["double", ...
     "double", ...
     "double", ...
     "string", ...
     "double", ...
     "double", ...
     "double", ...
     "string", ...
     "double", ...
     "double", ...
     "string", ...
     "double", ...
     "double", ...
     "string"];


%% =========================================================
% PARAMETERS AVAILABLE IN VIEWER
% ==========================================================

% x and y are stored too, but are not colour-field choices.
% velocity is derived only for visualisation:
% velocity = sqrt(u^2 + v^2)

cacheNames = ...
    ["x", ...
     "y", ...
     "rho", ...
     "drhodt", ...
     "pressure", ...
     "u", ...
     "v", ...
     "velocity", ...
     "dudt", ...
     "dvdt"];

allowedVariables = cacheNames(3:end);

parameterLabels = ...
    ["Density", ...
     "Density rate", ...
     "Pressure", ...
     "x velocity", ...
     "y velocity", ...
     "Velocity magnitude", ...
     "x acceleration", ...
     "y acceleration"];

colorLabels = ...
    ["Density (kg/m^3)", ...
     "Density rate (kg/m^3/s)", ...
     "Pressure (Pa)", ...
     "u (m/s)", ...
     "v (m/s)", ...
     "Velocity magnitude (m/s)", ...
     "du/dt (m/s^2)", ...
     "dv/dt (m/s^2)"];


if ~any(plotVariable == allowedVariables)

    error( ...
        "plotVariable must be one of: %s", ...
        strjoin(allowedVariables,", "));

end

fieldColumn = find(cacheNames == plotVariable,1);

initialVariableIndex = ...
    find(allowedVariables == plotVariable,1);


%% =========================================================
% READ FIRST FRAME
% ==========================================================

T0 = readtable(firstFilename,opts);

ID0 = T0.ID;

type0 = lower(strtrim(string(T0.Type)));

boundary = type0 == "boundary";
fluid = type0 == "fluid";

boundaryRows = find(boundary);
fluidRows = find(fluid);

numberOfParticles = height(T0);


if any(~(boundary | fluid))
    error("Invalid particle Type in the first frame.");
end


if ~isequal(ID0,(0:numberOfParticles-1)')
    error("Particle IDs are not sequential from 0 to N-1.");
end


%% =========================================================
% INITIALISE GPU / CPU CACHE
% ==========================================================

gpuEnabled = false;
gpuInfo = [];

numberOfCachedFields = numel(cacheNames);

cacheBytes = ...
    double(numberOfParticles) * ...
    double(numberOfCachedFields) * ...
    double(numberOfFrames) * 8;


if useGPU

    try

        gpuInfo = gpuDevice;
        gpuEnabled = true;

    catch gpuError

        warning( ...
            "GPU unavailable; using CPU memory instead.\n%s", ...
            gpuError.message);

    end

end


if gpuEnabled && ...
        cacheBytes > maximumGPUMemoryFraction*double(gpuInfo.AvailableMemory)

    warning( ...
        ['The frame cache needs %.2f GiB and exceeds the selected ' ...
         '%.0f%% GPU-memory limit. Using CPU memory instead.'], ...
        cacheBytes/1024^3, ...
        100*maximumGPUMemoryFraction);

    gpuEnabled = false;

end


if gpuEnabled

    frameGPU = gpuArray.zeros( ...
        numberOfParticles, ...
        numberOfCachedFields, ...
        numberOfFrames, ...
        'double');

    frameCPU = [];

    fprintf( ...
        "GPU cache enabled on %s (%.3f GiB).\n", ...
        gpuInfo.Name, ...
        cacheBytes/1024^3);

else

    frameCPU = zeros( ...
        numberOfParticles, ...
        numberOfCachedFields, ...
        numberOfFrames, ...
        'double');

    frameGPU = [];

    fprintf( ...
        "CPU cache enabled (%.3f GiB).\n", ...
        cacheBytes/1024^3);

end


%% =========================================================
% LOAD ALL FRAMES
% ==========================================================

% C++ line 3 contains:
% L2norm Pressure, KE, dr0, reported particle count
summaryValues = nan(numberOfFrames,4);

progressInterval = max(1,floor(numberOfFrames/20));


for k = 1:numberOfFrames

    filename = fullfile(files(k).folder,files(k).name);

    if k == 1
        T = T0;
    else
        T = readtable(filename,opts);
    end


    %% Validate particle identity

    if height(T) ~= numberOfParticles || ...
            ~isequal(T.ID,ID0)

        error( ...
            "Particle layout changed in %s", ...
            files(k).name);

    end


    if ~isequal( ...
            lower(strtrim(string(T.Type))), ...
            type0)

        error( ...
            "Particle Type changed in %s", ...
            files(k).name);

    end


    %% Derived velocity magnitude

    velocity = sqrt(T.u.^2 + T.v.^2);


    %% Cache the numerical data

    frameValues = ...
        [T.x, ...
         T.y, ...
         T.rho, ...
         T.drhodt, ...
         T.pressure, ...
         T.u, ...
         T.v, ...
         velocity, ...
         T.dudt, ...
         T.dvdt];


    if gpuEnabled
        frameGPU(:,:,k) = gpuArray(frameValues);
    else
        frameCPU(:,:,k) = frameValues;
    end


    %% Read C++ summary values from line 3

    summaryLine = readmatrix( ...
        filename, ...
        'Range','A3:D3');

    if numel(summaryLine) >= 4
        summaryValues(k,:) = summaryLine(1,1:4);
    end


    if mod(k,progressInterval) == 0 || ...
            k == numberOfFrames

        fprintf( ...
            "Loaded frame %d of %d.\n", ...
            k, ...
            numberOfFrames);

    end

end


clear T T0 frameValues velocity summaryLine;


%% =========================================================
% FIXED AXIS LIMITS OVER THE COMPLETE ANIMATION
% ==========================================================

if gpuEnabled

    allX = gather(reshape(frameGPU(:,1,:),[],1));
    allY = gather(reshape(frameGPU(:,2,:),[],1));

else

    allX = reshape(frameCPU(:,1,:),[],1);
    allY = reshape(frameCPU(:,2,:),[],1);

end


allX = allX(isfinite(allX));
allY = allY(isfinite(allY));

xMin = min(allX);
xMax = max(allX);

yMin = min(allY);
yMax = max(allY);

span = max(xMax-xMin,yMax-yMin);

if span <= 0 || ~isfinite(span)
    span = 1.0;
end

axisMargin = 0.03*span;

xLimits = [xMin-axisMargin, xMax+axisMargin];
yLimits = [yMin-axisMargin, yMax+axisMargin];

clear allX allY;


%% =========================================================
% INITIAL FRAME AND COLOR LIMIT
% ==========================================================

firstFrame = getFrame(1);

colorLabel = colorLabels(initialVariableIndex);

fieldLimitCache = ...
    nan(numberOfCachedFields,2);

[fieldMin,fieldMax] = ...
    getColorLimits(fieldColumn);


%% =========================================================
% FIGURE
% ==========================================================

fig = figure( ...
    'Color','w', ...
    'Position',[80 80 1250 700], ...
    'Name','Interactive EWCSPH Cylinder Viewer', ...
    'NumberTitle','off');


ax = axes( ...
    fig, ...
    'Position',[0.07 0.18 0.79 0.75]);

hold(ax,'on');


%% Continuous field on the existing polar grid
% C++ orders each ring by angle, with NTheta points per ring.
% Closing the angular seam does not add new simulation values.
nTheta = numel(boundaryRows);
nRings = numel(fluidRows)/nTheta;
if nTheta < 3 || nRings < 2 || nRings ~= round(nRings)
    error("Expected at least two complete fluid rings on the polar grid.");
end
surfaceRows = reshape(fluidRows,nTheta,nRings);
surfaceRows = [surfaceRows; surfaceRows(1,:)];
meshX = firstFrame(surfaceRows(:),1);
meshY = firstFrame(surfaceRows(:),2);
meshC = firstFrame(surfaceRows(:),fieldColumn);
meshSize = size(surfaceRows);
hFluid = surf(ax, ...
    reshape(meshX,meshSize), ...
    reshape(meshY,meshSize), ...
    zeros(meshSize), ...
    reshape(meshC,meshSize), ...
    'FaceColor','interp', 'EdgeColor','none');
view(ax,2);


%% Cylinder boundary

hBoundary = scatter( ...
    ax, ...
    firstFrame(boundary,1), ...
    firstFrame(boundary,2), ...
    42, ...
    [0 0 0], ...
    'filled');


%% Velocity arrows
arrowRows = fluidRows(1:arrowStride:end);
hVelocity = quiver(ax, ...
    firstFrame(arrowRows,1), ...
    firstFrame(arrowRows,2), ...
    firstFrame(arrowRows,6), ...
    firstFrame(arrowRows,7), ...
    arrowScale, 'Color',arrowColor, 'LineWidth',0.8);
if ~showVelocityArrows
    hVelocity.Visible = 'off';
end
% Keep selection attached to the field nodes and cylinder.
hVelocity.HitTest = 'off';
hVelocity.PickableParts = 'none';

axis(ax,'equal');

xlim(ax,xLimits);
ylim(ax,yLimits);

xlabel(ax,'x (m)');
ylabel(ax,'y (m)');

grid(ax,'off');
box(ax,'on');

colormap(ax,turbo);

c = colorbar(ax);
c.Label.String = colorLabel;

clim(ax,[fieldMin fieldMax]);

titleHandle = title(ax,'');


%% =========================================================
% NAVIGATION CONTROLS
% ==========================================================

firstButton = uicontrol( ...
    fig, ...
    'Style','pushbutton', ...
    'String','|<', ...
    'Units','normalized', ...
    'Position',[0.07 0.065 0.055 0.055]);


backButton = uicontrol( ...
    fig, ...
    'Style','pushbutton', ...
    'String','<', ...
    'Units','normalized', ...
    'Position',[0.13 0.065 0.055 0.055]);


playButton = uicontrol( ...
    fig, ...
    'Style','pushbutton', ...
    'String','Play', ...
    'Units','normalized', ...
    'Position',[0.19 0.065 0.075 0.055]);


nextButton = uicontrol( ...
    fig, ...
    'Style','pushbutton', ...
    'String','>', ...
    'Units','normalized', ...
    'Position',[0.27 0.065 0.055 0.055]);


lastButton = uicontrol( ...
    fig, ...
    'Style','pushbutton', ...
    'String','>|', ...
    'Units','normalized', ...
    'Position',[0.33 0.065 0.055 0.055]);


historyButton = uicontrol( ...
    fig, ...
    'Style','pushbutton', ...
    'String','History', ...
    'Units','normalized', ...
    'Position',[0.395 0.065 0.075 0.055]);


parameterText = uicontrol( ...
    fig, ...
    'Style','text', ...
    'String','Displayed parameter:', ...
    'BackgroundColor','w', ...
    'HorizontalAlignment','right', ...
    'Units','normalized', ...
    'Position',[0.49 0.018 0.13 0.035]);


parameterMenu = uicontrol( ...
    fig, ...
    'Style','popupmenu', ...
    'String',cellstr(parameterLabels), ...
    'Value',initialVariableIndex, ...
    'BackgroundColor','w', ...
    'Units','normalized', ...
    'Position',[0.625 0.018 0.235 0.043], ...
    'TooltipString', ...
    'Choose the parameter used to colour fluid particles');


if numberOfFrames > 1

    sliderStep = ...
        [1/(numberOfFrames-1), ...
         min(1,10/(numberOfFrames-1))];

else

    sliderStep = [1 1];

end


frameSlider = uicontrol( ...
    fig, ...
    'Style','slider', ...
    'Min',1, ...
    'Max',max(1,numberOfFrames), ...
    'Value',1, ...
    'SliderStep',sliderStep, ...
    'Units','normalized', ...
    'Position',[0.49 0.075 0.32 0.035]);


frameLabel = uicontrol( ...
    fig, ...
    'Style','text', ...
    'String','', ...
    'BackgroundColor','w', ...
    'HorizontalAlignment','left', ...
    'Units','normalized', ...
    'Position',[0.82 0.058 0.17 0.065]);


%% =========================================================
% STATE AND TIMER
% ==========================================================

playTimer = timer( ...
    'ExecutionMode','fixedSpacing', ...
    'Period',1/playbackFPS, ...
    'BusyMode','drop', ...
    'TimerFcn',@timerTick);


state.index = 1;
state.currentFrame = firstFrame;
state.selectedID = [];

guidata(fig,state);


%% Callbacks

firstButton.Callback = @(~,~)showFrame(1);

backButton.Callback = @(~,~)stepFrame(-1);

playButton.Callback = @togglePlay;

nextButton.Callback = @(~,~)stepFrame(1);

lastButton.Callback = @(~,~)showFrame(numberOfFrames);

historyButton.Callback = @showSelectedHistory;

parameterMenu.Callback = @parameterChanged;

frameSlider.Callback = @sliderMoved;

fig.WindowKeyPressFcn = @keyPressed;

fig.CloseRequestFcn = @closeViewer;


%% Particle data cursor

dataCursor = datacursormode(fig);

dataCursor.Enable = 'on';

dataCursor.UpdateFcn = @particleDataTip;


showFrame(1);


%% =========================================================
% NESTED FUNCTIONS
% ==========================================================

    function [fieldMinLocal,fieldMaxLocal] = ...
            getColorLimits(column)

        %% Manual limits

        if ~isempty(manualColorLimits)

            fieldMinLocal = manualColorLimits(1);
            fieldMaxLocal = manualColorLimits(2);

            return;

        end


        %% Use cached limits if already calculated

        if all(isfinite(fieldLimitCache(column,:)))

            fieldMinLocal = ...
                fieldLimitCache(column,1);

            fieldMaxLocal = ...
                fieldLimitCache(column,2);

            return;

        end


        %% Get ALL fluid values for this parameter

        if gpuEnabled

            values = reshape( ...
                frameGPU(fluid,column,:), ...
                [],1);

            values = ...
                sort(values(isfinite(values)));

        else

            values = reshape( ...
                frameCPU(fluid,column,:), ...
                [],1);

            values = ...
                sort(values(isfinite(values)));

        end


        numberOfValues = numel(values);


        %% Percentile limits

        if numberOfValues == 0

            fieldMinLocal = 0.0;
            fieldMaxLocal = 1.0;

        else

            lowerIndex = max( ...
                1, ...
                ceil( ...
                    colorPercentiles(1)/100 * ...
                    numberOfValues));

            upperIndex = min( ...
                numberOfValues, ...
                ceil( ...
                    colorPercentiles(2)/100 * ...
                    numberOfValues));


            if gpuEnabled

                fieldMinLocal = ...
                    gather(values(lowerIndex));

                fieldMaxLocal = ...
                    gather(values(upperIndex));

            else

                fieldMinLocal = ...
                    values(lowerIndex);

                fieldMaxLocal = ...
                    values(upperIndex);

            end

        end


        %% Velocity magnitude cannot be negative

        if cacheNames(column) == "velocity"

            fieldMinLocal = 0.0;

        end


        %% Prevent identical colour limits

        if fieldMinLocal == fieldMaxLocal

            fieldMaxLocal = ...
                fieldMinLocal + ...
                max(1.0,abs(fieldMinLocal)*0.01);

        end


        fieldLimitCache(column,:) = ...
            [fieldMinLocal fieldMaxLocal];

    end


    %% ---------------------------------------------------------
    % Return one cached frame to CPU
    % ----------------------------------------------------------

    function frame = getFrame(frameIndex)

        if gpuEnabled

            frame = ...
                gather(frameGPU(:,:,frameIndex));

        else

            frame = ...
                frameCPU(:,:,frameIndex);

        end

    end


    %% ---------------------------------------------------------
    % Display requested frame
    % ----------------------------------------------------------

    function showFrame(requestedIndex)

        requestedIndex = ...
            max( ...
                1, ...
                min( ...
                    numberOfFrames, ...
                    round(requestedIndex)));


        frame = ...
            getFrame(requestedIndex);


        %% Continuous field values (same mesh connectivity every frame)
        meshX = frame(surfaceRows(:),1);
        meshY = frame(surfaceRows(:),2);
        meshC = frame(surfaceRows(:),fieldColumn);
        meshC(~isfinite(meshC)) = NaN;
        set(hFluid, ...
            'XData',reshape(meshX,meshSize), ...
            'YData',reshape(meshY,meshSize), ...
            'CData',reshape(meshC,meshSize));


        %% Velocity arrow update
        arrowU = frame(arrowRows,6);
        arrowV = frame(arrowRows,7);
        invalidArrow = ~isfinite(arrowU) | ~isfinite(arrowV);
        arrowU(invalidArrow) = NaN;
        arrowV(invalidArrow) = NaN;
        set(hVelocity, ...
            'XData',frame(arrowRows,1), ...
            'YData',frame(arrowRows,2), ...
            'UData',arrowU, ...
            'VData',arrowV);

        %% Boundary positions

        set( ...
            hBoundary, ...
            'XData',frame(boundary,1), ...
            'YData',frame(boundary,2));


        %% C++ frame summary

        L2pressure = ...
            summaryValues(requestedIndex,1);

        KE = ...
            summaryValues(requestedIndex,2);


        titleHandle.String = sprintf( ...
            ['EWCSPH Cylinder, h/dp = %.2f, t = %.3f s' ...
             '   |   KE = %.4g   |   L_2(P) = %.4g'], ...
            hCoefficient, ...
            time(requestedIndex), ...
            KE, ...
            L2pressure);


        frameSlider.Value = ...
            requestedIndex;


        frameLabel.String = sprintf( ...
            'Frame %d / %d\nt = %.3f s', ...
            requestedIndex, ...
            numberOfFrames, ...
            time(requestedIndex));


        currentState = guidata(fig);

        currentState.index = ...
            requestedIndex;

        currentState.currentFrame = ...
            frame;

        guidata(fig,currentState);


        drawnow limitrate;

    end


    %% ---------------------------------------------------------
    % Step one frame
    % ----------------------------------------------------------

    function stepFrame(direction)

        stopPlayback;

        currentState = guidata(fig);

        showFrame( ...
            currentState.index + direction);

    end


    %% ---------------------------------------------------------
    % Play / pause
    % ----------------------------------------------------------

    function togglePlay(~,~)

        if strcmp(playTimer.Running,'off')

            playButton.String = 'Pause';

            start(playTimer);

        else

            stopPlayback;

        end

    end


    %% ---------------------------------------------------------
    % Timer
    % ----------------------------------------------------------

    function timerTick(~,~)

        if ~isvalid(fig)
            return;
        end

        currentState = guidata(fig);

        if currentState.index >= numberOfFrames

            stopPlayback;

        else

            showFrame( ...
                currentState.index + 1);

        end

    end


    %% ---------------------------------------------------------
    % Stop playback
    % ----------------------------------------------------------

    function stopPlayback

        if strcmp(playTimer.Running,'on')
            stop(playTimer);
        end

        if isvalid(playButton)
            playButton.String = 'Play';
        end

    end


    %% ---------------------------------------------------------
    % Slider
    % ----------------------------------------------------------

    function sliderMoved(source,~)

        stopPlayback;

        showFrame(source.Value);

    end


    %% ---------------------------------------------------------
    % Change displayed variable
    % ----------------------------------------------------------

    function parameterChanged(source,~)

        stopPlayback;

        selectedIndex = ...
            source.Value;

        plotVariable = ...
            allowedVariables(selectedIndex);

        fieldColumn = ...
            find(cacheNames == plotVariable,1);

        c.Label.String = ...
            colorLabels(selectedIndex);


        [newFieldMin,newFieldMax] = ...
            getColorLimits(fieldColumn);

        clim( ...
            ax, ...
            [newFieldMin newFieldMax]);


        currentState = guidata(fig);

        showFrame(currentState.index);

    end


    %% ---------------------------------------------------------
    % Keyboard controls
    % ----------------------------------------------------------

    function keyPressed(~,event)

        switch event.Key

            case 'space'

                togglePlay([],[]);

            case 'leftarrow'

                stepFrame(-1);

            case 'rightarrow'

                stepFrame(1);

            case 'home'

                stopPlayback;
                showFrame(1);

            case 'end'

                stopPlayback;
                showFrame(numberOfFrames);

        end

    end


    %% ---------------------------------------------------------
    % Particle data tip
    % ----------------------------------------------------------

    function output = particleDataTip(~,event)

        currentState = guidata(fig);

        dataIndex = event.DataIndex;

        frame = ...
            currentState.currentFrame;


        if isequal(event.Target,hFluid)

            if dataIndex < 1 || ...
                    dataIndex > numel(surfaceRows)

                output = ...
                    {'No fluid particle selected'};

                return;

            end

            row = ...
                surfaceRows(dataIndex);

            particleType = ...
                'Fluid';


        elseif isequal(event.Target,hBoundary)

            if dataIndex < 1 || ...
                    dataIndex > numel(boundaryRows)

                output = ...
                    {'No boundary particle selected'};

                return;

            end

            row = ...
                boundaryRows(dataIndex);

            particleType = ...
                'Boundary';


        else

            output = ...
                {'Unknown plotted object'};

            return;

        end


        currentState.selectedID = ...
            ID0(row);

        guidata(fig,currentState);


        output = { ...
            sprintf('Type: %s',particleType), ...
            sprintf('ID: %d',ID0(row)), ...
            sprintf('x: %.6g m',frame(row,1)), ...
            sprintf('y: %.6g m',frame(row,2)), ...
            sprintf('rho: %.6g kg/m^3',frame(row,3)), ...
            sprintf('drho/dt: %.6g kg/m^3/s',frame(row,4)), ...
            sprintf('pressure: %.6g Pa',frame(row,5)), ...
            sprintf('u: %.6g m/s',frame(row,6)), ...
            sprintf('v: %.6g m/s',frame(row,7)), ...
            sprintf('|V|: %.6g m/s',frame(row,8)), ...
            sprintf('du/dt: %.6g m/s^2',frame(row,9)), ...
            sprintf('dv/dt: %.6g m/s^2',frame(row,10))};

    end


    %% ---------------------------------------------------------
    % History of selected physical particle
    % ----------------------------------------------------------

    function showSelectedHistory(~,~)

        stopPlayback;

        currentState = guidata(fig);


        if isempty(currentState.selectedID)

            errordlg( ...
                'Pause and click a particle first.', ...
                'No particle selected');

            return;

        end


        selectedID = ...
            currentState.selectedID;

        row = ...
            find(ID0 == selectedID,1);


        % x, y, u, v, pressure, rho
        historyColumns = ...
            [1 2 6 7 5 3];

        history = zeros( ...
            numberOfFrames, ...
            numel(historyColumns));


        for columnIndex = 1:numel(historyColumns)

            if gpuEnabled

                history(:,columnIndex) = ...
                    gather(reshape( ...
                        frameGPU( ...
                            row, ...
                            historyColumns(columnIndex), ...
                            :), ...
                        [],1));

            else

                history(:,columnIndex) = ...
                    reshape( ...
                        frameCPU( ...
                            row, ...
                            historyColumns(columnIndex), ...
                            :), ...
                        [],1);

            end

        end


        historyFigure = figure( ...
            'Color','w', ...
            'Position',[130 60 1100 760]);


        tiledlayout( ...
            historyFigure, ...
            3,2, ...
            'TileSpacing','compact', ...
            'Padding','compact');


        nexttile;

        plot( ...
            time, ...
            history(:,1), ...
            'LineWidth',1.3);

        ylabel('x (m)');
        grid on;


        nexttile;

        plot( ...
            time, ...
            history(:,2), ...
            'LineWidth',1.3);

        ylabel('y (m)');
        grid on;


        nexttile;

        plot( ...
            time, ...
            history(:,3), ...
            'LineWidth',1.3);

        ylabel('u (m/s)');
        grid on;


        nexttile;

        plot( ...
            time, ...
            history(:,4), ...
            'LineWidth',1.3);

        ylabel('v (m/s)');
        grid on;


        nexttile;

        plot( ...
            time, ...
            history(:,5), ...
            'LineWidth',1.3);

        xlabel('Time (s)');
        ylabel('Pressure (Pa)');
        grid on;


        nexttile;

        plot( ...
            time, ...
            history(:,6), ...
            'LineWidth',1.3);

        xlabel('Time (s)');
        ylabel('\rho (kg/m^3)');
        grid on;


        sgtitle(sprintf( ...
            'C++ History for Particle ID %d, h/dp = %.2f', ...
            selectedID, ...
            hCoefficient));

    end


    %% ---------------------------------------------------------
    % Clean close
    % ----------------------------------------------------------

    function closeViewer(~,~)

        stopPlayback;

        if isvalid(playTimer)
            delete(playTimer);
        end

        delete(fig);

    end

end