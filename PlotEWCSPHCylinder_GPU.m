function PlotEWCSPHCylinder_Streamlines_GPU
% Continuous viewer with velocity, pressure and derived vorticity colormaps.
% Vorticity = dv/dx - du/dy; positive = counterclockwise (red).
% Negative = clockwise (blue). Derivatives are calculated on the CPU.
% Boundary derivative accuracy depends on the simulation boundary values.
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
    "EWCSPHCylinder_dr0_0.014142_Nr_31_Ntheta_22_Rr0_20.000000";

% Choose ONE h/dp case written by the C++ code:
% available in the current C++ code: 1.2, 1.5, 1.8, 2.0
hCoefficient = 2;

% Initial displayed variable:
% "rho", "drhodt", "pressure", "u", "v",
% "velocity", "dudt", or "dvdt"
plotVariable = "vorticity";

% Leave empty for percentile-based automatic limits.
manualColorLimits = [];

% Global limits over all fluid particles and all frames.
% Example: [5 95] ignores the lowest/highest 5% when choosing colours.
colorPercentiles = [5 95];

playbackFPS = 20;

% Saved MP4 settings.
videoFPS = 20;
videoQuality = 95;

showStreamlines = false;
streamGridSize = 220;
streamSeedCount = 25;

% Velocity arrows over the continuous field.
showVelocityArrows = false;
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
     "dvdt", ...
     "vorticity"];

allowedVariables = cacheNames(3:end);

parameterLabels = ...
    ["Density", ...
     "Density rate", ...
     "Pressure", ...
     "x velocity", ...
     "y velocity", ...
     "Velocity magnitude", ...
     "x acceleration", ...
     "y acceleration", ...
     "Vorticity"];

colorLabels = ...
    ["Density (kg/m^3)", ...
     "Density rate (kg/m^3/s)", ...
     "Pressure (Pa)", ...
     "u (m/s)", ...
     "v (m/s)", ...
     "Velocity magnitude (m/s)", ...
     "du/dt (m/s^2)", ...
     "dv/dt (m/s^2)", ...
     "Vorticity (s^{-1})"];


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
    vorticity = calculateVorticity(T.x,T.y,T.u,T.v);


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
         T.dvdt, ...
         vorticity];


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


%% Click-selection display
% Fluid particles remain visually hidden. Clicking the continuous surface
% selects the nearest actual SPH fluid particle. Clicking the cylinder
% selects the nearest boundary particle.
hSelectedParticle = plot( ...
    ax,NaN,NaN,'o', ...
    'MarkerSize',10, ...
    'LineWidth',1.6, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor','w', ...
    'LineStyle','none', ...
    'HitTest','off', ...
    'PickableParts','none');

hParticleInfo = text( ...
    ax,0.015,0.985,'', ...
    'Units','normalized', ...
    'VerticalAlignment','top', ...
    'HorizontalAlignment','left', ...
    'BackgroundColor','w', ...
    'EdgeColor',[0.25 0.25 0.25], ...
    'Margin',5, ...
    'Interpreter','none', ...
    'Visible','off', ...
    'HitTest','off', ...
    'PickableParts','none');


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

% Streamlines use instantaneous velocities on the fixed polar grid.
hStreamlines = gobjects(0);
streamCenter = [mean(firstFrame(boundary,1)),mean(firstFrame(boundary,2))];
radius = hypot(firstFrame(:,1)-streamCenter(1),firstFrame(:,2)-streamCenter(2));
wallRadius = mean(radius(boundary));
outerRadius = max(radius);
streamRadii = mean(reshape(radius,nTheta,[]),1);
angles = mod(atan2(firstFrame(1:nTheta,2)-streamCenter(2), ...
                  firstFrame(1:nTheta,1)-streamCenter(1)),2*pi);
[angles,angleOrder] = sort(angles);
angles = [angles(end)-2*pi; angles; angles(1)+2*pi];
[streamX,streamY] = meshgrid( ...
    linspace(streamCenter(1)-outerRadius,streamCenter(1)+outerRadius,streamGridSize), ...
    linspace(streamCenter(2)-outerRadius,streamCenter(2)+outerRadius,streamGridSize));
queryR = hypot(streamX-streamCenter(1),streamY-streamCenter(2));
queryTheta = mod(atan2(streamY-streamCenter(2),streamX-streamCenter(1)),2*pi);
streamMask = queryR <= wallRadius | queryR >= outerRadius;
seedRelativeY = linspace(-0.85*outerRadius,0.85*outerRadius,streamSeedCount);
seedX = streamCenter(1)-sqrt((0.95*outerRadius)^2-seedRelativeY.^2);
seedY = streamCenter(2)+seedRelativeY;
wakeX = linspace(1.5*wallRadius,0.6*outerRadius,5);
seedX = [seedX,streamCenter(1)+wakeX,streamCenter(1)+wakeX];
seedY = [seedY,streamCenter(2)+0.5*wallRadius*ones(1,5), ...
               streamCenter(2)-0.5*wallRadius*ones(1,5)];

axis(ax,'equal');

xlim(ax,xLimits);
ylim(ax,yLimits);

xlabel(ax,'x (m)');
ylabel(ax,'y (m)');

grid(ax,'off');
box(ax,'on');

setFieldColormap;

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


saveVideoButton = uicontrol( ...
    fig, ...
    'Style','pushbutton', ...
    'String','Save MP4', ...
    'Units','normalized', ...
    'Position',[0.395 0.015 0.075 0.043], ...
    'TooltipString','Save all loaded frames as an MP4 animation');


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

saveVideoButton.Callback = @saveAnimation;

parameterMenu.Callback = @parameterChanged;

frameSlider.Callback = @sliderMoved;

fig.WindowKeyPressFcn = @keyPressed;

fig.CloseRequestFcn = @closeViewer;


%% Particle selection by mouse click
% Do not use MATLAB data-cursor indexing on the interpolated surface.
% A surface click is converted to the nearest real SPH particle instead.
hFluid.HitTest = 'on';
hFluid.PickableParts = 'all';
hFluid.ButtonDownFcn = @selectNearestParticle;

hBoundary.HitTest = 'on';
hBoundary.PickableParts = 'all';
hBoundary.ButtonDownFcn = @selectNearestParticle;


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


        if cacheNames(column) == "vorticity"
            magnitude = max(abs([fieldMinLocal fieldMaxLocal]));
            if magnitude == 0 || ~isfinite(magnitude)
                magnitude = 1;
            end
            fieldMinLocal = -magnitude;
            fieldMaxLocal = magnitude;
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

        updateStreamlines(frame);

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

        % Keep the selected physical particle highlighted as frames change.
        if ~isempty(currentState.selectedID)
            selectedRow = find(ID0 == currentState.selectedID,1);
            if ~isempty(selectedRow)
                updateSelectedParticleDisplay(selectedRow,frame);
            end
        end


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
        setFieldColormap;


        [newFieldMin,newFieldMax] = ...
            getColorLimits(fieldColumn);

        clim( ...
            ax, ...
            [newFieldMin newFieldMax]);


        currentState = guidata(fig);

        showFrame(currentState.index);

    end


    %% ---------------------------------------------------------
    % Save animation as MP4
    % ----------------------------------------------------------

    function saveAnimation(~,~)

        stopPlayback;

        currentState = guidata(fig);
        originalFrame = currentState.index;

        % Use the currently selected displayed quantity in the filename.
        defaultName = sprintf( ...
            'EWCSPH_hdp_%g_%s.mp4', ...
            hCoefficient, ...
            char(plotVariable));

        [fileName,pathName] = uiputfile( ...
            '*.mp4', ...
            'Save animation as', ...
            fullfile(simulationFolder,defaultName));

        if isequal(fileName,0)
            return;
        end

        outputFile = fullfile(pathName,fileName);

        writer = VideoWriter(outputFile,'MPEG-4');
        writer.FrameRate = videoFPS;
        writer.Quality = videoQuality;

        % Save current interface state.
        controls = [ ...
            firstButton,backButton,playButton,nextButton,lastButton, ...
            historyButton,saveVideoButton,parameterText,parameterMenu, ...
            frameSlider,frameLabel];

        controlVisibility = get(controls,'Visible');
        oldAxesPosition = ax.Position;
        oldSelectedVisibility = hSelectedParticle.Visible;
        oldInfoVisibility = hParticleInfo.Visible;

        % Temporarily clear the selected particle so its marker/info box
        % does not reappear when showFrame updates each movie frame.
        movieState = currentState;
        movieState.selectedID = [];
        guidata(fig,movieState);

        try
            % Hide viewer controls and particle-selection annotation so the
            % exported movie contains only the scientific visualisation.
            set(controls,'Visible','off');
            hSelectedParticle.Visible = 'off';
            hParticleInfo.Visible = 'off';

            % Use more of the figure area for the exported animation.
            ax.Position = [0.07 0.08 0.80 0.86];

            open(writer);

            for frameIndex = 1:numberOfFrames
                showFrame(frameIndex);
                drawnow;

                movieFrame = getframe(fig);
                writeVideo(writer,movieFrame);
            end

            close(writer);

        catch exportError

            try
                close(writer);
            catch
            end

            % Restore the viewer before reporting the error.
            ax.Position = oldAxesPosition;
            for controlIndex = 1:numel(controls)
                controls(controlIndex).Visible = controlVisibility{controlIndex};
            end
            guidata(fig,currentState);
            hSelectedParticle.Visible = oldSelectedVisibility;
            hParticleInfo.Visible = oldInfoVisibility;
            showFrame(originalFrame);

            rethrow(exportError);
        end

        % Restore interactive viewer.
        ax.Position = oldAxesPosition;
        for controlIndex = 1:numel(controls)
            controls(controlIndex).Visible = controlVisibility{controlIndex};
        end
        guidata(fig,currentState);
        hSelectedParticle.Visible = oldSelectedVisibility;
        hParticleInfo.Visible = oldInfoVisibility;
        showFrame(originalFrame);

        fprintf('\nAnimation saved to:\n%s\n',outputFile);
        msgbox(sprintf('Animation saved to:\n%s',outputFile), ...
            'Animation saved');

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
    % Select nearest actual SPH particle from a mouse click
    % ----------------------------------------------------------

    function selectNearestParticle(source,~)

        % Left click only.
        if ~strcmp(fig.SelectionType,'normal')
            return;
        end

        stopPlayback;

        currentState = guidata(fig);
        frame = currentState.currentFrame;

        % Click location in axes coordinates.
        point = ax.CurrentPoint;
        xClick = point(1,1);
        yClick = point(1,2);

        % Select only from the physical particle set represented by the
        % clicked graphics object.
        if isequal(source,hBoundary)
            candidateRows = boundaryRows;
        else
            candidateRows = fluidRows;
        end

        dx = frame(candidateRows,1) - xClick;
        dy = frame(candidateRows,2) - yClick;

        [~,nearestIndex] = min(dx.^2 + dy.^2);
        row = candidateRows(nearestIndex);

        currentState.selectedID = ID0(row);
        guidata(fig,currentState);

        updateSelectedParticleDisplay(row,frame);

    end


    %% ---------------------------------------------------------
    % Update selected-particle marker and numerical information
    % ----------------------------------------------------------

    function updateSelectedParticleDisplay(row,frame)

        if boundary(row)
            particleType = 'Boundary';
        else
            particleType = 'Fluid';
        end

        % Highlight only the selected particle; the rest of the fluid
        % particles remain hidden so the field stays continuous.
        set( ...
            hSelectedParticle, ...
            'XData',frame(row,1), ...
            'YData',frame(row,2));

        hParticleInfo.String = sprintf( ...
            ['Selected %s particle\n' ...
             'ID = %d\n' ...
             'x = %.6g m\n' ...
             'y = %.6g m\n' ...
             'rho = %.6g kg/m^3\n' ...
             'drho/dt = %.6g kg/m^3/s\n' ...
             'pressure = %.6g Pa\n' ...
             'u = %.6g m/s\n' ...
             'v = %.6g m/s\n' ...
             '|V| = %.6g m/s\n' ...
             'du/dt = %.6g m/s^2\n' ...
             'dv/dt = %.6g m/s^2\n' ...
             'vorticity = %.6g s^-1'], ...
            particleType, ...
            ID0(row), ...
            frame(row,1), ...
            frame(row,2), ...
            frame(row,3), ...
            frame(row,4), ...
            frame(row,5), ...
            frame(row,6), ...
            frame(row,7), ...
            frame(row,8), ...
            frame(row,9), ...
            frame(row,10), ...
            frame(row,11));

        hParticleInfo.Visible = 'on';

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


    % Instantaneous streamlines; these are not time-dependent trajectories.
    function updateStreamlines(frame)
        delete(hStreamlines(isgraphics(hStreamlines)));
        hStreamlines = gobjects(0);
        if ~showStreamlines
            return;
        end
        U = reshape(frame(:,6),nTheta,[]);
        V = reshape(frame(:,7),nTheta,[]);
        U = U(angleOrder,:);
        V = V(angleOrder,:);
        U = [U(end,:);U;U(1,:)];
        V = [V(end,:);V;V(1,:)];
        % Polar interpolation avoids interpolating through the cylinder.
        Fu = griddedInterpolant({angles,streamRadii},U,'linear','none');
        Fv = griddedInterpolant({angles,streamRadii},V,'linear','none');
        gridU = Fu(queryTheta,queryR);
        gridV = Fv(queryTheta,queryR);
        gridU(streamMask) = NaN;
        gridV(streamMask) = NaN;
        paths = stream2(streamX,streamY,gridU,gridV,seedX,seedY,[0.15 2500]);
        for k = 1:numel(paths)
            points = paths{k};
            if size(points,1) < 2
                continue;
            end
            rr = hypot(points(:,1)-streamCenter(1),points(:,2)-streamCenter(2));
            invalid = rr <= wallRadius | rr >= outerRadius | any(~isfinite(points),2);
            points(invalid,:) = NaN;
            h = plot3(ax,points(:,1),points(:,2),1e-3*ones(size(points,1),1), ...
                'Color',[0.15 0.15 0.15],'LineWidth',0.8, ...
                'HitTest','off','PickableParts','none');
            hStreamlines(end+1,1) = h;
        end
    end

    %% Vorticity from derivatives on the nonuniform polar grid
    % Uses the boundary ring too. Radial end derivatives are one-sided.
    % Angular derivatives are periodic; no interpolation across the cylinder.
    function omega = calculateVorticity(xValues,yValues,uValues,vValues)
        nt = numel(boundaryRows);
        nr = numberOfParticles/nt;
        if nt < 3 || nr ~= round(nr) || nr < 3
            error("Vorticity requires a complete polar grid with at least three rings.");
        end
        X = reshape(xValues,nt,nr);
        Y = reshape(yValues,nt,nr);
        U = reshape(uValues,nt,nr);
        V = reshape(vValues,nt,nr);
        radius = hypot(X,Y);
        radii = mean(radius,1);
        angles = unwrap(atan2(Y(:,1),X(:,1)));
        angleStep = 2*pi/nt;
        if any(diff(radii) <= 0) || ...
                max(abs(radius-radii),[],'all') > 1e-5*max(radii) || ...
                max(abs(diff(angles)-angleStep)) > 1e-5
            error("Expected fixed concentric rings with uniformly spaced angles.");
        end
        dUdr = zeros(size(U));
        dVdr = zeros(size(V));
        for angleIndex = 1:nt
            dUdr(angleIndex,:) = gradient(U(angleIndex,:),radii);
            dVdr(angleIndex,:) = gradient(V(angleIndex,:),radii);
        end
        dUdTheta = (circshift(U,-1,1)-circshift(U,1,1))/(2*angleStep);
        dVdTheta = (circshift(V,-1,1)-circshift(V,1,1))/(2*angleStep);
        cosTheta = X./radius;
        sinTheta = Y./radius;
        % d/dx = cos(theta)*d/dr - sin(theta)/r*d/dtheta
        % d/dy = sin(theta)*d/dr + cos(theta)/r*d/dtheta
        omegaMesh = cosTheta.*dVdr - sinTheta./radius.*dVdTheta ...
            - sinTheta.*dUdr - cosTheta./radius.*dUdTheta;
        omega = omegaMesh(:);
    end

    function setFieldColormap
        if plotVariable == "vorticity"
            anchors = [0.10 0.25 0.80; 1 1 1; 0.80 0.10 0.10];
            colours = interp1([0 0.5 1],anchors,linspace(0,1,257));
            colormap(ax,colours);
        else
            colormap(ax,turbo);
        end
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