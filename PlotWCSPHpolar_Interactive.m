function PlotWCSPHpolar_Interactive
% Interactive viewer for the 17-column WCSPHpolar C++ CSV output.
% Displayed physical quantities are read directly from C++.

%% User settings
simulationFolder = ...
    "WCSPHpolar_dp_0.500000_h_1.000000_Nparticles_859_wendland";

plotVariable = "velocity";
colorLimits = [0 4];
playbackFPS = 20;

%% Find and sort C++ frames
files = dir(fullfile(simulationFolder,"WCSPHpolar_hdp_*_t_*.csv"));

if isempty(files)
    error("No C++ timestep CSV files found in:\n%s",simulationFolder);
end

numberOfFrames = numel(files);
time = nan(numberOfFrames,1);

for k = 1:numberOfFrames
    token = regexp(files(k).name, ...
        '_t_([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\.csv$', ...
        'tokens','once');

    if isempty(token)
        error("Could not read time from %s",files(k).name);
    end

    time(k) = str2double(token{1});
end

[time,order] = sort(time);
files = files(order);
firstFilename = fullfile(files(1).folder,files(1).name);

%% Exact 17-column C++ format
opts = delimitedTextImportOptions("NumVariables",17);
opts.DataLines = [7 Inf];
opts.Delimiter = ",";

opts.VariableNames = ...
    ["ID","r","z","empty1", ...
     "rho","drhodt","pressure","empty2", ...
     "u_r","u_z","velocity","empty3", ...
     "pairAr","pairAz","totalAr","totalAz","Type"];

opts.VariableTypes = ...
    ["double","double","double","string", ...
     "double","double","double","string", ...
     "double","double","double","string", ...
     "double","double","double","double","string"];

T0 = readtable(firstFilename,opts);

if ~ismember(plotVariable,string(T0.Properties.VariableNames))
    error("C++ output does not contain the column %s.",plotVariable);
end

%% Fixed particle identities
ID0 = T0.ID;
type0 = lower(strtrim(string(T0.Type)));
boundary = type0 == "boundary";
fluid = type0 == "fluid";
fluidRows = find(fluid);
numberOfParticles = height(T0);

if any(~(boundary | fluid))
    error("Invalid C++ Type column in the first frame.");
end

if ~isequal(ID0,(0:numberOfParticles-1)')
    error("C++ particle IDs are not sequential.");
end

metadata = readmatrix(firstFilename,'Range','A3:D3');
dp = metadata(3);
boundthick = 3*dp;

%% Display-only wall at the symmetry axis
leftWallR = (-0.5*dp:-dp:-boundthick+0.5*dp)';
leftWallZ = (min(T0.z):dp:max(T0.z))';
[leftRGrid,leftZGrid] = meshgrid(leftWallR,leftWallZ);
leftWallRPlot = leftRGrid(:);
leftWallZPlot = leftZGrid(:);

%% Field label
switch plotVariable
    case "rho",      colorLabel = "Density (kg/m^3)";
    case "drhodt",   colorLabel = "Density rate (kg/m^3/s)";
    case "pressure", colorLabel = "Pressure (Pa)";
    case "u_r",      colorLabel = "Radial velocity (m/s)";
    case "u_z",      colorLabel = "Vertical velocity (m/s)";
    case "velocity", colorLabel = "Velocity magnitude (m/s)";
    case "pairAr",   colorLabel = "Pair radial acceleration (m/s^2)";
    case "pairAz",   colorLabel = "Pair vertical acceleration (m/s^2)";
    case "totalAr",  colorLabel = "Total radial acceleration (m/s^2)";
    case "totalAz",  colorLabel = "Total vertical acceleration (m/s^2)";
    otherwise,        colorLabel = plotVariable;
end

%% Figure
fig = figure('Color','w','Position',[80 80 1250 650], ...
    'Name','Interactive Axisymmetric WCSPH Viewer', ...
    'NumberTitle','off');

ax = axes(fig,'Position',[0.07 0.17 0.80 0.76]);
hold(ax,'on');

field0 = T0.(char(plotVariable));

hFluid = scatter(ax,T0.r(fluid),T0.z(fluid),30,field0(fluid), ...
    'filled','MarkerEdgeColor','none');

hBoundary = scatter(ax, ...
    [T0.r(boundary);leftWallRPlot], ...
    [T0.z(boundary);leftWallZPlot], ...
    150,[0 0 0],'filled');

axis(ax,'equal');
xlim(ax,[min(leftWallRPlot)-dp,max(T0.r)+dp]);
ylim(ax,[min(T0.z)-dp,max(T0.z)+dp]);
xlabel(ax,'r (m)');
ylabel(ax,'z (m)');
grid(ax,'off');
box(ax,'off');
colormap(ax,turbo);
c = colorbar(ax);
c.Label.String = colorLabel;
caxis(ax,colorLimits);
titleHandle = title(ax,'');

%% Navigation buttons
firstButton = uicontrol(fig,'Style','pushbutton','String','|<', ...
    'Units','normalized','Position',[0.07 0.055 0.055 0.055], ...
    'TooltipString','First frame');

backButton = uicontrol(fig,'Style','pushbutton','String','<', ...
    'Units','normalized','Position',[0.13 0.055 0.055 0.055], ...
    'TooltipString','Previous frame');

playButton = uicontrol(fig,'Style','pushbutton','String','Play', ...
    'Units','normalized','Position',[0.19 0.055 0.075 0.055], ...
    'TooltipString','Play or pause');

nextButton = uicontrol(fig,'Style','pushbutton','String','>', ...
    'Units','normalized','Position',[0.27 0.055 0.055 0.055], ...
    'TooltipString','Next frame');

lastButton = uicontrol(fig,'Style','pushbutton','String','>|', ...
    'Units','normalized','Position',[0.33 0.055 0.055 0.055], ...
    'TooltipString','Last frame');

historyButton = uicontrol(fig,'Style','pushbutton','String','History', ...
    'Units','normalized','Position',[0.395 0.055 0.075 0.055], ...
    'TooltipString','History of last clicked particle');

if numberOfFrames > 1
    sliderStep = [1/(numberOfFrames-1),min(1,10/(numberOfFrames-1))];
else
    sliderStep = [1 1];
end

frameSlider = uicontrol(fig,'Style','slider', ...
    'Min',1,'Max',max(1,numberOfFrames),'Value',1, ...
    'SliderStep',sliderStep, ...
    'Units','normalized','Position',[0.49 0.065 0.32 0.035]);

frameLabel = uicontrol(fig,'Style','text','String','', ...
    'BackgroundColor','w','HorizontalAlignment','left', ...
    'Units','normalized','Position',[0.82 0.052 0.16 0.060]);

%% State and timer
playTimer = timer('ExecutionMode','fixedSpacing', ...
    'Period',1/playbackFPS,'BusyMode','drop','TimerFcn',@timerTick);

state.index = 1;
state.currentTable = T0;
state.selectedID = [];
state.fluidRows = fluidRows;
guidata(fig,state);

firstButton.Callback = @(~,~)showFrame(1);
backButton.Callback = @(~,~)stepFrame(-1);
playButton.Callback = @togglePlay;
nextButton.Callback = @(~,~)stepFrame(1);
lastButton.Callback = @(~,~)showFrame(numberOfFrames);
historyButton.Callback = @showSelectedHistory;
frameSlider.Callback = @sliderMoved;
fig.WindowKeyPressFcn = @keyPressed;
fig.CloseRequestFcn = @closeViewer;

dataCursor = datacursormode(fig);
dataCursor.Enable = 'on';
dataCursor.UpdateFcn = @particleDataTip;

showFrame(1);

%% Nested callbacks
    function showFrame(requestedIndex)
        requestedIndex = max(1,min(numberOfFrames,round(requestedIndex)));
        filename = fullfile(files(requestedIndex).folder,files(requestedIndex).name);
        T = readtable(filename,opts);

        if height(T) ~= numberOfParticles || ~isequal(T.ID,ID0)
            error("Particle layout changed in %s",files(requestedIndex).name);
        end

        if ~isequal(lower(strtrim(string(T.Type))),type0)
            error("Particle Type changed in %s",files(requestedIndex).name);
        end

        field = T.(char(plotVariable));
        rFluid = T.r(fluid);
        zFluid = T.z(fluid);
        fieldFluid = field(fluid);
        valid = isfinite(rFluid) & isfinite(zFluid) & isfinite(fieldFluid);
        rFluid(~valid) = NaN;
        zFluid(~valid) = NaN;
        fieldFluid(~valid) = NaN;

        set(hFluid,'XData',rFluid,'YData',zFluid,'CData',fieldFluid);
        set(hBoundary,'XData',[T.r(boundary);leftWallRPlot], ...
            'YData',[T.z(boundary);leftWallZPlot]);

        titleHandle.String = sprintf( ...
            'Axisymmetric SPH Radial Dam-Break (single side), t = %.3f s', ...
            time(requestedIndex));

        frameSlider.Value = requestedIndex;
        frameLabel.String = sprintf('Frame %d / %d\nt = %.3f s', ...
            requestedIndex,numberOfFrames,time(requestedIndex));

        currentState = guidata(fig);
        currentState.index = requestedIndex;
        currentState.currentTable = T;
        guidata(fig,currentState);
        drawnow;
    end

    function stepFrame(direction)
        stopPlayback;
        currentState = guidata(fig);
        showFrame(currentState.index+direction);
    end

    function togglePlay(~,~)
        if strcmp(playTimer.Running,'off')
            playButton.String = 'Pause';
            start(playTimer);
        else
            stopPlayback;
        end
    end

    function timerTick(~,~)
        if ~isvalid(fig), return; end
        currentState = guidata(fig);

        if currentState.index >= numberOfFrames
            stopPlayback;
        else
            showFrame(currentState.index+1);
        end
    end

    function stopPlayback
        if strcmp(playTimer.Running,'on'), stop(playTimer); end
        if isvalid(playButton), playButton.String = 'Play'; end
    end

    function sliderMoved(source,~)
        stopPlayback;
        showFrame(source.Value);
    end

    function keyPressed(~,event)
        switch event.Key
            case 'space',      togglePlay([],[]);
            case 'leftarrow',  stepFrame(-1);
            case 'rightarrow', stepFrame(1);
            case 'home',       stopPlayback; showFrame(1);
            case 'end',        stopPlayback; showFrame(numberOfFrames);
        end
    end

    function output = particleDataTip(~,event)
        currentState = guidata(fig);
        dataIndex = event.DataIndex;

        if dataIndex < 1 || dataIndex > numel(currentState.fluidRows)
            output = {'No fluid particle selected'};
            return;
        end

        row = currentState.fluidRows(dataIndex);
        T = currentState.currentTable;
        currentState.selectedID = T.ID(row);
        guidata(fig,currentState);

        output = { ...
            sprintf('ID: %d',T.ID(row)), ...
            sprintf('r: %.6g m',T.r(row)), ...
            sprintf('z: %.6g m',T.z(row)), ...
            sprintf('velocity: %.6g m/s',T.velocity(row)), ...
            sprintf('u_r: %.6g m/s',T.u_r(row)), ...
            sprintf('u_z: %.6g m/s',T.u_z(row)), ...
            sprintf('rho: %.6g kg/m^3',T.rho(row)), ...
            sprintf('pressure: %.6g Pa',T.pressure(row)), ...
            sprintf('totalAr: %.6g m/s^2',T.totalAr(row)), ...
            sprintf('totalAz: %.6g m/s^2',T.totalAz(row))};
    end

    function showSelectedHistory(~,~)
        stopPlayback;
        currentState = guidata(fig);

        if isempty(currentState.selectedID)
            errordlg('Pause and click a fluid particle first.', ...
                'No particle selected');
            return;
        end

        selectedID = currentState.selectedID;
        history = nan(numberOfFrames,6);

        for frameIndex = 1:numberOfFrames
            filename = fullfile(files(frameIndex).folder,files(frameIndex).name);
            T = readtable(filename,opts);
            row = find(T.ID == selectedID,1);
            history(frameIndex,:) = [T.r(row),T.z(row),T.u_r(row), ...
                T.u_z(row),T.pressure(row),T.totalAz(row)];
        end

        historyFigure = figure('Color','w','Position',[130 60 1100 760]);
        tiledlayout(historyFigure,3,2,'TileSpacing','compact','Padding','compact');
        nexttile; plot(time,history(:,1),'LineWidth',1.3); ylabel('r (m)'); grid on;
        nexttile; plot(time,history(:,2),'LineWidth',1.3); ylabel('z (m)'); grid on;
        nexttile; plot(time,history(:,3),'LineWidth',1.3); ylabel('u_r (m/s)'); grid on;
        nexttile; plot(time,history(:,4),'LineWidth',1.3); ylabel('u_z (m/s)'); grid on;
        nexttile; plot(time,history(:,5),'LineWidth',1.3); xlabel('Time (s)'); ylabel('Pressure (Pa)'); grid on;
        nexttile; plot(time,history(:,6),'LineWidth',1.3); yline(0,'k--'); xlabel('Time (s)'); ylabel('Total a_z (m/s^2)'); grid on;
        sgtitle(sprintf('C++ History for Particle ID %d',selectedID));
    end

    function closeViewer(~,~)
        stopPlayback;
        if isvalid(playTimer), delete(playTimer); end
        delete(fig);
    end
end