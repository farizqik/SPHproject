function PlotEWCSPHCylinder_int_GPU(simulationFolder)

%% =========================================================

% DATASET SELECTION

% ==========================================================

% MATLAB does not duplicate simulation parameters. All physical, numerical,

% geometry, case, and reference-state values are read from the C++ file

% SimulationConstants.csv.  If no folder is supplied, select one interactively.

if nargin < 1 || strlength(string(simulationFolder)) == 0

    selectedFolder = uigetdir(pwd,'Select EWCSPH simulation output folder');

    if isequal(selectedFolder,0)

        return;

    end

    simulationFolder = string(selectedFolder);

else

    simulationFolder = string(simulationFolder);

end

constantsFilename = fullfile(simulationFolder,'SimulationConstants.csv');

if ~isfile(constantsFilename)

    error(['SimulationConstants.csv was not found in:\n%s\n\n' ...
           'Run the metadata-enabled C++ solver first.'],simulationFolder);

end

constantCells = readcell(constantsFilename,'Delimiter',',');

if size(constantCells,2) < 2 || size(constantCells,1) < 2

    error('Invalid SimulationConstants.csv format in %s.',constantsFilename);

end

constantNames = strtrim(string(constantCells(2:end,1)));

constantValues = string(constantCells(2:end,2));

% Simulation/reference values: ALL read from C++ metadata.

referenceDensity = constantNumber('rho0');

referenceSpeed = constantNumber('inletVelocity');

referenceUx = constantNumber('inletVelocityX');

referenceUy = constantNumber('inletVelocityY');

flowAngleDeg = constantNumber('flowAngleDeg');

inletVelocity = referenceSpeed;

cylinderDiameter = constantNumber('Dcylinder');

kinematicViscosity = constantNumber('viscosity');

gammaEOS = constantNumber('gammaEOS');

eosB = constantNumber('B');

reynoldsNumber = constantNumber('Reynolds');

nTheta = round(constantNumber('NTheta'));

metadataNboundary = round(constantNumber('Nboundary'));

metadataNfluid = round(constantNumber('Nfluid'));

metadataNparticles = round(constantNumber('Nparticles'));

metadataHasCylinder = logical(round(constantNumber('hasCylinder')));

kernelName = constantText('kernel');

if nTheta < 1 || ~isfinite(nTheta)

    error('Invalid NTheta in SimulationConstants.csv.');

end

fprintf('\nLoaded C++ metadata from:\n%s\n',constantsFilename);

fprintf('Kernel: %s | NTheta: %d | flow angle: %.6g deg\n', ...
    kernelName,nTheta,flowAngleDeg);

%% =========================================================

% VISUALISATION / POST-PROCESSING SETTINGS ONLY

% ==========================================================

% These settings change only how already-computed C++ results are displayed.

plotVariable = "velocity";

manualColorLimits = [];

colorPercentiles = [5 95];

playbackFPS = 20;

videoFPS = 20;

videoQuality = 95;

showStreamlines = false;

% Denser streamline interpolation/seeding for both SPH and potential flow.
streamGridSize = 320;

streamSeedCount = 100;

showVelocityArrows = false;

arrowStride = 3;

arrowScale = 1.5;

arrowColor = [0 0 0];

useGPU = true;

maximumGPUMemoryFraction = 0.60;

if arrowStride < 1 || arrowStride ~= round(arrowStride)

    error("arrowStride must be a positive integer.");

end

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

%% Select h/dp from C++ output (no hard-coded MATLAB simulation setting)

availableH = unique(allH(isfinite(allH)));

if isempty(availableH)

    error("No valid h/dp values were found in the C++ timestep filenames.");

elseif numel(availableH) == 1

    hCoefficient = availableH(1);

else

    [selectionIndex,confirmed] = listdlg( ...
        'PromptString','Select C++ h/dp result to visualise:', ...
        'SelectionMode','single', ...
        'ListString',cellstr(compose('%.12g',availableH)));

    if ~confirmed

        return;

    end

    hCoefficient = availableH(selectionIndex);

end

hTolerance = 1e-10*max(1,abs(hCoefficient));

selected = abs(allH-hCoefficient) <= hTolerance;

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

opts = delimitedTextImportOptions("NumVariables",16);

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
     "velocity", ...
     "empty3", ...
     "dudt", ...
     "dvdt", ...
     "vorticity", ...
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
     "double", ...
     "string", ...
     "double", ...
     "double", ...
     "double", ...
     "string"];

%% =========================================================

% PARAMETERS AVAILABLE IN VIEWER

% ==========================================================

% x and y are stored too, but are not colour-field choices.

% All numerical fields below are read directly from the C++ output.

% MATLAB only displays/interpolates them for visualisation.

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

buffer = type0 == "buffer";

boundaryRows = find(boundary);

fluidRows = find(fluid);

bufferRows = find(buffer);

numberOfParticles = height(T0);

% Case identity comes from the C++ metadata; cross-check against CSV Type.

hasCylinder = metadataHasCylinder;

if hasCylinder ~= ~isempty(boundaryRows)

    error(['C++ metadata/particle-type mismatch: hasCylinder=%d but the ' ...
           'timestep CSV contains %d boundary particles.'], ...
           hasCylinder,numel(boundaryRows));

end

if metadataNboundary ~= numel(boundaryRows)

    error('Nboundary in SimulationConstants.csv does not match the timestep CSV.');

end

if metadataNfluid ~= numel(fluidRows)

    error('Nfluid in SimulationConstants.csv does not match the timestep CSV.');

end

if metadataNparticles ~= numberOfParticles

    error('Nparticles in SimulationConstants.csv does not match the timestep CSV.');

end

if any(~(boundary | fluid | buffer))

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

% C++ line 2 contains the summary labels.

% C++ line 3 contains the corresponding summary values.

% The updated cylinder solver writes the potential-flow-based L2 values
% directly in this summary row. MATLAB reads them instead of recomputing them.

summaryValues = nan(numberOfFrames,7);

% Updated C++ timestep summary:
%   1 L2norm Velocity
%   2 L2norm ux
%   3 L2norm vy
%   4 L2norm rho
%   5 KE
%   6 dr0
%   7 Nfluid
%
% Do not recompute the L2 errors in MATLAB.  The C++ solver now evaluates
% velocity against the analytical cylinder potential-flow field, so MATLAB
% reads those exact diagnostics from line 3 of each timestep file.

l2DensityHistory = nan(numberOfFrames,1);
l2UHistory = nan(numberOfFrames,1);
l2VHistory = nan(numberOfFrames,1);
l2VelocityHistory = nan(numberOfFrames,1);

progressInterval = max(1,floor(numberOfFrames/20));

for k = 1:numberOfFrames

    filename = fullfile(files(k).folder,files(k).name);

    if k == 1
        T = T0;
    else
        T = readtable(filename,opts);
    end

    %% Validate particle identity
    if height(T) ~= numberOfParticles || ~isequal(T.ID,ID0)
        error("Particle layout changed in %s",files(k).name);
    end

    if ~isequal(lower(strtrim(string(T.Type))),type0)
        error("Particle Type changed in %s",files(k).name);
    end

    %% Cache numerical fields supplied directly by C++
    frameValues = ...
        [T.x, ...
         T.y, ...
         T.rho, ...
         T.drhodt, ...
         T.pressure, ...
         T.u, ...
         T.v, ...
         T.velocity, ...
         T.dudt, ...
         T.dvdt, ...
         T.vorticity];

    if gpuEnabled
        frameGPU(:,:,k) = gpuArray(frameValues);
    else
        frameCPU(:,:,k) = frameValues;
    end

    %% Read the C++ summary values from line 3
    summaryLine = readmatrix(filename,'Range','A3:G3');

    if numel(summaryLine) >= 7
        summaryValues(k,:) = summaryLine(1,1:7);

        l2VelocityHistory(k) = summaryValues(k,1);
        l2UHistory(k) = summaryValues(k,2);
        l2VHistory(k) = summaryValues(k,3);
        l2DensityHistory(k) = summaryValues(k,4);
    end

    if mod(k,progressInterval) == 0 || k == numberOfFrames
        fprintf("Loaded frame %d of %d.\n",k,numberOfFrames);
    end
end

clear T T0 frameValues summaryLine;

%% =========================================================

% SUMMARY LABELS / L2 COLUMNS

% ==========================================================

summaryHeader = readcell(firstFilename,'Range','A2:G2');
summaryLabels = strtrim(string(summaryHeader(1,1:7)));
l2Columns = find(contains(lower(summaryLabels),"l2"));

if isempty(l2Columns)
    warning("No L2 column was found in C++ summary line 2.");
end

%% =========================================================
% C++ POTENTIAL-FLOW COMPARISON FILES
% ==========================================================
%
% Global comparison is small enough to load once.
% Local comparison can contain three rows per angle at every solver timestep,
% so it is read selectively by the comparison callbacks.

globalComparisonFilename = "";
localComparisonFilename = "";
globalComparison = table();
comparisonAvailable = false;

if metadataHasCylinder

    globalFiles = dir(fullfile( ...
        simulationFolder, ...
        "PotentialFlowGlobalComparison_hdp_*.csv"));

    localFiles = dir(fullfile( ...
        simulationFolder, ...
        "PotentialFlowLocalComparison_hdp_*.csv"));

    globalH = nan(numel(globalFiles),1);
    localH = nan(numel(localFiles),1);

    for kk = 1:numel(globalFiles)
        token = regexp( ...
            globalFiles(kk).name, ...
            ['PotentialFlowGlobalComparison_hdp_(' numberPattern ')\.csv$'], ...
            'tokens','once');
        if ~isempty(token)
            globalH(kk) = str2double(token{1});
        end
    end

    for kk = 1:numel(localFiles)
        token = regexp( ...
            localFiles(kk).name, ...
            ['PotentialFlowLocalComparison_hdp_(' numberPattern ')\.csv$'], ...
            'tokens','once');
        if ~isempty(token)
            localH(kk) = str2double(token{1});
        end
    end

    globalMatch = find(isfinite(globalH) & abs(globalH-hCoefficient) <= hTolerance,1);
    localMatch = find(isfinite(localH) & abs(localH-hCoefficient) <= hTolerance,1);

    if ~isempty(globalMatch) && ~isempty(localMatch)

        globalComparisonFilename = fullfile( ...
            globalFiles(globalMatch).folder, ...
            globalFiles(globalMatch).name);

        localComparisonFilename = fullfile( ...
            localFiles(localMatch).folder, ...
            localFiles(localMatch).name);

        globalComparison = readtable( ...
            globalComparisonFilename, ...
            'VariableNamingRule','preserve');

        requiredGlobalColumns = ...
            ["t","timestep", ...
             "FxPressureSPH","FyPressureSPH","FPressureSPH", ...
             "CDPressureSPH","CLPressureSPH","CFPressureSPH", ...
             "FxViscousSPH","FyViscousSPH","FViscousSPH", ...
             "CDViscousSPH","CLViscousSPH","CFViscousSPH", ...
             "FxSPH","FySPH","FSPH", ...
             "CDSPH","CLSPH","CFSPH", ...
             "FxExact","FyExact","FExact", ...
             "CDExact","CLExact","CFExact", ...
             "errorFx","errorFy","errorF", ...
             "errorCD","errorCL","errorCF"];

        missingGlobal = setdiff( ...
            requiredGlobalColumns, ...
            string(globalComparison.Properties.VariableNames));

        if isempty(missingGlobal)
            comparisonAvailable = true;
            fprintf("Loaded updated C++ global comparison: %s\n",globalComparisonFilename);
            fprintf("Using updated C++ local comparison: %s\n",localComparisonFilename);
        else
            warning( ...
                "Global comparison CSV is missing: %s", ...
                strjoin(missingGlobal,", "));
        end
    else
        warning( ...
            "Updated potential-flow comparison files were not found for h/dp = %.6g.", ...
            hCoefficient);
    end
end

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
    'Name','Interactive EWCSPH Viewer', ...
    'NumberTitle','off');

%% Tabs

tabGroup = uitabgroup( ...
    fig, ...
    'Units','normalized', ...
    'Position',[0 0 1 1]);

flowTab = uitab( ...
    tabGroup, ...
    'Title','Flow Field');

l2Tab = uitab( ...
    tabGroup, ...
    'Title','L2 vs Time');

globalComparisonTab = uitab( ...
    tabGroup, ...
    'Title','Global Comparison');

localComparisonTab = uitab( ...
    tabGroup, ...
    'Title','Local Comparison');

% =========================================================
% FLOW FIELD TAB -- SPH VS ANALYTICAL POTENTIAL FLOW
% =========================================================
% Both panels use the same selected parameter, the same spatial limits,
% and the same colour limits so colours have an identical meaning.

ax = axes( ...
    flowTab, ...
    'Position',[0.055 0.20 0.36 0.69]);

axExact = axes( ...
    flowTab, ...
    'Position',[0.545 0.20 0.36 0.69]);

hold(ax,'on');
hold(axExact,'on');

%% Continuous fields on the existing polar grid
% C++ orders each ring by angle, with NTheta points per ring.
% Closing the angular seam does not add new simulation values.

if ~isscalar(nTheta) || ~isfinite(nTheta) || nTheta < 3 || nTheta ~= round(nTheta)
    error("Invalid NTheta in SimulationConstants.csv: %g",nTheta);
end

nRings = numel(fluidRows)/nTheta;

if nRings < 2 || nRings ~= round(nRings)
    error(["The physical-fluid count (%d) is not an integer number of rings " ...
           "for NTheta = %d. Check the C++ particle counts and SimulationConstants.csv."], ...
           numel(fluidRows),nTheta);
end

baseSurfaceRows = reshape(fluidRows,nTheta,nRings);
firstFluidRingRows = baseSurfaceRows(:,1);
surfaceRows = [baseSurfaceRows; baseSurfaceRows(1,:)];
meshSize = size(surfaceRows);

meshX = firstFrame(surfaceRows(:),1);
meshY = firstFrame(surfaceRows(:),2);
meshC = firstFrame(surfaceRows(:),fieldColumn);

meshCExact = zeros(size(meshC));

hFluid = surf(ax, ...
    reshape(meshX,meshSize), ...
    reshape(meshY,meshSize), ...
    zeros(meshSize), ...
    reshape(meshC,meshSize), ...
    'FaceColor','interp','EdgeColor','none');

hExactFluid = surf(axExact, ...
    reshape(meshX,meshSize), ...
    reshape(meshY,meshSize), ...
    zeros(meshSize), ...
    reshape(meshCExact,meshSize), ...
    'FaceColor','interp','EdgeColor','none');

view(ax,2);
view(axExact,2);
hExactFluid.HitTest = 'off';
hExactFluid.PickableParts = 'none';

%% Optional cylinder boundary on both panels
if hasCylinder
    boundaryX = firstFrame(boundary,1);
    boundaryY = firstFrame(boundary,2);
else
    boundaryX = NaN;
    boundaryY = NaN;
end

hBoundary = scatter(ax,boundaryX,boundaryY,36,[0 0 0],'filled');
hExactBoundary = scatter(axExact,boundaryX,boundaryY,36,[0 0 0],'filled');
hExactBoundary.HitTest = 'off';
hExactBoundary.PickableParts = 'none';

%% Click-selection display -- SPH panel only
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
    arrowScale,'Color',arrowColor,'LineWidth',0.8);

exactU0 = zeros(size(arrowRows));
exactV0 = zeros(size(arrowRows));

hExactVelocity = quiver(axExact, ...
    firstFrame(arrowRows,1), ...
    firstFrame(arrowRows,2), ...
    exactU0,exactV0, ...
    arrowScale,'Color',arrowColor,'LineWidth',0.8);

if ~showVelocityArrows
    hVelocity.Visible = 'off';
    hExactVelocity.Visible = 'off';
end

hVelocity.HitTest = 'off';
hVelocity.PickableParts = 'none';
hExactVelocity.HitTest = 'off';
hExactVelocity.PickableParts = 'none';

%% Streamline geometry shared by both panels
hStreamlines = gobjects(0);
hExactStreamlines = gobjects(0);

if hasCylinder
    streamCenter = [mean(firstFrame(boundary,1)),mean(firstFrame(boundary,2))];
else
    streamCenter = [mean(firstFrame(firstFluidRingRows,1)), ...
                    mean(firstFrame(firstFluidRingRows,2))];
end

radius = hypot(firstFrame(:,1)-streamCenter(1),firstFrame(:,2)-streamCenter(2));

if hasCylinder
    innerRadius = mean(radius(boundary));
else
    innerRadius = min(radius(fluid));
end

wallRadius = innerRadius;
outerRadius = max(radius);
streamRadii = mean(reshape(radius,nTheta,[]),1);

angleRingRows = 1:nTheta;
angles = mod(atan2(firstFrame(angleRingRows,2)-streamCenter(2), ...
                  firstFrame(angleRingRows,1)-streamCenter(1)),2*pi);
[angles,angleOrder] = sort(angles);
angles = [angles(end)-2*pi; angles; angles(1)+2*pi];

[streamX,streamY] = meshgrid( ...
    linspace(streamCenter(1)-outerRadius,streamCenter(1)+outerRadius,streamGridSize), ...
    linspace(streamCenter(2)-outerRadius,streamCenter(2)+outerRadius,streamGridSize));

queryR = hypot(streamX-streamCenter(1),streamY-streamCenter(2));
queryTheta = mod(atan2(streamY-streamCenter(2),streamX-streamCenter(1)),2*pi);
streamMask = queryR <= innerRadius | queryR >= outerRadius;

flowUnit = [cosd(flowAngleDeg),sind(flowAngleDeg)];
normalUnit = [-flowUnit(2),flowUnit(1)];

% Dense upstream seeding across the cylinder diameter/domain.
seedCross = linspace(-0.92*outerRadius,0.92*outerRadius,streamSeedCount);
upstreamDistance = 0.96*outerRadius;
seedX = streamCenter(1)-upstreamDistance*flowUnit(1)+seedCross*normalUnit(1);
seedY = streamCenter(2)-upstreamDistance*flowUnit(2)+seedCross*normalUnit(2);

% Add several downstream rows of seeds to resolve the wake region in SPH.
if hasCylinder
    wakeDistance = linspace(1.25*innerRadius,0.70*outerRadius,10);
    wakeOffsets = [-0.75 -0.40 0 0.40 0.75]*innerRadius;
    for io = 1:numel(wakeOffsets)
        seedX = [seedX, streamCenter(1)+wakeDistance*flowUnit(1)+wakeOffsets(io)*normalUnit(1)]; %#ok<AGROW>
        seedY = [seedY, streamCenter(2)+wakeDistance*flowUnit(2)+wakeOffsets(io)*normalUnit(2)]; %#ok<AGROW>
    end
end

%% Axes formatting
for a = [ax axExact]
    axis(a,'equal');
    xlim(a,xLimits);
    ylim(a,yLimits);
    xlabel(a,'x (m)');
    ylabel(a,'y (m)');
    grid(a,'off');
    box(a,'on');
end

setFieldColormap;

c = colorbar(ax);
cExact = colorbar(axExact);
c.Label.String = colorLabel;
cExact.Label.String = colorLabel;

[fieldMin,fieldMax] = getComparisonColorLimits(fieldColumn);
clim(ax,[fieldMin fieldMax]);
clim(axExact,[fieldMin fieldMax]);

% Keep both spatial views synchronized during zoom/pan.
linkaxes([ax axExact],'xy');

titleHandle = title(ax,'SPH');
titleExactHandle = title(axExact,'Analytical potential flow');

%% =========================================================

% L2 ERROR VS TIME TAB -- FOUR SEPARATE PLOT BOXES

% ==========================================================

% Keep the four errors separate so each history has its own y-scale.

% Layout inside the same tab:

%   top-left     : L2(rho)

%   top-right    : L2(u)

%   bottom-left  : L2(v)

%   bottom-right : L2(velocity vector)

l2DensityPlot = l2DensityHistory;

l2UPlot = l2UHistory;

l2VPlot = l2VHistory;

l2VelocityPlot = l2VelocityHistory;

% A logarithmic axis cannot display exact zero. Keep the stored values

% unchanged and hide only zero/non-positive values from the display.

l2DensityPlot(l2DensityPlot <= 0) = NaN;

l2UPlot(l2UPlot <= 0) = NaN;

l2VPlot(l2VPlot <= 0) = NaN;

l2VelocityPlot(l2VelocityPlot <= 0) = NaN;

% ---------------------------------------------------------

% 1) Density L2 -- top left

% ---------------------------------------------------------

axL2Density = axes( ...
    l2Tab, ...
    'Position',[0.08 0.56 0.38 0.36]);

semilogy( ...
    axL2Density, ...
    time, ...
    l2DensityPlot, ...
    '-', ...
    'LineWidth',1.6);

grid(axL2Density,'on');

box(axL2Density,'on');

xlabel(axL2Density,'Time, t (s)');

ylabel(axL2Density,'Normalized L_2 error');

title(axL2Density,'Density: L_2(\rho)','Interpreter','tex');

hCurrentTimeDensity = xline( ...
    axL2Density, ...
    time(1), ...
    '--', ...
    'Current frame', ...
    'HandleVisibility','off');

% ---------------------------------------------------------

% 2) x-velocity L2 -- top right

% ---------------------------------------------------------

axL2U = axes( ...
    l2Tab, ...
    'Position',[0.56 0.56 0.38 0.36]);

semilogy( ...
    axL2U, ...
    time, ...
    l2UPlot, ...
    '-', ...
    'LineWidth',1.6);

grid(axL2U,'on');

box(axL2U,'on');

xlabel(axL2U,'Time, t (s)');

ylabel(axL2U,'Normalized L_2 error');

title(axL2U,'x velocity: L_2(u)','Interpreter','tex');

hCurrentTimeU = xline( ...
    axL2U, ...
    time(1), ...
    '--', ...
    'Current frame', ...
    'HandleVisibility','off');

% ---------------------------------------------------------

% 3) y-velocity L2 -- bottom left

% ---------------------------------------------------------

axL2V = axes( ...
    l2Tab, ...
    'Position',[0.08 0.10 0.38 0.36]);

semilogy( ...
    axL2V, ...
    time, ...
    l2VPlot, ...
    '-', ...
    'LineWidth',1.6);

grid(axL2V,'on');

box(axL2V,'on');

xlabel(axL2V,'Time, t (s)');

ylabel(axL2V,'Normalized L_2 error');

title(axL2V,'y velocity: L_2(v)','Interpreter','tex');

hCurrentTimeV = xline( ...
    axL2V, ...
    time(1), ...
    '--', ...
    'Current frame', ...
    'HandleVisibility','off');

% ---------------------------------------------------------

% 4) velocity-vector L2 -- bottom right

% ---------------------------------------------------------

axL2Velocity = axes( ...
    l2Tab, ...
    'Position',[0.56 0.10 0.38 0.36]);

semilogy( ...
    axL2Velocity, ...
    time, ...
    l2VelocityPlot, ...
    '-', ...
    'LineWidth',1.6);

grid(axL2Velocity,'on');

box(axL2Velocity,'on');

xlabel(axL2Velocity,'Time, t (s)');

ylabel(axL2Velocity,'Normalized L_2 error');

title(axL2Velocity,'Velocity vector: L_2(\bf{u})','Interpreter','tex');

hCurrentTimeVelocity = xline( ...
    axL2Velocity, ...
    time(1), ...
    '--', ...
    'Current frame', ...
    'HandleVisibility','off');


%% =========================================================
% GLOBAL COMPARISON TAB
% ==========================================================

globalGroupNames = ...
    ["F_x: pressure / viscous / total / exact / error", ...
     "F_y: pressure / viscous / total / exact / error", ...
     "|F|: pressure / viscous / total / exact / error", ...
     "C_D: pressure / viscous / total / exact / error", ...
     "C_L: pressure / viscous / total / exact / error", ...
     "C_F: pressure / viscous / total / exact / error"];

uicontrol( ...
    globalComparisonTab, ...
    'Style','text', ...
    'String','Global quantity:', ...
    'BackgroundColor','w', ...
    'HorizontalAlignment','right', ...
    'Units','normalized', ...
    'Position',[0.06 0.93 0.16 0.04]);

globalGroupMenu = uicontrol( ...
    globalComparisonTab, ...
    'Style','popupmenu', ...
    'String',cellstr(globalGroupNames), ...
    'Value',1, ...
    'BackgroundColor','w', ...
    'Units','normalized', ...
    'Position',[0.23 0.93 0.48 0.045]);

globalComparisonInfo = uicontrol( ...
    globalComparisonTab, ...
    'Style','text', ...
    'String','', ...
    'BackgroundColor','w', ...
    'HorizontalAlignment','left', ...
    'Units','normalized', ...
    'Position',[0.73 0.91 0.24 0.07]);

axGlobalValue = axes( ...
    globalComparisonTab, ...
    'Position',[0.07 0.14 0.42 0.70]);

axGlobalError = axes( ...
    globalComparisonTab, ...
    'Position',[0.55 0.14 0.40 0.70]);

hGlobalCurrentTimeValue = gobjects(1);
hGlobalCurrentTimeError = gobjects(1);

if ~comparisonAvailable
    globalComparisonInfo.String = ...
        'Updated C++ comparison CSVs were not found for the selected h/dp.';
end

%% =========================================================
% LOCAL COMPARISON TAB
% ==========================================================

localGroupNames = ...
    ["u velocity", ...
     "v velocity", ...
     "Radial velocity u_r", ...
     "Tangential velocity u_t", ...
     "Velocity magnitude", ...
     "C_u", ...
     "C_v", ...
     "C_{ur}", ...
     "C_{ut}", ...
     "C_V", ...
     "Pressure", ...
     "C_p", ...
     "Wall shear stress", ...
     "Pressure-force x", ...
     "Pressure-force y", ...
     "Pressure-force magnitude", ...
     "Viscous-force x", ...
     "Viscous-force y", ...
     "Viscous-force magnitude", ...
     "Total force x", ...
     "Total force y", ...
     "Total force magnitude", ...
     "Local force coefficient x", ...
     "Local force coefficient y", ...
     "Local resultant force coefficient", ...
     "Normalized velocity-vector error"];

uicontrol( ...
    localComparisonTab, ...
    'Style','text', ...
    'String','Plot against:', ...
    'BackgroundColor','w', ...
    'HorizontalAlignment','right', ...
    'Units','normalized', ...
    'Position',[0.03 0.935 0.10 0.035]);

localModeMenu = uicontrol( ...
    localComparisonTab, ...
    'Style','popupmenu', ...
    'String',{'Angle at current viewer time','Time at selected angle'}, ...
    'Value',1, ...
    'BackgroundColor','w', ...
    'Units','normalized', ...
    'Position',[0.14 0.935 0.22 0.04]);

uicontrol( ...
    localComparisonTab, ...
    'Style','text', ...
    'String','Parameter:', ...
    'BackgroundColor','w', ...
    'HorizontalAlignment','right', ...
    'Units','normalized', ...
    'Position',[0.37 0.935 0.08 0.035]);

localGroupMenu = uicontrol( ...
    localComparisonTab, ...
    'Style','popupmenu', ...
    'String',cellstr(localGroupNames), ...
    'Value',5, ...
    'BackgroundColor','w', ...
    'Units','normalized', ...
    'Position',[0.46 0.935 0.25 0.04]);

uicontrol( ...
    localComparisonTab, ...
    'Style','text', ...
    'String','Layer:', ...
    'BackgroundColor','w', ...
    'HorizontalAlignment','right', ...
    'Units','normalized', ...
    'Position',[0.72 0.935 0.055 0.035]);

localLayerMenu = uicontrol( ...
    localComparisonTab, ...
    'Style','popupmenu', ...
    'String',{'Boundary (0)','Fluid ring 1 (1)','Fluid ring 2 (2)'}, ...
    'Value',1, ...
    'BackgroundColor','w', ...
    'Units','normalized', ...
    'Position',[0.78 0.935 0.18 0.04]);

uicontrol( ...
    localComparisonTab, ...
    'Style','text', ...
    'String','Angle (deg):', ...
    'BackgroundColor','w', ...
    'HorizontalAlignment','right', ...
    'Units','normalized', ...
    'Position',[0.03 0.885 0.10 0.035]);

thetaDegrees = (0:nTheta-1)'*(360/nTheta);
thetaDegreeStrings = cellstr(compose('%.6g',thetaDegrees));

localThetaMenu = uicontrol( ...
    localComparisonTab, ...
    'Style','popupmenu', ...
    'String',thetaDegreeStrings, ...
    'Value',1, ...
    'BackgroundColor','w', ...
    'Units','normalized', ...
    'Position',[0.14 0.885 0.10 0.04]);

localRefreshButton = uicontrol( ...
    localComparisonTab, ...
    'Style','pushbutton', ...
    'String','Refresh', ...
    'Units','normalized', ...
    'Position',[0.26 0.885 0.10 0.04]);

localComparisonInfo = uicontrol( ...
    localComparisonTab, ...
    'Style','text', ...
    'String','', ...
    'BackgroundColor','w', ...
    'HorizontalAlignment','left', ...
    'Units','normalized', ...
    'Position',[0.38 0.87 0.58 0.065]);

axLocalValue = axes( ...
    localComparisonTab, ...
    'Position',[0.07 0.12 0.42 0.71]);

axLocalError = axes( ...
    localComparisonTab, ...
    'Position',[0.55 0.12 0.40 0.71]);

%% =========================================================

% NAVIGATION CONTROLS

% ==========================================================

firstButton = uicontrol( ...
    flowTab, ...
    'Style','pushbutton', ...
    'String','|<', ...
    'Units','normalized', ...
    'Position',[0.07 0.065 0.055 0.055]);

backButton = uicontrol( ...
    flowTab, ...
    'Style','pushbutton', ...
    'String','<', ...
    'Units','normalized', ...
    'Position',[0.13 0.065 0.055 0.055]);

playButton = uicontrol( ...
    flowTab, ...
    'Style','pushbutton', ...
    'String','Play', ...
    'Units','normalized', ...
    'Position',[0.19 0.065 0.075 0.055]);

nextButton = uicontrol( ...
    flowTab, ...
    'Style','pushbutton', ...
    'String','>', ...
    'Units','normalized', ...
    'Position',[0.27 0.065 0.055 0.055]);

lastButton = uicontrol( ...
    flowTab, ...
    'Style','pushbutton', ...
    'String','>|', ...
    'Units','normalized', ...
    'Position',[0.33 0.065 0.055 0.055]);

comparisonButton = uicontrol( ...
    flowTab, ...
    'Style','pushbutton', ...
    'String','Comparison', ...
    'Units','normalized', ...
    'Position',[0.395 0.065 0.085 0.055], ...
    'TooltipString','Open the updated local SPH-versus-potential-flow comparison tab');

historyButton = uicontrol( ...
    flowTab, ...
    'Style','pushbutton', ...
    'String','History', ...
    'Units','normalized', ...
    'Position',[0.07 0.015 0.075 0.043]);

saveVideoButton = uicontrol( ...
    flowTab, ...
    'Style','pushbutton', ...
    'String','Save MP4', ...
    'Units','normalized', ...
    'Position',[0.15 0.015 0.075 0.043], ...
    'TooltipString','Save all loaded frames as an MP4 animation');

parameterText = uicontrol( ...
    flowTab, ...
    'Style','text', ...
    'String','Displayed parameter:', ...
    'BackgroundColor','w', ...
    'HorizontalAlignment','right', ...
    'Units','normalized', ...
    'Position',[0.49 0.018 0.13 0.035]);

parameterMenu = uicontrol( ...
    flowTab, ...
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
    flowTab, ...
    'Style','slider', ...
    'Min',1, ...
    'Max',max(1,numberOfFrames), ...
    'Value',1, ...
    'SliderStep',sliderStep, ...
    'Units','normalized', ...
    'Position',[0.49 0.075 0.32 0.035]);

frameLabel = uicontrol( ...
    flowTab, ...
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

comparisonButton.Callback = @openComparisonTab;

if ~hasCylinder || ~comparisonAvailable

    comparisonButton.Enable = 'off';

    comparisonButton.TooltipString = ...
        'Updated cylinder comparison output is unavailable for this dataset.';

    globalGroupMenu.Enable = 'off';
    localModeMenu.Enable = 'off';
    localGroupMenu.Enable = 'off';
    localLayerMenu.Enable = 'off';
    localThetaMenu.Enable = 'off';
    localRefreshButton.Enable = 'off';

end

saveVideoButton.Callback = @saveAnimation;


globalGroupMenu.Callback = @updateGlobalComparisonPlot;
localModeMenu.Callback = @updateLocalComparisonPlot;
localGroupMenu.Callback = @updateLocalComparisonPlot;
localLayerMenu.Callback = @updateLocalComparisonPlot;
localThetaMenu.Callback = @updateLocalComparisonPlot;
localRefreshButton.Callback = @updateLocalComparisonPlot;

parameterMenu.Callback = @parameterChanged;

frameSlider.Callback = @sliderMoved;

fig.WindowKeyPressFcn = @keyPressed;

fig.CloseRequestFcn = @closeViewer;

tabGroup.SelectionChangedFcn = @comparisonTabChanged;

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

if comparisonAvailable
    updateGlobalComparisonPlot([],[]);
    updateLocalComparisonPlot([],[]);
end

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

        %% SPH and analytical potential-flow fields
        meshX = frame(surfaceRows(:),1);
        meshY = frame(surfaceRows(:),2);

        meshC = frame(surfaceRows(:),fieldColumn);
        meshC(~isfinite(meshC)) = NaN;

        exactField = potentialFlowField(frame,fieldColumn);
        meshCExact = exactField(surfaceRows(:));
        meshCExact(~isfinite(meshCExact)) = NaN;

        set(hFluid, ...
            'XData',reshape(meshX,meshSize), ...
            'YData',reshape(meshY,meshSize), ...
            'CData',reshape(meshC,meshSize));

        set(hExactFluid, ...
            'XData',reshape(meshX,meshSize), ...
            'YData',reshape(meshY,meshSize), ...
            'CData',reshape(meshCExact,meshSize));

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

        [exactArrowU,exactArrowV] = potentialFlowVelocity( ...
            frame(arrowRows,1),frame(arrowRows,2));

        set(hExactVelocity, ...
            'XData',frame(arrowRows,1), ...
            'YData',frame(arrowRows,2), ...
            'UData',exactArrowU, ...
            'VData',exactArrowV);

        updateStreamlines(frame);

        %% Boundary positions

        if hasCylinder
            set(hBoundary, ...
                'XData',frame(boundary,1), ...
                'YData',frame(boundary,2));
            set(hExactBoundary, ...
                'XData',frame(boundary,1), ...
                'YData',frame(boundary,2));
        else
            set(hBoundary,'XData',NaN,'YData',NaN);
            set(hExactBoundary,'XData',NaN,'YData',NaN);
        end

        %% Case-aware title

        currentL2rho = l2DensityHistory(requestedIndex);

        currentL2vel = l2VelocityHistory(requestedIndex);

        if hasCylinder
            titleHandle.String = sprintf( ...
                'SPH: %s, Re_D = %.0f, t = %.3f s | L_2 = %.3g', ...
                char(parameterLabels(fieldColumn-2)),reynoldsNumber, ...
                time(requestedIndex),currentL2vel);
            titleExactHandle.String = sprintf( ...
                'Potential flow: %s | angle = %.1f deg', ...
                char(parameterLabels(fieldColumn-2)),flowAngleDeg);
        else
            titleHandle.String = sprintf( ...
                'SPH: %s, t = %.3f s | L_2 = %.3g', ...
                char(parameterLabels(fieldColumn-2)),time(requestedIndex),currentL2vel);
            titleExactHandle.String = sprintf( ...
                'Uniform analytical reference: %s', ...
                char(parameterLabels(fieldColumn-2)));
        end

        frameSlider.Value = ...
            requestedIndex;

        frameLabel.String = sprintf( ...
            'Frame %d / %d\nt = %.3f s', ...
            requestedIndex, ...
            numberOfFrames, ...
            time(requestedIndex));

        % Keep the L2 time marker synchronized with the displayed frame.

        if isgraphics(hCurrentTimeDensity)

            hCurrentTimeDensity.Value = time(requestedIndex);

        end

        if isgraphics(hCurrentTimeU)

            hCurrentTimeU.Value = time(requestedIndex);

        end

        if isgraphics(hCurrentTimeV)

            hCurrentTimeV.Value = time(requestedIndex);

        end

        if isgraphics(hCurrentTimeVelocity)

            hCurrentTimeVelocity.Value = time(requestedIndex);

        end

        if isgraphics(hGlobalCurrentTimeValue)
            hGlobalCurrentTimeValue.Value = time(requestedIndex);
        end

        if isgraphics(hGlobalCurrentTimeError)
            hGlobalCurrentTimeError.Value = time(requestedIndex);
        end

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
        cExact.Label.String = ...
            colorLabels(selectedIndex);

        setFieldColormap;

        [newFieldMin,newFieldMax] = ...
            getComparisonColorLimits(fieldColumn);

        clim(ax,[newFieldMin newFieldMax]);
        clim(axExact,[newFieldMin newFieldMax]);

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
            comparisonButton,historyButton,saveVideoButton, ...
            parameterText,parameterMenu,frameSlider,frameLabel];

        controlVisibility = get(controls,'Visible');

        oldAxesPosition = ax.Position;
        oldExactAxesPosition = axExact.Position;

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

            ax.Position = [0.045 0.08 0.40 0.84];
            axExact.Position = [0.545 0.08 0.40 0.84];

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
            axExact.Position = oldExactAxesPosition;

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
        axExact.Position = oldExactAxesPosition;

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

        if isequal(source,hBoundary) && hasCylinder

            candidateRows = boundaryRows;

        else

            candidateRows = fluidRows;

        end

        if isempty(candidateRows)

            return;

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

    %% ---------------------------------------------------------
    % Updated C++ SPH-versus-potential-flow comparison
    %
    % Global quantities come from PotentialFlowGlobalComparison_hdp_*.csv.
    % Local angle/time quantities come from
    % PotentialFlowLocalComparison_hdp_*.csv.
    % MATLAB plots the C++ values; it does not recompute the force solution.
    % ----------------------------------------------------------

    %% ---------------------------------------------------------
    % Refresh comparison plots when their tab is opened
    % ----------------------------------------------------------

    function comparisonTabChanged(~,event)

        if ~comparisonAvailable
            return;
        end

        if isequal(event.NewValue,globalComparisonTab)
            updateGlobalComparisonPlot([],[]);
        elseif isequal(event.NewValue,localComparisonTab)
            updateLocalComparisonPlot([],[]);
        end

    end

    %% ---------------------------------------------------------
    % Open the comparison tab
    % ----------------------------------------------------------

    function openComparisonTab(~,~)

        if ~comparisonAvailable
            msgbox( ...
                'Updated C++ potential-flow comparison files are unavailable.', ...
                'Comparison unavailable');
            return;
        end

        stopPlayback;
        tabGroup.SelectedTab = localComparisonTab;
        updateLocalComparisonPlot([],[]);

    end

    %% ---------------------------------------------------------
    % Global comparison versus time
    % ----------------------------------------------------------

    function updateGlobalComparisonPlot(~,~)

        if ~comparisonAvailable
            cla(axGlobalValue);
            cla(axGlobalError);
            text(axGlobalValue,0.5,0.5, ...
                'Updated C++ global comparison data are unavailable.', ...
                'HorizontalAlignment','center');
            axis(axGlobalValue,'off');
            axis(axGlobalError,'off');
            return;
        end

        cla(axGlobalValue);
        cla(axGlobalError);
        axis(axGlobalValue,'on');
        axis(axGlobalError,'on');
        hold(axGlobalValue,'on');
        hold(axGlobalError,'on');

        groupIndex = globalGroupMenu.Value;

        switch groupIndex
            case 1
                valueNames = ["FxPressureSPH","FxViscousSPH","FxSPH","FxExact"];
                valueLabels = ["SPH pressure","SPH viscous","SPH total","Exact"];
                errorName = "errorFx";
                errorLabel = "SPH total - Exact";
                yText = 'F_x (N/m)';
                titleText = 'Global x-force';
            case 2
                valueNames = ["FyPressureSPH","FyViscousSPH","FySPH","FyExact"];
                valueLabels = ["SPH pressure","SPH viscous","SPH total","Exact"];
                errorName = "errorFy";
                errorLabel = "SPH total - Exact";
                yText = 'F_y (N/m)';
                titleText = 'Global y-force';
            case 3
                valueNames = ["FPressureSPH","FViscousSPH","FSPH","FExact"];
                valueLabels = ["SPH pressure","SPH viscous","SPH total","Exact"];
                errorName = "errorF";
                errorLabel = "SPH total - Exact";
                yText = '|F| (N/m)';
                titleText = 'Global resultant force';
            case 4
                valueNames = ["CDPressureSPH","CDViscousSPH","CDSPH","CDExact"];
                valueLabels = ["SPH pressure","SPH viscous","SPH total","Exact"];
                errorName = "errorCD";
                errorLabel = "SPH total - Exact";
                yText = 'C_D';
                titleText = 'Drag coefficient';
            case 5
                valueNames = ["CLPressureSPH","CLViscousSPH","CLSPH","CLExact"];
                valueLabels = ["SPH pressure","SPH viscous","SPH total","Exact"];
                errorName = "errorCL";
                errorLabel = "SPH total - Exact";
                yText = 'C_L';
                titleText = 'Lift coefficient';
            otherwise
                valueNames = ["CFPressureSPH","CFViscousSPH","CFSPH","CFExact"];
                valueLabels = ["SPH pressure","SPH viscous","SPH total","Exact"];
                errorName = "errorCF";
                errorLabel = "SPH total - Exact";
                yText = 'C_F';
                titleText = 'Resultant force coefficient';
        end

        for jj = 1:numel(valueNames)
            if ismember(valueNames(jj),string(globalComparison.Properties.VariableNames))
                plot( ...
                    axGlobalValue, ...
                    globalComparison.t, ...
                    globalComparison.(valueNames(jj)), ...
                    'LineWidth',1.25, ...
                    'DisplayName',valueLabels(jj));
            end
        end

        if ismember(errorName,string(globalComparison.Properties.VariableNames))
            plot( ...
                axGlobalError, ...
                globalComparison.t, ...
                globalComparison.(errorName), ...
                'LineWidth',1.25, ...
                'DisplayName',errorLabel);
        end

        yline(axGlobalError,0,':','HandleVisibility','off');

        currentState = guidata(fig);
        if isempty(currentState)
            currentViewerTime = time(1);
        else
            currentViewerTime = time(currentState.index);
        end

        hGlobalCurrentTimeValue = xline( ...
            axGlobalValue, ...
            currentViewerTime, ...
            ':', ...
            'Current viewer time', ...
            'HandleVisibility','off');

        hGlobalCurrentTimeError = xline( ...
            axGlobalError, ...
            currentViewerTime, ...
            ':', ...
            'Current viewer time', ...
            'HandleVisibility','off');

        xlabel(axGlobalValue,'Time, t (s)');
        ylabel(axGlobalValue,yText);
        title(axGlobalValue,sprintf( ...
            '%s: SPH vs exact',titleText));
        grid(axGlobalValue,'on');
        box(axGlobalValue,'on');
        legend(axGlobalValue,'Location','best');

        xlabel(axGlobalError,'Time, t (s)');
        ylabel(axGlobalError,['Error in ' yText]);
        title(axGlobalError,sprintf('%s: error',titleText));
        grid(axGlobalError,'on');
        box(axGlobalError,'on');
        legend(axGlobalError,'Location','best');

        sgtitle(globalComparisonTab,sprintf( ...
            'Global comparison | Re_D = %.4g | h/dp = %.4g', ...
            reynoldsNumber,hCoefficient));

        hold(axGlobalValue,'off');
        hold(axGlobalError,'off');

        globalComparisonInfo.String = sprintf( ...
            '%d solver timesteps | t = %.4g to %.4g s', ...
            height(globalComparison), ...
            globalComparison.t(1), ...
            globalComparison.t(end));

    end

    %% ---------------------------------------------------------
    % Local comparison: versus angle or versus time
    % ----------------------------------------------------------

    function updateLocalComparisonPlot(~,~)

        if ~comparisonAvailable
            cla(axLocalValue);
            cla(axLocalError);
            text(axLocalValue,0.5,0.5, ...
                'Updated C++ local comparison data are unavailable.', ...
                'HorizontalAlignment','center');
            axis(axLocalValue,'off');
            axis(axLocalError,'off');
            return;
        end

        groupIndex = localGroupMenu.Value;

        [names,labels,yText,titleText,boundaryOnly] = ...
            localComparisonSpec(groupIndex);

        if boundaryOnly
            localLayerMenu.Value = 1;
            localLayerMenu.Enable = 'off';
            layer = 0;
        else
            localLayerMenu.Enable = 'on';
            layer = localLayerMenu.Value - 1;
        end

        modeIndex = localModeMenu.Value;

        cla(axLocalValue);
        cla(axLocalError);
        axis(axLocalValue,'on');
        axis(axLocalError,'on');
        hold(axLocalValue,'on');
        hold(axLocalError,'on');

        if modeIndex == 1
            % Angle at the nearest C++ solver time to the current viewer frame.
            localThetaMenu.Enable = 'off';

            currentState = guidata(fig);
            if isempty(currentState)
                viewerTime = time(1);
            else
                viewerTime = time(currentState.index);
            end

            [~,globalIndex] = min(abs(globalComparison.t-viewerTime));
            comparisonTime = globalComparison.t(globalIndex);

            localTable = readLocalAngleSlice(globalIndex,layer,names);

            if isempty(localTable)
                localComparisonInfo.String = 'No matching local rows were found.';
                hold(axLocalValue,'off');
                hold(axLocalError,'off');
                return;
            end

            [xData,orderLocal] = sort(localTable.thetaDeg);
            localTable = localTable(orderLocal,:);
            xLabelText = '\theta (deg)';

            localComparisonInfo.String = sprintf( ...
                'Viewer t = %.6g s | C++ comparison t = %.6g s | layer = %d', ...
                viewerTime,comparisonTime,layer);

            titleSuffix = sprintf('at t = %.6g s, layer %d',comparisonTime,layer);

        else
            % Time history at one selected angle in degrees.
            localThetaMenu.Enable = 'on';

            thetaIndex = localThetaMenu.Value - 1;
            thetaDegRequested = thetaDegrees(localThetaMenu.Value);

            localTable = readLocalTimeHistory(layer,thetaIndex,names);

            if isempty(localTable)
                localComparisonInfo.String = 'No matching local time-history rows were found.';
                hold(axLocalValue,'off');
                hold(axLocalError,'off');
                return;
            end

            xData = localTable.t;
            xLabelText = 'Time, t (s)';

            thetaDegSelected = localTable.thetaDeg(1);

            localComparisonInfo.String = sprintf( ...
                '\theta = %.6g deg | layer = %d | %d timesteps', ...
                thetaDegSelected,layer,height(localTable));

            titleSuffix = sprintf( ...
                'at \theta = %.6g deg, layer %d', ...
                thetaDegSelected,layer);
        end

        % ------------------------------------------------------
        % Split the selected parameter into value and error plots.
        % ------------------------------------------------------
        valueCount = 2;

        if groupIndex == 13 || (groupIndex >= 17 && groupIndex <= 19) || groupIndex == 26
            valueCount = 1;
        end

        if groupIndex >= 14 && groupIndex <= 16
            valueCount = 2;
        end

        % Value plot
        if groupIndex == 13
            % Free-slip potential flow: exact wall shear = 0.
            plot(axLocalValue,xData,localTable.tauWall, ...
                'LineWidth',1.25,'DisplayName','SPH');
            plot(axLocalValue,xData,zeros(size(xData)), ...
                '--','LineWidth',1.1,'DisplayName','Exact = 0');

            plot(axLocalError,xData,localTable.tauWall, ...
                'LineWidth',1.25,'DisplayName','SPH - Exact');

        elseif groupIndex >= 17 && groupIndex <= 19
            % Potential flow has zero viscous wall force.
            sphName = names(1);
            plot(axLocalValue,xData,localTable.(sphName), ...
                'LineWidth',1.25,'DisplayName',labels(1));
            plot(axLocalValue,xData,zeros(size(xData)), ...
                '--','LineWidth',1.1,'DisplayName','Exact = 0');

            plot(axLocalError,xData,localTable.(sphName), ...
                'LineWidth',1.25,'DisplayName','SPH - Exact');

        elseif groupIndex >= 14 && groupIndex <= 16
            % Pressure-force exact values are the potential-flow forces.
            sphName = names(1);
            exactName = names(2);
            plot(axLocalValue,xData,localTable.(sphName), ...
                'LineWidth',1.25,'DisplayName',labels(1));
            plot(axLocalValue,xData,localTable.(exactName), ...
                '--','LineWidth',1.25,'DisplayName',labels(2));

            errorData = localTable.(sphName) - localTable.(exactName);
            plot(axLocalError,xData,errorData, ...
                'LineWidth',1.25,'DisplayName','SPH - Exact');

        elseif groupIndex == 26
            % This quantity is already an error norm; there is no separate
            % SPH/exact scalar pair for it.
            plot(axLocalValue,xData,localTable.(names(1)), ...
                'LineWidth',1.25,'DisplayName',labels(1));
            plot(axLocalError,xData,localTable.(names(1)), ...
                'LineWidth',1.25,'DisplayName',labels(1));

        else
            % Standard case: first two series are SPH and exact, remaining
            % series are error measures written by C++.
            for jj = 1:min(valueCount,numel(names))
                if ismember(names(jj),string(localTable.Properties.VariableNames))
                    plot(axLocalValue,xData,localTable.(names(jj)), ...
                        'LineWidth',1.25,'DisplayName',labels(jj));
                end
            end

            for jj = valueCount+1:numel(names)
                if ismember(names(jj),string(localTable.Properties.VariableNames))
                    plot(axLocalError,xData,localTable.(names(jj)), ...
                        'LineWidth',1.25,'DisplayName',labels(jj));
                end
            end
        end

        yline(axLocalError,0,':','HandleVisibility','off');

        xlabel(axLocalValue,xLabelText);
        ylabel(axLocalValue,yText);
        title(axLocalValue,sprintf('%s: SPH vs exact',titleText));
        grid(axLocalValue,'on');
        box(axLocalValue,'on');
        legend(axLocalValue,'Location','best');

        xlabel(axLocalError,xLabelText);
        ylabel(axLocalError,['Error in ' yText]);
        title(axLocalError,sprintf('%s: error',titleText));
        grid(axLocalError,'on');
        box(axLocalError,'on');
        legend(axLocalError,'Location','best');

        sgtitle(localComparisonTab,sprintf('%s %s',titleText,titleSuffix));

        hold(axLocalValue,'off');
        hold(axLocalError,'off');

    end

    %% ---------------------------------------------------------
    % Definition of every local comparison group
    % ----------------------------------------------------------

    function [names,labels,yText,titleText,boundaryOnly] = ...
            localComparisonSpec(groupIndex)

        boundaryOnly = false;

        switch groupIndex
            case 1
                names = ["uSPH","uExact","errorU"];
                labels = ["SPH","Exact","Error"];
                yText = 'u (m/s)';
                titleText = 'x velocity';
            case 2
                names = ["vSPH","vExact","errorV"];
                labels = ["SPH","Exact","Error"];
                yText = 'v (m/s)';
                titleText = 'y velocity';
            case 3
                names = ["urSPH","urExact","errorUr"];
                labels = ["SPH","Exact","Error"];
                yText = 'u_r (m/s)';
                titleText = 'Radial velocity';
            case 4
                names = ["utSPH","utExact","errorUt"];
                labels = ["SPH","Exact","Error"];
                yText = 'u_t (m/s)';
                titleText = 'Tangential velocity';
            case 5
                names = ["velocitySPH","velocityExact", ...
                         "errorVelocityMagnitude", ...
                         "absErrorVelocityMagnitude", ...
                         "errorVelocityVector"];
                labels = ["SPH magnitude","Exact magnitude", ...
                          "Signed magnitude error", ...
                          "Absolute magnitude error", ...
                          "Vector error magnitude"];
                yText = 'Velocity (m/s)';
                titleText = 'Velocity magnitude and errors';
            case 6
                names = ["CuSPH","CuExact","errorCu"];
                labels = ["SPH","Exact","Error"];
                yText = 'C_u';
                titleText = 'x-velocity coefficient';
            case 7
                names = ["CvSPH","CvExact","errorCv"];
                labels = ["SPH","Exact","Error"];
                yText = 'C_v';
                titleText = 'y-velocity coefficient';
            case 8
                names = ["CurSPH","CurExact","errorCur"];
                labels = ["SPH","Exact","Error"];
                yText = 'C_{ur}';
                titleText = 'Radial-velocity coefficient';
            case 9
                names = ["CutSPH","CutExact","errorCut"];
                labels = ["SPH","Exact","Error"];
                yText = 'C_{ut}';
                titleText = 'Tangential-velocity coefficient';
            case 10
                names = ["CVelocitySPH","CVelocityExact","errorCVelocity"];
                labels = ["SPH","Exact","Error"];
                yText = 'C_V';
                titleText = 'Velocity-magnitude coefficient';
            case 11
                names = ["pSPH","pExact","errorPressure"];
                labels = ["SPH","Exact","Error"];
                yText = 'Pressure (Pa)';
                titleText = 'Cylinder surface pressure';
                boundaryOnly = true;
            case 12
                names = ["CpSPH","CpExact","errorCp"];
                labels = ["SPH","Exact","Error"];
                yText = 'C_p';
                titleText = 'Pressure coefficient';
                boundaryOnly = true;
            case 13
                names = "tauWall";
                labels = "SPH wall shear";
                yText = '\tau_w (Pa)';
                titleText = 'Wall shear stress';
                boundaryOnly = true;
            case 14
                names = ["dFxPressureSPH","dFxExact"];
                labels = ["SPH pressure","Exact"];
                yText = 'dF_{x,p} (N/m)';
                titleText = 'Local x pressure force';
                boundaryOnly = true;
            case 15
                names = ["dFyPressureSPH","dFyExact"];
                labels = ["SPH pressure","Exact"];
                yText = 'dF_{y,p} (N/m)';
                titleText = 'Local y pressure force';
                boundaryOnly = true;
            case 16
                names = ["dFPressureSPH","dFExact"];
                labels = ["SPH pressure","Exact"];
                yText = 'dF_p (N/m)';
                titleText = 'Local pressure-force magnitude';
                boundaryOnly = true;
            case 17
                names = "dFxViscousSPH";
                labels = "SPH viscous";
                yText = 'dF_{x,\mu} (N/m)';
                titleText = 'Local x viscous force';
                boundaryOnly = true;
            case 18
                names = "dFyViscousSPH";
                labels = "SPH viscous";
                yText = 'dF_{y,\mu} (N/m)';
                titleText = 'Local y viscous force';
                boundaryOnly = true;
            case 19
                names = "dFViscousSPH";
                labels = "SPH viscous";
                yText = 'dF_\mu (N/m)';
                titleText = 'Local viscous-force magnitude';
                boundaryOnly = true;
            case 20
                names = ["dFxSPH","dFxExact","errorFx"];
                labels = ["SPH total","Exact","Error"];
                yText = 'dF_x (N/m)';
                titleText = 'Local total x force';
                boundaryOnly = true;
            case 21
                names = ["dFySPH","dFyExact","errorFy"];
                labels = ["SPH total","Exact","Error"];
                yText = 'dF_y (N/m)';
                titleText = 'Local total y force';
                boundaryOnly = true;
            case 22
                names = ["dFSPH","dFExact", ...
                         "errorForceMagnitude", ...
                         "absErrorForceMagnitude", ...
                         "errorForceVector"];
                labels = ["SPH magnitude","Exact magnitude", ...
                          "Signed magnitude error", ...
                          "Absolute magnitude error", ...
                          "Vector error magnitude"];
                yText = 'dF (N/m)';
                titleText = 'Local total-force magnitude and errors';
                boundaryOnly = true;
            case 23
                names = ["CfxLocalSPH","CfxLocalExact","errorCfx"];
                labels = ["SPH","Exact","Error"];
                yText = 'C_{fx,local}';
                titleText = 'Local x-force coefficient';
                boundaryOnly = true;
            case 24
                names = ["CfyLocalSPH","CfyLocalExact","errorCfy"];
                labels = ["SPH","Exact","Error"];
                yText = 'C_{fy,local}';
                titleText = 'Local y-force coefficient';
                boundaryOnly = true;
            case 25
                names = ["CFLocalSPH","CFLocalExact","errorCF"];
                labels = ["SPH","Exact","Error"];
                yText = 'C_{F,local}';
                titleText = 'Local resultant-force coefficient';
                boundaryOnly = true;
            otherwise
                names = "errorVelocityVectorNormalized";
                labels = "Normalized vector error";
                yText = '||u-u_{exact}|| / U_\infty';
                titleText = 'Normalized velocity-vector error';
        end

        names = string(names);
        labels = string(labels);

    end

    %% ---------------------------------------------------------
    % Read one solver-timestep block from the large local CSV
    % ----------------------------------------------------------

    function localTable = readLocalAngleSlice(globalIndex,layer,names)

        rowsPerSolverStep = 3*nTheta;
        firstDataRow = 2 + (globalIndex-1)*rowsPerSolverStep;
        lastDataRow = firstDataRow + rowsPerSolverStep - 1;

        selectedNames = unique( ...
            ["t","thetaIndex","layer","thetaDeg",string(names)], ...
            'stable');

        localOpts = detectImportOptions( ...
            localComparisonFilename, ...
            'VariableNamingRule','preserve');

        localOpts.SelectedVariableNames = cellstr(selectedNames);
        localOpts.DataLines = [firstDataRow lastDataRow];
        localOpts.ExtraColumnsRule = 'ignore';

        localTable = readtable( ...
            localComparisonFilename, ...
            localOpts);

        localTable = localTable(localTable.layer == layer,:);

    end

    %% ---------------------------------------------------------
    % Read one angle/layer history from the large local CSV
    % ----------------------------------------------------------

    function localTable = readLocalTimeHistory(layer,thetaIndex,names)

        selectedNames = unique( ...
            ["t","thetaIndex","layer","thetaDeg",string(names)], ...
            'stable');

        % Read only the columns required for the selected graph.
        % This scans the large local CSV but avoids loading the many unused
        % comparison columns into memory.
        localOpts = detectImportOptions( ...
            localComparisonFilename, ...
            'VariableNamingRule','preserve');

        localOpts.SelectedVariableNames = cellstr(selectedNames);
        localOpts.ExtraColumnsRule = 'ignore';

        localTable = readtable( ...
            localComparisonFilename, ...
            localOpts);

        keep = ...
            localTable.layer == layer & ...
            localTable.thetaIndex == thetaIndex;

        localTable = localTable(keep,:);

        if ~isempty(localTable)
            [~,orderHistory] = sort(localTable.t);
            localTable = localTable(orderHistory,:);
        end

    end

    % Instantaneous streamlines; these are not time-dependent trajectories.

    function updateStreamlines(frame)
        delete(hStreamlines(isgraphics(hStreamlines)));
        delete(hExactStreamlines(isgraphics(hExactStreamlines)));
        hStreamlines = gobjects(0);
        hExactStreamlines = gobjects(0);

        if ~showStreamlines
            return;
        end

        % ---------------- SPH streamline field ----------------
        U = reshape(frame(:,6),nTheta,[]);
        V = reshape(frame(:,7),nTheta,[]);
        U = U(angleOrder,:);
        V = V(angleOrder,:);
        U = [U(end,:);U;U(1,:)];
        V = [V(end,:);V;V(1,:)];

        Fu = griddedInterpolant({angles,streamRadii},U,'linear','none');
        Fv = griddedInterpolant({angles,streamRadii},V,'linear','none');
        gridU = Fu(queryTheta,queryR);
        gridV = Fv(queryTheta,queryR);
        gridU(streamMask) = NaN;
        gridV(streamMask) = NaN;

        paths = stream2(streamX,streamY,gridU,gridV,seedX,seedY,[0.10 4000]);
        hStreamlines = drawStreamlinePaths(ax,paths);

        % -------- analytical potential-flow streamline field --------
        [gridUExact,gridVExact] = potentialFlowVelocity(streamX,streamY);
        gridUExact(streamMask) = NaN;
        gridVExact(streamMask) = NaN;

        exactPaths = stream2( ...
            streamX,streamY,gridUExact,gridVExact,seedX,seedY,[0.10 4000]);
        hExactStreamlines = drawStreamlinePaths(axExact,exactPaths);
    end

    function handles = drawStreamlinePaths(targetAxes,paths)
        handles = gobjects(0);
        for kk = 1:numel(paths)
            points = paths{kk};
            if size(points,1) < 2
                continue;
            end

            rr = hypot(points(:,1)-streamCenter(1),points(:,2)-streamCenter(2));
            invalid = rr <= wallRadius | rr >= outerRadius | any(~isfinite(points),2);
            points(invalid,:) = NaN;

            hLine = plot3(targetAxes,points(:,1),points(:,2), ...
                1e-3*ones(size(points,1),1), ...
                'Color',[0.12 0.12 0.12],'LineWidth',0.65, ...
                'HitTest','off','PickableParts','none');
            handles(end+1,1) = hLine; %#ok<AGROW>
        end
    end

    %% Analytical potential-flow fields used only for visual comparison.
    function [uExact,vExact] = potentialFlowVelocity(xq,yq)
        if ~hasCylinder
            uExact = referenceUx*ones(size(xq));
            vExact = referenceUy*ones(size(yq));
            return;
        end

        dx = xq-streamCenter(1);
        dy = yq-streamCenter(2);
        rq = hypot(dx,dy);
        thq = atan2(dy,dx);
        beta = thq-deg2rad(flowAngleDeg);

        rqSafe = max(rq,innerRadius);
        a2r2 = (innerRadius^2)./(rqSafe.^2);

        ur = referenceSpeed.*(1-a2r2).*cos(beta);
        ut = -referenceSpeed.*(1+a2r2).*sin(beta);

        uExact = ur.*cos(thq)-ut.*sin(thq);
        vExact = ur.*sin(thq)+ut.*cos(thq);

        if hasCylinder
            inside = rq < innerRadius;
            uExact(inside) = NaN;
            vExact(inside) = NaN;
        end
    end

    function exactField = potentialFlowField(frame,column)
        xq = frame(:,1);
        yq = frame(:,2);
        [uExact,vExact] = potentialFlowVelocity(xq,yq);
        speedExact = hypot(uExact,vExact);

        % Bernoulli gauge pressure with p_inf = 0.
        pExact = 0.5*referenceDensity.*(referenceSpeed^2-speedExact.^2);

        % Same Tait EOS used by C++ initialization.
        eosArgument = max(1+pExact/eosB,eps);
        rhoExact = referenceDensity.*eosArgument.^(1/gammaEOS);

        exactField = zeros(size(xq));
        switch cacheNames(column)
            case "rho"
                exactField = rhoExact;
            case "drhodt"
                exactField(:) = 0;
            case "pressure"
                exactField = pExact;
            case "u"
                exactField = uExact;
            case "v"
                exactField = vExact;
            case "velocity"
                exactField = speedExact;
            case "dudt"
                exactField(:) = 0;
            case "dvdt"
                exactField(:) = 0;
            case "vorticity"
                exactField(:) = 0;
            otherwise
                exactField(:) = NaN;
        end
    end

    function [fieldMinLocal,fieldMaxLocal] = getComparisonColorLimits(column)
        if ~isempty(manualColorLimits)
            fieldMinLocal = manualColorLimits(1);
            fieldMaxLocal = manualColorLimits(2);
            return;
        end

        [fieldMinLocal,fieldMaxLocal] = getColorLimits(column);
        exactValues = potentialFlowField(firstFrame,column);
        exactValues = exactValues(fluid & isfinite(exactValues));

        if ~isempty(exactValues)
            exactValues = sort(exactValues);
            nExact = numel(exactValues);
            lo = max(1,ceil(colorPercentiles(1)/100*nExact));
            hi = min(nExact,ceil(colorPercentiles(2)/100*nExact));
            fieldMinLocal = min(fieldMinLocal,exactValues(lo));
            fieldMaxLocal = max(fieldMaxLocal,exactValues(hi));
        end

        if cacheNames(column) == "velocity"
            fieldMinLocal = 0;
        end

        if fieldMinLocal == fieldMaxLocal
            delta = max(1,abs(fieldMinLocal)*0.01);
            fieldMinLocal = fieldMinLocal-delta;
            fieldMaxLocal = fieldMaxLocal+delta;
        end
    end

    function setFieldColormap
        if plotVariable == "vorticity"
            anchors = [0.10 0.25 0.80;1 1 1;0.80 0.10 0.10];
            colours = interp1([0 0.5 1],anchors,linspace(0,1,257));
            colormap(ax,colours);
            colormap(axExact,colours);
        else
            colormap(ax,turbo);
            colormap(axExact,turbo);
        end
    end

    %% ---------------------------------------------------------

    % Read a numeric C++ constant from SimulationConstants.csv

    % ----------------------------------------------------------

    function value = constantNumber(name)

        idx = find(strcmpi(constantNames,string(name)),1);

        if isempty(idx)

            error("Required C++ constant '%s' is missing from %s.",name,constantsFilename);

        end

        value = str2double(constantValues(idx));

        if ~isscalar(value) || isnan(value)

            error("C++ constant '%s' is not numeric in %s.",name,constantsFilename);

        end

    end

    %% ---------------------------------------------------------

    % Read a text C++ constant from SimulationConstants.csv

    % ----------------------------------------------------------

    function value = constantText(name)

        idx = find(strcmpi(constantNames,string(name)),1);

        if isempty(idx)

            error("Required C++ constant '%s' is missing from %s.",name,constantsFilename);

        end

        value = strtrim(constantValues(idx));

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
