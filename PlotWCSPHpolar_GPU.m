clear;
clc;
close all;

%% =========================================================
% USER SETTINGS
% ==========================================================

simulationFolder = ...
    "WCSPHpolar_dp_0.250000_h_0.500000_Nparticles_2528_wendland";

plotVariable = "pressure";
trackedParticleID = 1048;      
saveVideo = true;
saveParticleHistory = true;
colorLimits = [];              

fluidMarkerSize = 30;
boundaryMarkerSize = 150;
videoFrameRate = 60;
animationFPS = 60;

useGPU = true;
maximumGPUMemoryFraction = 0.60;

%% =========================================================
% FIND AND SORT C++ OUTPUT FILES
% ==========================================================

files = dir(fullfile(simulationFolder,"WCSPHpolar_hdp_*_t_*.csv"));

if isempty(files)
    error("No C++ timestep CSV files were found in:\n%s",simulationFolder);
end

numberOfFrames = numel(files);
time = nan(numberOfFrames,1);

for k = 1:numberOfFrames
    token = regexp(files(k).name, ...
        '_t_([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\.csv$', ...
        'tokens','once');

    if isempty(token)
        error("Could not read time from filename: %s",files(k).name);
    end

    time(k) = str2double(token{1});
end

[time,order] = sort(time);
files = files(order);

%% =========================================================
% EXACT 17-COLUMN C++ FORMAT
% ==========================================================

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

cacheNames = ...
    ["r","z","rho","drhodt","pressure","u_r","u_z", ...
     "velocity","pairAr","pairAz","totalAr","totalAz"];

allowedVariables = cacheNames(3:end);

if ~any(plotVariable == allowedVariables)
    error("plotVariable must be one of: %s",strjoin(allowedVariables,", "));
end

fieldColumn = find(cacheNames == plotVariable,1);
firstFilename = fullfile(files(1).folder,files(1).name);
T0 = readtable(firstFilename,opts);

metadata0 = readmatrix(firstFilename,'Range','A3:D3');
hdpData = readmatrix(firstFilename,'Range','A1:B1');

if numel(metadata0) < 4 || numel(hdpData) < 2
    error("The C++ metadata rows are incomplete in %s",files(1).name);
end

dp = metadata0(3);
NfluidMetadata = round(metadata0(4));
hdp = hdpData(2);
h = hdp*dp;

ID0 = T0.ID;
type0 = lower(strtrim(string(T0.Type)));
boundary = type0 == "boundary";
fluid = type0 == "fluid";

if any(~(boundary | fluid))
    badRows = find(~(boundary | fluid),1);
    error("Invalid Type entry in first CSV at data row %d.",badRows);
end

Nparticles = height(T0);
Nboundary = sum(boundary);
Nfluid = sum(fluid);

if Nfluid ~= NfluidMetadata
    error("CSV Type column gives %d fluid particles, but metadata gives %d.", ...
        Nfluid,NfluidMetadata);
end

if ~isequal(ID0,(0:Nparticles-1)')
    error("Particle IDs must be sequential from 0 to %d.",Nparticles-1);
end

%% =========================================================
% INITIALISE GPU AND ALLOCATE FRAME CACHE
% ==========================================================

gpuEnabled = false;
gpuInfo = [];

if useGPU
    try
        gpuInfo = gpuDevice;
        gpuEnabled = true;
    catch gpuError
        warning("GPU unavailable; using a CPU frame cache instead.\n%s", ...
            gpuError.message);
    end
end

numberOfCachedFields = numel(cacheNames);
cacheBytes = double(Nparticles)*double(numberOfCachedFields)* ...
    double(numberOfFrames)*8;

if gpuEnabled && ...
        cacheBytes > maximumGPUMemoryFraction*double(gpuInfo.AvailableMemory)
    warning("The complete frame cache requires %.2f GiB, which exceeds the selected %.0f%% GPU-memory limit. Using CPU memory instead.", ...
        cacheBytes/1024^3,100*maximumGPUMemoryFraction);
    gpuEnabled = false;
end

frameGPU = [];
frameCPU = [];

if gpuEnabled
    frameGPU = gpuArray.zeros( ...
        Nparticles,numberOfCachedFields,numberOfFrames,'double');
    fprintf("GPU cache enabled on %s (%.2f GiB).\n", ...
        gpuInfo.Name,cacheBytes/1024^3);
else
    frameCPU = zeros( ...
        Nparticles,numberOfCachedFields,numberOfFrames,'double');
    fprintf("CPU cache enabled (%.2f GiB).\n",cacheBytes/1024^3);
end

%% =========================================================
% READ EACH CSV ONCE AND CACHE ITS NUMERICAL COLUMNS
% ==========================================================

KE = nan(numberOfFrames,1);
progressInterval = max(1,floor(numberOfFrames/20));

for k = 1:numberOfFrames
    filename = fullfile(files(k).folder,files(k).name);

    if k == 1
        T = T0;
        metadata = metadata0;
    else
        T = readtable(filename,opts);
        metadata = readmatrix(filename,'Range','A3:D3');
    end

    if height(T) ~= Nparticles || ~isequal(T.ID,ID0)
        error("Particle layout changed at t = %.6f s.",time(k));
    end

    % if ~isequal(lower(strtrim(string(T.Type))),type0)
    %     error("Particle Type changed at t = %.6f s.",time(k));
    % end

    frameValues = [T.r,T.z,T.rho,T.drhodt,T.pressure,T.u_r,T.u_z, ...
        T.velocity,T.pairAr,T.pairAz,T.totalAr,T.totalAz];

    if gpuEnabled
        frameGPU(:,:,k) = gpuArray(frameValues);
    else
        frameCPU(:,:,k) = frameValues;
    end

    KE(k) = metadata(2);

    if mod(k,progressInterval) == 0 || k == numberOfFrames
        fprintf("Loaded frame %d of %d.\n",k,numberOfFrames);
    end
end

clear T frameValues;

fprintf("\nC++ OUTPUT INFORMATION\n");
fprintf("--------------------------------------\n");
fprintf("Files              = %d\n",numberOfFrames);
fprintf("Initial time       = %.6f s\n",time(1));
fprintf("Final time         = %.6f s\n",time(end));
fprintf("dp                 = %.6f m\n",dp);
fprintf("h                  = %.6f m\n",h);
fprintf("Boundary particles = %d\n",Nboundary);
fprintf("Fluid particles    = %d\n",Nfluid);
fprintf("Total particles    = %d\n",Nparticles);

%% =========================================================
% GLOBAL COLOUR LIMITS
% ==========================================================

if isempty(colorLimits)
    if gpuEnabled
        allFieldValues = reshape(frameGPU(fluid,fieldColumn,:),[],1);
        allFieldValues = allFieldValues(isfinite(allFieldValues));
        allFieldValues = sort(allFieldValues);
        numberOfValues = numel(allFieldValues);

        if numberOfValues == 0
            error("No finite fluid values exist for %s.",plotVariable);
        end

        lowerIndex = max(1,round(0.05*numberOfValues));
        upperIndex = min(numberOfValues,round(0.95*numberOfValues));
        fieldMin = gather(allFieldValues(lowerIndex));
        fieldMax = gather(allFieldValues(upperIndex));
        clear allFieldValues;
    else
        allFieldValues = reshape(frameCPU(fluid,fieldColumn,:),[],1);
        allFieldValues = sort(allFieldValues(isfinite(allFieldValues)));
        numberOfValues = numel(allFieldValues);

        if numberOfValues == 0
            error("No finite fluid values exist for %s.",plotVariable);
        end

        lowerIndex = max(1,round(0.05*numberOfValues));
        upperIndex = min(numberOfValues,round(0.95*numberOfValues));
        fieldMin = allFieldValues(lowerIndex);
        fieldMax = allFieldValues(upperIndex);
        clear allFieldValues;
    end

    if plotVariable == "velocity"
        fieldMin = 0.0;
    end

    if fieldMin == fieldMax
        fieldMax = fieldMin + 1.0;
    end
else
    fieldMin = colorLimits(1);
    fieldMax = colorLimits(2);
end

%% =========================================================
% FIELD LABEL
% ==========================================================

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
end

%% =========================================================
% KINETIC ENERGY
% ==========================================================

figure('Color','w','Position',[150 150 800 550]);
plot(time,KE,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Kinetic energy (J)');
title('Kinetic Energy vs Time');
grid on;
box on;

%% =========================================================
% SINGLE-SIDE PLOT GEOMETRY
% ==========================================================

if gpuEnabled
    firstFrame = gather(frameGPU(:,:,1));
else
    firstFrame = frameCPU(:,:,1);
end

r0 = firstFrame(:,1);
z0 = firstFrame(:,2);
field0 = firstFrame(:,fieldColumn);

boundthick = 3*dp;
leftWallR = (-0.5*dp:-dp:-boundthick+0.5*dp)';
leftWallZ = (min(z0):dp:max(z0))';
[leftWallRGrid,leftWallZGrid] = meshgrid(leftWallR,leftWallZ);
leftWallRPlot = leftWallRGrid(:);
leftWallZPlot = leftWallZGrid(:);

xBoundaryPlot = [r0(boundary);leftWallRPlot];
zBoundaryPlot = [z0(boundary);leftWallZPlot];

xmin = min(leftWallRPlot)-dp;
xmax = max(r0)+dp;
ymin = min(z0)-dp;
ymax = max(z0)+dp;

%% =========================================================
% FIGURE AND VIDEO
% ==========================================================

fig = figure('Color','w','Position',[100 100 1200 550]);
hold on;

hFluid = scatter(r0(fluid),z0(fluid),fluidMarkerSize,field0(fluid), ...
    'filled','MarkerEdgeColor','none');
hBoundary = scatter(xBoundaryPlot,zBoundaryPlot,boundaryMarkerSize, ...
    [0 0 0],'filled');

dcm = datacursormode(fig);
dcm.Enable = 'on';

axis equal;
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('r (m)');
ylabel('z (m)');
grid off;
box off;
colormap(turbo);
c = colorbar;
c.Label.String = colorLabel;
caxis([fieldMin fieldMax]);
title(sprintf( ...
    'Axisymmetric SPH Radial Dam-Break, t = %.3f s',time(1)));

if isprop(hFluid,'DataTipTemplate')
    hFluid.DataTipTemplate.DataTipRows(end+1) = ...
        dataTipTextRow('Particle ID',ID0(fluid));
end

[~,simulationName] = fileparts(char(simulationFolder));
videoName = sprintf('%s_%s_single_side.mp4', ...
    simulationName,char(plotVariable));

if saveVideo
    video = VideoWriter(videoName,'MPEG-4');
    video.FrameRate = videoFrameRate;
    open(video);
end

%% =========================================================
% TRACKED PARTICLE STORAGE
% ==========================================================

trackingEnabled = ~isempty(trackedParticleID);

if trackingEnabled
    trackedRow = find(ID0 == trackedParticleID,1);
    if isempty(trackedRow)
        error("Particle ID %d does not exist.",trackedParticleID);
    end

    trackR = nan(numberOfFrames,1);
    trackZ = nan(numberOfFrames,1);
    trackUr = nan(numberOfFrames,1);
    trackUz = nan(numberOfFrames,1);
    trackRho = nan(numberOfFrames,1);
    trackPressure = nan(numberOfFrames,1);
    trackPairAr = nan(numberOfFrames,1);
    trackPairAz = nan(numberOfFrames,1);
    trackTotalAr = nan(numberOfFrames,1);
    trackTotalAz = nan(numberOfFrames,1);
end

%% =========================================================
% ANIMATION
% ==========================================================

for k = 1:numberOfFrames
    if gpuEnabled
        frame = gather(frameGPU(:,:,k));
    else
        frame = frameCPU(:,:,k);
    end

    rFluidPlot = frame(fluid,1);
    zFluidPlot = frame(fluid,2);
    fieldFluidPlot = frame(fluid,fieldColumn);
    valid = isfinite(rFluidPlot) & isfinite(zFluidPlot) & ...
        isfinite(fieldFluidPlot);

    rFluidPlot(~valid) = NaN;
    zFluidPlot(~valid) = NaN;
    fieldFluidPlot(~valid) = NaN;

    set(hFluid,'XData',rFluidPlot,'YData',zFluidPlot, ...
        'CData',fieldFluidPlot);
    set(hBoundary,'XData',[frame(boundary,1);leftWallRPlot], ...
        'YData',[frame(boundary,2);leftWallZPlot]);

    title(sprintf( ...
        'Axisymmetric SPH Radial Dam-Break, t = %.3f s',time(k)));

    if trackingEnabled
        trackR(k) = frame(trackedRow,1);
        trackZ(k) = frame(trackedRow,2);
        trackRho(k) = frame(trackedRow,3);
        trackPressure(k) = frame(trackedRow,5);
        trackUr(k) = frame(trackedRow,6);
        trackUz(k) = frame(trackedRow,7);
        trackPairAr(k) = frame(trackedRow,9);
        trackPairAz(k) = frame(trackedRow,10);
        trackTotalAr(k) = frame(trackedRow,11);
        trackTotalAz(k) = frame(trackedRow,12);
    end

    if saveVideo
        drawnow;
        writeVideo(video,getframe(fig));
    else
        drawnow;
        pause(1/animationFPS);
    end
end

if saveVideo
    close(video);
    fprintf("Video saved as %s\n",videoName);
end

%% =========================================================
% TRACKED PARTICLE OUTPUT
% ==========================================================

if trackingEnabled
    particleID = repmat(trackedParticleID,numberOfFrames,1);
    particleHistory = table( ...
        time,particleID,trackR,trackZ,trackUr,trackUz,trackRho, ...
        trackPressure,trackPairAr,trackPairAz,trackTotalAr,trackTotalAz, ...
        'VariableNames', ...
        ["time","ID","r","z","u_r","u_z","rho","pressure", ...
         "pairAr","pairAz","totalAr","totalAz"]);

    if saveParticleHistory
        historyFilename = fullfile(simulationFolder, ...
            sprintf('particle_%d_history.csv',trackedParticleID));
        writetable(particleHistory,historyFilename);
        fprintf("Particle history saved as %s\n",historyFilename);
    end

    figure('Color','w','Position',[150 80 1100 750]);
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
    nexttile; plot(time,trackR,'LineWidth',1.3); ylabel('r (m)'); grid on;
    nexttile; plot(time,trackZ,'LineWidth',1.3); ylabel('z (m)'); grid on;
    nexttile; plot(time,trackUr,'LineWidth',1.3); ylabel('u_r (m/s)'); grid on;
    nexttile; plot(time,trackUz,'LineWidth',1.3); ylabel('u_z (m/s)'); grid on;
    nexttile; plot(time,trackPressure,'LineWidth',1.3); ...
        xlabel('Time (s)'); ylabel('Pressure (Pa)'); grid on;
    nexttile; plot(time,trackTotalAz,'LineWidth',1.3); yline(0,'k--'); ...
        xlabel('Time (s)'); ylabel('Total a_z (m/s^2)'); grid on;
    sgtitle(sprintf('C++ History for Particle ID %d',trackedParticleID));
end
