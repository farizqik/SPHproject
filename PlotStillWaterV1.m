clear;
clc;
close all;

%% =========================================================
% USER SETTINGS
% ==========================================================

simulationFolder = ...
    "Stillwater_dp_0.500000_h_1.000000_Nparticles_672_wendland";

% plotVariable = "pressure";
% plotVariable = "density";
 plotVariable = "velocity";

saveVideo = true;

videoName = simulationFolder + ".mp4";


%% =========================================================
% READ PARAMETERS FROM FOLDER NAME
% ==========================================================

folderName = string(simulationFolder);

tokens = regexp( ...
    folderName, ...
    'Stillwater_dp_([0-9.]+)_h_([0-9.]+)_Nparticles_([0-9]+)_([A-Za-z]+)$', ...
    'tokens', ...
    'once');

if isempty(tokens)
    error("Folder name does not match the C++ format.");
end

dp = str2double(tokens{1});
h = str2double(tokens{2});
NparticlesFolder = str2double(tokens{3});
kernel = string(tokens{4});

hdp = h/dp;


%% =========================================================
% GEOMETRY
% Same as C++
% ==========================================================

tanklength = 25.0;
tankheight = 10.0;

waterlength = 15.0;

freeboard = 2.0;
waterheight = tankheight - freeboard;

boundthick = 3*dp;


%% =========================================================
% FIND OUTPUT CSV FILES
% ==========================================================

pattern = sprintf( ...
    "Stillwater_hdp_%.6f_t_*.csv", ...
    hdp);

files = dir( ...
    fullfile(simulationFolder, pattern));

if isempty(files)
    error("No CSV files found.");
end


%% =========================================================
% EXTRACT TIME FROM FILENAMES
% ==========================================================

time = zeros(length(files),1);

for n = 1:length(files)

    timeToken = regexp( ...
        files(n).name, ...
        '_t_([0-9]+\.[0-9]+)', ...
        'tokens', ...
        'once');

    if isempty(timeToken)

        error( ...
            "Could not extract time from %s", ...
            files(n).name);

    end

    time(n) = str2double(timeToken{1});

end


%% Sort according to physical time

[time, order] = sort(time);

files = files(order);


fprintf("\n");
fprintf("OUTPUT FILES\n");
fprintf("--------------------------------------\n");

fprintf("Timesteps    = %d\n", length(files));
fprintf("Initial time = %.6f s\n", time(1));
fprintf("Final time   = %.6f s\n", time(end));


%% =========================================================
% CSV IMPORT SETTINGS
%
% C++ CSV structure:
%
% line 1 = h/dp
% line 2 = metadata header
% line 3 = metadata values
% line 4 = blank
% line 5 = time
% line 6 = particle variable names
% line 7 = first particle
%
% Columns:
%
% 1  ID
% 2  x
% 3  y
% 4  empty
% 5  rho
% 6  drhodt
% 7  pressure
% 8  empty
% 9  u
% 10 v
% 11 empty
% 12 dudt
% 13 dvdt
% 14 type
% ==========================================================

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
     "type"];

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
% READ FIRST TIMESTEP
% ==========================================================

filename = fullfile( ...
    files(1).folder, ...
    files(1).name);

T = readtable(filename, opts);


%% =========================================================
% READ INITIAL PARTICLE DATA
% ==========================================================

ID0 = T.ID;

x0 = T.x;
y0 = T.y;

rhoInitial = T.rho;
drhodtInitial = T.drhodt;
pressureInitial = T.pressure;

uInitial = T.u;
vInitial = T.v;

type0 = lower(strtrim(string(T.type)));


%% =========================================================
% PARTICLE IDENTITY
%
% Particle identity now comes DIRECTLY from C++.
% No x/y criterion.
% No geometry reconstruction.
% ==========================================================

boundary = type0 == "boundary";
fluid = type0 == "fluid";


%% =========================================================
% PARTICLE COUNTS
% ==========================================================

Nparticles = height(T);

Nboundary = sum(boundary);
Nfluid = sum(fluid);


fprintf("\n");
fprintf("PARTICLE INFORMATION\n");
fprintf("--------------------------------------\n");

fprintf("Boundary particles = %d\n", Nboundary);
fprintf("Fluid particles    = %d\n", Nfluid);
fprintf("Total particles    = %d\n", Nparticles);


%% =========================================================
% CHECK AGAINST FOLDER NAME
% ==========================================================

if Nparticles ~= NparticlesFolder

    error( ...
        "CSV contains %d particles but folder says %d.", ...
        Nparticles, ...
        NparticlesFolder);

end


%% =========================================================
% CHECK PARTICLE IDs
%
% Expected C++ IDs:
%
% 0, 1, 2, ..., Nparticles-1
% ==========================================================

expectedID = (0:Nparticles-1)';

if ~isequal(ID0, expectedID)

    warning( ...
        "Particle IDs are not sequential from 0 to %d.", ...
        Nparticles-1);

end


%% =========================================================
% DISPLAY ID RANGES
% ==========================================================

boundaryIDs = ID0(boundary);
fluidIDs = ID0(fluid);

fprintf("\n");
fprintf("PARTICLE ID RANGES\n");
fprintf("--------------------------------------\n");

fprintf( ...
    "Boundary IDs : %d - %d\n", ...
    min(boundaryIDs), ...
    max(boundaryIDs));

fprintf( ...
    "Fluid IDs    : %d - %d\n", ...
    min(fluidIDs), ...
    max(fluidIDs));


%% =========================================================
% INITIAL PARTICLE CHECK
% ==========================================================

firstBoundaryRow = find(boundary,1,'first');
lastBoundaryRow  = find(boundary,1,'last');

firstFluidRow = find(fluid,1,'first');


fprintf("\n");
fprintf("INITIAL PARTICLE CHECK\n");
fprintf("--------------------------------------\n");

fprintf( ...
    "First boundary:\nID = %d, x = %.3f, y = %.3f\n", ...
    ID0(firstBoundaryRow), ...
    x0(firstBoundaryRow), ...
    y0(firstBoundaryRow));

fprintf( ...
    "\nLast boundary:\nID = %d, x = %.3f, y = %.3f\n", ...
    ID0(lastBoundaryRow), ...
    x0(lastBoundaryRow), ...
    y0(lastBoundaryRow));

fprintf( ...
    "\nFirst fluid:\nID = %d, x = %.3f, y = %.3f\n", ...
    ID0(firstFluidRow), ...
    x0(firstFluidRow), ...
    y0(firstFluidRow));


%% =========================================================
% AXIS LIMITS
%
% Based on actual initial coordinates from C++
% ==========================================================

xmin = min(x0) - dp;
xmax = max(x0) + dp;

ymin = min(y0) - dp;
ymax = max(y0) + dp;


%% =========================================================
% INITIAL FIELD
% ==========================================================

switch plotVariable

    case "pressure"

        field0 = pressureInitial;

        colorLabel = ...
            "Pressure (Pa)";


    case "density"

        field0 = rhoInitial;

        colorLabel = ...
            "Density (kg/m^3)";


    case "velocity"

        field0 = sqrt( ...
            uInitial.^2 + ...
            vInitial.^2);

        colorLabel = ...
            "Velocity magnitude (m/s)";


    otherwise

        error( ...
            "Use pressure, density, or velocity.");

end


%% =========================================================
% COLOUR LIMITS
% ==========================================================

rho0_ref = 1000.0;
g = 9.81;


switch plotVariable

    case "pressure"

        fieldMin = 0.0;

        fieldMax = ...
            rho0_ref * g * waterheight;


    case "density"

        initialFluidField = field0(fluid);

        initialFluidField = ...
            initialFluidField( ...
            isfinite(initialFluidField));

        fieldMin = min(initialFluidField);
        fieldMax = max(initialFluidField);

        if fieldMin == fieldMax

            fieldMin = fieldMin - 1;
            fieldMax = fieldMax + 1;

        end


    case "velocity"

        fieldMin = 0.0;

        fieldMax = ...
            sqrt(g*waterheight);

end


fprintf("\n");
fprintf("COLOUR RANGE\n");
fprintf("--------------------------------------\n");

fprintf("Minimum = %.6f\n", fieldMin);
fprintf("Maximum = %.6f\n", fieldMax);


%% =========================================================
% IDEAL HYDROSTATIC PRESSURE
% ==========================================================

if plotVariable == "pressure"

    pressureHydrostatic = ...
        rho0_ref * g * ...
        (waterheight - y0(fluid));

    pressureHydrostatic = ...
        max(pressureHydrostatic,0.0);


    figHydro = figure( ...
        'Color','w', ...
        'Position',[1050 100 900 700]);

    hold on;


    %% Fluid

    scatter( ...
        x0(fluid), ...
        y0(fluid), ...
        25, ...
        pressureHydrostatic, ...
        'filled', ...
        'MarkerEdgeColor','none');


    %% Boundary

    scatter( ...
        x0(boundary), ...
        y0(boundary), ...
        80, ...
        [0 0 0], ...
        'filled');


    axis equal;

    xlim([xmin xmax]);
    ylim([ymin ymax]);

    xlabel('x (m)');
    ylabel('y (m)');

    grid off;
    box off;

    colormap(turbo);

    cHydro = colorbar;

    cHydro.Label.String = ...
        'Ideal hydrostatic pressure (Pa)';

    clim([fieldMin fieldMax]);

    title( ...
        'Ideal Hydrostatic Pressure');

end

%% =========================================================
% KINETIC ENERGY VS TIME
% ==========================================================

KE = zeros(length(files),1);

for n = 1:length(files)

    filename = fullfile( ...
        files(n).folder, ...
        files(n).name);

    % C++ CSV:
    % line 2 = L2norm Pressure, KE, dp, Nparticles
    % line 3 = numerical values
    %
    % Read line 3, columns A:D
    metadata = readmatrix( ...
        filename, ...
        'Range','A3:D3');

    % Column 2 contains kinetic energy
    KE(n) = metadata(2);

end


%% =========================================================
% PLOT KINETIC ENERGY VS TIME
% ==========================================================

figure( ...
    'Color','w', ...
    'Position',[200 200 800 600]);

plot( ...
    time, ...
    KE, ...
    'LineWidth',1.5);

xlabel('Time (s)');
ylabel('Kinetic Energy (J)');

title('Kinetic Energy vs Time');

grid on;
box on;

set(gca, ...
    'FontSize',12, ...
    'LineWidth',1);

fprintf("\nKINETIC ENERGY\n");
fprintf("--------------------------------------\n");
fprintf("Initial KE = %.6e J\n", KE(1));
fprintf("Maximum KE = %.6e J\n", max(KE));
fprintf("Final KE   = %.6e J\n", KE(end));



%% =========================================================
% CREATE SIMULATION FIGURE
% ==========================================================

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 900 700]);

hold on;


%% =========================================================
% FLUID PARTICLES
% ==========================================================

hFluid = scatter( ...
    x0(fluid), ...
    y0(fluid), ...
    25, ...
    field0(fluid), ...
    'filled', ...
    'MarkerEdgeColor','none');


%% =========================================================
% BOUNDARY PARTICLES
% ==========================================================

hBoundary = scatter( ...
    x0(boundary), ...
    y0(boundary), ...
    80, ...
    [0 0 0], ...
    'filled');


%% =========================================================
% FIGURE SETTINGS
% ==========================================================

axis equal;

xlim([xmin xmax]);
ylim([ymin ymax]);

xlabel('x (m)');
ylabel('y (m)');

grid off;
box off;

colormap(turbo);

c = colorbar;

c.Label.String = ...
    colorLabel;

clim([fieldMin fieldMax]);

title(sprintf( ...
    'SPH Still Water, t = %.3f s', ...
    time(1)));

set(gca, ...
    'FontSize',12, ...
    'LineWidth',1);


%% =========================================================
% VIDEO
% ==========================================================

if saveVideo

    video = VideoWriter( ...
        videoName, ...
        'MPEG-4');

    video.FrameRate = 5;

    open(video);

end


%% =========================================================
% ANIMATION LOOP
% ==========================================================

for n = 1:length(files)

    %% -----------------------------------------------------
    % Read timestep
    % ------------------------------------------------------

    filename = fullfile( ...
        files(n).folder, ...
        files(n).name);

    T = readtable(filename, opts);


    %% -----------------------------------------------------
    % Check particle number
    % ------------------------------------------------------

    if height(T) ~= Nparticles

        error( ...
            "Particle count changed at t = %.6f.", ...
            time(n));

    end


    %% -----------------------------------------------------
    % Read variables
    % ------------------------------------------------------

    ID = T.ID;

    x = T.x;
    y = T.y;

    rho = T.rho;
    pressure = T.pressure;

    u = T.u;
    v = T.v;

    type = ...
        lower(strtrim(string(T.type)));


    %% -----------------------------------------------------
    % Check ID consistency
    %
    % Particle identity must remain unchanged.
    % ------------------------------------------------------

    if ~isequal(ID, ID0)

        error( ...
            "Particle ID ordering changed at t = %.6f.", ...
            time(n));

    end


    %% -----------------------------------------------------
    % Check particle types
    % ------------------------------------------------------

    if ~isequal(type, type0)

        error( ...
            "Particle type changed at t = %.6f.", ...
            time(n));

    end


    %% -----------------------------------------------------
    % Select plotted field
    % ------------------------------------------------------

    switch plotVariable

        case "pressure"

            field = pressure;


        case "density"

            field = rho;


        case "velocity"

            field = ...
                sqrt(u.^2 + v.^2);

    end


    %% -----------------------------------------------------
    % Check NaN / Inf
    % ------------------------------------------------------

    valid = ...
        isfinite(x) & ...
        isfinite(y) & ...
        isfinite(field);


    if any(~valid)

        fprintf( ...
            "Warning: %d invalid particles at t = %.4f s\n", ...
            sum(~valid), ...
            time(n));

    end


    xplot = x;
    yplot = y;
    fieldplot = field;

    xplot(~valid) = NaN;
    yplot(~valid) = NaN;
    fieldplot(~valid) = NaN;


    %% -----------------------------------------------------
    % UPDATE FLUID PARTICLES
    %
    % Particle identity comes directly from C++ "type".
    % It does NOT depend on x or y.
    % ------------------------------------------------------

    set( ...
        hFluid, ...
        'XData', xplot(fluid), ...
        'YData', yplot(fluid), ...
        'CData', fieldplot(fluid));


    %% -----------------------------------------------------
    % UPDATE BOUNDARY PARTICLES
    % ------------------------------------------------------

    set( ...
        hBoundary, ...
        'XData', xplot(boundary), ...
        'YData', yplot(boundary));


    %% -----------------------------------------------------
    % TITLE
    % ------------------------------------------------------

    title(sprintf( ...
        'SPH Still Water, t = %.3f s', ...
        time(n)));

    drawnow;

    pause(0.0001);


    %% -----------------------------------------------------
    % VIDEO FRAME
    % ------------------------------------------------------

    if saveVideo

        frame = getframe(fig);

        writeVideo( ...
            video, ...
            frame);

    end

end


%% =========================================================
% CLOSE VIDEO
% ==========================================================

if saveVideo

    close(video);

    fprintf( ...
        "\nVideo saved as %s\n", ...
        videoName);

end