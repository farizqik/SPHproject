clear;
clc;
close all;

%% =========================================================
%  USER SETTINGS
% ==========================================================

simulationFolder = ...
    "Stillwater_dx_0.500000_h_1.000000_Ncol_50_Nrow_50_wendland";

folderName = string(simulationFolder);

tokens = regexp( ...
    folderName, ...
    'Stillwater_dx_([0-9.]+)_h_([0-9.]+)_Ncol_([0-9]+)_Nrow_([0-9]+)_([A-Za-z]+)$', ...
    'tokens', ...
    'once');

if isempty(tokens)

    error( ...
        "Folder name does not match expected format:\n%s", ...
        simulationFolder);

end


%% =========================================================
% SIMULATION PARAMETERS FROM FOLDER
% ==========================================================

dx = str2double(tokens{1});

h = str2double(tokens{2});

Ncol = str2double(tokens{3});

Nrow = str2double(tokens{4});

kernel = string(tokens{5});


%% Assuming square initial particle spacing

dy = dx;


%% Calculate dx/h automatically

dxh = dx/h;
dyh = dy/h;

%% Choose what to plot

% Available:
% "pressure"
% "density"
% "velocity"

plotVariable = "pressure";
% plotVariable = "density";
% plotVariable = "velocity";

% Save animation?
saveVideo = true;

videoName = string(simulationFolder)+ ".mp4" ;

%% =========================================================
%  DOMAIN SETTINGS
% ==========================================================

Nboundary = 3;

fluidXmin = 0.0;

fluidXmax = ...
    (Ncol - 2*Nboundary - 1) * dx;

fluidYmin = 0.0;

fluidYmax = ...
    (Nrow - Nboundary - 1) * dy;

fprintf("Fluid domain:\n");
fprintf("x = %.2f to %.2f m\n", ...
    fluidXmin, fluidXmax);

fprintf("y = %.2f to %.2f m\n", ...
    fluidYmin, fluidYmax);

%% =========================================================
%  FIND OUTPUT FILES
% ==========================================================

pattern = sprintf( ...
    "Stillwater_dxh_%.6f_t_*.csv", ...
    dxh);

files = dir(fullfile(simulationFolder, pattern));

if isempty(files)
    error("No CSV files found.");
end

%% =========================================================
%  EXTRACT TIME FROM FILE NAMES
% ==========================================================

time = zeros(length(files),1);

for n = 1:length(files)

    tokens = regexp( ...
        files(n).name, ...
        '_t_([0-9]+\.[0-9]+)?', ...
        'tokens', ...
        'once');

    if isempty(tokens)

        error( ...
            "Could not extract time from file: %s", ...
            files(n).name);

    end

    time(n) = str2double(tokens{1});

end

%% Sort according to physical time

[time, order] = sort(time);

files = files(order);

fprintf("\n");
fprintf("Number of timesteps found = %d\n", ...
    length(files));

fprintf("Initial time = %.6f s\n", ...
    time(1));

fprintf("Final time   = %.6f s\n", ...
    time(end));

%% =========================================================
%  READ FIRST TIMESTEP
% ==========================================================

filename = fullfile( ...
    files(1).folder, ...
    files(1).name);

data = readmatrix( ...
    filename, ...
    'NumHeaderLines', 7);

%% =========================================================
%  COLUMN MAPPING
% ==========================================================
%
% C++ output:
%
% column 1  = x
% column 2  = y
% column 3  = empty
% column 4  = rho
% column 5  = drhodt
% column 6  = pressure
% column 7  = empty
% column 8  = u
% column 9  = v
% column 10 = empty
% column 11 = dudt
% column 12 = dvdt

x0 = data(:,1);
y0 = data(:,2);

rhoInitial = data(:,4);
pressureInitial = data(:,6);


uInitial = data(:,8);
vInitial = data(:,9);

%% =========================================================
%  IDENTIFY BOUNDARY AND FLUID PARTICLES
% ==========================================================

tol = 1e-4 * dx  ;

boundary = ...
    x0 < fluidXmin - tol | ...
    x0 > fluidXmax + tol | ...
    y0 < fluidYmin - tol;



fluid = ~boundary;

fprintf("\n");
fprintf("Number of fluid particles    = %d\n", ...
    sum(fluid));

fprintf("Number of boundary particles = %d\n", ...
    sum(boundary));

%% =========================================================
%  AXIS LIMITS
% ==========================================================

xmin = min(x0) - dx;
xmax = max(x0) + dx;

ymin = min(y0) - dy;
ymax = max(y0) + dy;

%% =========================================================
%  INITIAL FIELD
% ==========================================================

switch plotVariable

    case "pressure"

        field0 = pressureInitial;
        colorLabel = "Pressure (Pa)";

    case "density"

        field0 = rhoInitial;
        colorLabel = "Density (kg/m^3)";

    case "velocity"

        field0 = sqrt( ...
            uInitial.^2 + ...
            vInitial.^2);

        colorLabel = ...
            "Velocity magnitude (m/s)";

    otherwise

        error( ...
            "Unknown plotVariable. Use pressure, density or velocity.");

end

%% =========================================================
%  COLOUR LIMITS
%
% Use only INITIAL fluid field.
%
% This prevents later numerical pressure spikes from
% destroying the colour contrast.
% ==========================================================

initialFluidField = field0(fluid);

initialFluidField = ...
    initialFluidField( ...
    isfinite(initialFluidField));

if isempty(initialFluidField)
    error("Initial field contains no finite values.");
end

% fieldMin = min(initialFluidField);
% fieldMax = max(initialFluidField);

rho0_ref = 1000.0;
g = 9.81;

H = fluidYmax;

switch plotVariable

    case "pressure"

        fieldMin = 0.0;
        fieldMax = rho0_ref * g * H;

    case "density"

        initialFluidField = field0(fluid);
        initialFluidField = initialFluidField(isfinite(initialFluidField));

        fieldMin = min(initialFluidField);
        fieldMax = max(initialFluidField);

        if fieldMin == fieldMax
            fieldMin = fieldMin - 1;
            fieldMax = fieldMax + 1;
        end

    case "velocity"

        fieldMin = 0.0;

        % choose a useful fixed upper limit
        fieldMax = sqrt(g * H);

end

fprintf("\nColour range\n");
fprintf("----------------------\n");
fprintf("minimum = %.4f\n", fieldMin);
fprintf("maximum = %.4f\n", fieldMax);

%% Avoid zero colour range

if fieldMin == fieldMax
    fieldMin = fieldMin - 1;
    fieldMax = fieldMax + 1;
end

fprintf("\n");
fprintf("Initial colour minimum = %g\n", ...
    fieldMin);

fprintf("Initial colour maximum = %g\n", ...
    fieldMax);

%% Additional pressure diagnostic

if plotVariable == "pressure"

    fprintf("\n");
    fprintf( ...
        "Initial minimum fluid pressure = %.4f Pa\n", ...
        min(pressureInitial(fluid)));

    fprintf( ...
        "Initial maximum fluid pressure = %.4f Pa\n", ...
        max(pressureInitial(fluid)));

end

%% =========================================================
%  IDEAL HYDROSTATIC PRESSURE FIGURE
% ==========================================================

if plotVariable == "pressure"

    % Fluid height
    H = fluidYmax;

    % Ideal hydrostatic pressure at each initial fluid particle
    pressureHydrostatic = ...
        rho0_ref * g * (H - y0(fluid));

    % Do not allow negative pressure above free surface
    pressureHydrostatic = ...
        max(pressureHydrostatic, 0.0);

    % Create ideal hydrostatic pressure figure
    figHydro = figure( ...
        'Color', 'w', ...
        'Position', [1050 100 900 700]);

    scatter( ...
        x0(fluid), ...
        y0(fluid), ...
        25, ...
        pressureHydrostatic, ...
        'filled', ...
        'MarkerEdgeColor', 'none');

    hold on;

    scatter( ...
        x0(boundary), ...
        y0(boundary), ...
        90, ...
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

    set(gca, ...
        'FontSize', 12, ...
        'LineWidth', 1);

end

%% =========================================================
%  CREATE RED COLORMAP
%
% Low value  = light red
% High value = dark red
% ==========================================================

nColor = 256;

lightRed = [1.00, 0.80, 0.80];
darkRed  = [0.30, 0.00, 0.00];

r = linspace( ...
    lightRed(1), ...
    darkRed(1), ...
    nColor)';

g = linspace( ...
    lightRed(2), ...
    darkRed(2), ...
    nColor)';

b = linspace( ...
    lightRed(3), ...
    darkRed(3), ...
    nColor)';

redMap = [r g b];

%% =========================================================
%  CREATE FIGURE
% ==========================================================

fig = figure( ...
    'Color', 'w', ...
    'Position', [100 100 900 700]);

hold on;

%% =========================================================
%  FLUID PARTICLES
% ==========================================================

hFluid = scatter( ...
    x0(fluid), ...
    y0(fluid), ...
    25, ...
    field0(fluid), ...
    'filled', ...
    'MarkerEdgeColor','none');

%% =========================================================
%  BOUNDARY PARTICLES
% ==========================================================

hBoundary = scatter( ...
    x0(boundary), ...
    y0(boundary), ...
    90, ...
    [0 0 0], ...
    'filled');

%% =========================================================
%  FIGURE FORMATTING
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
c.Label.String = colorLabel;

clim([fieldMin fieldMax]);

title(sprintf( ...
    'SPH Still Water, t = %.3f s', ...
    time(1)));

set(gca, ...
    'FontSize', 12, ...
    'LineWidth', 1);

%% =========================================================
%  OPTIONAL VIDEO
% ==========================================================

if saveVideo

    video = VideoWriter( ...
        videoName, ...
        'MPEG-4');

    video.FrameRate = 5;

    open(video);

end

%% =========================================================
%  ANIMATION LOOP
% ==========================================================

for n = 1:length(files)

    %% -----------------------------------------------------
    % Read timestep
    % ------------------------------------------------------

    filename = fullfile( ...
        files(n).folder, ...
        files(n).name);

    data = readmatrix( ...
        filename, ...
        'NumHeaderLines', 7);

    %% -----------------------------------------------------
    % Read variables
    % ------------------------------------------------------

    x = data(:,1);
    y = data(:,2);

    rho = data(:,4);
    pressure = data(:,6);

    u = data(:,8);
    v = data(:,9);

    %% -----------------------------------------------------
    % Choose plotted field
    % ------------------------------------------------------

    switch plotVariable

        case "pressure"

            field = pressure;

        case "density"

            field = rho;

        case "velocity"

            field = sqrt( ...
                u.^2 + ...
                v.^2);

    end

    %% -----------------------------------------------------
    % Check NaN and Inf
    % ------------------------------------------------------

    valid = ...
        isfinite(x) & ...
        isfinite(y) & ...
        isfinite(field);

    invalid = ~valid;

    if any(invalid)

        fprintf( ...
            "Warning: %d invalid particles at t = %.4f s\n", ...
            sum(invalid), ...
            time(n));

    end

    %% Hide invalid particles

    xplot = x;
    yplot = y;
    fieldplot = field;

    xplot(~valid) = NaN;
    yplot(~valid) = NaN;
    fieldplot(~valid) = NaN;

    %% -----------------------------------------------------
    % Update fluid particles
    % ------------------------------------------------------

    set( ...
        hFluid, ...
        'XData', xplot(fluid), ...
        'YData', yplot(fluid), ...
        'CData', fieldplot(fluid));

    %% -----------------------------------------------------
    % Update boundary particle positions
    % ------------------------------------------------------

    set( ...
        hBoundary, ...
        'XData', xplot(boundary), ...
        'YData', yplot(boundary));

    %% -----------------------------------------------------
    % Update title
    % ------------------------------------------------------

    title(sprintf( ...
        'SPH Still Water, t = %.3f s', ...
        time(n)));

    drawnow;
    pause(0.2)

    %% -----------------------------------------------------
    % Save video frame
    % ------------------------------------------------------

    if saveVideo

        frame = getframe(fig);

        writeVideo( ...
            video, ...
            frame);

    end

end

%% =========================================================
%  CLOSE VIDEO
% ==========================================================

if saveVideo

    close(video);

    fprintf("\n");
    fprintf( ...
        "Video saved as %s\n", ...
        videoName);



end