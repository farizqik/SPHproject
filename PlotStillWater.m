clear;
clc;
close all;

%% =========================================================
%  USER SETTINGS
% ==========================================================

kernel = "wendland";
dxh = 0.5;

% Available:
% "pressure"
% "density"
% "velocity"

plotVariable = "pressure";

% Save animation?
saveVideo = true;

videoName = "SPH_still_water.mp4";

%% =========================================================
%  DOMAIN SETTINGS
% ==========================================================

dx = 0.5;
dy = 0.5;

Ncol = 50;
Nrow = 20;

Nboundary = 1;

% Fluid domain:
%
% C++ coordinates:
%
% x = (i - Nboundary)*dx
% y = (j - Nboundary)*dy
%
% Left, right and bottom:
% 3 boundary layers
%
% Top:
% free surface

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
    "Stillwater_dxh_%.6f_t_*_%s.csv", ...
    dxh, kernel);

files = dir(pattern);

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
        '_t_([0-9]+\.[0-9]+)_', ...
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

boundary = ...
    x0 < fluidXmin | ...
    x0 > fluidXmax | ...
    y0 < fluidYmin;

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

fieldMin = min(initialFluidField);
fieldMax = max(initialFluidField);

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

    video.FrameRate = 30;

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