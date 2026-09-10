%% Replot digitized data from Fig. 11: motion of a collapsing water-column front
% Data were manually digitized from the supplied raster image and are
% approximate. x = t*sqrt(2g/L), y = z/L.

clear; close all; clc;

%% Your WCSPH data
simulationFolder = ...
    "WCSPHpolar_dp_0.250000_h_0.500000_Nparticles_2928_wendland";

cartesianSimulationFolder = ...
    "WCSPH_dp_0.250000_h_0.500000_Nparticles_3264_wendland";

L = 8.0;                   % Initial water-column length (m)
tankLength = 4*L;          % Outer-wall position (m)
g = 9.81;                  % Gravitational acceleration (m/s^2)
frontLayerFactor = 1.0;    % Use the lowest 1.0*dp fluid layer

files = dir(fullfile(simulationFolder,"WCSPHpolar_hdp_*_t_*.csv"));

if isempty(files)
    error("No C++ timestep CSV files found in:\n%s",simulationFolder);
end

numberOfFrames = numel(files);
simulationTime = nan(numberOfFrames,1);

for k = 1:numberOfFrames
    token = regexp(files(k).name, ...
        '_t_([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\.csv$', ...
        'tokens','once');

    if isempty(token)
        error("Could not read time from %s",files(k).name);
    end

    simulationTime(k) = str2double(token{1});
end

[simulationTime,order] = sort(simulationTime);
files = files(order);

firstFilename = fullfile(files(1).folder,files(1).name);
metadata = readmatrix(firstFilename,'Range','A3:D3');
dp = metadata(3);

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

frontPosition = nan(numberOfFrames,1);

for k = 1:numberOfFrames
    filename = fullfile(files(k).folder,files(k).name);
    T = readtable(filename,opts);

    fluid = lower(strtrim(string(T.Type))) == "fluid";
    validFluid = fluid & isfinite(T.r) & isfinite(T.z);

    if ~any(validFluid)
        continue;
    end

    minimumFluidZ = min(T.z(validFluid));
    frontLayer = validFluid & ...
        T.z <= minimumFluidZ + frontLayerFactor*dp;

    if any(frontLayer)
        % Convert the outermost particle-centre position to the approximate
        % water-front interface position.
        frontPosition(k) = max(T.r(frontLayer)) + 0.5*dp;
    end
end

% The reference graph describes the freely advancing front. Stop using the
% numerical result once it reaches this model's outer wall and reflects.
firstWallHit = find(frontPosition >= tankLength - 0.5*dp,1,'first');
if ~isempty(firstWallHit) && firstWallHit < numberOfFrames
    frontPosition(firstWallHit+1:end) = NaN;
end

yourX = simulationTime*sqrt(2*g/L);
yourY = frontPosition/L;

% Keep only valid data within the published comparison range.
keep = isfinite(yourX) & isfinite(yourY) & yourX <= 3.25;
yourX = yourX(keep);
yourY = yourY(keep);

%% Cartesian WCSPH data
cartesianFiles = dir(fullfile( ...
    cartesianSimulationFolder,"WCSPH_hdp_*_t_*.csv"));

if isempty(cartesianFiles)
    error("No Cartesian C++ timestep CSV files found in:\n%s", ...
        cartesianSimulationFolder);
end

numberOfCartesianFrames = numel(cartesianFiles);
cartesianTime = nan(numberOfCartesianFrames,1);

for k = 1:numberOfCartesianFrames
    token = regexp(cartesianFiles(k).name, ...
        '_t_([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\.csv$', ...
        'tokens','once');

    if isempty(token)
        error("Could not read time from %s",cartesianFiles(k).name);
    end

    cartesianTime(k) = str2double(token{1});
end

[cartesianTime,order] = sort(cartesianTime);
cartesianFiles = cartesianFiles(order);

firstCartesianFilename = fullfile( ...
    cartesianFiles(1).folder,cartesianFiles(1).name);
cartesianMetadata = readmatrix(firstCartesianFilename,'Range','A3:D3');
cartesianDp = cartesianMetadata(3);

cartesianOpts = delimitedTextImportOptions("NumVariables",14);
cartesianOpts.DataLines = [7 Inf];
cartesianOpts.Delimiter = ",";
cartesianOpts.VariableNames = ...
    ["ID","x","y","empty1", ...
     "rho","drhodt","pressure","empty2", ...
     "u","v","empty3","dudt","dvdt","Type"];
cartesianOpts.VariableTypes = ...
    ["double","double","double","string", ...
     "double","double","double","string", ...
     "double","double","string","double","double","string"];

cartesianFrontPosition = nan(numberOfCartesianFrames,1);

for k = 1:numberOfCartesianFrames
    filename = fullfile(cartesianFiles(k).folder,cartesianFiles(k).name);
    T = readtable(filename,cartesianOpts);

    fluid = lower(strtrim(string(T.Type))) == "fluid";
    validFluid = fluid & isfinite(T.x) & isfinite(T.y);

    if ~any(validFluid)
        continue;
    end

    minimumFluidY = min(T.y(validFluid));
    frontLayer = validFluid & ...
        T.y <= minimumFluidY + frontLayerFactor*cartesianDp;

    if any(frontLayer)
        cartesianFrontPosition(k) = ...
            max(T.x(frontLayer)) + 0.5*cartesianDp;
    end
end

% Exclude the Cartesian solution after its first outer-wall impact.
firstCartesianWallHit = find( ...
    cartesianFrontPosition >= tankLength - 0.5*cartesianDp,1,'first');

if ~isempty(firstCartesianWallHit) && ...
        firstCartesianWallHit < numberOfCartesianFrames
    cartesianFrontPosition(firstCartesianWallHit+1:end) = NaN;
end

cartesianX = cartesianTime*sqrt(2*g/L);
cartesianY = cartesianFrontPosition/L;

keep = isfinite(cartesianX) & isfinite(cartesianY) & cartesianX <= 3.25;
cartesianX = cartesianX(keep);
cartesianY = cartesianY(keep);

% Experimental data: Martin & Moyce, initial aspect ratio 1.125
MM1125 = [ ...
    0.37 1.11
    0.78 1.23
    1.54 1.88
    1.92 2.32
    2.24 2.77
    2.60 3.21
    2.94 3.65];

% Experimental data: Martin & Moyce, initial aspect ratio 2.25
MM225 = [ ...
    0.39 1.11
    0.82 1.21
    1.18 1.42
    1.42 1.66
    1.59 1.88
    1.83 2.10
    1.97 2.32
    2.20 2.55
    2.30 2.77
    2.52 2.99
    2.66 3.21
    2.85 3.42
    2.98 3.65
    3.14 3.88];

% Experimental data: Koshizuka et al.
Koshizuka = [ ...
    0.00 1.00
    0.38 1.11
    0.75 1.25
    1.15 1.50
    1.52 1.88
    1.94 2.23
    2.33 2.59
    2.73 2.98
    3.12 3.60];

% Digitized calculated curves. Using explicit points avoids presenting a
% fitted equation as though it were reported by the original authors.
SOLA_VOF = [ ...
    0.00 1.00; 0.10 1.00; 0.25 1.04; 0.50 1.17; 0.75 1.33; 1.00 1.52;
    1.25 1.75; 1.50 2.00; 1.75 2.29; 2.00 2.61; 2.25 2.97;
    2.50 3.37; 2.70 3.72; 2.84 4.00];

MPS = [ ...
    0.00 1.00; 0.10 1.00; 0.25 1.03; 0.50 1.13; 0.75 1.29; 1.00 1.48;
    1.25 1.70; 1.50 1.96; 1.75 2.25; 2.00 2.57; 2.25 2.93;
    2.50 3.32; 2.70 3.68; 2.86 4.00];

% Smooth the digitized curve points without changing the experimental data.
xCurve = linspace(0,2.80,300);
ySOLA = interp1(SOLA_VOF(:,1),SOLA_VOF(:,2),xCurve,'pchip');
yMPS  = interp1(MPS(:,1),MPS(:,2),xCurve,'pchip');

figure('Color','w','Position',[200 120 570 650]);
hold on;

h1 = plot(MM1125(:,1),MM1125(:,2),'o','MarkerSize',5.5, ...
    'MarkerFaceColor',[0.20 0.20 0.20],'MarkerEdgeColor',[0.20 0.20 0.20], ...
    'LineStyle','none');
h2 = plot(MM225(:,1),MM225(:,2),'o','MarkerSize',5.5, ...
    'MarkerFaceColor','none','MarkerEdgeColor',[0.30 0.30 0.30], ...
    'LineStyle','none');
h3 = plot(Koshizuka(:,1),Koshizuka(:,2),'s','MarkerSize',5.5, ...
    'MarkerFaceColor','none','MarkerEdgeColor',[0.30 0.30 0.30], ...
    'LineStyle','none');
h4 = plot(xCurve,ySOLA,'--','Color',[0.30 0.30 0.30],'LineWidth',1.2);
h5 = plot(xCurve,yMPS,'-','Color',[0.10 0.10 0.10],'LineWidth',1.3);
h6 = plot(yourX,yourY,'-','Color',[0.85 0.10 0.10], ...
    'LineWidth',0.8,'Marker','o','MarkerSize',4, ...
    'MarkerIndices',1:max(1,round(numel(yourX)/20)):numel(yourX), ...
    'MarkerFaceColor','none');
markerInterval = max(1,round(numel(cartesianX)/20));

h7 = plot(cartesianX,cartesianY,'s', ...
    'Color',[0.10 0.35 0.85], ...
    'LineStyle','none', ...
    'MarkerSize',6, ...
    'MarkerIndices',1:markerInterval:numel(cartesianX), ...
    'MarkerFaceColor','none');
xlim([0 3.25]); ylim([1 4]);
xticks(0:0.5:3.0); yticks(1:0.5:4.0);
xlabel('$t\sqrt{2g/L}$','Interpreter','latex','FontSize',13);
ylabel('$z/L$','Interpreter','latex','FontSize',13);
grid on; box on;
set(gca,'FontName','Times New Roman','FontSize',11, ...
    'GridColor',[0.72 0.72 0.72],'GridAlpha',0.55,'LineWidth',0.8);
axis square;

legend([h1 h2 h3 h4 h5 h7 h6], ...
    {'Exp. (Martin & Moyce, 1.125)', ...
     'Exp. (Martin & Moyce, 2.25)', ...
     'Exp. (Koshizuka et al.)', ...
     'Cal. (SOLA-VOF)', ...
     'Cal. (MPS)', ...
     'Present Cartesian WCSPH', ...
     'Present Axisymmetric WCSPH'}, ...
    'Location','northwest','Box','off','FontSize',9);

title('Motion of the leading edge in collapse of a water column', ...
    'FontWeight','normal','FontSize',12);

% Optional export:
% exportgraphics(gcf,'leading_edge_collapse_replot.png','Resolution',300);
