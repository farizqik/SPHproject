%% Post-process polar SPH partition-of-unity and continuity output
% Compatible with the CSV layout written by Pasted code (2)(2).cpp.

clear;
clc;
close all;

%% Plot every polar SPH CSV file in the current folder
% Expected filename form: 2DPartitionpolar_hdp_1.200000.csv
%
% The script creates:
% 1. One r-z figure for PGauss, PCubic, and PWedn for every h/dp
% 2. One r-z figure for all six gradient components for every h/dp
% 3. One r-z figure for all three continuity results for every h/dp
% 4. Two L2-error summary graphs versus h/dp
% 5. PNG files and a summary CSV saved in MATLAB_plots_polar

filePattern = '2DPartitionpolar_hdp_*.csv';
outputFolder = 'MATLAB_plots_polar';
saveFigures = true;

if saveFigures && ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

files = dir(filePattern);

if isempty(files)
    error('No files matching "%s" were found in the current folder.', ...
        filePattern);
end

nFiles = numel(files);

hOverDpValues = zeros(nFiles,1);
L2Gauss = zeros(nFiles,1);
L2Cubic = zeros(nFiles,1);
L2Wedn = zeros(nFiles,1);
L2normdrhodtGauss = zeros(nFiles,1);
L2normdrhodtCubic = zeros(nFiles,1);
L2normdrhodtWedn  = zeros(nFiles,1);


fprintf('Found %d CSV files.\n',nFiles);

for n = 1:nFiles

    filename = files(n).name;
    fprintf('Reading %s\n',filename);

    %% Read metadata
    % C++ layout:
    % row 1: h/dp; row 2: L2 headings; row 3: L2 values;
    % row 4: blank; row 5: data headings; row 6 onward: particle data.
    allCells = readcell(filename,'Delimiter',',');
    if size(allCells,1) < 6 || size(allCells,2) < 20
        error('File %s does not match the C++ polar CSV layout.',filename);
    end

    hOverDp = allCells{1,2};
    hOverDpValues(n) = hOverDp;

    % All six L2 norms are written by C++ on row 3.
    L2Gauss(n) = allCells{3,1};
    L2Cubic(n) = allCells{3,2};
    L2Wedn(n) = allCells{3,3};
    L2normdrhodtGauss(n) = allCells{3,4};
    L2normdrhodtCubic(n) = allCells{3,5};
    L2normdrhodtWedn(n)  = allCells{3,6};


    %% Read numerical data
    % Row 5 contains the variable names, so numerical data begin at row 6.
    data = readmatrix(filename, ...
        'Delimiter',',', ...
        'NumHeaderLines',5);

    % Remove empty rows and separator rows.
    data = data(~all(isnan(data),2),:);

    if size(data,2) < 20
        error(['File %s has only %d columns. The script expects the ' ...
               '20-column polar CSV format produced by the C++ code.'], ...
               filename,size(data,2));
    end

    %% Extract variables according to the C++ CSV layout
    id = data(:,1);
    r = data(:,2);
    z = data(:,3);

    PGauss = data(:,5);
    PCubic = data(:,6);
    PWedn = data(:,7);

    dPGaussR = data(:,9);
    dPGaussZ = data(:,10);

    dPCubicR = data(:,12);
    dPCubicZ = data(:,13);

    dPWednR = data(:,15);
    dPWednZ = data(:,16);

    drhodtGauss = data(:,18);
    drhodtCubic = data(:,19);
    drhodtWedn = data(:,20);

    %% Reconstruct the 2-D structured grid
    rUnique = unique(r);
    zUnique = unique(z);

    Nr = numel(rUnique);
    Nz = numel(zUnique);

    if Nr*Nz ~= numel(r)
        error(['The r-z coordinates in %s do not form a complete ' ...
               'structured grid.'],filename);
    end

    % C++ loops over z outside and r inside: r changes fastest.
    R = reshape(r,Nr,Nz);
    Z = reshape(z,Nr,Nz);

    PGaussGrid = reshape(PGauss,Nr,Nz);
    PCubicGrid = reshape(PCubic,Nr,Nz);
    PWednGrid = reshape(PWedn,Nr,Nz);

    dPGaussRGrid = reshape(dPGaussR,Nr,Nz);
    dPGaussZGrid = reshape(dPGaussZ,Nr,Nz);
    dPCubicRGrid = reshape(dPCubicR,Nr,Nz);
    dPCubicZGrid = reshape(dPCubicZ,Nr,Nz);
    dPWednRGrid = reshape(dPWednR,Nr,Nz);
    dPWednZGrid = reshape(dPWednZ,Nr,Nz);

    drhodtGaussGrid = reshape(drhodtGauss,Nr,Nz);
    drhodtCubicGrid = reshape(drhodtCubic,Nr,Nz);
    drhodtWednGrid = reshape(drhodtWedn,Nr,Nz);

    hdpText = strrep(sprintf('%.6f',hOverDp),'.','p');

    %% =====================================================
    % FIGURE 1: Partition of unity
    % ======================================================
    figP = figure( ...
        'Name',sprintf('Partition of unity, h/dp = %.6f',hOverDp), ...
        'Position',[100 100 1500 470]);

    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    pMin = min([PGauss;PCubic;PWedn]);
    pMax = max([PGauss;PCubic;PWedn]);

    nexttile;
    contourf(R,Z,PGaussGrid,30,'LineColor','none');
    axis equal tight;
    clim([pMin pMax]);
    colorbar;
    xlabel('r');
    ylabel('z');
    title('P_{Gauss}');

    nexttile;
    contourf(R,Z,PCubicGrid,30,'LineColor','none');
    axis equal tight;
    clim([pMin pMax]);
    colorbar;
    xlabel('r');
    ylabel('z');
    title('P_{Cubic}');

    nexttile;
    contourf(R,Z,PWednGrid,30,'LineColor','none');
    axis equal tight;
    clim([pMin pMax]);
    colorbar;
    xlabel('r');
    ylabel('z');
    title('P_{Wendland}');

    sgtitle(sprintf('Polar partition of unity, h/dp = %.6f',hOverDp));

    if saveFigures
        exportgraphics(figP, ...
            fullfile(outputFolder, ...
            ['PartitionUnity_hdp_' hdpText '.png']), ...
            'Resolution',300);
    end

    %% =====================================================
    % FIGURE 2: Gradient components
    % ======================================================
    figGrad = figure( ...
        'Name',sprintf('Kernel gradients, h/dp = %.6f',hOverDp), ...
        'Position',[100 50 800 900]);

    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    % Use one symmetric colour range for all gradient components.
    gradMaximum = max(abs([ ...
        dPGaussR;dPGaussZ; ...
        dPCubicR;dPCubicZ; ...
        dPWednR;dPWednZ]));

    if gradMaximum == 0
        gradMaximum = 1;
    end

    gradientLimits = [-gradMaximum gradMaximum];

    nexttile;
    contourf(R,Z,dPGaussRGrid,30,'LineColor','none');
    axis equal tight;
    clim(gradientLimits);
    colorbar;
    xlabel('r');
    ylabel('z');
    title('dP_{Gauss}/dr');

    nexttile;
    contourf(R,Z,dPGaussZGrid,30,'LineColor','none');
    axis equal tight;
    clim(gradientLimits);
    colorbar;
    xlabel('r');
    ylabel('z');
    title('dP_{Gauss}/dz');

    nexttile;
    contourf(R,Z,dPCubicRGrid,30,'LineColor','none');
    axis equal tight;
    clim(gradientLimits);
    colorbar;
    xlabel('r');
    ylabel('z');
    title('dP_{Cubic}/dr');

    nexttile;
    contourf(R,Z,dPCubicZGrid,30,'LineColor','none');
    axis equal tight;
    clim(gradientLimits);
    colorbar;
    xlabel('r');
    ylabel('z');
    title('dP_{Cubic}/dz');

    nexttile;
    contourf(R,Z,dPWednRGrid,30,'LineColor','none');
    axis equal tight;
    clim(gradientLimits);
    colorbar;
    xlabel('r');
    ylabel('z');
    title('dP_{Wendland}/dr');

    nexttile;
    contourf(R,Z,dPWednZGrid,30,'LineColor','none');
    axis equal tight;
    clim(gradientLimits);
    colorbar;
    xlabel('r');
    ylabel('z');
    title('dP_{Wendland}/dz');

    sgtitle(sprintf('Polar kernel-gradient sums, h/dp = %.6f',hOverDp));

    if saveFigures
        exportgraphics(figGrad, ...
            fullfile(outputFolder, ...
            ['GradientComponents_hdp_' hdpText '.png']), ...
            'Resolution',300);
    end

    %% =====================================================
    % FIGURE 3: Continuity drho/dt
    % ======================================================
    figRho = figure( ...
        'Name',sprintf('Continuity, h/dp = %.6f',hOverDp), ...
        'Position',[100 100 1500 470]);

    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    rhoMin = min([drhodtGauss;drhodtCubic;drhodtWedn]);
    rhoMax = max([drhodtGauss;drhodtCubic;drhodtWedn]);

    nexttile;
    contourf(R,Z,drhodtGaussGrid,30,'LineColor','none');
    axis equal tight;
    clim([rhoMin rhoMax]);
    colorbar;
    xlabel('r');
    ylabel('z');
    title('(d\rho/dt)_{Gauss}');

    nexttile;
    contourf(R,Z,drhodtCubicGrid,30,'LineColor','none');
    axis equal tight;
    clim([rhoMin rhoMax]);
    colorbar;
    xlabel('r');
    ylabel('z');
    title('(d\rho/dt)_{Cubic}');

    nexttile;
    contourf(R,Z,drhodtWednGrid,30,'LineColor','none');
    axis equal tight;
    clim([rhoMin rhoMax]);
    colorbar;
    xlabel('r');
    ylabel('z');
    title('(d\rho/dt)_{Wendland}');

    sgtitle(sprintf('Polar SPH continuity, h/dp = %.6f',hOverDp));

    if saveFigures
        exportgraphics(figRho, ...
            fullfile(outputFolder, ...
            ['Continuity_hdp_' hdpText '.png']), ...
            'Resolution',300);
    end
end

%% =========================================================
% FIGURE 4: Partition-of-unity L2 norm versus h/dp
% ==========================================================

[hOverDpValues,order] = sort(hOverDpValues);

L2Gauss = L2Gauss(order);
L2Cubic = L2Cubic(order);
L2Wedn = L2Wedn(order);
L2normdrhodtGauss = L2normdrhodtGauss(order);
L2normdrhodtCubic = L2normdrhodtCubic(order);
L2normdrhodtWedn  = L2normdrhodtWedn(order);

figL2 = figure( ...
    'Name','Partition-of-unity L2 norm versus h/dp', ...
    'Position',[200 150 900 600]);

plot(hOverDpValues,L2Gauss,'-o','LineWidth',1.5,'MarkerSize',7);
hold on;
plot(hOverDpValues,L2Cubic,'-s','LineWidth',1.5,'MarkerSize',7);
plot(hOverDpValues,L2Wedn,'-^','LineWidth',1.5,'MarkerSize',7);
hold off;

grid on;
box on;
xlabel('h/dp');
ylabel('L_2 norm');
title('Partition of unity L2 error versus h/dp');
legend('Gaussian','Cubic spline','Wendland C2', ...
    'Location','best');

if saveFigures
    exportgraphics(figL2, ...
        fullfile(outputFolder,'L2norm_vs_h_over_dp.png'), ...
        'Resolution',300);
end


%% =========================================================
% FIGURE 5: L2 norm of drho/dt versus h/dp
% ==========================================================

figL2rho = figure( ...
    'Name','L2 norm (d\rho/dt) versus h/dp', ...
    'Position',[200 150 900 600]);

plot(hOverDpValues,L2normdrhodtGauss,'-o','LineWidth',1.5,'MarkerSize',7);
hold on;
plot(hOverDpValues,L2normdrhodtCubic,'-s','LineWidth',1.5,'MarkerSize',7);
plot(hOverDpValues,L2normdrhodtWedn,'-^','LineWidth',1.5,'MarkerSize',7);
hold off;

grid on;
box on;
xlabel('h/dp');
ylabel('L_2 norm of (d\rho/dt)');
title('Polar continuity L2 error versus h/dp');
legend('Gaussian','Cubic spline','Wendland C2','Location','best');

if saveFigures
    exportgraphics(figL2rho, ...
        fullfile(outputFolder,'L2norm_drhodt_vs_h_over_dp.png'), ...
        'Resolution',300);
end


%% Save and display the summary values
summaryTable = table( ...
    hOverDpValues,L2Gauss,L2Cubic,L2Wedn, ...
    L2normdrhodtGauss,L2normdrhodtCubic,L2normdrhodtWedn, ...
    'VariableNames', ...
    {'h_over_dp', ...
     'L2_P_Gaussian','L2_P_Cubic','L2_P_Wendland', ...
     'L2_DrhoDt_Gaussian','L2_DrhoDt_Cubic','L2_DrhoDt_Wendland'});

disp(summaryTable);

if saveFigures
    writetable(summaryTable, ...
        fullfile(outputFolder,'L2norm_summary.csv'));
end

fprintf('\nAll plots have been created.\n');
fprintf('Results are saved in the folder: %s\n',outputFolder);



% 
% output_path = fullfile(pwd,'plot_all_SPH_2D.m');
% fid = fopen(output_path,'w','n','UTF-8');
% if fid == -1
%     error('Cannot open %s for writing.', output_path);
% end
% fprintf(fid, '%s', matlab_code);
% fclose(fid);
% fprintf('Created: %s\n', output_path);
% 
% fprintf('Created: %s\n', char(output_path));