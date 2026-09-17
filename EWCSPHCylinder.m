%% Post-process EWCSPH cylinder consistency tests

clear;
clc;
close all;


%% =========================================================
% USER SETTINGS
% ==========================================================

simulationFolder = ...
    "EWCSPHCylinder_dr0_0.001414_Nr_100_Ntheta_222_Rr0_40.000000";

filePattern = "particles_hdp_*.csv";

outputFolder = fullfile( ...
    simulationFolder, ...
    "MATLAB_plots");

saveFigures = true;

% Color scale percentile
% 95 = stronger clipping
% 98 = recommended
% 99 = weaker clipping
colorPercentile = 98;


if saveFigures && ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end


%% =========================================================
% FIND CSV FILES
% ==========================================================

files = dir( ...
    fullfile(simulationFolder,filePattern));

if isempty(files)

    error( ...
        'No files matching "%s" found.', ...
        filePattern);

end


nFiles = numel(files);


%% =========================================================
% STORAGE FOR L2 NORMS
% ==========================================================

hcoefValues = zeros(nFiles,1);

% Partition of unity
L2PGauss = zeros(nFiles,1);
L2PCubic = zeros(nFiles,1);
L2PWedn  = zeros(nFiles,1);

% Gradient
L2GradGauss = zeros(nFiles,1);
L2GradCubic = zeros(nFiles,1);
L2GradWedn  = zeros(nFiles,1);

% Continuity
L2RhoGauss = zeros(nFiles,1);
L2RhoCubic = zeros(nFiles,1);
L2RhoWedn  = zeros(nFiles,1);


fprintf("\nFound %d CSV files.\n",nFiles);


%% =========================================================
% CSV IMPORT SETTINGS
%
% Line 17 onward = particle data
% ==========================================================

opts = delimitedTextImportOptions("NumVariables",23);

opts.DataLines = [17 Inf];
opts.Delimiter = ",";

opts.VariableNames = ...
    ["ID", ...
     "r", ...
     "theta", ...
     "x", ...
     "y", ...
     "Type", ...
     "empty1", ...
     "PGauss", ...
     "PCubic", ...
     "PWedn", ...
     "empty2", ...
     "dPGaussX", ...
     "dPGaussY", ...
     "empty3", ...
     "dPCubicX", ...
     "dPCubicY", ...
     "empty4", ...
     "dPWednX", ...
     "dPWednY", ...
     "empty5", ...
     "drhodtGauss", ...
     "drhodtCubic", ...
     "drhodtWedn"];

opts.VariableTypes = ...
    ["double", ...
     "double", ...
     "double", ...
     "double", ...
     "double", ...
     "string", ...
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
     "string", ...
     "double", ...
     "double", ...
     "string", ...
     "double", ...
     "double", ...
     "double"];


%% =========================================================
% LOOP THROUGH h COEFFICIENTS
% ==========================================================

for n = 1:nFiles

    filename = fullfile( ...
        files(n).folder, ...
        files(n).name);

    fprintf("Reading %s\n",files(n).name);


    %% -----------------------------------------------------
    % Extract h coefficient
    % ------------------------------------------------------

    token = regexp( ...
        files(n).name, ...
        '_hdp_([0-9.]+)\.csv$', ...
        'tokens', ...
        'once');

    hcoef = str2double(token{1});

    hcoefValues(n) = hcoef;


    %% -----------------------------------------------------
    % Read C++ L2 data
    % ------------------------------------------------------

    L2data = readmatrix( ...
        filename, ...
        'Range','A14:F14');

    L2PGauss(n) = L2data(1);
    L2PCubic(n) = L2data(2);
    L2PWedn(n)  = L2data(3);

    L2RhoGauss(n) = L2data(4);
    L2RhoCubic(n) = L2data(5);
    L2RhoWedn(n)  = L2data(6);


    %% -----------------------------------------------------
    % Read particle table
    % ------------------------------------------------------

    T = readtable(filename,opts);


    x = T.x;
    y = T.y;
    r = T.r;

    type = lower(strtrim(string(T.Type)));

    boundary = type == "boundary";
    fluid = type == "fluid";

    %% =========================================================
    % PARTICLE DISTRIBUTION
    % Plot only once because particle positions do not change
    % with h coefficient
    % ==========================================================

    if n == 1

        figParticles = figure( ...
            'Color','w', ...
            'Position',[200 100 750 700]);

        % Fluid particles
        plot( ...
            x(fluid), ...
            y(fluid), ...
            '.', ...
            'MarkerSize',3);

        hold on;

        % Cylinder boundary particles
        plot( ...
            x(boundary), ...
            y(boundary), ...
            'o', ...
            'MarkerSize',4, ...
            'LineWidth',1);

        axis equal;
        grid on;
        box on;

        xlabel('x (m)');
        ylabel('y (m)');

        title('Polar Particle Distribution');

        legend( ...
            'Fluid', ...
            'Cylinder boundary', ...
            'Location','best');


        %% Save particle plot

        if saveFigures

            exportgraphics( ...
                figParticles, ...
                fullfile( ...
                outputFolder, ...
                'ParticleDistribution.png'), ...
                'Resolution',300);

        end

    end


    %% -----------------------------------------------------
    % Partition of unity
    % ------------------------------------------------------

    PGauss = T.PGauss;
    PCubic = T.PCubic;
    PWedn  = T.PWedn;


    %% -----------------------------------------------------
    % Derivatives
    % ------------------------------------------------------

    dPGaussX = T.dPGaussX;
    dPGaussY = T.dPGaussY;

    dPCubicX = T.dPCubicX;
    dPCubicY = T.dPCubicY;

    dPWednX = T.dPWednX;
    dPWednY = T.dPWednY;


    %% Gradient magnitude

    gradGauss = sqrt( ...
        dPGaussX.^2 + ...
        dPGaussY.^2);

    gradCubic = sqrt( ...
        dPCubicX.^2 + ...
        dPCubicY.^2);

    gradWedn = sqrt( ...
        dPWednX.^2 + ...
        dPWednY.^2);


    %% L2 gradient error

    L2GradGauss(n) = ...
        sqrt(mean(gradGauss.^2));

    L2GradCubic(n) = ...
        sqrt(mean(gradCubic.^2));

    L2GradWedn(n) = ...
        sqrt(mean(gradWedn.^2));


    %% -----------------------------------------------------
    % Continuity
    % ------------------------------------------------------

    drhodtGauss = T.drhodtGauss;
    drhodtCubic = T.drhodtCubic;
    drhodtWedn  = T.drhodtWedn;


    %% =====================================================
    % RECONSTRUCT POLAR GRID
    % ======================================================

    rUnique = unique(r);

    Nrings = numel(rUnique);

    Ntheta = sum(boundary);


    if Nrings*Ntheta ~= height(T)

        error( ...
            "Particle grid cannot be reshaped.");

    end


    X = reshape(x,Ntheta,Nrings);
    Y = reshape(y,Ntheta,Nrings);


    %% Partition

    PGaussGrid = reshape( ...
        PGauss,Ntheta,Nrings);

    PCubicGrid = reshape( ...
        PCubic,Ntheta,Nrings);

    PWednGrid = reshape( ...
        PWedn,Ntheta,Nrings);


    %% Gradient

    gradGaussGrid = reshape( ...
        gradGauss,Ntheta,Nrings);

    gradCubicGrid = reshape( ...
        gradCubic,Ntheta,Nrings);

    gradWednGrid = reshape( ...
        gradWedn,Ntheta,Nrings);


    %% Continuity

    rhoGaussGrid = reshape( ...
        drhodtGauss,Ntheta,Nrings);

    rhoCubicGrid = reshape( ...
        drhodtCubic,Ntheta,Nrings);

    rhoWednGrid = reshape( ...
        drhodtWedn,Ntheta,Nrings);


    %% -----------------------------------------------------
    % Close angular seam
    % ------------------------------------------------------

    X = [X;X(1,:)];
    Y = [Y;Y(1,:)];

    PGaussGrid = [PGaussGrid;PGaussGrid(1,:)];
    PCubicGrid = [PCubicGrid;PCubicGrid(1,:)];
    PWednGrid  = [PWednGrid;PWednGrid(1,:)];

    gradGaussGrid = [gradGaussGrid;gradGaussGrid(1,:)];
    gradCubicGrid = [gradCubicGrid;gradCubicGrid(1,:)];
    gradWednGrid  = [gradWednGrid;gradWednGrid(1,:)];

    rhoGaussGrid = [rhoGaussGrid;rhoGaussGrid(1,:)];
    rhoCubicGrid = [rhoCubicGrid;rhoCubicGrid(1,:)];
    rhoWednGrid  = [rhoWednGrid;rhoWednGrid(1,:)];


    %% =====================================================
    % FIGURE 1: PARTITION OF UNITY
    % ======================================================

    figP = figure( ...
        'Color','w', ...
        'Position',[50 100 1500 470]);

    tiledlayout(1,3, ...
        'TileSpacing','compact', ...
        'Padding','compact');


    % ------------------------------------------------------
    % Percentile color scale centered around exact P = 1
    % ------------------------------------------------------

    PAll = [PGauss;PCubic;PWedn];

    Perror = abs(PAll - 1.0);

    Perror = sort(Perror);

    indexP = ceil( ...
        colorPercentile/100 * numel(Perror));

    indexP = max(1,min(indexP,numel(Perror)));

    deltaP = Perror(indexP);

    if deltaP == 0
        deltaP = 0.01;
    end

    pLimits = [1-deltaP,1+deltaP];


    %% Gaussian

    nexttile;

    contourf( ...
        X,Y,PGaussGrid,30, ...
        'LineColor','none');

    addContourDataTips(gca,X,Y,PGaussGrid);

    axis equal tight;
    clim(pLimits);
    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');

    title('P_{Gaussian}');


    %% Cubic

    nexttile;

    contourf( ...
        X,Y,PCubicGrid,30, ...
        'LineColor','none');

    addContourDataTips(gca,X,Y,PCubicGrid);

    axis equal tight;
    clim(pLimits);
    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');

    title('P_{Cubic}');


    %% Wendland

    nexttile;

    contourf( ...
        X,Y,PWednGrid,30, ...
        'LineColor','none');

    addContourDataTips(gca,X,Y,PWednGrid);

    axis equal tight;
    clim(pLimits);
    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');

    title('P_{Wendland}');


    sgtitle(sprintf( ...
        'Partition of Unity, h coefficient = %.2f', ...
        hcoef));


    %% =====================================================
    % FIGURE 2: GRADIENT MAGNITUDE
    %
    % Exact = 0
    % ======================================================

    figGrad = figure( ...
        'Color','w', ...
        'Position',[50 100 1500 470]);

    tiledlayout(1,3, ...
        'TileSpacing','compact', ...
        'Padding','compact');


    % ------------------------------------------------------
    % Percentile gradient limit
    % ------------------------------------------------------

    gradAll = [ ...
        gradGauss; ...
        gradCubic; ...
        gradWedn];

    gradAll = sort(gradAll);

    indexGrad = ceil( ...
        colorPercentile/100 * numel(gradAll));

    indexGrad = ...
        max(1,min(indexGrad,numel(gradAll)));

    gradMax = gradAll(indexGrad);

    if gradMax == 0
        gradMax = 1;
    end


    %% Gaussian

    nexttile;

    contourf( ...
        X,Y,gradGaussGrid,30, ...
        'LineColor','none');

    addContourDataTips(gca,X,Y,gradGaussGrid);

    axis equal tight;

    clim([0 gradMax]);

    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');

    title('|∇P|_{Gaussian}');


    %% Cubic

    nexttile;

    contourf( ...
        X,Y,gradCubicGrid,30, ...
        'LineColor','none');

    addContourDataTips(gca,X,Y,gradCubicGrid);

    axis equal tight;

    clim([0 gradMax]);

    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');

    title('|∇P|_{Cubic}');


    %% Wendland

    nexttile;

    contourf( ...
        X,Y,gradWednGrid,30, ...
        'LineColor','none');

    addContourDataTips(gca,X,Y,gradWednGrid);

    axis equal tight;

    clim([0 gradMax]);

    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');

    title('|∇P|_{Wendland}');


    sgtitle(sprintf( ...
        'Kernel Gradient Error, h coefficient = %.2f', ...
        hcoef));


    %% =====================================================
    % FIGURE 3: CONTINUITY
    % ======================================================

    figRho = figure( ...
        'Color','w', ...
        'Position',[50 100 1500 470]);

    tiledlayout(1,3, ...
        'TileSpacing','compact', ...
        'Padding','compact');


    % ------------------------------------------------------
    % Symmetric percentile color scale
    % ------------------------------------------------------

    rhoAll = abs([ ...
        drhodtGauss; ...
        drhodtCubic; ...
        drhodtWedn]);

    rhoAll = sort(rhoAll);

    indexRho = ceil( ...
        colorPercentile/100 * numel(rhoAll));

    indexRho = ...
        max(1,min(indexRho,numel(rhoAll)));

    rhoMax = rhoAll(indexRho);

    if rhoMax == 0
        rhoMax = 1;
    end

    rhoLimits = [-rhoMax rhoMax];


    %% Gaussian

    nexttile;

    contourf( ...
        X,Y,rhoGaussGrid,30, ...
        'LineColor','none');

    addContourDataTips(gca,X,Y,rhoGaussGrid);

    axis equal tight;

    clim(rhoLimits);

    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');

    title('(d\rho/dt)_{Gaussian}');


    %% Cubic

    nexttile;

    contourf( ...
        X,Y,rhoCubicGrid,30, ...
        'LineColor','none');

    addContourDataTips(gca,X,Y,rhoCubicGrid);

    axis equal tight;

    clim(rhoLimits);

    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');

    title('(d\rho/dt)_{Cubic}');


    %% Wendland

    nexttile;

    contourf( ...
        X,Y,rhoWednGrid,30, ...
        'LineColor','none');

    addContourDataTips(gca,X,Y,rhoWednGrid);

    axis equal tight;

    clim(rhoLimits);

    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');

    title('(d\rho/dt)_{Wendland}');


    sgtitle(sprintf( ...
        'Continuity Equation, h coefficient = %.2f', ...
        hcoef));


    %% =====================================================
    % SAVE SPATIAL FIGURES
    % ======================================================

    hText = strrep( ...
        sprintf('%.2f',hcoef), ...
        '.', ...
        'p');


    if saveFigures

        exportgraphics( ...
            figP, ...
            fullfile( ...
                outputFolder, ...
                "PartitionUnity_" + hText + ".png"), ...
            'Resolution',300);


        exportgraphics( ...
            figGrad, ...
            fullfile( ...
                outputFolder, ...
                "Gradient_" + hText + ".png"), ...
            'Resolution',300);


        exportgraphics( ...
            figRho, ...
            fullfile( ...
                outputFolder, ...
                "Continuity_" + hText + ".png"), ...
            'Resolution',300);

    end

end


%% =========================================================
% SORT BY h COEFFICIENT
% ==========================================================

[hcoefValues,order] = sort(hcoefValues);

L2PGauss = L2PGauss(order);
L2PCubic = L2PCubic(order);
L2PWedn  = L2PWedn(order);

L2GradGauss = L2GradGauss(order);
L2GradCubic = L2GradCubic(order);
L2GradWedn  = L2GradWedn(order);

L2RhoGauss = L2RhoGauss(order);
L2RhoCubic = L2RhoCubic(order);
L2RhoWedn  = L2RhoWedn(order);


%% =========================================================
% FIGURE 4: PARTITION L2
% ==========================================================

figL2P = figure( ...
    'Color','w', ...
    'Position',[200 150 850 600]);

plot( ...
    hcoefValues,L2PGauss, ...
    '-o','LineWidth',1.5,'MarkerSize',7);

hold on;

plot( ...
    hcoefValues,L2PCubic, ...
    '-s','LineWidth',1.5,'MarkerSize',7);

plot( ...
    hcoefValues,L2PWedn, ...
    '-^','LineWidth',1.5,'MarkerSize',7);

grid on;
box on;

xlabel('h coefficient');
ylabel('L_2 norm');

title('Partition of Unity L_2 Error');

legend( ...
    'Gaussian', ...
    'Cubic spline', ...
    'Wendland C2', ...
    'Location','best');


%% =========================================================
% FIGURE 5: GRADIENT L2
% ==========================================================

figL2Grad = figure( ...
    'Color','w', ...
    'Position',[200 150 850 600]);

plot( ...
    hcoefValues,L2GradGauss, ...
    '-o','LineWidth',1.5,'MarkerSize',7);

hold on;

plot( ...
    hcoefValues,L2GradCubic, ...
    '-s','LineWidth',1.5,'MarkerSize',7);

plot( ...
    hcoefValues,L2GradWedn, ...
    '-^','LineWidth',1.5,'MarkerSize',7);

grid on;
box on;

xlabel('h coefficient');
ylabel('L_2 norm of |∇P|');

title('Kernel Gradient L_2 Error');

legend( ...
    'Gaussian', ...
    'Cubic spline', ...
    'Wendland C2', ...
    'Location','best');


%% =========================================================
% FIGURE 6: CONTINUITY L2
% ==========================================================

figL2Rho = figure( ...
    'Color','w', ...
    'Position',[200 150 850 600]);

plot( ...
    hcoefValues,L2RhoGauss, ...
    '-o','LineWidth',1.5,'MarkerSize',7);

hold on;

plot( ...
    hcoefValues,L2RhoCubic, ...
    '-s','LineWidth',1.5,'MarkerSize',7);

plot( ...
    hcoefValues,L2RhoWedn, ...
    '-^','LineWidth',1.5,'MarkerSize',7);

grid on;
box on;

xlabel('h coefficient');

ylabel( ...
    'L_2 norm of d\rho/dt');

title( ...
    'Continuity Equation L_2 Error');

legend( ...
    'Gaussian', ...
    'Cubic spline', ...
    'Wendland C2', ...
    'Location','best');


%% =========================================================
% SAVE L2 FIGURES
% ==========================================================

if saveFigures

    exportgraphics( ...
        figL2P, ...
        fullfile( ...
            outputFolder, ...
            'L2_PartitionUnity.png'), ...
        'Resolution',300);

    exportgraphics( ...
        figL2Grad, ...
        fullfile( ...
            outputFolder, ...
            'L2_Gradient.png'), ...
        'Resolution',300);

    exportgraphics( ...
        figL2Rho, ...
        fullfile( ...
            outputFolder, ...
            'L2_Continuity.png'), ...
        'Resolution',300);

end


%% =========================================================
% SUMMARY TABLE
% ==========================================================

summaryTable = table( ...
    hcoefValues, ...
    L2PGauss, ...
    L2PCubic, ...
    L2PWedn, ...
    L2GradGauss, ...
    L2GradCubic, ...
    L2GradWedn, ...
    L2RhoGauss, ...
    L2RhoCubic, ...
    L2RhoWedn, ...
    'VariableNames', ...
    {'h_coefficient', ...
     'L2_P_Gaussian', ...
     'L2_P_Cubic', ...
     'L2_P_Wendland', ...
     'L2_Grad_Gaussian', ...
     'L2_Grad_Cubic', ...
     'L2_Grad_Wendland', ...
     'L2_DrhoDt_Gaussian', ...
     'L2_DrhoDt_Cubic', ...
     'L2_DrhoDt_Wendland'});


disp(summaryTable);


if saveFigures

    writetable( ...
        summaryTable, ...
        fullfile( ...
            outputFolder, ...
            'L2_summary.csv'));

end


fprintf("\nAll plots created.\n");
fprintf("Saved in:\n%s\n",outputFolder);

%% Continuous value inspection on the filled contour
function addContourDataTips(ax,X,Y,field)
    % Remove the duplicated angular seam before building the interpolant.
    X = X(1:end-1,:);
    Y = Y(1:end-1,:);
    field = field(1:end-1,:);
    valid = isfinite(X) & isfinite(Y) & isfinite(field);
    xy = [X(valid),Y(valid)];
    values = field(valid);
    [xy,indices] = unique(xy,'rows');
    values = values(indices);
    F = scatteredInterpolant(xy(:,1),xy(:,2),values,'linear','none');
    radii = hypot(X(:),Y(:));
    data = struct('F',F,'rmin',min(radii),'rmax',max(radii));
    contours = findobj(ax,'Type','contour');
    set(contours,'UserData',data,'HitTest','on','PickableParts','all', ...
        'ButtonDownFcn',@showContinuousValue);

    % Keep scroll zoom and pan. A plain click runs the contour callback.
    ax.Interactions = [zoomInteraction panInteraction];
    axtoolbar(ax,{'zoomin','zoomout','pan','restoreview'});
end

function showContinuousValue(contourHandle,~)
    ax = ancestor(contourHandle,'axes');
    fig = ancestor(ax,'figure');
    if ~strcmp(get(fig,'SelectionType'),'normal')
        return;
    end
    position = get(ax,'CurrentPoint');
    xq = position(1,1);
    yq = position(1,2);
    data = get(contourHandle,'UserData');
    rq = hypot(xq,yq);
    if rq < data.rmin || rq > data.rmax
        return;
    end
    value = data.F(xq,yq);
    if ~isfinite(value)
        return;
    end

    % The anchor is not drawn: there are no particle dots on the contour.
    existingTip = findall(ax,'Tag','ContinuousValueAnchor');

    if ~isempty(existingTip)
        delete(existingTip);
        delete(findall(ax,'Tag','ContinuousValueText'));
        return;
    end
    anchor = line(ax,xq,yq,'LineStyle','none','Marker','none', ...
        'HitTest','off','PickableParts','none', ...
        'Tag','ContinuousValueAnchor');
    if isprop(anchor,'DataTipTemplate')
        anchor.DataTipTemplate.DataTipRows = [ ...
            dataTipTextRow('x (m)','XData','%.8g'); ...
            dataTipTextRow('y (m)','YData','%.8g'); ...
            dataTipTextRow('Interpolated value',value,'%.10g')];
        anchor.DataTipTemplate.Interpreter = 'none';
        datatip(anchor,xq,yq);
    else
        % Compatibility fallback for releases without line data templates.
        delete(findall(ax,'Tag','ContinuousValueText'));
        text(ax,xq,yq,sprintf('x = %.8g\ny = %.8g\nInterpolated value = %.10g', ...
            xq,yq,value),'BackgroundColor','white','EdgeColor','black', ...
            'Margin',4,'VerticalAlignment','bottom','Interpreter','none', ...
            'HitTest','off','PickableParts','none','Tag','ContinuousValueText');
    end
end
