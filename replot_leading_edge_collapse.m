%% Replot digitized data from Fig. 11: motion of a collapsing water-column front
% Data were manually digitized from the supplied raster image and are
% approximate. x = t*sqrt(2g/L), y = z/L.

clear; close all; clc;

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

xlim([0 3.25]); ylim([1 4]);
xticks(0:0.5:3.0); yticks(1:0.5:4.0);
xlabel('$t\sqrt{2g/L}$','Interpreter','latex','FontSize',13);
ylabel('$z/L$','Interpreter','latex','FontSize',13);
grid on; box on;
set(gca,'FontName','Times New Roman','FontSize',11, ...
    'GridColor',[0.72 0.72 0.72],'GridAlpha',0.55,'LineWidth',0.8);
axis square;

legend([h1 h2 h3 h4 h5], ...
    {'Exp. (Martin & Moyce, 1.125)', ...
     'Exp. (Martin & Moyce, 2.25)', ...
     'Exp. (Koshizuka et al.)', ...
     'Cal. (SOLA-VOF)', ...
     'Cal. (MPS)'}, ...
    'Location','northwest','Box','off','FontSize',9);

title('Motion of the leading edge in collapse of a water column', ...
    'FontWeight','normal','FontSize',12);

% Optional export:
% exportgraphics(gcf,'leading_edge_collapse_replot.png','Resolution',300);
