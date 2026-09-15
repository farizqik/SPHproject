clear;
clc;
close all;

%% Folder

simulationFolder = ...
    "EWCSPHCylinder_dr0_0.001414_Nr_20_Ntheta_222.144147_Rr0_40.000000_wendland";

filename = fullfile(simulationFolder, "particles.csv");

if ~isfile(filename)
    error("File not found: %s", filename);
end


%% Read particle data

opts = delimitedTextImportOptions("NumVariables",6);

opts.DataLines = [14 Inf];
opts.Delimiter = ",";

opts.VariableNames = ...
    ["ID","r","theta","x","y","Type"];

opts.VariableTypes = ...
    ["double","double","double","double","double","string"];

T = readtable(filename, opts);


%% Separate boundary and fluid

boundary = lower(T.Type) == "boundary";
fluid    = lower(T.Type) == "fluid";


%% Plot

figure;

plot(T.x(fluid), T.y(fluid), '.', ...
    'MarkerSize',10);

hold on;

plot(T.x(boundary), T.y(boundary), 'o', ...
    'MarkerSize',12);

axis equal;
grid on;
box on;

xlabel('x (m)');
ylabel('y (m)');

legend('Fluid','Cylinder boundary');

title('Polar Particle Distribution');