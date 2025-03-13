sim_title = "THz-f50-20-5b1";
outer_folder = "Jeremy/two_freq/";
data_folder = outer_folder + sim_title + "/DATA/";
plots_folder = data_folder + "plots/";
% mkdir(plots_folder)

% Scaling for spectra
orders =                    1e-5;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Constants;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load from Input File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = data_folder + 'sup-input.json';
json_str = fileread(fname);
json_arr = jsondecode(json_str);


pulse{1}.lambda0 = json_arr.laser.pulses(1).lambda0;
pulse{1}.omg0 = 2 * pi * const.cl / pulse{1}.lambda0;

try
    pulse{2}.lambda0 = json_arr.laser.pulses(2).lambda0;
    pulse{2}.omg0 = 2 * pi * const.cl / pulse{2}.lambda0;
catch 
    plot_spectral_energy =      false;
end

clearvars fname json_str;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load plotting class and preferences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path to SupData.m class
addpath("/gpfs/fs1/home/aellio18/Desktop/SUPER/super/analysis")

% Plotting defaults
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Electric Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% E-field at each z step --------------------------------------------------
EF = SupData(data_folder + 'EF-re.h5',true);
    
mkdir(plots_folder + "EF/")
figure; 

clearvars iz;
size_EF = size(EF.data);

for iz = 1:size_EF(3)-1
    imagesc(EF.axes(2).data, EF.axes(1).data, ...
        abs(squeeze(EF.data(:,:,iz)))); axis xy; hold off;
    colorbar; colormap(parula); 
    set(gca, 'ColorScale', 'log');
    a = colorbar;
    a.Label.String = "Electric Field";
    xlabel(EF.axes(2).name + EF.axes(2).units); 
    ylabel(EF.axes(1).name + EF.axes(1).units);
    title("z = " + num2str(EF.axes(3).data(iz)) + EF.axes(3).units);
    filename = plots_folder + "EF/" + sim_title +...
        "_" + num2str(iz) + "_EF.jpg";
    
    saveas(gcf,filename);
    clearvars filename;
end

clearvars size_EF iz;
