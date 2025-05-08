%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Keys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outer_folder = "../../super-jax/output/";
outer_folder = "/pscratch/sd/a/aellio/WORKING/";
% sim_title = "Jeremy-f33_STP_sol_filt_Nt/";
% sim_title = "Tanner-f85-01b_gaus_streak2/";
sim_title = "Tanner-f10_ff_v10025c/";
data_folder = outer_folder + sim_title;
plots_folder = data_folder + "plots/";

saveplots = false;
brute = true;   % true: perform brute-force calc on emission angle; 
                % false: use interpolation

if ~exist(plots_folder, 'dir')
    mkdir(plots_folder)
end

filename = data_folder + "E_full.h5";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("/global/u1/a/aellio/super_analysis");

Constants;

% Plotting defaults
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',18);
set(groot, 'DefaultLineLineWidth', 1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load grid quantities.
grid.r = h5read(filename, "/r");
grid.t = h5read(filename, "/t");
grid.z_skip = h5read(filename, "/z");
grid.dt = grid.t(2) - grid.t(1);
grid.Nt = length(grid.t);
grid.Nr = length(grid.r);

% Load full Electric field and electron density (large arrays).
EF = h5read(filename, "/__xarray_dataarray_variable__"); % complex.
filename = data_folder + "ne_full.h5";
neF = h5read(filename, "/__xarray_dataarray_variable__");
clearvars filename;

% Load smaller arrays.
filename = data_folder + "data_line.h5";
grid.z = h5read(filename, "/z");
grid.dz = grid.z(2) - grid.z(1); 
U = h5read(filename, "/U");
EOA = h5read(filename, "/E");
neOA = h5read(filename, "/ne");
streak = h5read(filename, "/streak");

% Load grid.omg
grid.omg = 2*pi/(max(grid.t)-min(grid.t)+grid.dt)*linspace(-grid.Nt/2,grid.Nt/2-1,grid.Nt); %frequency in rad/s 
grid.domg = grid.omg(2) - grid.omg(1);                                                      %frequency step in rad/s
grid.THz = grid.omg / (2 * pi * 1e12);

% Load grid.r
grid.dr = zeros(grid.Nr,1);
grid.dr(1) = 0.5*(grid.r(2)+grid.r(1));
for ir = 2:grid.Nr-1
    grid.dr(ir) = 0.5*(grid.r(ir+1)-grid.r(ir-1));
end
grid.dr(grid.Nr) = 0.5*(grid.r(grid.Nr) + grid.dr(grid.Nr-1) - grid.r(grid.Nr-1));
clearvars ir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate THz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Calc_THz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Emission Angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Calc_Emission_Angle;