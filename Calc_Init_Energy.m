%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the initial energy for multiple simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path to SupData.m class
addpath("/gpfs/fs1/home/aellio18/Desktop/SUPER/super/analysis")

% Simulations
% outer_folder = "h-runs/new/";
% sim_titles = ["h-ffa-24", "h-ffa-26", "h-ffh-24", ...
%     "h-ga-24", "h-ga-26", "h-gh-24"];
% 
% calc_energy(sim_titles, outer_folder);
% clearvars sim_titles outer_folder
% 
outer_folder = "h-runs_sub/";
sim_titles = ["h-ffa-24", "h-ffa-25", "h-ffa-26", "h-ffh-24", ...
    "h-ffh-25", "h-ffh-26"];

calc_energy(sim_titles, outer_folder);
clearvars sim_titles outer_folder

outer_folder = "h-runs_super/";
sim_titles = ["h-ffa-24", "h-ffa-25", "h-ffa-26", "h-ffh-24", ...
    "h-ffh-25", "h-ffh-26"];

calc_energy(sim_titles, outer_folder);
clearvars sim_titles outer_folder

% outer_folder = "h-runs/";
% sim_titles = ["h-ffa-25","h-ffh-25", "h-ffh-26","h-ga-25", "h-gh-25", "h-gh-26"];

% outer_folder = "i-runs/";
% sim_titles = ["i_640_sub","i_640_super", "i_800_sub","i_800_super"];
% 
% calc_energy(sim_titles, outer_folder);
% clearvars sim_titles outer_folder

function [] = calc_energy(sim_titles, outer_folder)
    init_U_arr.sim = sim_titles;
    init_U_arr.energy = zeros(length(sim_titles),1);
    
    for i=1:length(sim_titles)
        sim_name = sim_titles(i);
        data_folder = outer_folder + sim_name + "/DATA/";
        plots_folder = data_folder + "plots/";
    
        U = SupData(data_folder + 'U.h5');
        init_U_arr.energy(i) = U.data(1);
    end
    
    save(outer_folder + 'energy_constraint','init_U_arr');
end


