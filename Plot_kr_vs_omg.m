%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prelim quantities for transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("/gpfs/fs1/home/aellio18/Desktop/JUPPE")

ell = 0; % orbital angular momentum of pulse --> Bessel of order 0
grid.rmax = json_arr.grid.xmax(2); % Get maximum R
grid.nr = json_arr.grid.nx(2); % Get number of grid points in r, Nr

grid.nt = json_arr.grid.nt;

[grid.r, grid.kperp] = DHTarrays(grid.nr,grid.rmax,ell);
[grid.T_R2K,grid.T_K2R] = DHTtransformMatrix(grid.nr,grid.rmax,ell);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Perform FFT for E(r,t,z) --> E(r,omg,z) 
% Perform Hankel transform for E(r,omg,z) --> E(k_perp,omg,z)
% -------------------------------------------------------------------------
EAp = grid.nt * grid.dt * fftshift(ifft(EF.data(:,:,1),[],2),2);
EAp = grid.T_R2K*EAp;

% Get index of the focal position, zf.
zf = json_arr.laser.pulses(1).zf;
[~,zf_index] = min(abs(zf - Spectrum.axes(2).data));
zf_index = ceil(zf_index / json_arr.diag.fullskip);

% Transformed E-field at focus.
EAp_zf = grid.nt * grid.dt * fftshift(ifft(EF.data(:,:,zf_index),[],2),2);
EAp_zf = grid.T_R2K*EAp_zf;

% Transformed final E-field.
EAp_final = grid.nt * grid.dt * fftshift(ifft(Efinal.data,[],2),2);
EAp_final = grid.T_R2K*EAp_final;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,3,1);
imagesc(grid.omg_norm(i_zero:end), grid.kperp(1:end/5), abs(EAp(1:end/5,i_zero:end))); 
axis xy; 
maxValue = max(abs(EAp(:)));
colorbar; 
colormap(flipud(hot));
set(gca, 'ColorScale','log'); 
xlabel("$\omega [\omega_0]$"); 
ylabel("$k_r$ [m$^{-1}$]"); title('Unfiltered: Initial')

subplot(2,3,2);
imagesc(grid.omg_norm(i_zero:end), grid.kperp(1:end/5), abs(EAp_zf(1:end/5,i_zero:end))); 
axis xy; 
maxValue = max(abs(EAp_zf(:)));
colorbar; 
colormap(flipud(hot)); 
set(gca, 'ColorScale','log'); 
xlabel("$\omega [\omega_0]$");
ylabel("$k_r$ [m$^{-1}$]"); title('Unfiltered: At Focus')
 
subplot(2,3,3);
imagesc(grid.omg_norm(i_zero:end), grid.kperp(1:end/5), abs(EAp_final(1:end/5,i_zero:end))); 
axis xy; 
maxValue = max(abs(EAp_final(:)));
colorbar; 
colormap(flipud(hot));
set(gca, 'ColorScale','log'); 
xlabel("$\omega [\omega_0]$");
ylabel("$k_r$ [m$^{-1}$]"); title('Unfiltered: Final')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter to get THz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EAp = EAp .* full_filter;
EAp_zf = EAp_zf .* full_filter;
EAp_final = EAp_final .* full_filter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,4);
imagesc(grid.THz(i_zero:i_temp), grid.kperp(1:end/20), abs(EAp(1:end/20,i_zero:i_temp))); 
axis xy; 
maxValue = max(abs(EAp(:)));
caxis([0,0.75*maxValue]);
colorbar; 
colormap(flipud(hot));
set(gca, 'ColorScale','linear'); 
xlabel("f [THz]"); 
ylabel("$k_r$ [m$^{-1}$]"); title('THz: Initial')

subplot(2,3,5);
imagesc(grid.THz(i_zero:i_temp), grid.kperp(1:end/20), abs(EAp_zf(1:end/20,i_zero:i_temp))); 
axis xy; 
maxValue = max(abs(EAp_zf(:)));
caxis([0,0.75*maxValue]);
colorbar; 
colormap(flipud(hot));
set(gca, 'ColorScale','linear'); 
xlabel("f [THz]"); 
ylabel("$k_r$ [m$^{-1}$]"); title('THz: At Focus')

subplot(2,3,6);
imagesc(grid.THz(i_zero:i_temp), grid.kperp(1:end/20), abs(EAp_final(1:end/20,i_zero:i_temp))); 
axis xy; 
maxValue = max(abs(EAp_final(:)));
caxis([0,0.75*maxValue]);
colorbar; 
colormap(flipud(hot));
set(gca, 'ColorScale','linear'); 
xlabel("f [THz]"); 
ylabel("$k_r$ [m$^{-1}$]"); title('THz: Final')

filename = plots_folder + sim_title + "_kr_vs_omg.jpg";
saveas(gcf,filename);
clearvars filename;