%% Cleanup

%MODIFY this code to simplify (when generating data, etc.)
clear all; close all; clc;

%% Describe the setup for the forward models

% Bounds of the domain, and number of points in the discretization (last)
x_disc = [-90 90 90];
y_disc = [-45 45 90];

%Specify boundary types and boundary values (x / y constant head
%boundaries, no flux in z)
bdry_types = [1; 1; 1; 1; 0; 0];
bdry_vals = zeros(6,1);
bdrys = struct('types',bdry_types,'vals',bdry_vals);

% Locations of all wells (pumping and observation)
well_locs = [...
    -20 20; ...
    0 20; ...
    20 20; ...
    -20 0; ...
    0 0; ...
    20 0; ...
    -20 -20; ...
    0 -20; ...
    20 -20; ...
    ];

% List defining each observation. Columns are:
% pumping angular frequency, pumping well, Q_max, and observation well
% NOTE: The model will run the fastest if this list is sorted by angular
% frequency and pumping well (the first two columns)
%
V = 0.005;
P = [10 100 1000];

test_list = [];
for i = 1:1:numel(P)
    test_list = [test_list; ...
        (2*pi)/P(i) 5 V*pi/P(i) 1; ...
        (2*pi)/P(i) 5 V*pi/P(i) 2; ...
        (2*pi)/P(i) 5 V*pi/P(i) 3; ...
        (2*pi)/P(i) 5 V*pi/P(i) 4; ...
        (2*pi)/P(i) 5 V*pi/P(i) 6; ...
        (2*pi)/P(i) 5 V*pi/P(i) 7; ...
        (2*pi)/P(i) 5 V*pi/P(i) 8; ...
        (2*pi)/P(i) 5 V*pi/P(i) 9; ...
    ];
end

%% Create the "official" input files needed by the steady periodic model.

%This step should automatically create the input files needed to run all of
%the steady periodic models. 

[domain] = equigrid_setup(x_disc,y_disc);

[experiment] = OHT_create_inputs(well_locs,test_list,domain);

%% Calculate counts (used often for plotting and other functionality)

% Calculates number of observations, size of grid, and sets up a
% meshgridded set of points to do plotting of results.

num_omegas = size(experiment,1);
num_obs = size(test_list,1);
num_x = numel(domain.x) - 1;
num_y = numel(domain.y) - 1;
num_cells = num_x*num_y;

%% Setup grid of cell centers for plotting and geostatistical setups
[coords, grid] = plaid_cellcenter_coord(domain);

%% Visualize inputs

% Each time the angular frequency changes in "test_list", a new "omega
% group" is created. This is a set of models that get run together, which
% saves time.
%
% This block of code is just to make sure that all of the inputs for the
% steady periodic model have been created correctly. It will show:
% 1) In the command window, the current "omega group"
% 2) In Figure 1, the observation weighting for the current observation
% 3) In figure 2, the flux weighting for the given pumping test
%
% This will loop through all observations as you hit <Enter>, showing the
% model inputs for each individual model run that generates an observation.

for i = 1:1:num_omegas
    disp(['Omega = ', num2str(experiment(i).omega)]);
    num_testobs = size(experiment(i).tests,1);
    for j = 1:1:num_testobs
        pump_loc = experiment(i).tests(j,1);
        obs_loc = experiment(i).tests(j,2);
        figure(1);
        obsmap = reshape(full(experiment(i).obs(:,obs_loc)),num_y,num_x);
        om = pcolor(grid{1},grid{2},obsmap);
        set(om,'LineStyle','none')
        colorbar
        title(['Omega group ', num2str(i), ', Observation weights, test/obs pair #', num2str(j)]);
        figure(2)
        pumpmap = reshape(full(experiment(i).stims(:,pump_loc)),num_y,num_x);
        pm = pcolor(grid{1},grid{2},pumpmap);
        set(pm,'LineStyle','none')
        colorbar
        title(['Omega group ', num2str(i), ', Pump weights, test/obs pair #', num2str(j)]);
        pause
    end
end

%% For synthetic problem - true parameter field statistics

%Creation of true parameter field. With given constants, the parameter
%field will be homogeneous. In order to change this, simply change the
%var_lnK and var_lnSS values.

mean_lnK = -9.2;
var_lnK = 0;

mean_lnSs = -9.2;
var_lnSs = 0;

corr_x = 50;
corr_y = 10;

distmat_row = dimdist(coords(1,:),coords);
corr_row = exp(-(...
    (distmat_row(:,:,1)./corr_x).^2 + ...
    (distmat_row(:,:,2)./corr_y).^2).^.5);

randn('state',0)
[corr_relz] = toepmat_vector_math(corr_row,'r',[],2,[num_y num_x]);

lnK_true = corr_relz(:,1).*var_lnK.^.5 + mean_lnK;
lnSs_true = corr_relz(:,2).*var_lnSs.^.5 + mean_lnSs;
params_true = [lnK_true; lnSs_true];

lnK_true_grid = reshape(lnK_true,num_y,num_x);
lnSs_true_grid = reshape(lnSs_true,num_y,num_x);

lnK_range = [min(lnK_true) max(lnK_true)];
lnSs_range = [min(lnSs_true) max(lnSs_true)];

figure(1)
subplot(1,2,1)
pc1 = pcolor(grid{1},grid{2},lnK_true_grid);
set(pc1,'LineStyle','none')
colorbar
% caxis(lnK_range)
subplot(1,2,2)
pc2 = pcolor(grid{1},grid{2},lnSs_true_grid);
set(pc2,'LineStyle','none')
colorbar
% caxis(lnSs_range);

%% Perform all model runs to generate data

% Steady_periodic_multimodel_run performs the model runs to generate all
% observations. 
% IMPORTANT:
% The first output is simulated observations, output as
% follows: [A_obs(1); A_obs(2); A_obs(3) ... A_obs(numobs); B_obs(1) ...
% B_obs(numobs)]
%
% NOTE: B in this case is the imaginary part of the phasor. Note that this
% is not the same as the coefficient in front of the sine (which would be
% negative B)
%
% The second output is the full phasor field for each pumping test.
%
% H is the sensitivity matrix. Each row represents a different observation,
% organized the same way as sim_obs. Each column represents a different
% parameter sensitivity - the first set of columns is with respect to each
% K value, and the second set of columns is with respect to each Ss value.
% K and Ss grid cells are stepped through in sequential order by y, then x,
% then z, as would be done by meshgrid.
h_function = @(params) OHT_run_distribKSs(params, ...
    domain,bdrys,experiment,1);

h_homog_function = @(KSs) h_function([KSs(1)*ones(num_cells,1); KSs(2)*ones(num_cells,1)]);

Phi_function = @(params) OHT_run_distribKSs(params, ...
    domain,bdrys,experiment,2);

H_tilde_func = @(params) OHT_run_distribKSs(params, ...
    domain,bdrys,experiment,3);

%Timing of forward model running
tic
sim_obs = h_function(params_true);
toc

%Timing of forward model running, returning full phasor field
tic
Phi_true = Phi_function(params_true);
toc

%Timing of adjoint sensitivity calculations
tic
H_adj = H_tilde_func(params_true);
toc

%% Calculate sensitivities using FD approach

H_FD = zeros(2*num_obs,2*num_cells);

%Time the calculation of the sensitivity matrix using the FD approach
tic
params_base = params_true;
h_base = h_function(params_base);
for j = 1:1:(2*num_cells)
    del = 0.1;
    params_mod = params_true;
    params_mod(j) = params_mod(j) + del;
    h_mod = h_function(params_mod);
    H_FD(:,j) = (h_mod - h_base)./del;
    disp(['Parameter #', num2str(j)])
end
toc

%% Cross-checks - visualize results and sensitivities

H_adj_AK = H_adj((1:num_obs),(1:num_cells));
H_adj_BK = H_adj((num_obs+1):(2*num_obs),(1:num_cells));

H_adj_ASs = H_adj((1:num_obs),((num_cells+1):(2*num_cells)));
H_adj_BSs = H_adj((num_obs+1):(2*num_obs),((num_cells+1):(2*num_cells)));

H_FD_AK = H_FD((1:num_obs),(1:num_cells));
H_FD_BK = H_FD((num_obs+1):(2*num_obs),(1:num_cells));

H_FD_ASs = H_FD((1:num_obs),((num_cells+1):(2*num_cells)));
H_FD_BSs = H_FD((num_obs+1):(2*num_obs),((num_cells+1):(2*num_cells)));

A_obs = sim_obs(1:num_obs);
B_obs = sim_obs((num_obs+1):(2*num_obs));
amp = (A_obs.^2 + B_obs.^2).^.5;
phase = atan2(-B_obs,A_obs);

for i = 1:1:num_obs
    ampK_adj_sens = (A_obs(i)./amp(i)).*H_adj_AK(i,:) + (B_obs(i)./amp(i)).*H_adj_BK(i,:);
    ampK_adj_sensmap = reshape(ampK_adj_sens,num_y,num_x);
    ampSs_adj_sens = (A_obs(i)./amp(i)).*H_adj_ASs(i,:) + (B_obs(i)./amp(i)).*H_adj_BSs(i,:);
    ampSs_adj_sensmap = reshape(ampSs_adj_sens,num_y,num_x);

    ampK_FD_sens = (A_obs(i)./amp(i)).*H_FD_AK(i,:) + (B_obs(i)./amp(i)).*H_FD_BK(i,:);
    ampK_FD_sensmap = reshape(ampK_FD_sens,num_y,num_x);
    ampSs_FD_sens = (A_obs(i)./amp(i)).*H_FD_ASs(i,:) + (B_obs(i)./amp(i)).*H_FD_BSs(i,:);
    ampSs_FD_sensmap = reshape(ampSs_FD_sens,num_y,num_x);
    
    ampK_range = [min(ampK_FD_sens) max(ampK_FD_sens)];
    ampSs_range = [min(ampSs_FD_sens) max(ampSs_FD_sens)];

    figure(3)
    subplot(1,3,1)
    pc1 = pcolor(grid{1},grid{2},ampK_adj_sensmap);
    set(pc1,'LineStyle','none')
    set(gca,'FontSize',14)
    axis equal
    axis([x_disc(1) x_disc(2) y_disc(1) y_disc(2)])
    title('Adjoint Sensitivity Visualization: d(amp)/d(lnK)')
    colorbar
    caxis(ampK_range)
    subplot(1,3,2)
    pc2 = pcolor(grid{1},grid{2},ampK_FD_sensmap);
    set(pc2,'LineStyle','none')
    set(gca,'FontSize',14)
    axis equal
    axis([x_disc(1) x_disc(2) y_disc(1) y_disc(2)])
    title('FD Sensitivity Visualization: d(amp)/d(lnK)')
    colorbar
    caxis(ampK_range)
    subplot(1,3,3)
    pc3 = pcolor(grid{1},grid{2},(ampK_FD_sensmap - ampK_adj_sensmap));
    set(pc3,'LineStyle','none')
    set(gca,'FontSize',14)
    axis equal
    axis([x_disc(1) x_disc(2) y_disc(1) y_disc(2)])
    title('FD - Adjoint difference: d(amp)/d(lnK)')
    colorbar
    caxis(ampK_range)
    set(3,'Position',[0 0 2000 300])

    figure(4)
    subplot(1,3,1)
    pc1 = pcolor(grid{1},grid{2},ampSs_adj_sensmap);
    set(pc1,'LineStyle','none')
    set(gca,'FontSize',14)
    axis equal
    axis([x_disc(1) x_disc(2) y_disc(1) y_disc(2)])
    title('Adjoint Sensitivity Visualization: d(amp)/d(lnSs)')
    colorbar
    caxis(ampSs_range)
    subplot(1,3,2)
    pc2 = pcolor(grid{1},grid{2},ampSs_FD_sensmap);
    set(pc2,'LineStyle','none')
    set(gca,'FontSize',14)
    axis equal
    axis([x_disc(1) x_disc(2) y_disc(1) y_disc(2)])
    title('FD Sensitivity Visualization: d(amp)/d(lnSs)')
    colorbar
    caxis(ampSs_range)
    subplot(1,3,3)
    pc3 = pcolor(grid{1},grid{2},(ampSs_FD_sensmap - ampSs_adj_sensmap));
    set(gca,'FontSize',14)
    set(pc3,'LineStyle','none')
    axis equal
    axis([x_disc(1) x_disc(2) y_disc(1) y_disc(2)])
    title('FD - Adjoint difference: d(amp)/d(lnSs)')
    colorbar
    caxis(ampSs_range)
    set(4,'Position',[0 350 2000 300])
    
    pause
end

