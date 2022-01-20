%% Cleanup

%MODIFY this code to simplify (when generating data, etc.)
clear all; close all; clc;

%% Describe the setup for the forward models

%Specify domain
% [domain] = equigrid_setup(x_disc,y_disc);
domain = struct('x',[],'y',[],'z',[]);
domain.x = [-50:2:50];
domain.y = [-50:2:50];
domain.z = [0 1];

xmin = min(domain.x); xmax = max(domain.x);
ymin = min(domain.y); ymax = max(domain.y);


%Specify boundary types and boundary values (x / y constant head
%boundaries, no flux in z)
bdry_types = [1; 1; 1; 1; 0; 0];
bdry_vals = zeros(6,1);
bdry_L = 1e-5*ones(6,1);
bdrys = struct('types',bdry_types,'vals',bdry_vals,'leaks',bdry_L);

% Locations of all wells (pumping and observation)
well_locs = [...
    -20 -20; ...
    -20 0; ...
    -20 20; ...
    0 -20; ...
    0 0; ...
    0 20; ...
    20 -20; ...
    20 0; ...
    20 20 ...
   ];
num_wells = size(well_locs,1);

% List defining each observation. Columns are:
% pumping angular frequency, pumping well, Q_max, and observation well
% NOTE: The model will run the fastest if this list is sorted by angular
% frequency and pumping well (the first two columns)
%
V = 0.01;
P = [10 50 100 200 400 800 1600];
% P = [1600];
% P = [10];

test_list = [];
for i = 1:1:numel(P)
    for j = 1:1: num_wells
        for k = (j+1):1:num_wells
            test_list = [...
                test_list; ...
                (2*pi)/P(i) j V*pi/P(i) k ...
                ];
        end
    end
end

%% Create the "official" input files needed by the steady periodic model.

%This step should automatically create the input files needed to run all of
%the steady periodic models. 

[experiment] = OHT_create_inputs(well_locs,test_list,domain);

%% Calculate counts (used often for plotting and other functionality)

% Calculates number of observations, size of grid, and sets up a
% meshgridded set of points to do plotting of results.

num_omegas = size(experiment,1);
numobs = size(test_list,1);
num_x = numel(domain.x) - 1;
num_y = numel(domain.y) - 1;
num_cells = num_x*num_y;

%% Setup grid of cell centers for plotting and geostatistical setups
[coords, cgrid] = plaid_cellcenter_coord(domain);

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

figure(1);
for i = 1:1:num_omegas
    disp(['Period = ', num2str(2*pi./experiment(i).omega)]);
    num_testobs = size(experiment(i).tests,1);
    for j = 1:1:num_testobs
        pump_loc = experiment(i).tests(j,1);
        obs_loc = experiment(i).tests(j,2);
        subplot(1,2,1)
        obsmap = reshape(full(experiment(i).obs(:,obs_loc)),num_y,num_x);
        om = pcolor(cgrid{1},cgrid{2},obsmap);
        set(om,'LineStyle','none')
        colorbar
        title(['Period group ', num2str(i), ', Observation weights, test/obs pair #', num2str(j)]);
        axis equal
        axis([xmin xmax ymin ymax])
        subplot(1,2,2)
        pumpmap = reshape(full(experiment(i).stims(:,pump_loc)),num_y,num_x);
        pm = pcolor(cgrid{1},cgrid{2},pumpmap);
        set(pm,'LineStyle','none')
        colorbar
        title(['Period group ', num2str(i), ', Pump weights, test/obs pair #', num2str(j)]);
        axis equal
        axis([xmin xmax ymin ymax])
        pause(0.1)
    end
end

%% For synthetic problem - true parameter field statistics

%Checkerboard case
xl_check = 20; x_offset = 10;
yl_check = 20; y_offset = 10;
check_pattern = sign(sin(pi*(cgrid{1}-x_offset)/xl_check)).* ...
    sign(sin(pi*(cgrid{2}-y_offset)/yl_check));

lnK_mean = -9.2; lnK_jump = 1;
lnSs_mean = -11.2; lnSs_jump = 0.05;

lnK_true_grid = lnK_mean + check_pattern.*lnK_jump;
lnSs_true_grid = lnSs_mean + check_pattern.*lnSs_jump;

lnK_true = reshape(lnK_true_grid,num_cells,1);
lnSs_true = reshape(lnSs_true_grid,num_cells,1);

% 
% 
% %Geostatistical Case
% mean_lnK = -9.2;
% var_lnK = 4;
% 
% mean_lnSs = -11.2;
% var_lnSs = 0.1;
% 
% corr_x = 20;
% corr_y = 20;
% 
% distmat_row = dimdist(coords(1,:),coords);
% corr_row = exp(-(...
%     (distmat_row(:,:,1)./corr_x).^2 + ...
%     (distmat_row(:,:,2)./corr_y).^2).^.5);
% 
% randn('state',0)
% [corr_relz] = toepmat_vector_math(corr_row,'r',[],2,[num_y num_x]);
% 
% lnK_true = corr_relz(:,1).*var_lnK.^.5 + mean_lnK;
% lnSs_true = corr_relz(:,2).*var_lnSs.^.5 + mean_lnSs;
% lnK_true_grid = reshape(lnK_true,num_y,num_x);
% lnSs_true_grid = reshape(lnSs_true,num_y,num_x);

params_true = [lnK_true; lnSs_true];

lnK_range = [min(lnK_true) max(lnK_true)];
lnSs_range = [min(lnSs_true) max(lnSs_true)];

figure(2)
set(2,'Position',[100 100 1200 400])
subplot(1,2,1)
pc1 = pcolor(cgrid{1},cgrid{2},lnK_true_grid);
set(pc1,'LineStyle','none')
title('ln(K) field')
axis equal
axis([xmin xmax ymin ymax])
colorbar
caxis(lnK_range)
subplot(1,2,2)
pc2 = pcolor(cgrid{1},cgrid{2},lnSs_true_grid);
set(pc2,'LineStyle','none')
title('ln(Ss field)')
axis equal
axis([xmin xmax ymin ymax])
colorbar
caxis(lnSs_range);

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
% TODO: May change this so that it is the given pumping test for each
% observation. Also, need to add description of structure.
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

tic
sim_obs = h_function(params_true);
toc

tic
Phi_true = Phi_function(params_true);
toc

tic

params_init = [-9*ones(num_cells,1); -9*ones(num_cells,1)];

tic
H_adj = H_tilde_func(params_init);
toc

tic
Phi_init = Phi_function(params_init);
toc

%% Cross-checks - visualize results and sensitivities

%These cross-checks are used simply to make sure that the model appears to
%be running correctly. We calculate the amplitude and phase first (which
%are more intuitive than A and B), and then look at how these vary with
%space.
%
%The second set of plots looks at the sensitivity structure of the signal
%amplitude to K and Ss, which should look like "blobs" centered around the
%pumping and observation well for each observation.

A_obs = sim_obs(1:numobs);
B_obs = sim_obs((numobs+1):(2*numobs));
amp = (A_obs.^2 + B_obs.^2).^.5;
phase = atan2(-B_obs,A_obs);

A_full = Phi_true(1:num_cells,:);
B_full = Phi_true((num_cells+1):(2*num_cells),:);
amp_full = (A_full.^2 + B_full.^2).^.5;
phase_full = atan2(-B_full,A_full);

H_AK = H_adj((1:numobs),(1:num_cells));
H_BK = H_adj((numobs+1):(2*numobs),(1:num_cells));

H_ASs = H_adj((1:numobs),((num_cells+1):(2*num_cells)));
H_BSs = H_adj((numobs+1):(2*numobs),((num_cells+1):(2*num_cells)));


%Check all phasor fields produced through numerical model running
figure(3)
num_totalstims = size(Phi_init,2);
for i = 1:1:num_totalstims
    amp_field = reshape(log(amp_full(:,i)),num_y,num_x);
    phase_field = reshape(phase_full(:,i),num_y,num_x);
    subplot(1,2,1);
    pc3 = pcolor(cgrid{1},cgrid{2},amp_field);
    set(pc3,'LineStyle','none')
    title(['Pumping test ', num2str(i), ', ln(Amplitude) field'])
    axis equal
    axis([xmin xmax ymin ymax])
    subplot(1,2,2);
    pc4 = pcolor(cgrid{1},cgrid{2},phase_field);
    set(pc4,'LineStyle','none')
    title(['Pumping test ', num2str(i), ', Phase field'])
    axis equal
    axis([xmin xmax ymin ymax])
    pause(0.1)
end

figure(4)
for i = 1:1:numobs
    ampK_sens = (A_obs(i)./amp(i)).*H_AK(i,:) + (B_obs(i)./amp(i)).*H_BK(i,:);
    ampK_sensmap = reshape(ampK_sens,num_y,num_x);
    subplot(1,2,1)
    pc1 = pcolor(cgrid{1},cgrid{2},ampK_sensmap);
    set(pc1,'LineStyle','none')
    title(['Obs. ', num2str(i), ', Amplitude sensitivity to K'])
    axis equal
    axis([xmin xmax ymin ymax])
    ampSs_sens = (A_obs(i)./amp(i)).*H_ASs(i,:) + (B_obs(i)./amp(i)).*H_BSs(i,:);
    ampSs_sensmap = reshape(ampSs_sens,num_y,num_x);
    sensmapSs = reshape(ampSs_sensmap,num_y,num_x);
    subplot(1,2,2)
    pc2 = pcolor(cgrid{1},cgrid{2},sensmapSs);
    set(pc2,'LineStyle','none')
    title(['Obs. ', num2str(i), ', Amplitude sensitivity to Ss'])
    axis equal
    axis([xmin xmax ymin ymax])
    pause(0.1)
end

%% Prior setup

var_lnK_est = 4;
var_lnSs_est = 0.1;

corr_x_est = 20;
corr_y_est = 20;

distmat_row = dimdist(coords(1,:),coords);
corr_row = exp(-(...
    (distmat_row(:,:,1)./corr_x_est).^2 + ...
    (distmat_row(:,:,2)./corr_y_est).^2).^.5);


QK_row = var_lnK_est.*corr_row;
QSs_row = var_lnSs_est.*corr_row;

Qprod_func = @(invec) covar_product_K_Ss(QK_row,QSs_row,invec,num_x,num_y);

%% Inversion

lnK_homog_guess = -9;
lnSs_homog_guess = -11;

%TODO: Use raw phasor data to generate noisy time-series data, then
%re-extract using LSQ process

data_errorcov_guess = 1e-8;

params_init = [lnK_homog_guess*ones(num_cells,1); lnSs_homog_guess*ones(num_cells,1)];
beta_init = [lnK_homog_guess; lnSs_homog_guess];
X = [ones(num_cells,1) zeros(num_cells,1); zeros(num_cells,1) ones(num_cells,1)];
R = data_errorcov_guess.*eye(2*numobs);
y = sim_obs;

tic;
[params_best, beta_best, H_local, negloglike] = ql_geostat_inv(y,params_init,beta_init,X,R,Qprod_func,h_function,H_tilde_func);
inv_time = toc;

negloglike_func = @(s,beta) negloglike_eval(y,X,s,beta,Q,R,h_function);

%% Plotting of true vs. results

figure(5)
subplot(1,2,1)
plot(well_locs(:,1),well_locs(:,2),'ok')
hold on
pc1 = pcolor(cgrid{1},cgrid{2},reshape(params_best(1:num_cells),num_y,num_x));
plot(well_locs(:,1),well_locs(:,2),'ok')
set(pc1,'LineStyle','none')
hold off
title('Estimated ln(K[m/s])')
xlabel('x (m)')
ylabel('y (m)')
colorbar
caxis(lnK_range)
axis equal
axis([xmin xmax ymin ymax])
set(gca,'FontSize',16)


subplot(1,2,2)
hold on
pc1 = pcolor(cgrid{1},cgrid{2},reshape(params_best((num_cells+1):(2*num_cells)),num_y,num_x));
plot(well_locs(:,1),well_locs(:,2),'ok')
set(pc1,'LineStyle','none')
hold off
title('Estimated ln(Ss [1/m])')
xlabel('x (m)')
ylabel('y (m)')
colorbar
caxis(lnSs_range)
axis equal
axis([xmin xmax ymin ymax])
set(gca,'FontSize',16)

set(5,'Position',[100 100 1200 400])

figure(6)
subplot(1,2,1)
hold on
pc2 = pcolor(cgrid{1},cgrid{2},reshape(params_true(1:num_cells),num_y,num_x));
plot(well_locs(:,1),well_locs(:,2),'ok')
set(pc2,'LineStyle','none')
hold off
title('True ln(K[m/s])')
xlabel('x (m)')
ylabel('y (m)')
colorbar
caxis(lnK_range)
axis equal
axis([xmin xmax ymin ymax])
set(gca,'FontSize',16)

subplot(1,2,2)
hold on
pc2 = pcolor(cgrid{1},cgrid{2},reshape(params_true((num_cells+1):(2*num_cells)),num_y,num_x));
plot(well_locs(:,1),well_locs(:,2),'ok')
set(pc2,'LineStyle','none')
hold off
title('True ln(Ss[1/m])')
xlabel('x (m)')
ylabel('y (m)')
colorbar
caxis(lnSs_range)
axis equal
axis([xmin xmax ymin ymax])
set(gca,'FontSize',16)

set(6,'Position',[100 100 1200 400])

%% Save output

save('inverse_test_2D_try2_out.mat')

