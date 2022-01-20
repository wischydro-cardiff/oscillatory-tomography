%% Cleanup

clear all; close all; clc;

%% Describe the setup for the forward models

% Bounds of the domain, and number of points in the discretization (last)
x_disc = [-300 300 300];
y_disc = [-300 300 300];

% Locations of all wells (pumping and observation)
well_locs = [...
    0 0; ...
    0 30; ...
    60 0; ...
    0 -90; ...
    -120 0; ...
    ];

%Colors/symbols associated with each well number
format_types = {'*k','ob','^b','or','^r','og','^g','ok','^k'};

%Specify boundary types and boundary values (x / y constant head
%boundaries, no flux in z)
bdry_types = [1; 1; 1; 1; 0; 0];
bdry_vals = zeros(6,1);
bdry_L = 3e-4*ones(6,1);
bdrys = struct('types',bdry_types,'vals',bdry_vals,'leaks',bdry_L);

% List defining each observation. Columns are:
% pumping angular frequency, pumping well, Q_max, and observation well
% NOTE: The model will run the fastest if this list is sorted by angular
% frequency and pumping well (the first two columns)
%
periods = logspace(1,4,20);
test_list = [];
Qmax = -1e-3;
for ind = 1:1:numel(periods)
    test_list = [test_list; ...
        (2*pi)/periods(ind) 1 Qmax 2; ...
        (2*pi)/periods(ind) 1 Qmax 3; ...
        (2*pi)/periods(ind) 1 Qmax 4; ...
        (2*pi)/periods(ind) 1 Qmax 5; ...
    ];
end

% Volumes used in various pumping tests
volumes = abs(Qmax).*periods/pi;

%% Create the "official" input files needed by the steady periodic model.

%This step will automatically create the input files needed to run all of
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
xmin = min(domain.x); xmax = max(domain.x);
ymin = min(domain.y); ymax = max(domain.y);


%% Setup grid of cell centers for plotting and geostatistical setups
[coords, cgrid] = plaid_cellcenter_coord(domain);

%% For synthetic problem - true parameter field statistics

%Setup a homogeneous or heterogeneous parameter field with the given mean
%ln(T) / ln(S) values and the associated variances and correlation lengths.

mean_lnT = log(3e-4);
var_lnT = 0; %1;

mean_lnS = log(1e-5);
var_lnS = 0; %1;

corr_x = 55;
corr_y = 2;

distmat_row = dimdist(coords(1,:),coords);
corr_row = exp(-(...
    (distmat_row(:,:,1)./corr_x).^2 + ...
    (distmat_row(:,:,2)./corr_y).^2).^.5);

randn('state',0)
[corr_relz] = toepmat_vector_math(corr_row,'r',[],2,[num_y num_x]);

lnK_true = corr_relz(:,1).*var_lnT.^.5 + mean_lnT;
lnSs_true = corr_relz(:,2).*var_lnS.^.5 + mean_lnS;
params_true = [lnK_true; lnSs_true];

lnK_true_grid = reshape(lnK_true,num_y,num_x);
lnSs_true_grid = reshape(lnSs_true,num_y,num_x);

figure(1)
subplot(1,2,1)
pc1 = pcolor(cgrid{1},cgrid{2},log10(exp(lnK_true_grid)));
set(pc1,'LineStyle','none')
colorbar
% caxis(K_range)
subplot(1,2,2)
pc2 = pcolor(cgrid{1},cgrid{2},log10(exp(lnSs_true_grid)));
set(pc2,'LineStyle','none')
colorbar
% caxis(Ss_range);

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
sim_Phi = Phi_function(params_true);
toc


%%

A_full = sim_Phi(1:num_cells,:);
B_full = sim_Phi((num_cells+1):(2*num_cells),:);
amp_full = (A_full.^2 + B_full.^2).^.5;
phase_full = atan2(-B_full,A_full);

%Check all phasor fields produced through numerical model running
figure(3)
num_totalstims = size(sim_Phi,2);
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
    pause
end

%% 

%Extract data from the model runs, convert into amplitudes and phases of
%drawdown in order to perform cross-comparison with analytical models
%developed by Black and Kipp

A_obs = -sim_obs(1:num_obs);
B_obs = -sim_obs((num_obs+1):(2*num_obs));
amp = (A_obs.^2 + B_obs.^2).^.5;
phase = atan2(-B_obs,A_obs);

%Store the data generated using the numerical model in a structure
%synth_data
synth_data = zeros(num_obs,8);

for ind = 1:1:num_obs
    synth_data(ind,1) = test_list(ind,4); %Well number
    synth_data(ind,2) = norm(well_locs(test_list(ind,2),:) - well_locs(test_list(ind,4),:)); %Radius
    synth_data(ind,3) = 2*pi./test_list(ind,1); %period
    synth_data(ind,4) = abs(test_list(ind,3)); %Qmax
    synth_data(ind,5) = amp(ind); %Amplitude
    if phase(ind) < 0
        synth_data(ind,6) = phase(ind) + 2*pi; %Phase
    else
        synth_data(ind,6) = phase(ind);
    end
    synth_data(ind,7) = A_obs(ind); %Phasor real part
    synth_data(ind,8) = B_obs(ind); %Phasor imaginary part
end

%% Compare amplitude / phase graphs for various distances and frequencies

%Store the solution from an analytical model (the phase / amplitude
%developed in Black and Kipp, as corrected by Cardiff

analyt_soln = zeros(num_obs,8);

for ind = 1:1:num_obs
    r = synth_data(ind,2);
    om = 2.*pi./synth_data(ind,3);
    Q_max = synth_data(ind,4);
    S = exp(mean_lnS);
    T = exp(mean_lnT);
    analyt_soln(ind,1) = synth_data(ind,1);
    analyt_soln(ind,2) = r;
    analyt_soln(ind,3) = 2*pi./om;
    analyt_soln(ind,4) = Q_max;
    u = ((om.*S.*r.^2)./(2.*T)).^.5;
    phasor_2D_cardiff = Q_max./(2.*pi.*T).*besselk(0,(u+1i*u));
    amp_2D_cardiff = Q_max./(2.*pi.*T).*abs(besselk(0,(u+1i*u)));
    phase_2D_cardiff = -angle(besselk(0,(u+1i*u)));
    if phase_2D_cardiff < 0
        phase_2D_cardiff = phase_2D_cardiff + 2*pi;
    end
    analyt_soln(ind,5) = amp_2D_cardiff;
    analyt_soln(ind,6) = phase_2D_cardiff;
    analyt_soln(ind,7) = real(phasor_2D_cardiff);
    analyt_soln(ind,8) = imag(phasor_2D_cardiff);
end

%% Analytical / Numerical comparison

%Produce plots showing the amplitude and phase given by the numerical model
%and the analytical model. Calculate simple metrics such as the mean
%absolute difference in order to quantify any errors.
figurenum = 2;
figureh = figure(figurenum); clf
set(figurenum,'Position',[100 100 1600 600])
subplot(1,2,1)
rad_list = unique(synth_data(:,2));
num_rad = numel(rad_list);
for j = 1:1:num_rad
    crad = rad_list(j);
    crad_ind = synth_data(:,2) == crad;
    hold on
    plot(synth_data(crad_ind,3),synth_data(crad_ind,5),'^k');
    plot(analyt_soln(crad_ind,3),analyt_soln(crad_ind,5),'k');
    hold off
end
grid on
set(gca,'XScale','log','YScale','log')
set(gca,'FontSize',16)
xlabel('Period (s)')
ylabel('Amplitude (m)')

subplot(1,2,2)
rad_list = unique(synth_data(:,2));
num_rad = numel(rad_list);
for j = 1:1:num_rad
    crad = rad_list(j);
    crad_ind = synth_data(:,2) == crad;
    hold on
    plot(synth_data(crad_ind,3),synth_data(crad_ind,6),'^k');
    plot(analyt_soln(crad_ind,3),analyt_soln(crad_ind,6),'k');
    hold off
end
grid on
set(gca,'XScale','log','YScale','linear')
set(gca,'FontSize',16)
xlabel('Period (s)')
ylabel('Phase Delay (radians)')
legend('Numerical Results','Analytical Solution')

%Add annotations in positions that look good
annotation(figureh,'textbox',...
    [0.128197311965536 0.322333333333337 0.0354266357421875 0.0383333333333333],'String',{'120m'},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FontName','Helvetica');

annotation(figureh,'textbox',...
    [0.130005478469442 0.424000000000003 0.030560302734375 0.0383333333333333],'String',{'90m'},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FontName','Helvetica');

annotation(figureh,'textbox',...
    [0.128130478469442 0.537333333333336 0.030560302734375 0.0383333333333333],'String',{'60m'},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FontName','Helvetica');

annotation(figureh,'textbox',...
    [0.132505478469442 0.665666666666669 0.030560302734375 0.0383333333333333],'String',{'30m'},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FontName','Helvetica');

annotation(figureh,'textbox',...
    [0.700630478469442 0.249000000000004 0.030560302734375 0.0383333333333333],'String',{'30m'},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FontName','Helvetica');

annotation(figureh,'textbox',...
    [0.720072311965536 0.497333333333339 0.0354266357421875 0.0383333333333333],'String',{'120m'},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FontName','Helvetica');

annotation(figureh,'textbox',...
    [0.706880478469442 0.439000000000005 0.030560302734375 0.0383333333333333],'String',{'90m'},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FontName','Helvetica');

annotation(figureh,'textbox',...
    [0.701880478469442 0.350666666666671 0.030560302734375 0.0383333333333333],'String',{'60m'},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FontName','Helvetica');

amp_abs_error = abs((analyt_soln(:,5) - synth_data(:,5))./analyt_soln(:,5));
mean_amp_abs_error = mean(amp_abs_error);

phase_abs_error = abs((analyt_soln(:,6) - synth_data(:,6)))
%Correction for any phase wrapping issues
phase_abs_error = (phase_abs_error).*(phase_abs_error < pi) + ...
    (2*pi - phase_abs_error).*(phase_abs_error >= pi);
mean_phase_abs_error = mean(phase_abs_error)




%% Process synthetic data results using homogeneous assumption

%These steps can be used to investigate how much error would be made (in
%terms of parameter estimation) by using the results of the analytical
%model with associated discretization errors to determine effective T / S
%values using the analytical model. Based on approach presented by
%Rasmussen for inverting amplitude and phase data


results = zeros(num_obs,5);
for j = 1:1:num_obs

    well_num = synth_data(j,1);
    r = synth_data(j,2);
    om = 2*pi./synth_data(j,3);
    Q_max = synth_data(j,4);
    amp_obs = synth_data(j,5);
    phase_obs = synth_data(j,6);
    period = 2*pi./(om);
    
    c = [-0.12665 2.8642 -0.47779 0.16586 -0.076402 0.03089];
    ras_sum = 0;
    for k = 0:1:5
        ras_sum = ras_sum + c(k+1)*(log(phase_obs)).^k;
    end
    
    diff_est = om*r.^2./(exp(ras_sum));
    u_est = ((om.*r.^2)./(2.*diff_est)).^.5;
    
    T_est = Q_max./(2.*pi.*amp_obs).*abs(besselk(0,(u_est+1i*u_est)));
    S_est = T_est./diff_est;
    results(j,:) = [well_num period T_est S_est diff_est];
end

%% Plotting results

%Plotting of parameter estimation results when using the phase / amplitude
%analysis developed by Rasmussen.

xscale = 'log'; yscale = 'log';
per_range = [min(results(:,2)) max(results(:,2))];
%Minimum trusted signal amplitude for plotting
amp_tol = 1e-4;
K_central_vec = exp(mean_lnT)*ones(2,1);
K_max_vec = exp(max(lnK_true))*ones(2,1);
K_min_vec = exp(min(lnK_true))*ones(2,1);
Ss_central_vec = exp(mean_lnS)*ones(2,1);
Ss_max_vec = exp(max(lnSs_true))*ones(2,1);
Ss_min_vec = exp(min(lnSs_true))*ones(2,1);
D_central_vec = exp(mean_lnT - mean_lnS)*ones(2,1);
D_max_vec = exp(max(lnK_true - lnSs_true))*ones(2,1);
D_min_vec = exp(min(lnK_true - lnSs_true))*ones(2,1);


format_list = cell(num_obs,1);
for j = 1:1:num_obs
    format_list{j} = format_types{results(j,1)};
end

figure(4); clf;
hold on
for j = 1:1:num_obs
    if synth_data(j,5) > amp_tol
        plot(results(j,2),results(j,5),format_list{j})
    end
end
plot(per_range,D_central_vec,'g')
plot(per_range,D_max_vec,'r')
plot(per_range,D_min_vec,'b')
set(gca,'XScale',xscale,'YScale',yscale)
hold off
xlabel('Period (s)')
ylabel('D (m^2/s)')
title('Effective Diffusivity')

figure(5); clf;
hold on
for j = 1:1:num_obs
    if synth_data(j,5) > amp_tol
        plot(results(j,2),results(j,3),format_list{j})
    end
end
plot(per_range,K_central_vec,'g')
plot(per_range,K_max_vec,'r')
plot(per_range,K_min_vec,'b')
set(gca,'XScale',xscale,'YScale',yscale)
hold off
xlabel('Period (s)')
ylabel('T (m^2/s)')
title('Effective Transmissivity')

figure(6); clf;
hold on
for j = 1:1:num_obs
    if synth_data(j,5) > amp_tol
        plot(results(j,2),results(j,4),format_list{j})
    end
end
plot(per_range,Ss_central_vec,'g')
plot(per_range,Ss_max_vec,'r')
plot(per_range,Ss_min_vec,'b')
set(gca,'XScale',xscale,'YScale',yscale)
hold off
xlabel('Period (s)')
ylabel('S (-)')
title('Effective Storativity')

