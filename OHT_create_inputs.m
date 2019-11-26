function [experiment] = OHT_create_inputs(well_locs,test_list,domain)
   
%OHT_create_inputs: This function takes simple user inputs -
%the well locations, a matrix defining the list of tests, and
%the domain / discretization information, and creates the inputs
% required by the OHT_run_distribKSs.
%
%[experiment] =
%    OHT_create_inputs(well_locs,test_list,domain)
%
%where:
%   OUTPUTS:
%      -experiment is a structure array containing 1 element per testing
%      frequency. The fields within the structure array are: omega, tests,
%      stims, and obs
%   INPUTS:
%       -well_locs (numwells x dim) is a matrix containing the location of each
%       well on a row. dim may be 2 or 3. numwells is the number of wells.
%       -test_list is a (numobs x 4) or a (numobs x 6) vector. (numobs x 4)
%       vectors are used for monopole tests, and (numobs x 6) are used for
%       dipole tests. The columns containing the following information:
%          Column 1: Angular frequency of pumping (radians / time)
%          Column 2: Well performing pumping (numbered according to entries in
%          well_locs)
%          Column 3: Peak flow rate of pumping (L^3/T), representing the
%          flowrate at the peak / trough of a sinusoid
%          Column 4: (If dipole testing) For second well, well performing pumping
%          (numbered according to entries in well_locs)
%          Column 5: (If dipole testing) For second well, peak flow rate of
%          pumping (L^3/T), representing the flowrate at the peak / trough of
%          a sinusoid
%          Final Column: Observation well (numbered according to entries in
%          well_locs)
%       -domain is a structure containing the fields x, y, and z, which are the
%       locations of the model grid cell boundaries.
%
% Note: In the columns of "test_list" the values for peak flow rates will
% generally be given as real numbers, representing an oscillating source
% with zero phase offset relative to a cosine. However, if desired, complex
% numbers may be supplied to the columns of test_list that represent peak
% flow rates in order to produce oscillating sources with different phase
% offsets.
%
% Code by Michael Cardiff
% 10/2014, Last Updated: 1/2016

%% Map stimulations to the given grid

if ~isfield(domain,'x')
    error('X boundary information must be supplied in domain_disc.x');
else
    xb = domain.x;
end

if ~isfield(domain,'y')
    error('Y boundary information must be supplied in domain_disc.y');
else
    yb = domain.y;
end

if isfield(domain,'z')
    zb = domain.z;
else
    zb = [0 1];
end

xc = (xb(2:end) + xb(1:(end-1)))./2;
yc = (yb(2:end) + yb(1:(end-1)))./2;
zc = (zb(2:end) + zb(1:(end-1)))./2;

num_z = numel(zc);

if (num_z > 1)
    if (size(well_locs,2) < 3)
        error('well_locs must be given in 3D (points or intervals) if z is discretized')
    else
        well_weights = grid_idw(well_locs,xc,yc,zc)';
    end
else
    if (size(well_locs,2) ~= 2)
        error('well_locs must be given in 2D if z has no discretization');
    else
        well_weights = grid_idw(well_locs,xc,yc)';
    end
end

%% setup a different structure for each group of omegas

%Pre-process to check how many times the frequency changes
freq_count = 1;
num_totalobs = size(test_list,1);
for t = 2:1:num_totalobs;
    if test_list(t,1) ~= test_list(t-1,1)
        freq_count = freq_count + 1;
    end
end
omega_cell = cell(freq_count,1);
obs_list_cell = cell(freq_count,1);
obsweight_list_cell = cell(freq_count,1);
inflows_list_cell = cell(freq_count,1);

num_test_cols = size(test_list,2);

if num_test_cols == 4 %monopole testing
    first_ingroup = 1;
    for f = 1:1:freq_count
        omega_cell{f} = test_list(first_ingroup,1);
        %Find the last entry with the same omega
        if (first_ingroup + 1) > num_totalobs
            last_ingroup = first_ingroup;
        else
            last_ingroup = first_ingroup + 1;
            found_mismatch = 0;
            while (last_ingroup <= num_totalobs) && (found_mismatch == 0)
                if test_list(last_ingroup,1) == test_list(first_ingroup,1)
                    last_ingroup = last_ingroup + 1;
                else
                    found_mismatch = 1;
                end
            end
            last_ingroup = last_ingroup - 1;
        end
        [group_obslocs,~,obslocs_map] = unique(test_list(first_ingroup:last_ingroup,4));
        [group_pumptests,~,pumptests_map] = unique(test_list(first_ingroup:last_ingroup,[2 3]),'rows');
        group_obsweight = well_weights(:,[group_obslocs]);
        group_inflows = well_weights(:,[group_pumptests(:,1)]);
        for i = 1:1:size(group_pumptests,1)
            group_inflows(:,i) = group_inflows(:,i).*group_pumptests(i,2);
        end
        group_obs_list = [pumptests_map obslocs_map];
        obs_list_cell{f} = group_obs_list;
        obsweight_list_cell{f} = group_obsweight;
        inflows_list_cell{f} = group_inflows;
        first_ingroup = last_ingroup + 1;
    end
elseif num_test_cols == 6 %dipole testing
    first_ingroup = 1;
    for f = 1:1:freq_count
        omega_cell{f} = test_list(first_ingroup,1);
        %Find the last entry with the same omega
        if (first_ingroup + 1) > num_totalobs
            last_ingroup = first_ingroup;
        else
            last_ingroup = first_ingroup + 1;
            found_mismatch = 0;
            while (last_ingroup <= num_totalobs) && (found_mismatch == 0)
                if test_list(last_ingroup,1) == test_list(first_ingroup,1)
                    last_ingroup = last_ingroup + 1;
                else
                    found_mismatch = 1;
                end
            end
            last_ingroup = last_ingroup - 1;
        end
        [group_obslocs,~,obslocs_map] = unique(test_list(first_ingroup:last_ingroup,6));
        [group_pumptests,~,pumptests_map] = unique(test_list(first_ingroup:last_ingroup,[2 3 4 5]),'rows');
        group_obsweight = well_weights(:,[group_obslocs]);
        group_inflows1 = well_weights(:,[group_pumptests(:,1)]);
        group_inflows2 = well_weights(:,[group_pumptests(:,3)]);
        group_inflows = zeros(size(group_inflows));
        for i = 1:1:size(group_pumptests,1)
            group_inflows(:,i) = group_inflows1(:,i).*group_pumptests(i,2) + ...
                group_inflows2(:,i).*group_pumptests(i,4);
        end
        group_obs_list = [pumptests_map obslocs_map];
        obs_list_cell{f} = group_obs_list;
        obsweight_list_cell{f} = group_obsweight;
        inflows_list_cell{f} = group_inflows;
        first_ingroup = last_ingroup + 1;
    end
else
    error('test_list must have either 4 columns (monopole tests) or 6 columns (dipole tests)');
end
    
experiment = struct('omega',omega_cell,...
    'tests',obs_list_cell,...
    'stims',inflows_list_cell, ...
    'obs',obsweight_list_cell);

