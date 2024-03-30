% channel data with k-wave
clear; clc; close all;
dir_ = pwd;
dir_kwave = [dir_ '/k-wave-toolbox-version-1/k-Wave'];
addpath(dir_kwave);

% simulation settings
DATA_CAST       = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
%% Define acoustic parameters & transducers
% acoustic parameters
Acoustic.fc = 6.5e6;
Acoustic.C = 1540;
Acoustic.density = 1000;
Acoustic.lambda = Acoustic.C / Acoustic.fc;

% transducer
Transducer.num_element = 128;
Transducer.ele_pitch = 430e-6;
Transducer.ele_width = 420e-6;
Transducer.ele_height = 5e-3;
Transducer.ele_kerf = Transducer.ele_pitch - Transducer.ele_width;
Transducer.elevation_focus = 25e-3;
Transducer.transmit_fnum = 3;

%% Define domain
size_x = 72e-3; % depth
size_y = 55e-3; % lateral
size_z = 80e-3; % elevational

offset = 0.5*[size_x, 0, 0];

points_per_wavelength = 4;
delta_x = Acoustic.C / (points_per_wavelength * Acoustic.fc);
delta_y = delta_x;
delta_z = delta_x;
delta_ = [delta_x, delta_y, delta_z];

Nx = round(size_x/delta_x);
Ny = round(size_y/delta_y);
Nz = round(size_z/delta_z);

kgrid = kWaveGrid(Nx, delta_x, Ny, delta_y, Nz, delta_z);

% temporal grid
Acoustic.fs = 4*Acoustic.fc;
round_trip_time = 2*size_x / Acoustic.C;

kgrid.dt = 1/Acoustic.fs;
kgrid.t_array = 0:kgrid.dt:round_trip_time;
kgrid.Nt = numel(kgrid.t_array);

PML_Size = round(4*Acoustic.lambda/delta_z);  % size of the PML layers in grid points
%% Define medium
medium.sound_speed = Acoustic.C * ones(Nx, Ny, Nz);
medium.density = Acoustic.density * ones(Nx, Ny, Nz);

mask_ = zeros(Nx, Ny, Nz);
point_target_pos = [20e-3, 0, 0; ...% x, y, z (depth, lateral, elevation)
                    30e-3, 0, 0]; 
% point_target_pos = point_target_pos - offset;
point_target_spl = round(point_target_pos ./ repmat(delta_, size(point_target_pos,1), 1));
point_target_spl(:,2) = point_target_spl(:,2) + round(0.5*Ny);
point_target_spl(:,3) = point_target_spl(:,3) + round(0.5*Nz);
for p_idx = 1:size(point_target_spl,1)
    tmp = point_target_spl(p_idx, :);
    mask_(tmp(1), tmp(2), tmp(3)) = 1;
end

medium.density(logical(mask_)) = 2*Acoustic.density;
                
%% Define input signal
source_strength = 1e6; % [MPa]
num_cycles = 5;

input_signal = toneBurst(Acoustic.fs, Acoustic.fc, num_cycles);
% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength ./ (Acoustic.C * Acoustic.density)) .* input_signal;

%% Define k-wave transducer
transducer.number_elements = Transducer.num_element;
transducer.element_width = round(Transducer.ele_width/delta_y);
transducer.element_length = round(Transducer.ele_height/delta_z); % element height
transducer.element_spacing = round(Transducer.ele_kerf/delta_y);
transducer.radius = inf;

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, ... % depth
                            Ny/2 - transducer_width/2, ... % lateral
                            Nz/2 - transducer.element_length/2]); % elevation

% properties used to derive the beamforming delays
transducer.sound_speed = Acoustic.C;                    % sound speed [m/s]
transducer.focus_distance = 30e-3;              % focus distance [m]
transducer.elevation_focus_distance = 25e-3;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
number_active_elements = round(transducer.focus_distance/Transducer.transmit_fnum/Transducer.ele_pitch);
transducer.active_elements = ones(transducer.number_elements, 1);

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

%%
% set the input settings
input_args = {...
    'PMLInside', false, 'PMLSize', PML_Size, ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};

sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});

