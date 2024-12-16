%% Comprehensive Solar Plant Simulation with SCADA-Based Cyber Attack (Averaged over 1000 Runs)

%% Function to Add Generator Information
% This function dynamically adds a generator to the system's generator matrix.
% It assigns key parameters such as active power (Pg), reactive power (Qg),
% and limits for power generation, among others.
function mpc = add_generator(mpc, bus, Pg, Qg, Qmax, Qmin, Vg, Pmax, Pmin)
    new_gen = zeros(1, size(mpc.gen, 2)); % Initialize new generator row with zeros
    new_gen(1) = bus;   % Bus number where the generator is connected
    new_gen(2) = Pg;    % Active power generation
    new_gen(3) = Qg;    % Reactive power generation
    new_gen(4) = Qmax;  % Maximum reactive power
    new_gen(5) = Qmin;  % Minimum reactive power
    new_gen(6) = Vg;    % Voltage magnitude setpoint
    new_gen(7) = 100;   % Base MVA of the machine
    new_gen(8) = 1;     % Machine status (1 = online, 0 = offline)
    new_gen(9) = Pmax;  % Maximum active power output
    new_gen(10) = Pmin; % Minimum active power output
    mpc.gen = [mpc.gen; new_gen]; % Append the new generator to the generator matrix
end

%% Initialization
clear; close all; clc;

% Load MATPOWER test case
mpc = loadcase('case141');

% Add buses for solar plant arrays
solar_buses = 119:128; % Create 10 solar buses
num_panels_per_bus = 1e5; % Each bus represents 100,000 solar panels

for i = 1:length(solar_buses)
    if ~ismember(solar_buses(i), mpc.bus(:, 1))
        % Append a new bus with default parameters if it doesn't already exist
        mpc.bus = [mpc.bus; solar_buses(i), 1, 0, 0, 0, 0, 1, 1.0, 0, 345, 1, 1.1, 0.9];
    end
end

% Add solar generators to the newly created buses
panel_output =  0.00007; % Base panel output in MW
for i = 1:length(solar_buses)
    % Total output of 100,000 panels per bus
    total_output_per_bus = panel_output * num_panels_per_bus;
    mpc = add_generator(mpc, solar_buses(i), total_output_per_bus, 0, total_output_per_bus, 0, 1.0, total_output_per_bus, 0);
end

% Add cost data for new generators
new_gencost = [2, 0, 0, 3, 0.02, 2, 0]; % Quadratic cost model parameters
for i = 1:length(solar_buses)
    mpc.gencost = [mpc.gencost; new_gencost]; % Append cost data for each solar generator
end

%% Simulation Parameters
n_runs = 500; % Number of simulation runs
hours = 1:24; % Time period for the simulation (24 hours)
solar_profile = [0, 0, 10, 30, 60, 100, 120, 150, 200, 300, 400, 380, 350, 300, 200, 100, 50, 0]; 

% Convert solar profile to MW by considering random fluctuations for real-world simulation
% Each panel generates between 50% and 150% of its nominal output
random_factors = @(n) 0.5 + rand(1, n); % Random scaling factors (50% to 150%)
array_profile = solar_profile * panel_output * num_panels_per_bus * length(solar_buses) * 1e6 / 1e3; % Convert to MW

% Identify generator indices corresponding to solar buses
[~, solar_gen_indices] = ismember(solar_buses, mpc.gen(:, 1));

% Define daylight hours (6 AM to 6 PM)
daylight_hours = 6:18;  % Hours between 6 AM and 6 PM (inclusive)

% Set array profile to 0 outside daylight hours
full_array_profile = zeros(1, 24); % Default to zero
full_array_profile(daylight_hours) = array_profile(daylight_hours); % Apply profile only during daylight hours

n_steps = length(full_array_profile); % Number of simulation steps (24 hours)

% Initialize accumulators for averaging
power_log_original_avg = zeros(1, n_steps);
power_log_attacked_avg = zeros(1, n_steps);

%% Simulation Parameters for Power Flow
mpopt = mpoption('pf.tol', 1e-4, 'verbose', 0); % Increase tolerance and disable verbosity for speed

%% Run Simulations
for run = 1:n_runs
    disp("Run#" + run);
    % Generate a random attack probability for each run
    attack_probability = 0.1 + 0.7 * rand(); % Randomly vary between 0.2 and 0.8

    % Temporary logs for each run
    power_log_original = zeros(1, n_steps);
    power_log_attacked = zeros(1, n_steps);

    for t = 1:n_steps
        % Apply randomized solar generation profile (only during daylight hours)
        if t >= daylight_hours(1) && t <= daylight_hours(end)
            % Generate random power values for solar generators only
            random_power_values = full_array_profile(t) .* (0.5 + rand(length(solar_gen_indices), 1));
            mpc.gen(solar_gen_indices, 2) = random_power_values;
        else
            mpc.gen(solar_gen_indices, 2) = 0; % No solar generation outside daylight hours
        end

        % Run power flow for normal operations
        try
            results_original = runpf(mpc, mpopt);
            power_log_original(t) = sum(results_original.gen(:, 2)); % Aggregate generator outputs
        catch
            power_log_original(t) = 0; % Assign default value in case of failure
        end

        % Simulate cyber-attack by randomly turning generator on and off
        attack_flag = rand() < attack_probability; 

        if attack_flag
            % Ensure some generation remains online
            mpc.gen(solar_gen_indices, 2) = zeros(size(solar_gen_indices)); % Turn off solar generators
        end

        % Run power flow for attacked scenario
        try
            results_attacked = runpf(mpc, mpopt);
            power_log_attacked(t) = sum(results_attacked.gen(:, 2)); % Aggregate generator outputs
        catch
            power_log_attacked(t) = 0; % Assign default value in case of failure
        end
    end

    % Accumulate results for averaging
    power_log_original_avg = power_log_original_avg + power_log_original;
    power_log_attacked_avg = power_log_attacked_avg + power_log_attacked;
end

% Calculate averages
power_log_original_avg = power_log_original_avg / n_runs;
power_log_attacked_avg = power_log_attacked_avg / n_runs;

%% Plotting SCADA Simulation Results (Combined)
figure;

plot(hours, power_log_original_avg, '-b', 'LineWidth', 2, 'DisplayName', 'Normal Output');
hold on;
plot(hours, power_log_attacked_avg, '-r', 'LineWidth', 2, 'DisplayName', 'Attacked Output');
hold off;

xlabel('Time (Hours)');
ylabel('MW');
title('Averaged Solar Generation Outputs (Normal vs. Attacked)');
legend('show');
grid on;
ylim([0, max([power_log_original_avg, power_log_attacked_avg])]); % Ensure y-axis starts at 0
