%%% Simulates the propagation of a plane wave through a 1D homogeneous medium
%%% Analyzes the accuracy of wavesim and PSTD (with varying time steps)
%%% Gerwin Osnabrugge 2015

%close all;
addpath('..');

%% options for grid (gopt) and for simulation (sopt) 
PPW=4; %points per wavelength = lambda/h
sopt.lambda = 1; %in mu %lambda_0 = 1; %wavelength in vacuum (in um)
sopt.callback_interval = 25;
sopt.gpu_enabled = false; % only cuda cards supported

dt_relative_range = (1/2).^(0.5:0.5:3); % PSTD timestep

mopt.lambda = sopt.lambda;
mopt.pixel_size = sopt.lambda/PPW;
mopt.boundary_widths = [0, 25*PPW];
mopt.boundary_strength = 0.2;
mopt.boundary_type = 'PML3';
N = [1 round(50*PPW)]; % size of medium (in pixels)
n1 = 1; %refractive index medium

%% define a plane wave source and homogeneous medium
source = sparse(N(1), N(2));
source_pos = 1;%:PPW;%3*PPW;
source(:,source_pos) = 1; % plane wave source
sample = SampleMedium(n1*ones(N), mopt);

%% reserve space for output data
relative_error = zeros(1, length(dt_relative_range)+1);
iterations_per_wavelength = zeros(1, length(dt_relative_range)+1);

%% Run wavesim simulations
sim = wavesim(sample, sopt);
iterations_per_wavelength(1) = sim.iterations_per_cycle;
[E_wavesim, state] = exec(sim, source);

%% calculate exact solution analytically
k0 = n1*2*pi/sopt.lambda;
E_theory=homogeneous_medium_analytic_solution(2*pi/sopt.lambda*n1, mopt.pixel_size, sim.x_range-(source_pos-1)*mopt.pixel_size);

% calculate relative error
difference=E_wavesim(round(N(1)/2),:)-E_theory;
relative_error(1)=mean2(abs(difference).^2) / mean2(abs(E_theory).^2);

%% Run PSTD simulations
for k = 1:length(dt_relative_range)
    sopt.dt_relative = dt_relative_range(k);
    sim_PSTD = PSTD(sample, sopt);
    iterations_per_wavelength(k+1) = sim_PSTD.iterations_per_cycle;
    [E_PSTD, state] = exec(sim_PSTD, source);
    difference=E_PSTD(round(N(1)/2),:)-E_theory;
    relative_error(k+1)=mean2(abs(difference).^2) / mean2(abs(E_theory).^2);
end

%% Plot relative errors
figure(1); clf; 
loglog(iterations_per_wavelength(1), relative_error(1), 's', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
hold on;
loglog(iterations_per_wavelength(2:end),  relative_error(2:end), '+r','MarkerSize', 10, 'LineWidth', 2.0);
grid on;
legend('Modified Born','PSTD', 'Location', 'NorthWest');
title('Accuracy in 1D-homogeneous medium','FontSize',14);
xlabel('Iterations per wavelength','FontSize',14);
ylabel('Relative error','FontSize',14);
set(gca,'FontSize',14);
xlim([0.1, iterations_per_wavelength(end)*2]);
ylim([1E-12, 10]);
