%%% Simulates the wave propagation of a point source in a random medium
%%% Gerwin Osnabrugge 2015

clear all; close all;
addpath('..');

%% options for grid (gopt) and for simulation (sopt) 
PPW=4; %points per wavelength = lambda/h
sopt.lambda = 1; %in mu %lambda_0 = 1; %wavelength in vacuum (in um)
sopt.energy_threshold = 1E-16;
sopt.callback_interval = 10;
sopt.max_cycles = 400;

mopt.lambda = sopt.lambda;
mopt.pixel_size = sopt.lambda/PPW;
mopt.boundary_widths = [0, 0, 0]; %periodic boundaries
mopt.boundary_strength = 0;
mopt.boundary_type = 'PML3';
N = [32*PPW 32*PPW, 32*PPW]; % size of medium (in pixels)

%% Construct random medium
% real refractive index
n0 = 1.3;          % mean
n_var = 0.1;     % variance

% imaginary refractive index
a0 = 0.05;       % mean
a_var = 0.02;    % variance

% randomly generate complex refractive index map
n_sample = 1.0*(n0 + n_var * randn(N)) + 1.0i*(a0 + a_var * randn(N));

% low pass filter to remove sharp edges
n_fft = fftn(n_sample);

W = @(n) [zeros(1,n*3/8), ones(1,n/4), zeros(1,n*3/8)];
window = bsxfun(@times, W(N(2))' * W(N(1)), reshape(W(N(3)), [1,1,N(3)]));
n_sample = ifftn(n_fft.*fftshift(window));
n_sample = max(real(n_sample), 1.0) + 1.0i * max(imag(n_sample), 0.0);

% construct sample object
sample = SampleMedium(n_sample, mopt); 

%% define a point source at the medium center
source = zeros(N(1), N(2), N(3));
source(end/2,end/2,end/2) = 1; % point source

%% wavesim simulation
sim = wavesim(sample, sopt);
iterations_per_wavelength(1) = sim.iterations_per_cycle;
E = exec(sim, source);

%% plot resulting field amplitude
figure(1); clf;

x = (-N(2)/2 + 1 : N(2)/2) /PPW;
y = (-N(1)/2 + 1 : N(1)/2) /PPW;
z = (-N(3)/2 + 1 : N(3)/2) /PPW;

slice(x,y,z, log(abs(E)), [0 8 16], 17, 17);
set(findobj(gca,'Type','Surface'),'EdgeColor','none')
axis square;
xlabel('x (\lambda)','FontSize',16);
ylabel('y (\lambda)','FontSize',16);
zlabel('z (\lambda)','FontSize',16);
h = colorbar;
set(get(h,'Title'),'String','log|E|','FontSize',18,'FontName','Times New Roman');
set(gca,'FontSize',14);
xlim([0 16]);

