%%% simulates phase conjugation experiment in random medium
%%% Experiment 1: point source at focus, used to record forward propagating scattered light
%
%%% Experiment 2: phase conjugated propagation of result of Experiment 1
% --> forms sharp focus again
%%% Gerwin Osnabrugge 2015

clear all; close all;
addpath('..');

%% options for grid (gopt) and for simulation (sopt) 
PPW=4; %points per wavelength = lambda/h
sopt.lambda = 1; %in mu %lambda_0 = 1; %wavelength in vacuum (in um)
sopt.energy_threshold = 1E-16;
sopt.callback_interval = 10;
sopt.max_cycles = 6000;

mopt.lambda = sopt.lambda;
mopt.pixel_size = sopt.lambda/PPW;
mopt.boundary_widths = [24*PPW, 24*PPW]; %periodic boundaries
mopt.boundary_strength = 0.2;
mopt.boundary_type = 'PML3';
N = [128*PPW 128*PPW]; % size of medium (in pixels)

%% Construct random medium
% real refractive index
n0 = 1.3;          % mean
n_var = 0.1;     % variance

% randomly generate complex refractive index map
n_sample = 1.0*(n0 + n_var * randn(N));

% low pass filter to remove sharp edges
n_fft = fft2(n_sample);
cutoff = 1/16; % 1 is no cutoff frequency
window = [zeros(1,N(2)*(1/2-cutoff)), ones(1,N(2)*(2*cutoff)), zeros(1,N(2)*(1/2-cutoff))]' * ...
    	 [zeros(1,N(1)*(1/2-cutoff)), ones(1,N(1)*(2*cutoff)), zeros(1,N(1)*(1/2-cutoff))];
n_sample = real(ifft2(n_fft.*fftshift(window)));

% construct sample object
sample = SampleMedium(n_sample, mopt); 

%% Experiment 1: recording phase
source1 = sparse(N(1), N(2));
source1(end*3/4,end/2) = 1; % point source
sim = wavesim(sample, sopt);
E_recording = exec(sim, source1);

%% Experiment 2: playback phase
source2 = sparse(N(1), N(2));
source2(1,:) = conj(E_recording(1,:)); % conjugate result at top boundary
sim = wavesim(sample, sopt);
E_playback = exec(sim, source2);

%% plot resulting field amplitude
fig = figure(1); clf;
set(fig,'Position',get(fig,'Position')+[0,0, 600, 0]);

x = (-N(2)/2 + 1 : N(2)/2) /PPW;
y = (-N(1)/2 + 1 : N(1)/2) /PPW;

% plot recording phase
subplot(1,2,1);
imagesc(x,y,log(abs(E_recording)));
title('Recording Phase','FontSize',16);
axis square;
xlabel('x (\lambda)','FontSize',16);
ylabel('y (\lambda)','FontSize',16);
h = colorbar;
set(get(h,'Title'),'String','log|E|','FontSize',18,'FontName','Times New Roman');
set(gca,'FontSize',14);

% plot playback phase
subplot(1,2,2);
imagesc(x,y,log(abs(E_playback)));
axis square;
title('Playback Phase','FontSize',16);
xlabel('x (\lambda)','FontSize',16);
ylabel('y (\lambda)','FontSize',16);
h = colorbar;
set(get(h,'Title'),'String','log|E|','FontSize',18,'FontName','Times New Roman');
set(gca,'FontSize',14);

