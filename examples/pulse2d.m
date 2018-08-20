%%% Simulate the pulse propagation of a point source with a specified bandwidth 
%%% in a 2D disordered medium

clear all; close all;
addpath('..');

%% options for grid (gopt) and for simulation (sopt) 
PPW=4;                          %points per wavelength = lambda/h
sopt.lambda = 1;                %in mu %lambda_0 = 1; %wavelength in vacuum (in um)
sopt.energy_threshold = 1E-6;
sopt.callback_interval = 1000;
sopt.singlePrecision = true;

mopt.lambda = sopt.lambda;
mopt.pixel_size = sopt.lambda/PPW;
mopt.boundary_widths = [30*PPW, 30*PPW]; % absorbing boundaries
mopt.boundary_strength = 0.2;
mopt.boundary_type = 'PML3';
N = [64*PPW 64*PPW]; % size of medium (in pixels)

%% Construct random medium
% real refractive index
n0 = 1;          % mean
n_var = 0.02;    % variance

% imaginary refractive index
a0 = 0;       % mean
a_var = 0;    % variance

% randomly generate complex refractive index map
n_sample = 1.0*(n0 + n_var * randn(N)) + 1.0i*(a0 + a_var * randn(N));

% construct sample object
sample = SampleMedium(n_sample, mopt); 

%% define a point source with spectrum (800-1200nm) at the medium center
% source frequencies
lambda = 1;                                 % center wavelength (in um)
Nfreq = 40;                                 % number of independent frequencies simulated
lambda_set = linspace(0.8,1.2,Nfreq)*lambda;% source bandwidth

source = zeros(N(1), N(2),Nfreq);
source(end/2,end/2,:) = gausswin(Nfreq);    % point source with Gaussian spectrum

%% wavesim simulation for all source frequencies
% Preallocate data for fields
E_set = zeros(N(1), N(2), Nfreq);

% run simulations
for f = 1:Nfreq
    sopt.lambda = lambda_set(f);
    sim = wavesim(sample, sopt);
    E_set(:,:,f) = exec(sim, source(:,:,f));
end

%% field in time domain
% Inverse Fourier transform E(frequency) -> E(time)
E_time = ifftshift(ifft(E_set,Nfreq,3),3);

% set axes
x = (-N(2)/2 + 1 : N(2)/2) /PPW;
y = (-N(1)/2 + 1 : N(1)/2) /PPW;
f = 3e8*2*pi./(lambda_set*10^-6);                  % frequency: f = ck0 [in Hz]
t_set = (-Nfreq/2+1:Nfreq/2)* 1/( 2*(f(1)-f(2)) );

% only store signal where time is positive
E_time = E_time(:,:,t_set >= 0);
t_set = t_set(t_set >= 0) * 10^15;       % in fs

% plot all fields as funcion of time
figure(2); clf;

for t = 1:Nfreq/2
    figure(2); imagesc(x,y,abs(E_time(:,:,t)).^2);
    title(['t = ',num2str(t_set(t)),' fs']);
    xlabel('x / \lambda','FontSize',16);
    ylabel('y / \lambda','FontSize',16);
    h = colorbar;
    set(get(h,'Title'),'String','log|E|','FontSize',18,'FontName','Times New Roman');
    set(gca,'FontSize',16);
    pause(0.1);    
end