clear 
close all
clc
filepath = fileparts(mfilename('fullpath'));
addpath(genpath(filepath));
%% code font settings
%%%% Set "Arial" as the Default font
set(0,'defaultAxesFontSize',14);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultTextFontSize',14);
set(0,'defaultTextFontName','Arial');

set(0,'defaultUipanelFontName','Arial');
set(0,'defaultUicontrolFontName','Arial');
%%
Nt     = 100000;  % Num. of sample
dt     = 0.01;   % time step for numerical integration, unit: msec
time   = linspace(0, Nt-1, Nt)*dt; % time vector; unit : msec
%%%%% parameter settings
I_list =  0.19:0.002:0.25;
tau    =  20;
a      = -0.3;
b      =  1.4;

X0     = [0, 0]; % initial value of state variables
                 % X0(1): membrane potential, v
                 % X0(2): recovery variable,  w
%%%%% parameter settings
%% Solve differential equation
freqs  = zeros(size(I_list));
X      = zeros(Nt, length(X0), length(I_list));

h = waitbar(0,'running');
for j = 1:length(I_list)
    X(1,:,j) = X0;
    I        = I_list(j);
    for i = 2:Nt
        X_now  = X(i-1,:,j);
        %%%%% Numerical integral scheme with 4th order Runge Kutta method
        X(i,:,j) = runge_kutta(X_now, dt, @FitzHughNagumo, I, tau, a, b);
    end
    
    pks      = findpeaks(X(:,1,j));
    Npks     = sum(pks>0.9*max(pks));
    T        = (Nt * dt) * 10^-3;
    freqs(j) = Npks/T; 
    
    %%%%%% Progress bar
    if mod(floor(j/length(I_list)*100), 5) == 0
        waitbar(j/length(I_list), h, ['Progress...', num2str(floor(j/length(I_list)*100)) , '%'])
    end
end

close(h)
%%
I_plt  = [0.19, 0.20, 0.21, 0.22, 0.23, 0.24];

fig = figure(1);
figure_setting(60, 70, fig);

sfh1 = subplot(2,4,2:3);
plot(I_list, freqs, '-o', 'linewidth', 1);
xlabel('\it I_{ext}')
ylabel('frequency (Hz)')
sfh1.Position = sfh1.Position + [0, 0, 0, -0.05];

xlim([I_list(1), I_list(end)])
ylim([0, 45])

for i = 1:length(I_plt)
    sfh2 = subplot(4,3,6+i);
    plot(time, X(:,1, I_list == I_plt(i)), 'LineWidth', 2);
    xlabel('time (ms)')
    ylabel(' \it V')
    xticks(0:200:Nt*dt)
    xlim([0, Nt*dt])
    title(['\it I_{ext} = ', num2str(I_plt(i))])

    sfh2.Position = sfh2.Position + [0, 0, 0, -0.02];
end
%%
fname = [filepath, filesep, 'figures', filesep, 'ex3', filesep, 'result'];
figure_save(fig, fname)