clear 
close all
clc
filepath = fileparts(mfilename('fullpath'));
addpath(genpath(filepath));
%% code font settings
%%%% Set "Arial" as the Default font
set(0,'defaultAxesFontSize',16);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultTextFontSize',16);
set(0,'defaultTextFontName','Arial');

set(0,'defaultUipanelFontName','Arial');
set(0,'defaultUicontrolFontName','Arial');
%%
Nt     = 100000;  % Num. of sample
dt     = 0.01;   % time step for numerical integration, unit: msec
time   = linspace(0, Nt-1, Nt)*dt; % time vector; unit : msec
%%%%% parameter settings
I_list =  linspace(0.19, 0.25, Nt);
tau    =  20;
a      = -0.3;
b      =  1.4;

X0     = [0, 0]; % initial value of state variables
                 % X0(1): membrane potential, v
                 % X0(2): recovery variable,  w
%%%%% parameter settings
%% Solve differential equation
X      = zeros(Nt, length(X0));
X(1,:) = X0;

for i = 2:Nt
    X_now  = X(i-1,:);
    I      = I_list(i);
    %%%%% Numerical integral scheme with 4th order Runge Kutta method
    X(i,:) = runge_kutta(X_now, dt, @FitzHughNagumo, I, tau, a, b);
end
%%
fig = figure(1);
figure_setting(50, 30, fig);

%%%%%%%
sfh1 = subplot(1,2,1,'parent', fig);
plot(X(:,1), X(:,2), 'r', 'LineWidth', 1);
xlabel('membrane potential \it V')
ylabel('recovery variable \it w')
title('phase space')
axis square
sfh1.Position = sfh1.Position + [0.02, 0.07, -0.05, -0.05,];

sfh2 = subplot(2,3,3,'parent', fig);
plot(time, X(:,1), 'LineWidth', 2);
ylabel(' \it V')
xticks(0:200:Nt*dt)
xticklabels('')
sfh2.Position = sfh2.Position + [-.15, -0.13, 0.2, 0];

sfh3 = subplot(2,3,6,'parent', fig);
plot(time, I_list, 'k');
xlabel('time (ms)')
ylabel(' \it I_{ext}')
xticks(0:200:Nt*dt)
sfh3.Position = sfh3.Position + [-.15, 0.24, 0.2, 0.00];
sfh3.Position(3) = sfh2.Position(3);
sfh3.Position(4) = 0.1;
%%
fname = [filepath, filesep, 'figures', filesep, 'ex2', filesep, 'result'];
figure_save(fig, fname)