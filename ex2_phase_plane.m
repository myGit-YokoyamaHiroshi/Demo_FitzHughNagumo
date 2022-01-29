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
Nt     = 50000;  % Num. of sample
fs     = 100;    % sampling frequency 
dt     = 1/fs;   % time step for numerical integration

time   = linspace(0, Nt-1, Nt)/fs; % time vector; unit : seconds
%%%%% parameter settings
I      =  0.23;  % 
tau    =  20;
a      = -0.3;
b      =  1.4;
X0     = [0, 0]; % initial value of state variables
                 % X0(1): membrane potential, v
                 % X0(2): recovery variable,  w
%%%%% parameter settings
%% Solve differential equation
X_trj      = zeros(Nt, length(X0));% Set inital value
X_trj(1,:) = X0;

for i = 2:Nt
    X_now  = X_trj(i-1,:);
    %%%%% Numerical integral scheme with 4th order Runge Kutta method
    X_trj(i,:) = runge_kutta(X_now, dt, @FitzHughNagumo, I, tau, a, b);
end
%% Get Nullcline
[v_null, w_null] = get_nullcline_FizHughNagumo(I, tau, a, b, -1.2, 1.2);
%% Calculate vector field
v = linspace(-1.2, 1.2, 60);
w = linspace(-1.2, 1.2, 60);
%%%% Make 2-D grid space
[V, W] = meshgrid(v,w);

dV = zeros(size(V));
dW = zeros(size(W));
for i = 1:numel(V)
    %%%% Calculate vector at the grid point [V(i), W(i)] 
    Fprime = FitzHughNagumo([V(i), W(i)], I, tau, a, b);
    dV(i)  = Fprime(1);
    dW(i)  = Fprime(2);
end
norm = sqrt(dV.^2 + dW.^2);
dV   = dV./norm;
dW   = dW./norm;
%% Determine the stability of each equilibrium points
%%%%% Calculate equilibrium points
[v_eq, w_eq] = salve_equilibria_FitzHughNagumo(I, tau, a, b);
eqpt_labels  = cell(1, length(v_eq));
for i = 1:length(v_eq)
    X = [v_eq(i), w_eq(i)];
    %%%%%%% Get jacobian matrix at equilibrium point [v_eq(i), w_eq(i)]
    J = jacobian_matrix_FitzHughNagumo(X, I, tau, a, b);
    [eigvec, eigvalue] = eig(J);
    eigvalue = diag(eigvalue);
    
    %%%%%%% Determine its stability
    if all(imag(eigvalue)==0)
        if all(real(eigvalue)>0)
            eqpt_labels{i} = 'Unstable node';
        elseif all(real(eigvalue)<0)
            eqpt_labels{i} = 'Stable node';
        else
            eqpt_labels{i} = 'Saddle node';
        end
    else
        if all(real(eigvalue)<0)
            eqpt_labels{i} = 'Stable focus';
        elseif any(real(eigvalue)>0)
            eqpt_labels{i} = 'Unstable focus';
        else
            eqpt_labels{i} = 'Center (Hopf)';
        end
    end
end
%% plot results
fig = figure(1);
figure_setting(60, 60, fig)

sfh1 = subplot(2,1,1,'parent', fig);
plot(time, X_trj(:,1), 'LineWidth', 3);
hold on
plot(time, X_trj(:,2), 'LineWidth', 3);
hold off

xlabel('time (s)')
ylabel('v, w')
lgnd = legend({'membrane potential \it v', 'recovery variable \it w'}, 'location', 'northeastoutside');

%%% plot trajectory
sfh2 = subplot(2,1,2,'parent', fig);
plot(X_trj(:,1), X_trj(:,2), 'r', 'LineWidth', 3);
trj_labels = {'trajectory'};

hold on
%%% plot vector fields
quiver(V, W, dV, dW, 'color', [0.6, 0.6, 0.6],'HandleVisibility','off'); 

%%% plot nullcline
plot(v_null, w_null, 'k-', 'LineWidth', 1,'HandleVisibility','off'); 

for i = 1:length(v_eq)
    %%% plot equilibrium point
    h_eq(i) = scatter(v_eq(i), w_eq(i), 100,'filled');
end
xlabel('membrane potential \it v')
ylabel('recovery variable \it w')
legend([trj_labels, eqpt_labels], 'location', 'northeastoutside')
axis tight equal;

title('phase plane')
xlim([-1.2, 1.2])
ylim([-1.0, 1.0])
alpha(0.8)
sfh2.Position = sfh2.Position - [0.05, 0, 0, 0];

fname = [filepath, filesep, 'figures', filesep, 'ex2', filesep, 'stability'];
figure_save(fig, fname)
hold off
