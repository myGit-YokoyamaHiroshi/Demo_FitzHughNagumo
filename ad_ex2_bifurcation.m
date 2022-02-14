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
dt     = 0.01;   % time step for numerical integration, unit: msec
time   = linspace(0, Nt-1, Nt)*dt; % time vector; unit : msec
%%%%% parameter settings
I      =  linspace(0, .5, 200);
tau    =  20;
a      = -0.3;
b      =  1.4;
% I      =  linspace(0, 4, 100);  % 
% tau    =  20;
% a      = -0.6;
% b      =  0.5;

X0     = [0, 0]; % initial value of state variables
                 % X0(1): membrane potential, v
                 % X0(2): recovery variable,  w
%%%%% parameter settings
%%
color_list  = turbo(6);
eqpt_labels = {};
eqpt_idx    = [];
V_list      = [];
I_list      = [];


cnt         = 1;
for i = 1:length(I)
    %%%%% Calculate equilibrium points
    [v_eq, w_eq] = solve_equilibria_FitzHughNagumo(I(i), tau, a, b);
    for j = 1:length(v_eq)
        X = [v_eq(j), w_eq(j)];
        %%%%%%% Get jacobian matrix at equilibrium point [v_eq(i), w_eq(i)]
        J = jacobian_matrix_FitzHughNagumo(X, I, tau, a, b);
        [eigvec, eigvalue] = eig(J);
        eigvalue    = diag(eigvalue);
        
        V_list(cnt) = v_eq(j);
        I_list(cnt) = I(i);
        %%%%%%% Determine its stability
        if all(imag(eigvalue)==0)
            if all(real(eigvalue)>0)
                labels        = 'Unstable node';
                eqpt_idx(cnt) = 1;
            elseif all(real(eigvalue)<0)
                labels        = 'Stable node';
                eqpt_idx(cnt) = 2;
            else
                labels        = 'Saddle node';
                eqpt_idx(cnt) = 3;
            end
        else
            if all(real(eigvalue)<0)
                labels        = 'Stable focus';
                eqpt_idx(cnt) = 4;
            elseif any(real(eigvalue)>0)
                labels        = 'Unstable focus';
                eqpt_idx(cnt) = 5;
            else
                labels        = 'Center (Hopf)';
                eqpt_idx(cnt) = 6;
            end
        end
        
        if cnt == 1
            eqpt_labels = {labels};
        else
            eqpt_labels = [eqpt_labels, {labels}];
        end
        cnt = cnt + 1;
    end
end
%% Determine regions where the system has a limit cycle
periodic = zeros(size(V_list));
for i = 1:length(I)
    idx = find(I_list==I(i));
    tmp = zeros(size(idx));
    for j = 1:length(idx)
        if contains(eqpt_labels{idx(j)}, 'Stable') %&& length(idx)~=1
            tmp(j) = 1;
        end
    end

    if sum(tmp)==0
        periodic(i) = 1;
    end
end

Iperi   = I(periodic==1);
Vmaxmin = zeros(length(Iperi), 2);
for i = 1:length(Iperi)
    X      = zeros(Nt, length(X0));
    X(1,:) = X0;
    
    for j = 2:Nt
        X_now  = X(j-1,:);
        %%%%% Numerical integral scheme with 4th order Runge Kutta method
        X(j,:) = runge_kutta(X_now, dt, @FitzHughNagumo, Iperi(i), tau, a, b);
    end
    
    Vmaxmin(i,:) = [min(X(:,1)), max(X(:,1))];
end

  
%%
fig = figure(1);
figure_setting(40, 20, fig)
%%%%%% Show the region where the system has a limit cycle
plot(Iperi, Vmaxmin, 'b', 'linewidth', 4,'HandleVisibility','off')
hold on
X  =[Iperi, fliplr(Iperi)];
Y = [Vmaxmin(:,1).', fliplr(Vmaxmin(:,2).')];
h = fill(X, Y, 'k','DisplayName','limit cycle');
set(h,'facealpha',.1)

%%%%%% Show the I-V curve 
for i = 1:6
    if sum(eqpt_idx==i) ~=0
        idx = find(eqpt_idx==i);
        scatter(I_list(idx), V_list(idx), ...
            'MarkerEdgeColor', color_list(i,:),...
            'MarkerFaceColor', color_list(i,:),...
            'DisplayName', eqpt_labels{idx(1)})

        legend('-DynamicLegend');
    end
end


xlabel('parameter \it I')
ylabel('Membrane potential V_*')
legend('location', 'southeast')

fname = [filepath, filesep, 'figures', filesep, 'ad_ex2', filesep, 'bifurcation'];
figure_save(fig, fname)
%%
