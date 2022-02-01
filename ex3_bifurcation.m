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
I      =  linspace(0, 1, 100);
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
%%
fig = figure(1);
figure_setting(40, 20, fig)
hold on
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

fname = [filepath, filesep, 'figures', filesep, 'ex3', filesep, 'bifurcation'];
figure_save(fig, fname)
%%
