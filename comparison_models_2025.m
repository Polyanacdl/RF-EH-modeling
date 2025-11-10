clc
close all

%% Code implementation for the calculations used to validate the baseline model

% Authors: [Lacerda, P. C.; Mariano, A. A.; Brante, G.]
% Date: [nov - 11 - 2025]

%% --- Parameters ---
C = 50e-3;     % Supercapacitor value [F]
Resr = 160e-3; % Equivalent series resistance [Ω]
V = 2;         % Target voltage [V]
Vact = 1.2;    % Boost converter turn-on voltage [V] 
Vmin = 1.02;   % Boost converter shutdown voltage [V]

%% --- Measurement data (Ref. [14] Table IV) ---
v1 = 1.290; teh1 = 200; % Sample 2 (50 cm)
v2 = 1.850; teh2 = 450; % Sample 3 (50 cm)

%% --- DC power (Eq. 4.8) ---
Pdc1 = (C * v1^2) / (2 * teh1); 
Pdc2 = (C * v2^2) / (2 * teh2); 
Pdc_avg = (Pdc1 + Pdc2) / 2;    

%% --- Time constant (Eq. 4.9) ---
taoc1 = -teh1 / log(1 - v1 / V); 
taoc2 = -teh2 / log(1 - v2 / V); 
taoc_avg = (taoc1 + taoc2) / 2;  

%% --- Initial EH phase ---
teh_ini = - taoc_avg * log(1 - 1.999 / V); % Time until Vc = 1.999 ≈ 2 V (V_target) (Eq. 4.9)
teh_ini_v = 0 : 0.001 : teh_ini;           
vc = V*(1 - exp(-teh_ini_v / taoc_avg)); % Charging curve

%% --- Initial TX phase ---
ttx_ini = 20.9869;               % Ref. [14], Table V
ttx_ini_v = 0 : 0.001 : ttx_ini; 
iout = 1.8e-3;                   % Avg. output current during TX phase [A] — Ref. [14], Table III
R2 = (V - Vmin) / iout - (ttx_ini / C) - Resr; % Equivalent resistance R2 (Eq. 4.12)
vd = V - (iout * (Resr + R2) + (iout * ttx_ini_v) ./ C);  

%% --- Operational cycles of EH/TX phases ---
teh_op = 54.74;  ttx_op = 2.46;  % Measured [s] - Ref. [14], Table VI
teh_op_v = 0 : 0.001 : teh_op;   
ttx_op_v = 0 : 0.001 : ttx_op;   

taoc_op = -teh_op/log(1 - ((Vact - vd(end))/(V - Vmin))); % (Eq. 4.15)
vc1 = (V - Vmin) * (1 - exp(-teh_op_v ./ taoc_op)) + vd(end); % (Eq. 4.17)
R2_new = (vc1(end) - vd(end)) / iout - (ttx_op / C) - Resr + (V * (1 -exp(-ttx_op / taoc_op)))/iout;  % (Eq. 4.19)
vd1 = vc1(end) - (iout * (Resr + R2_new) + (iout * ttx_op_v) / C) + V * (1 - exp(-ttx_op_v / taoc_op)); % (Eq. 4.19)

%% --- Composite time and voltage vectors ---
V1 = [vc vd vc1 vd1]; 
t1 = 0 : 0.001 : (teh_ini + ttx_ini + teh_op + ttx_op + 0.002);  
V2 = [V1 vc1 vd1 vc1 vd1 vc1 vd1 vc1 vd1]; 
t2 = 0 : 0.001 : (t1(end) + 4 * teh_op + 4 * ttx_op + 0.008); 

%% --- Plot: Initial EH phase ---
figure(1)
plot(teh_ini_v, vc, 'LineWidth', 2)   
grid on
xlabel('Time [s]'); ylabel('Voltage [V]')
title('Initial EH Phase Voltage Curve')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

%% --- Plot: Initial TX phase ---
figure(2)
plot([ttx_ini_v(1) ttx_ini_v],[V vd], 'LineWidth', 2)
grid on
xlabel('Time [s]'); ylabel('Voltage [V]')
title('Initial TX Phase Voltage Curve')
set(gca,'FontSize',14,'FontWeight','bold')

%% --- Plot: EH/TX cycles with measurements ---
v_points = [1.29 1.85 2 1.7698 1.02 1.2 1.0819 1.02];
t_points = [200 450 teh_ini teh_ini teh_ini+20.9869 t1(end)-ttx_op t1(end)-ttx_op t1(end)];

figure(3)
plot (t2, V2, 'LineWidth', 2)
grid on
hold on
plot(t_points, v_points,'*', 'LineWidth', 2)
legend('Proposed model', 'Measurements')
xlabel('Time [s]'); ylabel('Voltage [V]')
title('Initial EH/TX Phases and EH/TX Cycles of Voltage Curve')
set(gca,'FontSize',14,'FontWeight','bold')

%% --- Comparison with Ref. [9] model ---
vmax = 2;            % Máximum operation voltage [V]
vmin = 0;            % Minimum operation voltage [V]
cmin = 50e-3;        % Minimum capacitance [F]
pleak = 10e-6;       % Leakage power [W]
pcons = pleak;       % Consumption power [W]
Pharv = 488.2e-6;    % Harvested power
pmuh = 0.85;         % PMU efficiency
p_harv = Pharv*pmuh; % Effective harvested power [w]

% Charging 
tcharge_ini = -(vmax^2 / pcons) * cmin * log((vmax - (vmax * p_harv) / pcons) / (vmin - (vmax * p_harv) / pcons));
tcharge_ini_v = 0:0.001 : tcharge_ini;

vc_calc = ((vmax * p_harv) / pcons) * (1 - exp( -tcharge_ini / ((vmax^2 / pcons) * cmin)));
vc_calc_v = ((vmax * p_harv) / pcons) * (1 - exp( -tcharge_ini_v ./ ((vmax^2 / pcons) * cmin)));

% Plot comparison
figure(4)
plot(teh_ini_v,vc,'LineWidth', 2)
grid
hold on
plot([teh1 teh2],[v1 v2],'*','LineWidth', 2)
plot(tcharge_ini_v,vc_calc_v,'--','LineWidth', 2)
legend('Proposed model','Measurements','Model from [9]')
xlabel('Time [s]'); ylabel('Voltage [V]')
title('Initial EH Phase Voltage Curve')
set(gca,'FontSize',14,'FontWeight','bold')

% Discharging
tdisch_ini = ttx_ini;          
pmul = 0.92;                   % PMU efficiency
pcons2 = 5.9e-3/pmul +  pleak; % Consumed power during discharge - Ref. Leemput (eq. 7)
vd_calc = (vmax * (exp(-tdisch_ini/((vmax^2/pcons2)*cmin)))); 

tdisch_ini_v = ttx_ini_v;  
vd_calc_v = (vmax * (exp(-tdisch_ini_v./((vmax^2/pcons2)*cmin))));
t3 = teh_ini(end):0.001:teh_ini(end)+tdisch_ini; 
v3 = vd_calc_v;  

% Plot comparison
figure(5)
plot([ttx_ini_v(1) ttx_ini_v],[V vd],'-','LineWidth', 2)
grid on
hold on
plot(tdisch_ini_v,vd_calc_v,'--','LineWidth', 2)
set(gca,'FontSize',14,'FontWeight','bold')
legend('Proposed model','Model from [9]')
title('Initial TX Phase Voltage Curve')

% Cycles 
tcharge_op = teh_op;  % Measurement
pcons3 = 150e-06;    % Defined parameter related to consumed power during charging

% Charging 
vc_calc_op = ((1.2*p_harv)/pcons3*(1-exp(-tcharge_op/((1.2^2/pcons3)*cmin)))+vd_calc*(exp(-tcharge_op/((1.2^2/pcons3)*cmin))));
tcharge_op_v = 0:0.001:tcharge_op; 
vc_calc_op_v = ((1.2*p_harv)/pcons3*(1-exp(-tcharge_op_v./((1.2^2/pcons3)*cmin)))+vd_calc*(exp(-tcharge_op_v./((1.2^2/pcons3)*cmin))));

% Composite time and voltage vectors
t4 = t3(end):0.001:t3(end)+tcharge_op; 
v4 = vc_calc_op_v; 

% Discharging
tdisch_op = ttx_op; 
pcons4 = 5.9e-3/pmul +  pleak; 
vd_calc_op = ((1.2*p_harv)/pcons4*(1-exp(-tdisch_op/((1.2^2/pcons4)*cmin)))+vc_calc_op*(exp(-tdisch_op/((1.2^2/pcons4)*cmin))));

tdisch_op_v = 0:0.001:tdisch_op; 
vd_calc_op_v = ((1.2*p_harv)/pcons4*(1-exp(-tdisch_op_v./((1.2^2/pcons4)*cmin)))+vc_calc_op*(exp(-tdisch_op_v./((1.2^2/pcons4)*cmin))));

t4 = t3(end):0.001:t3(end)+5*tcharge_op+5*tdisch_op+0.009; 
v4 = [vc_calc_op_v vd_calc_op_v vc_calc_op_v vd_calc_op_v vc_calc_op_v vd_calc_op_v vc_calc_op_v vd_calc_op_v vc_calc_op_v vd_calc_op_v];

% Plot comparison
figure(6)
plot(t2,V2,'-','LineWidth', 2)
grid on
hold on
plot(t_points, v_points,'*','LineWidth', 2)
plot(tcharge_ini_v,vc_calc_v,'--', 'LineWidth', 2)
plot(t3,v3, '--', 'LineWidth', 2)
plot(t4,v4, '--', 'LineWidth', 2)
legend('Proposed Model','Measurements','Model from [9]')
set(gca,'FontSize',14,'FontWeight','bold')
title('Initial EH/TX Phases and EH/TX Cycles of Voltage Curve')
