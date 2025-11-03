clc
close all

% Code implementation for the calculations used to validate Model 1


C = 50e-3;     % Supercapacitor value [F]
Resr = 160e-3; % Equivalent series resistance of the supercapacitor [ohms]
V = 2;         % Target voltage [V]
Vact = 1.2;    % Turn-on voltage required to activate the boost converter [V] 
Vmin = 1.02;   % Minimum voltage at which the boost converter shuts down [V]

v1 = 1.290; % Measured voltage (sample 2), distance = 50 cm - Ref. [14], Table IV
v2 = 1.850; % Measured voltage (sample 3), distance = 50 cm - Ref. [14], Table IV

teh1 = 200; % Measured time (sample 2), distance = 50 cm - Ref. Ref. [14], Table IV
teh2 = 450; % Measured time (sample 3), distance = 50 cm - Ref. Ref. [14], Table IV

Pdc1 = (C * v1^2) / (2 * teh1); % DC power calculation (Eq. 4.8) using sample 2 values
Pdc2 = (C * v2^2) / (2 * teh2); % DC power calculation (Eq. 4.8) using sample 3 values
Pdc_avg = (Pdc1 + Pdc2) / 2;    % Average DC power

taoc1 = -teh1 / log(1 - v1 / V); % Time constant taoc calculation (Eq. 4.9) using sample 2 values
taoc2 = -teh2 / log(1 - v2 / V); % Time constant taoc calculation (Eq. 4.9) using sample 3 values
taoc_avg = (taoc1 + taoc2) / 2;  % Average taoc calculation

teh_ini = - taoc_avg * log(1 - 1.999 / V); % Calculation of the initial EH phase time (Eq. 4.9), for Vc = 1.999 ≈ 2 V (V_target)
teh_ini_v = 0 : 0.001 : teh_ini;           % Time vector for the initial EH phase

vc = V*(1 - exp(-teh_ini_v / taoc_avg)); % Voltage curve of the EH initial phase (Eq. 4.9)

ttx_ini = 20.9869;               % Time of initial TX phase - average value from measurements - Ref. [14], Table V
ttx_ini_v = 0 : 0.001 : ttx_ini; % Time vector for the initial TX phase 

iout = 1.8e-3; % Average output current during TX phase [A] — Ref. [14], Table III

R2 = (V - Vmin) / iout - (ttx_ini / C) - Resr; % Equivalent resistance R2 (Eq. 4.12)

vd = V - (iout * (Resr + R2) + (iout * ttx_ini_v) ./ C); % Voltage curve of the initial TX phase (Eq. 4.12)

teh_op = 54.74;  % Measured time of the operational TEH phase [s] - Ref. [14], Table VI
ttx_op = 2.46;   % Measured time of the operational TTX phase [s] - Ref. [14], Table VI
teh_op_v = 0 : 0.001 : teh_op;   % Time vector for the operational TEH phase
ttx_op_v = 0 : 0.001 : ttx_op;   % Time vector for the operational TTX phase

vc1 = (V - Vmin) * (1 - exp(-teh_op_v ./ taoc_avg)) + vd(end); % Voltage curve of the first operational EH cycle (Eq. 4.17)

R2_new = (vc1(end) - vd(end)) / iout - (ttx_op / C) - Resr + (V * (1 -exp(-ttx_op / taoc_avg)))/iout;  % New equivalent resistance R2 (Eq. 4.19)

vd1 = vc1(end) - (iout * (Resr + R2_new) + (iout * ttx_op_v) / C) + V * (1 - exp(-ttx_op_v / taoc_avg)); % Voltage curve of the first operational TX cycle (Eq. 4.19)

V1 = [vc vd vc1 vd1]; % Voltage vector with initial EH/TX cycle and EH/TX cycle 1 

t1 = 0 : 0.001 : (teh_ini + ttx_ini + teh_op + ttx_op + 0.002); % Time vector with initial EH/TX cycle and EH/TX cycle 1 

V2 = [V1 vc1 vd1 vc1 vd1 vc1 vd1 vc1 vd1]; % Voltage vector with initial EH/TX cycle  and 5 EH/TX cycles 

t2 = 0 : 0.001 : (t1(end) + 4 * teh_op + 4 * ttx_op + 0.008); % Time vector with initial EH/TX cycle  and 5 EH/TX cycles 

% Result of initial EH phase
figure(1)
plot(teh_ini_v, vc, 'LineWidth', 2)   
grid on
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Time [s]')
ylabel('Voltage [V]')
title('Initial EH Phase Voltage Curve')

% Result of initial TX phase
figure(2)
plot([ttx_ini_v(1) ttx_ini_v],[V vd], 'LineWidth', 2)
grid on
set(gca,'FontSize',14,'FontWeight','bold')
xlabel('Time [s]')
ylabel('Voltage [V]')
title('Initial TX Phase Voltage Curve')

% Voltage and time vectors of measurement data
v_points = [1.29 1.85 2 1.7698 1.02 1.2 1.0819 1.02];
t_points = [200 450 teh_ini teh_ini teh_ini+20.9869 t1(end)-ttx_op t1(end)-ttx_op t1(end)];

% result with EH/TX cycles and measurements
figure(3)
plot (t2, V2, 'LineWidth', 2)
grid on
hold on
plot(t_points, v_points,'*', 'LineWidth', 2)
set(gca,'FontSize',14,'FontWeight','bold')
xlabel('Time [s]')
ylabel('Voltage [V]')
title('Initial EH/TX Phases and EH/TX Cycles of Voltage Curve')

%%% Ref. [9] model calculations %%%
vmax = 2;      % Máximum operation voltage [V]
vmin = 0;      % Minimum operation voltage [V]
cmin = 50e-3;  % Minimum capacitance [F]
pleak = 10e-6; % Leakage power [W]
pcons = pleak; % Consumption power [W]
Pharv = 488.2e-6;  % Harvested power
pmuh = 0.85;   % PMU efficiency
p_harv = Pharv*pmuh; % Effective harvested power [w]

% Charging time calculation [s] - initial EH phase
tcharge_ini = -(vmax^2 / pcons) * cmin * log((vmax - (vmax * p_harv) / pcons) / (vmin - (vmax * p_harv) / pcons))
% Time vector of the EH phase
tcharge_ini_v = 0:0.001 : tcharge_ini;

% Charging voltage calculation [V] - initial EH phase
vc_calc = ((vmax * p_harv) / pcons) * (1 - exp( -tcharge_ini / ((vmax^2 / pcons) * cmin)));
% Voltage vector of the EH phase
vc_calc_v = ((vmax * p_harv) / pcons) * (1 - exp( -tcharge_ini_v ./ ((vmax^2 / pcons) * cmin)));

% Result of initial EH phase 
figure(4)
plot(teh_ini_v,vc,'LineWidth', 2)
grid
hold on
plot([teh1 teh2],[v1 v2],'*','LineWidth', 2)
plot(tcharge_ini_v,vc_calc_v,'--','LineWidth', 2)
legend('Proposed model','Measurements','Model from [9]')
set(gca,'FontSize',14,'FontWeight','bold')
xlabel('Time [s]')
ylabel('Voltage [V]')
title('Initial EH Phase Voltage Curve')
