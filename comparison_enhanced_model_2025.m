clc
close all

%% Code implementation for the calculations used to validate the enhanced model

% Authors: [Lacerda, P. C.; Mariano, A. A.; Brante, G.]
% Date: [nov - 11 - 2025]

%% --- Parameters ---
C = 50e-3;      % Supercapacitor value [F]
Resr = 160e-3;  % Equivalent series resistance [Ω]
Vtarget = 1.8;  % Target voltage [V]

%% --- RF-DC converter data (P21XXCSR-EVB datasheet pg9 Vref = 1.2 v band 3) ---

% RF power reference points (dBm)
Prf_dBm1 = [-12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];

% Efficiency of RF-to-DC conversion
n1 = [0.005 0.08 0.22 0.33 0.44 0.50 0.54 0.56 0.56 0.54 0.525 0.51 0.54 0.58 0.605 0.62 0.635 0.63 0.62 0.60 0.575 0.55 0.52 0.48 0.45 0.42 0.38 0.35 0.325];

% RF power reference points (W)
Prf_W1 = 10.^((Prf_dBm1 - 30) / 10); 

% DC power reference points (W)
Pdc_W1 = Prf_W1.*n1;

% Plots: Typical performance graph - Powerharvester Efficiency vs. RF input power (dBm)
figure(1) 
plot(Prf_dBm1,n1,'*', 'MarkerSize', 8, 'LineWidth', 2) 
grid on
xlabel('P_R_F [dBm]'); ylabel('Efficiency \eta')
title('Typical Performance Graph - Powerharvester Efficiency vs. RFin [dBm]')
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 

figure(2)
semilogx(Prf_W1,n1,'*', 'MarkerSize', 8, 'LineWidth', 2) 
grid on
xlabel('P_R_F [W]'); ylabel('Efficiency \eta')
title('Typical Performance Graph - Powerharvester Efficiency vs. RFin [W]')
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 

%% --- Elena's model application ---
M1 = 0.012938483042989; % saturation of the output DC power 
a1 = 90;                % non-linear charging rate concerning the RF power
b1 = 1e-3;              % minimum voltage required for current to start flowing through diode

omega1 = 1/(1+exp(a1*b1));                   % constant that ensures a zero input/zero output response
phi1 = M1./(1+exp(-a1*(Prf_W1-b1)));         % logistic function related to the received RF power
Pdc_calc_W1 = (phi1-(M1*omega1))/(1-omega1); % harvested power by the rectenna

% Note: same equations of Elena's model, but different description
% omega = 1/(1+exp(a*b));
% phi = 1./(1+exp(-a*(Prf_W - b)));
% Pdc_calc_W  = M*(phi - omega)/(1-omega);

%% --- Results from measurements Vref = 1.2 V ---

% Measured time points
t1 = 0:10:430;  

% Measured voltage points
vc1 = [0.333 0.9 1.14 1.25 1.37 1.45 1.57 1.62 1.72 1.74 1.75 1.78 1.77 1.77 1.77 1.77 1.75 1.79 1.79 1.75 1.79...
    1.77 1.78 1.76 1.76 1.77 1.77 1.75 1.76 1.74 1.74 1.76 1.74 1.74 1.74 1.75 1.75 1.75 1.75 1.75 1.74 1.74...
    1.75 1.76];

% Adjustment of the Pdc value
Pdc_avg1 = 0.64e-3 % average DC power 

% Find Prf corresponding to Pdc_avg  
Prf_ref1 = interp1(Pdc_calc_W1, Prf_W1, Pdc_avg1, 'linear') 

% Time of the EH phase teh
teh1 = (C*Vtarget^2)/(2*Pdc_avg1)
teh1_v = 0:0.01:teh1; 

% Constant time taoc 
taoc1 = -teh1/log((Vtarget-0.99*Vtarget)/Vtarget)

% Variable equivalent resistance
R1_1 = taoc1/C - Resr

% Voltage curve of the EH phase
Vc1 = Vtarget*(1-exp(-teh1_v./taoc1));

% Plot: Elena's model vs Measurement points for Prf-to-Pdc conversion
figure(3)
semilogx(Prf_W1,Pdc_W1,'*','MarkerSize', 8, 'LineWidth', 2) 
grid on
hold on
semilogx(Prf_W1,Pdc_calc_W1,'-','LineWidth', 2) 
semilogx(Prf_ref1,Pdc_avg1,'*','MarkerSize', 8, 'LineWidth', 2) 
xlabel('P_R_F [W]'); ylabel('P_D_C [W]')
legend('measurements','simulation', 'Avg. P_D_C [W]')
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 
title('Nonlinear model and measurement points of RF-DC conversion')

% Plot: Comparison Measurements vs Model for Supercapacitor Charging Behavior at Vref = 1.2 V
figure(4)
plot(t1,vc1,'*','MarkerSize', 8, 'LineWidth', 2) 
grid on
hold on
plot(teh1_v,Vc1,'LineWidth', 2)
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 
xlim([0 450])  
ylim([0 1.9])
legend('measurements','simulation')
xlabel('Time [s]', 'FontSize', 14, 'FontWeight', 'bold') 
ylabel('Voltage [V]', 'FontSize', 14, 'FontWeight', 'bold') 
title('EH Phase - Comparison measurements and simulation - V_r_e_f = 1.2 [V]', 'FontSize', 14, 'FontWeight', 'bold') 

%% --- RF-DC converter data (P21XXCSR-EVB datasheet pg10 Vref = 0.9 v band 3) ---

% RF power reference points (dBm)
Prf_dBm2 = [-14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];

% Efficiency of RF-to-DC conversion
n2 = [0.015 0.08 0.20 0.32 0.39 0.455 0.485 0.515 0.516 0.50 0.48 0.475 0.50 0.53 0.55 0.57 0.575 0.574 0.57 0.55 0.535 0.51 0.48 0.45 0.42 0.39 0.36 0.34 0.30 0.275 0.25];

% RF power reference points (W)
Prf_W2 = 10.^((Prf_dBm2 - 30) / 10); 

% DC power reference points (W)
Pdc_W2 = Prf_W2.*n2;

% Plot: Typical performance graph - Powerharvester Efficiency vs. RF input power(dBm)
figure(5)
plot(Prf_dBm2,n2,'*','MarkerSize', 8, 'LineWidth', 2)
grid on
xlabel('P_R_F [dBm]'); ylabel('Efficiency \eta')
title('Typical Performance Graph - Powerharvester Efficiency vs. RFin [dBm]')
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 

figure(6)
semilogx(Prf_W2,n2,'*', 'MarkerSize', 8, 'LineWidth', 2)
grid on
xlabel('P_R_F [W]'); ylabel('Efficiency \eta')
title('Typical Performance Graph - Powerharvester Efficiency vs. RFin [W]')
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 

%% --- Elena's model application ---
M2 = 0.010;             % saturation of the output DC power 
a2 = 88;                % non-linear charging rate concerning the RF power
b2 = 1e-3;              % minimum voltage required for current to start flowing through diode

omega2 = 1/(1+exp(a2*b2));                   % constant that ensures a zero input/zero output response
phi2 = M2./(1+exp(-a2*(Prf_W2-b2)));         % logistic function related to the received RF power
Pdc_calc_W2 = (phi2-(M2*omega2))/(1-omega2); % harvested power by the rectenna

%% --- Results from measurements Vref = 0.9 V ---

% Measured time points
t2 = 0:10:430; 

% Measured voltage points
vc2 = [0.3 0.9 1.01 1.09 1.18 1.23 1.29 1.37 1.44 1.49 1.54 1.60 1.64 1.69 1.72 1.75 1.76 1.77 1.78 1.78 1.77...
     1.78 1.79 1.79 1.78 1.77 1.79 1.78 1.77 1.77 1.78 1.76 1.77 1.76 1.76 1.75 1.74 1.74 1.74 1.74 1.75 1.76...
     1.75 1.76];

% Adjustment of the Pdc value
Pdc_avg2 = 0.48e-3 % average DC power 

% Find Prf corresponding to Pdc_avg  
Prf_ref2 = interp1(Pdc_calc_W2, Prf_W2, Pdc_avg2, 'linear') 

% Time of the EH phase teh
teh2 = (C*Vtarget^2)/(2*Pdc_avg2)
teh2_v = 0:0.01:teh2; 

% Constant time taoc
taoc2 = -teh2/log((Vtarget-0.99*Vtarget)/Vtarget)

% Variable equivalent resistance
R1_2 = taoc2/C - Resr

% Voltage curve of the EH phase calculation
Vc2 = Vtarget*(1-exp(-teh2_v./taoc2));

% Plot: Elena's model vs Measurement points for Prf-to-Pdc conversion
figure(7)
semilogx(Prf_W2,Pdc_W2,'*','MarkerSize', 8, 'LineWidth', 2) 
grid on
hold on
semilogx(Prf_W2,Pdc_calc_W2,'-','LineWidth', 2) 
semilogx(Prf_ref2,Pdc_avg2,'*','MarkerSize', 8, 'LineWidth', 2) 
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 
xlabel('P_R_F [W]'); ylabel('P_D_C [W]')
legend('measurements','simulation', 'Avg. P_D_C [W]')
title('Nonlinear model and measurement points of RF-DC conversion')

% Plot: Comparison Measurements vs Model for Supercapacitor Charging Behavior at Vref = 0.9 V
figure(8)
plot(t2,vc2,'*','MarkerSize', 8, 'LineWidth', 2) 
grid on
hold on
plot(teh2_v,Vc2,'LineWidth', 2)
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 
xlim([0 450])  
ylim([0 1.9])
legend('measurements','simulation')
xlabel('Time [s]', 'FontSize', 14, 'FontWeight', 'bold') 
ylabel('Voltage [V]', 'FontSize', 14, 'FontWeight', 'bold') 
title('EH Phase - Comparison measurements and simulation - V_r_e_f = 0.9 [V]', 'FontSize', 14, 'FontWeight', 'bold') 

%% --- RF-DC converter data (P21XXCSR-EVB datasheet pg11 Vref = 0.7 v band 3) ---

% RF power reference points (dBm)
Prf_dBm3 = [-15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];

% Efficiency of RF-to-DC conversion
n3 = [0.025 0.11 0.21 0.32 0.38 0.44 0.46 0.475 0.47 0.465 0.45 0.445 0.47 0.50 0.52 0.53 0.535 0.53 0.525 0.515 0.49 0.47 0.445 0.42 0.385 0.36 0.33 0.315 0.28 0.26 0.23 0.22];

% RF power reference points (W)
Prf_W3 = 10.^((Prf_dBm3 - 30) / 10); 

% DC power reference points (W)
Pdc_W3 = Prf_W3.*n3;

% Plot: Typical performance graph - Powerharvester Efficiency vs. RF input power(dBm)
figure(9)
plot(Prf_dBm3,n3,'*', 'MarkerSize', 8, 'LineWidth', 2)
grid on
xlabel('P_R_F [dBm]'); ylabel('Efficiency \eta')
title('Typical Performance Graph - Powerharvester Efficiency vs. RFin [dBm]')
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 

figure(10)
semilogx(Prf_W3,n3,'*', 'MarkerSize', 8, 'LineWidth', 2)
grid on
xlabel('P_R_F [dBm]'); ylabel('Efficiency \eta')
title('Typical Performance Graph - Powerharvester Efficiency vs. RFin [dBm]')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

%% --- Elena's model application ---
M3 = 0.009;             % saturation of the output DC power 
a3 = 85;                % non-linear charging rate concerning the RF power
b3 = 1e-3;              % minimum voltage required for current to start flowing through diode

omega3 = 1/(1+exp(a3*b3));                   % constant that ensures a zero input/zero output response
phi3 = M3./(1+exp(-a3*(Prf_W3-b3)));         % logistic function related to the received RF power
Pdc_calc_W3 = (phi3-(M3*omega3))/(1-omega3); % harvested power by the rectenna

%% --- Results from measurements Vref = 0.7 V ---

% Measured time points
t3 = 0:10:430; 

% Measured voltage points
vc3 = [0.3 0.82 0.93 1.00 1.08 1.16 1.24 1.30 1.35 1.40 1.44 1.47 1.48 1.52 1.54 1.56 1.57 1.59 1.61 1.62 1.64...
     1.65 1.68 1.72 1.73 1.75 1.76 1.76 1.77 1.76 1.77 1.76 1.76 1.76 1.75 1.75 1.75 1.74 1.73 1.73 1.72 1.72...
     1.72 1.72];

% Adjustment of the Pdc value
Pdc_avg3 = 0.39e-3  % average DC power 

% Find Prf corresponding to Pdc_avg  
Prf_ref3 = interp1(Pdc_calc_W3, Prf_W3, Pdc_avg3, 'linear') 

% Time of the EH phase teh 
teh3 = (C*Vtarget^2)/(2*Pdc_avg3)
teh3_v = 0:0.01:teh3; % create vector

% Constant time taoc 
taoc3 = -teh3/log((Vtarget-0.99*Vtarget)/Vtarget)

% Variable equivalent resistance
R1_3 = taoc3/C - Resr

% Voltage curve of the EH phase 
Vc3 = Vtarget*(1-exp(-teh3_v./taoc3));

% Plot: Elena's model vs Measurement points for Prf-to-Pdc conversion
figure(11)
semilogx(Prf_W3,Pdc_W3,'*','MarkerSize', 8, 'LineWidth', 2) 
grid on
hold on
semilogx(Prf_W3,Pdc_calc_W3,'-','LineWidth', 2) 
semilogx(Prf_ref3,Pdc_avg3,'*','MarkerSize', 8, 'LineWidth', 2) 
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 
xlabel('P_R_F [W]'); ylabel('P_D_C [W]')
legend('measurements','simulation', 'Avg. P_D_C [W]')
title('Nonlinear model and measurement points of RF-DC conversion')

% Plot: Comparison Measurements vs Model for Supercapacitor Charging Behavior at Vref = 0.7 V
figure(12)
plot(t3,vc3,'*','MarkerSize', 8, 'LineWidth', 2) 
grid on
hold on
plot(teh3_v,Vc3,'LineWidth', 2)
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 
xlim([0 450])  
ylim([0 1.9])
legend('measurements','simulation')
xlabel('Time [s]', 'FontSize', 14, 'FontWeight', 'bold') 
ylabel('Voltage [V]', 'FontSize', 14, 'FontWeight', 'bold') 
title('EH Phase - Comparison measurements and simulation - V_r_e_f = 0.7 [V]', 'FontSize', 14, 'FontWeight', 'bold') 

%% --- Results from measurements of the TX phase for distance(transmitter -> RF-EH receiver) = 15 cm Vref(boost) = 1.2 V payload(packet length) = 30 B ptx (BLE TX load) = 5 dBm ---

% Measured time points
t_tx = 0:0.010:1.13; 

% Measured voltage points
v_tx = [1.1815 1.1648 1.1753 1.1648 1.1666 1.1753 1.1857 1.1857 1.1724 1.1531 1.1504 1.1515 1.1334 1.1334 1.1187 1.1300...
    1.1222 1.1191 1.1125 1.0958 1.1121 1.1125 1.1180 1.1153 1.1125 1.1021 1.1230 1.1070 1.1228 1.1155 1.1069 1.1125 1.0922...
    1.0956 1.0976 1.1003 1.0812 1.0916 1.0770 1.0646 1.0615 1.0707 1.0641 1.0812 1.0812 1.0769 1.0882 1.0901 1.0793 1.0858...
    1.0812 1.0712 1.0603 1.0713 1.0526 1.0603 1.0498 1.0499 1.0498 1.0496 1.0469 1.0177 1.0455 1.0394 1.0466 1.0498 1.0394...
    1.0403 1.0358 1.0373 1.0409 1.0474 1.0289 1.0262 1.0185 1.0200 1.0247 1.0220 1.0080 1.0147 0.9930 0.9976 0.9771 1.0036...    
    1.0025 0.9900 0.9976 1.0001 0.9923 1.0105 0.9877 1.0051 1.0136 1.0141 0.9777 0.9976 0.9871 0.9767 0.9775 0.9698 0.9767...
    0.9662 0.9558 0.9662 0.9662 0.9743 0.9677 0.9704 0.9696 0.9679 0.9748 0.9708 0.9735 0.9662];

% Input parameters for TX phase
Vc = v_tx(1);
Vin = 1.2;
Vout = 3.3;
n = 0.80;
deltaI = 0;
Ilos = 0;
taoc = 12.5635;
Pdc_ref = 1.4e-3;
tact = 111e-3;
tsleep = 100e-3;
Iact = 7.6e-3;
Isleep = 1.6e-3;
ttx = tact+tsleep;

% Calculate the input current Iin_act related to the current Iact
Iin_act = Iact*(Vout/Vin)*(1/n)+(Iact*deltaI)/2 + Ilos;

% Calculate the input current Iin_sleep related to the current Isleep
Iin_sleep = Isleep*(Vout/Vin)*(1/n)+(Isleep*deltaI)/2 + Ilos;

% Packet transmission  - Automatic Mode
% Definition of the time for each phase
tact_v = 0:0.001:tact;
tsleep_v = 0:0.001:tsleep;
 
% Initialize the vectors
Vd_v = [];
ttx_v = 0:0.001:0;  % Total TX time
 
% Number of packets to be sent
num_packets = 5;
 
% Initialize the initial voltage value
Vc_current = Vc;

for k = 1:num_packets
     % Active mode (tact)
     Vd_act = Vc_current - Iin_act*(Resr) - (Iin_act*(tact_v))/C + Vtarget*(1-exp(-tact_v/taoc));
     
     % Update the voltage after active mode
      Vc_current = Vd_act(end);
 
     % Sleep/Idle Mode (tsleep)
     Vd_sleep = Vc_current - Iin_sleep*(Resr) - (Iin_sleep*(tsleep_v))/C + Vtarget*(1-exp(-tsleep_v/taoc));
     
     % Update the voltage after sleep/idle mode
     Vc_current = Vd_sleep(end);
     
     % Store the values by correcting the overlap
     if k == 1
         Vd_v = [Vd_act Vd_sleep(2:end)];  % First packet, take everything
     else
         Vd_v = [Vd_v Vd_act(2:end) Vd_sleep(2:end)];  % Remove duplication
     end
     
     % Update the total TX time
     ttx_v = 0:0.001:(k*ttx);
end

% Plot: TX phase - Comparison measurements and simulation - Vref = 1.2 [V] 
figure(13)
plot(t_tx,v_tx,'*','MarkerSize', 4, 'LineWidth', 2) 
grid on
hold on
plot(ttx_v,Vd_v,'-','LineWidth', 2)
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 
legend('measurements','simulation')
xlabel('Time [s]', 'FontSize', 14, 'FontWeight', 'bold') 
ylabel('Voltage [V]', 'FontSize', 14, 'FontWeight', 'bold') 
title('TX Phase - Comparison measurements and simulation - V_r_e_f = 1.2 [V] and payload length = 30 B', 'FontSize', 14, 'FontWeight', 'bold') 

% New EH phase
teh_new = 0:0.001:300e-3;
Vc_new = (1.8-0.9593)*(1-exp(-teh_new/taoc))+0.9593;
teh_new_shift = 0.955:0.001:300e-3+0.955;

% Plot: TX phase - Comparison measurements and simulation - Vref = 1.2 [V] 
figure(14)
plot(t_tx,v_tx,'*','MarkerSize', 4, 'LineWidth', 2) 
grid on
hold on
plot([0 ttx_v(1:955)],[Vd_v(1) Vd_v(1:955)],'-','LineWidth', 2) 
plot(teh_new_shift,Vc_new,'--','LineWidth', 2)
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 
legend('measurements','simulation - TX phase','simulation - EH phase')
xlabel('Time [s]', 'FontSize', 14, 'FontWeight', 'bold') 
ylabel('Voltage [V]', 'FontSize', 14, 'FontWeight', 'bold') 
title('TX Phase - Comparison measurements and simulation - V_r_e_f = 1.2 [V] and payload length = 30 B', 'FontSize', 14, 'FontWeight', 'bold') 

teh_recharge = 3.856; % Measured average recharging time

%% --- Results from measurements of the TX phase for distance(transmitter -> RF-EH receiver) = 15 cm Vref(boost) = 1.2 V payload(packet length) = 64 B Ptx (BLE TX load) = 5 dBm ---

% Measured time points
t_tx = 0:0.005:1.015; 

% Measured voltage points
v_tx = [1.1786 1.1544 1.1514 1.1544 1.1471 1.1660 1.1431 1.1520 1.1343 1.1515 1.1303 1.1303 1.1413 1.1246 1.1291 1.1222...
    1.1303 1.1182 1.1182 1.1182 1.1062 1.1171 1.1022 1.1062 1.1062 1.1062 1.0953 1.1125 1.0941 1.1182 1.1090 1.1062 1.1062...
    1.1062 1.1182 1.1073 1.1182 1.1154 1.1062 1.1137 1.1176 1.1245 1.1182 1.1141 1.0941 1.1062 1.0910 1.0948 1.0941 1.0941...
    1.0862 1.0941 1.0917 1.0865 1.0835 1.0882 1.0821 1.0821 1.0794 1.0916 1.0788 1.0700 1.0821 1.0700 1.0700 1.0674 1.0363...
    1.0623 1.0579 1.0640 1.0589 1.0778 1.0606 1.0675 1.0700 1.0812 1.0700 1.0812 1.0778 1.0700 1.0700 1.0778 1.0691 1.0761...
    1.0700 1.0743 1.0821 1.0674 1.0579 1.0579 1.0520 1.0579 1.0579 1.0579 1.0579 1.0501 1.0459 1.0459 1.0467 1.0415 1.0360...
    1.0338 1.0338 1.0338 1.0400 1.0338 1.0218 1.0266 1.0218 1.0139 1.0218 1.0218 1.0218 1.0293 1.0338 1.0311 1.0259 1.0228...
    1.0155 1.0218 1.0338 1.0218 1.0218 1.0297 1.0327 1.0338 1.0218 1.0338 1.0241 1.0338 1.0177 1.0108 1.0160 1.0206 1.0097...
    0.9976 1.0189 1.0017 1.0097 0.9919 1.0097 1.0069 0.9901 1.0006 1.0057 0.9976 0.9856 1.0087 0.9930 0.9856 0.9856 0.9856...
    0.9735 0.9735 0.9735 0.9782 0.9834 0.9826 0.9817 0.9843 0.9735 0.9856 0.9856 0.9834 0.9946 0.9856 0.9748 0.9856 0.9852...
    0.9856 0.9756 0.9856 0.9856 0.9749 0.9735 0.9622 0.9735 0.9635 0.9646 0.9615 0.9721 0.9615 0.9732 0.9663 0.9635 0.9646...
    0.9810 0.9629 0.9735 0.9735 0.9758 0.9735 0.9735 0.9735 0.9735 0.9669 0.9735 0.9784 0.9836 0.9735 0.9735 0.9735 0.9789...
    0.9619];

% Input parameters for TX phase
Vc = v_tx(1);
Vin = 1.2;
Vout = 3.3;
n = 0.80;
deltaI = 0;
Ilos = 0;
taoc = 12.5635;
Pdc_ref = 1.4e-3;
tact = 118e-3;
tsleep = 100e-3;
Iact = 8e-3;
Isleep = 1.6e-3;
ttx = tact+tsleep;

% Calculate the input current Iin_act related to the current Iact
Iin_act = Iact*(Vout/Vin)*(1/n)+(Iact*deltaI)/2 + Ilos;

% Calculate the input current Iin_sleep related to the current Isleep
Iin_sleep = Isleep*(Vout/Vin)*(1/n)+(Isleep*deltaI)/2 + Ilos;

% Packet transmission  - Automatic Mode
% Definition of the time for each phase
tact_v = 0:0.001:tact;
tsleep_v = 0:0.001:tsleep;
 
% Initialize the vectors
Vd_v = [];
ttx_v = 0:0.001:0;  % Total TX time
 
% Number of packets to be sent
num_packets = 5;
 
% Initialize the initial voltage value
Vc_current = Vc;

for k = 1:num_packets
     % Active mode (tact)
     Vd_act = Vc_current - Iin_act*(Resr) - (Iin_act*(tact_v))/C + Vtarget*(1-exp(-tact_v/taoc));
     
     % Update the voltage after active mode
      Vc_current = Vd_act(end);
 
     % Sleep/Idle Mode (tsleep)
     Vd_sleep = Vc_current - Iin_sleep*(Resr) - (Iin_sleep*(tsleep_v))/C + Vtarget*(1-exp(-tsleep_v/taoc));
     
     % Update the voltage after sleep/idle mode
     Vc_current = Vd_sleep(end);
     
     % Store the values by correcting the overlap
     if k == 1
         Vd_v = [Vd_act Vd_sleep(2:end)];  % First packet, take everything
     else
         Vd_v = [Vd_v Vd_act(2:end) Vd_sleep(2:end)];  % Remove duplication
     end
     
     % Update the total TX time
     ttx_v = 0:0.001:(k*ttx);
end

% Plot: TX phase - Comparison measurements and simulation - Vref = 1.2 [V]
figure(15)
plot(t_tx,v_tx,'*','MarkerSize', 4, 'LineWidth', 2) 
grid on
hold on
plot(ttx_v,Vd_v,'-','LineWidth', 2)
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 
legend('measurements','simulation')
xlabel('Time [s]', 'FontSize', 14, 'FontWeight', 'bold') 
ylabel('Voltage [V]', 'FontSize', 14, 'FontWeight', 'bold') 
title('TX Phase - Comparison measurements and simulation - V_r_e_f = 1.2 [V] and payload length = 64 B', 'FontSize', 14, 'FontWeight', 'bold') 

% New EH phase
teh_new = 0:0.001:300e-3;
Vc_new = (1.8-0.9645)*(1-exp(-teh_new/taoc))+0.9645;
teh_new_shift = 0.895:0.001:300e-3+0.895;

% Plot: TX phase - Comparison measurements and simulation - Vref = 1.2 [V]
figure(16)
plot(t_tx,v_tx,'*','MarkerSize', 4, 'LineWidth', 2) 
grid on
hold on
plot([0 ttx_v(1:895)],[Vd_v(1) Vd_v(1:895)],'-','LineWidth', 2) 
plot(teh_new_shift,Vc_new,'--','LineWidth', 2)
set(gca, 'FontSize', 14, 'FontWeight', 'bold') 
legend('measurements','simulation - TX phase','simulation - EH phase')
xlabel('Time [s]', 'FontSize', 14, 'FontWeight', 'bold') 
ylabel('Voltage [V]', 'FontSize', 14, 'FontWeight', 'bold') 
title('TX Phase - Comparison measurements and simulation - V_r_e_f = 1.2 [V] and payload length = 64 B', 'FontSize', 14, 'FontWeight', 'bold') 

teh_recharge = 3.720; % Measured average recharging time

%% --- Results from measurements of the TX phase for distance(transmitter -> RF-EH receiver) = 15 cm Vref(boost) = 1.2 V payload(packet length) = 128 B ptx (BLE TX load) = 5 dBm ---

% Measured time points
t_tx = 0:0.010:0.880; 

% Measured voltage points
v_tx = [1.2182 1.2032 1.2076 1.2008 1.1868 1.1764 1.1790 1.1710 1.1630 1.1489 1.1588 1.1549...
    1.1391 1.1354 1.1381 1.1292 1.1303 1.1308 1.1354 1.1308 1.1308 1.1336 1.1388 1.1378 1.1378...
    1.1358 1.1360 1.1219 1.1283 1.1081 1.1201 1.1072 1.1052 1.1031 1.1067 1.0983 1.0986 1.0864...
    1.0686 1.0825 1.0805 1.0786 1.0844 1.0904 1.0826 1.0913 1.0842 1.0830 1.0910 1.0745 1.0665...
    1.0658 1.0556 1.0521 1.0574 1.0514 1.0485 1.0475 1.0432 1.0411 1.0263 1.0022 1.0033 1.0313...
    1.0263 1.0287 1.0263 1.0263 1.0221 1.0322 1.0343 1.0263 1.0101 1.0102 1.0098 1.0125 1.0065...
    0.9957 0.9865 0.9836 0.9907 1.0022 1.0015 0.9941 0.9958 0.9963 1.0022 1.0072 0.9970];


% Input parameters for TX phase
Vc = v_tx(1);
Vin = 1.2;
Vout = 3.3;
n = 0.80;
deltaI = 0;
Ilos = 0;
taoc = 12.5635;
Pdc_ref = 1.4e-3;
tact = 128e-3;
tsleep = 100e-3;
Iact = 8.9e-3;
Isleep = 1.6e-3;
ttx = tact+tsleep;

% Calculate the input current Iin_act related to the current Iact
Iin_act = Iact*(Vout/Vin)*(1/n)+(Iact*deltaI)/2 + Ilos;

% Calculate the input current Iin_sleep related to the current Isleep
Iin_sleep = Isleep*(Vout/Vin)*(1/n)+(Isleep*deltaI)/2 + Ilos;

% Packet transmission  - Automatic Mode
% Definition of the time for each phase
tact_v = 0:0.001:tact;
tsleep_v = 0:0.001:tsleep;
 
% Initialize the vectors
Vd_v = [];
ttx_v = 0:0.001:0;  % Total TX time
 
% Number of packets to be sent
num_packets = 4;
 
% Initialize the initial voltage value
Vc_current = Vc;

for k = 1:num_packets
     % Active mode (tact)
     Vd_act = Vc_current - Iin_act*(Resr) - (Iin_act*(tact_v))/C + Vtarget*(1-exp(-tact_v/taoc));
     
     % Update the voltage after active mode
      Vc_current = Vd_act(end);
 
     % Sleep/Idle Mode (tsleep)
     Vd_sleep = Vc_current - Iin_sleep*(Resr) - (Iin_sleep*(tsleep_v))/C + Vtarget*(1-exp(-tsleep_v/taoc));
     
     % Update the voltage after sleep/idle mode
     Vc_current = Vd_sleep(end);
     
     % Store the values by correcting the overlap
     if k == 1
         Vd_v = [Vd_act Vd_sleep(2:end)];  % First packet, take everything
     else
         Vd_v = [Vd_v Vd_act(2:end) Vd_sleep(2:end)];  % Remove duplication
     end
     
     % Update the total TX time
     ttx_v = 0:0.001:(k*ttx);
end

% Plot: TX phase - Comparison measurements and simulation - Vref = 1.2 [V] 
figure(17)
plot(t_tx,v_tx,'*','MarkerSize', 4, 'LineWidth', 2) 
grid on
hold on
plot(ttx_v,Vd_v,'-','LineWidth', 2)
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
legend('measurements','simulation')
xlabel('Time [s]', 'FontSize', 14, 'FontWeight', 'bold') 
ylabel('Voltage [V]', 'FontSize', 14, 'FontWeight', 'bold') 
title('TX Phase - Comparison measurements and simulation - V_r_e_f = 1.2 [V] and payload length = 128 B', 'FontSize', 14, 'FontWeight', 'bold') 

% New EH phase
teh_new = 0:0.001:300e-3;
Vc_new = (1.8-0.9866)*(1-exp(-teh_new/taoc))+0.9866;
teh_new_shift = 0.769:0.001:300e-3+0.769;

% Plot: TX phase - Comparison measurements and simulation - Vref = 1.2 [V]
figure(18)
plot(t_tx,v_tx,'*','MarkerSize', 4, 'LineWidth', 2) 
grid on
hold on
plot([0 ttx_v(1:769)],[Vd_v(1) Vd_v(1:769)],'-','LineWidth', 2) 
plot(teh_new_shift,Vc_new,'--','LineWidth', 2)
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
legend('measurements','simulation - TX phase','simulation - EH phase')
xlabel('Time [s]', 'FontSize', 14, 'FontWeight', 'bold') 
ylabel('Voltage [V]', 'FontSize', 14, 'FontWeight', 'bold') 
title('TX Phase - Comparison measurements and simulation - V_r_e_f = 1.2 [V] and payload length = 128 B', 'FontSize', 14, 'FontWeight', 'bold') 

teh_recharge = 4.211; % Measured average recharging time


