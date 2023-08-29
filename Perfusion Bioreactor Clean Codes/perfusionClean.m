clear all; close all; clc
perfusionBioreactorData = readtable("perfusionBioreactorData");
t_exp = perfusionBioreactorData.t; % XData
Xv_exp = perfusionBioreactorData.Xv;
Xd_exp = perfusionBioreactorData.Xd;
% Xt_exp = perfusionBioreactorData.Xt;
% Xl_exp = perfusionBioreactorData.Xl;
Glc_exp = perfusionBioreactorData.Glc;
Glu_exp = perfusionBioreactorData.Glu;
Gln_exp = perfusionBioreactorData.Gln;
% Ami_exp = perfusionBioreactorData.Ami;
Lac_exp = perfusionBioreactorData.Lac;
Amm_exp = perfusionBioreactorData.Amm;
mAb_exp = perfusionBioreactorData.mAb;
c_exp = [Xv_exp,Xd_exp,Glc_exp,Glu_exp,Gln_exp,Lac_exp,Amm_exp,mAb_exp]; % YData

T_exp = zeros(2*length(t_exp)-1, 1); % initialize new matrix with zeros
T_exp(1:2:end) = t_exp; % copy original matrix into odd-numbered elements of new matrix

% compute means and insert into even-numbered elements of new matrix
for i = 2:2:length(T_exp)
    T_exp(i) = (T_exp(i-1) + T_exp(i+1))/2;
end

C_exp = zeros(2*size(c_exp,1)-1, size(c_exp,2)); % initialize new matrix with zeros
C_exp(1:2:end,:) = c_exp; % copy original matrix into odd-numbered elements of new matrix

% compute means and insert into even-numbered elements of new matrix
for i = 2:2:size(C_exp,1)
    C_exp(i,:) = (C_exp(i-1,:) + C_exp(i+1,:))/2;
end

% Define simulation span
tinit = t_exp(1);
tfin = t_exp(end);

% noOfIntervals = 30;
% delta_t = (tfin-tinit)/noOfIntervals;
% % tVector = (tinit:tfin)';
% tVector = linspace(tinit,tfin,noOfIntervals+1)';

delta_t = 1;
noOfIntervals = (tfin-tinit)/delta_t;
tVector = linspace(tinit,tfin,noOfIntervals+1)';

% Define the initial conditions
% Viable cells, Dead cells, Lysed cells, Glucose, Glutamate, Glutamine, Other amino acids, Lactate, Ammonia, mAb product
cinit = [1.75,0.00,4.94,2.46,0.00,0.00,0.00,0.00]; % Initial concentrations

SimTime = tinit;%zeros(0,1);%[]
SimState = cinit;%zeros(0,3);%[]

% Define simulation options
optionsSim = []; %odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.1); % []

% LSQCURVEFIT implements two different algorithms: trust-region-reflective and levenberg-Marquardt.
% Choose one via the option Algorithm: for instance, to choose Levenberg-Marquardt, set
optionsOpt = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');%[];

noOfUnknownParameter = 15;

Theta = zeros(noOfIntervals*noOfUnknownParameter,1);%zeros(0,noOfIntervals);

theta0 = zeros(noOfIntervals*noOfUnknownParameter,1);

theta0(:,1) = 10; % mumax
theta0(:,2) = 10; % Ksglc
theta0(:,3) = 10; % Ksglu
theta0(:,4) = 10; % Kilac
theta0(:,5) = 10; % Kiamm
theta0(:,6) = 10; % qXd
theta0(:,7) = 10; % qglc
theta0(:,8) = 10; % qglu
theta0(:,9) = 10; % qgln1
theta0(:,10) = 10; % qgln2
theta0(:,11) = 10; % qgln3
theta0(:,12) = 10; % qlac1
theta0(:,13) = 10; % qlac2
theta0(:,14) = 10; % qamm
theta0(:,15) = 10; % qmAb
%theta0(:,16) = 10; % Ksgln

density = 1; % 1 g/mL = 1000 kg/m3
% Define the model parameters
Fbleed = [0.00;0.00;0.00;0.00;0.00;0.00;0.00;0.00;0.00;0.00;202.30;436.70;280.70;205.70;150.50;119.40;161.00;145.75;120.75;104.00;60.40;49.20;0.00;0.60;0.00;0.60;0.00;0.00;0.00;0.00;0.00;0.00;0.00;0.00]/density;
Xv_target = 100; % Viable cells target density, [E6 Cells/mL]
cfeed = [0.00,0.00,9.030312,2.5269,0.00,0.00,0.00,0.00]; % Feed concentrations
V = 2500; % [mL]
CSPR = 0.05e-6; % 0.05 nL/cell/day as given in exp data
PR = 2; % 2 vvd as given in exp data
%proportioninHarvest = [0,0,1,1,1,1,1,1,1,1]; % For lysed
proportioninHarvest = [0,0,1,1,1,1,1,1]; % For total and lysed

LB = [];%-1e-6*ones(1,size(theta0,2));%-1e-6*ones(noOfIntervals,size(theta0,2));
UB = [];%1e-6*ones(1,size(theta0,2));%1e6*ones(noOfIntervals,size(theta0,2));
for i = 1:noOfIntervals
    theta = lsqcurvefit(@(theta,t)perfusionBioreactorSimulator(theta,t,Fbleed(i),Xv_target,cfeed,V,CSPR,PR,proportioninHarvest,cinit,optionsSim),theta0(i,:),t_exp,c_exp,LB,UB,optionsOpt); %t_exp(i:i+1),c_exp(i:i+1,:), LB(i,:),UB(i,:)
    %     theta = [Theta,theta];
    %     Theta =
    % Use theta for simulation:
    theta0(i+1,:) = theta;
    tSpan = [tVector(i),tVector(i+1)];
    c_sim = perfusionBioreactorSimulator(theta,tSpan,Fbleed(i),Xv_target,cfeed,V,CSPR,PR,proportioninHarvest,cinit,optionsSim);
    SimTime = [SimTime;tSpan(2)];
    SimState = [SimState;c_sim(end,:)];
    cinit = c_sim(end,:);
    Theta(i,:) = theta;%[Theta;theta];
    %loopCounter = loopCounter+1;
    %     LB = -theta;
    %     UB = theta;
end

time = SimTime;
tTheta = time(1:end-1);
Xv = SimState(:,1);
Xd = SimState(:,2);
Glc = SimState(:,3);
Glu = SimState(:,4);
Gln = SimState(:,5);
Lac = SimState(:,6);
Amm = SimState(:,7);
mAb = SimState(:,8);

figure(1)
subplot(2,4,1)
plot(t_exp,Xv_exp,'p',time,Xv); %xlim([tinit tfin]) %axis([0 t(end) 0 inf])
xlabel('Time (d)');
ylabel('Viable Cell Density (E6 Cells/mL)');
legend('Exp','Sim','Location','best')
subplot(2,4,2)
plot(t_exp,Xd_exp,'p',time,Xd); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Dead Cell Density (E6 Cells/mL)');
legend('Exp','Sim','Location','best')
subplot(2,4,3)
plot(t_exp,Glc_exp,'p',time,Glc); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Glucose Concentration (g/L)');
legend('Exp','Sim','Location','best')
subplot(2,4,4)
plot(t_exp,Glu_exp,'p',time,Glu); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Glutamate Concentration (g/L)');
legend('Exp','Sim','Location','best')
subplot(2,4,5)
plot(t_exp,Gln_exp,'p',time,Gln); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Glutamine Concentration (g/L)');
legend('Exp','Sim','Location','best')
subplot(2,4,6)
plot(t_exp,Lac_exp,'p',time,Lac); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Lactate Concentration (g/L)');
legend('Exp','Sim','Location','best')
subplot(2,4,7)
plot(t_exp,Amm_exp,'p',time,Amm); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Ammonia Concentration (g/L)');
legend('Exp','Sim','Location','best')
subplot(2,4,8) 
plot(t_exp,mAb_exp,'p',time,mAb); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('mAb Concentration (g/L)');
legend('Exp','Sim','Location','best')

figure(2)
subplot(4,3,1)
stairs(tTheta,Theta(:,1)); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('mu (1/d)');
grid on
subplot(4,3,2)
stairs(tTheta,Theta(:,6)); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('qXd (1/d)');
grid on
subplot(4,3,3)
stairs(tTheta,Theta(:,7)); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('qglc ((g/L)_{glc}/(E6 Cells/mL))');
grid on
subplot(4,3,4)
stairs(tTheta,Theta(:,8)); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('qglu ((g/L)_{glu}/(E6 Cells/mL))');
grid on
subplot(4,3,5)
stairs(tTheta,Theta(:,9)); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('qgln1 ((g/L)_{gln}/(g/L)_{glu})');
grid on
subplot(4,3,6)
stairs(tTheta,Theta(:,10)); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('qgln2 ((g/L)_{gln}/(E6 Cells/mL))');
grid on
subplot(4,3,7)
stairs(tTheta,Theta(:,11)); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('qgln3 (1/d)');
grid on
subplot(4,3,8)
stairs(tTheta,Theta(:,12)); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('qlac1 ((g/L)_{lac}/(g/L)_{glc})');
grid on
subplot(4,3,9)
stairs(tTheta,Theta(:,13)); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('qlac2 ((g/L)_{lac}/(E6 Cells/mL))');
grid on
subplot(4,3,10)
stairs(tTheta,Theta(:,14)); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('qamm ((g/L)_{amm}/(g/L)_{gln})');
grid on
subplot(4,3,11)
stairs(tTheta,Theta(:,15)); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('qmAb ((g/L)_{mAb}/(E6 Cells/mL))');
grid on
subplot(4,3,12)
stairs(tTheta,Fbleed); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Fbleed (g/d)');
grid on

% Back propagation to investigate cell metabolic kinetic rate equation:
% mu is YDATA and c(3), c(4), c(6) and c(7) are XDATA.

thetamu0 = [0.5;1.0;2.5;10.5;5.5];
optionsOptmu = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');%[];
LBmu = [];
UBmu = [];

cmu = [SimState(1:end-1,3),SimState(1:end-1,4),SimState(1:end-1,6),SimState(1:end-1,7)];%,SimState(1:end-1,5)]; 

%muHandle = @(thetamu)mumaxmu*(c(3)/(Ksglcmu+c(3)))*(c(4)/(Ksglumu+c(4)))*(Kilacmu/(Kilacmu+c(6)))*(Kiammmu/(Kiammmu+c(7)));
%muHandle = @(thetamu,cmu)thetamu(1)*(cmu(1)/(thetamu(2)+cmu(1)))*(cmu(2)/(thetamu(3)+cmu(2)))*(thetamu(4)/(thetamu(4)+cmu(3)))*(thetamu(5)/(thetamu(5)+cmu(4)));

%thetamu = lsqcurvefit(muHandle,thetamu0,cmu,Theta(:,1),LBmu,UBmu,optionsOptmu); %t_exp(i:i+1),c_exp(i:i+1,:), LB(i,:),UB(i,:)
thetamu = lsqcurvefit(@(thetamu,cmu)muFunction(thetamu,cmu),thetamu0,cmu,Theta(:,1),LBmu,UBmu,optionsOptmu); %t_exp(i:i+1),c_exp(i:i+1,:), LB(i,:),UB(i,:)

function c_sim = perfusionBioreactorSimulator(theta,t,Fbleed,Xv_target,cfeed,V,CSPR,PR,proportioninHarvest,cinit,optionsSim)

% Using the two lines below instead threw an error: %Error using lsqcurvefit %Function value and YDATA sizes are not equal. %Error in perfusionBioreactorModelSimulator (line 328) % theta = lsqcurvefit(@(theta,t)perfusionBioreactorSimulator(theta,t,Fbleed(i),Xv_target,cfeed,V,CSPR,PR,proportioninHarvest,cinit,optionsSim),theta0(i,:),t_exp,c_exp,LB,UB,optionsOpt);
% function [t_sim,c_sim] = perfusionBioreactorSimulator(theta,t,Fbleed,Xv_target,cfeed,V,CSPR,PR,proportioninHarvest,cinit,optionsSim)
% [t_sim,c_sim] = ode45(@(t,c,theta)perfusionBioreactor(t,c,theta,Fbleed,Xv_target,cfeed,V,CSPR,PR,proportioninHarvest),t,cinit,optionsSim,theta);

[~,c_sim] = ode45(@(t,c,theta)perfusionBioreactor(t,c,theta,Fbleed,Xv_target,cfeed,V,CSPR,PR,proportioninHarvest),t,cinit,optionsSim,theta);
end

%function dcdt = perfusionBioreactor(~,c,Fbleed,mumax,Ksglc,Ksglu,Ksgln,Ksami,Kmlac,Kmamm,qXd,qXl,qglc,qglu,qgln1,qgln2,qgln3,qami,qlac1,qlac2,qamm,qmAb,Xv_target,cfeed,V,CSPR,PR,proportioninHarvest)
function dcdt = perfusionBioreactor(~,c,theta,Fbleed,Xv_target,cfeed,V,CSPR,PR,proportioninHarvest)
%theta = zeros(noOfIntervals,1);
mumax  = theta(1);
Ksglc = theta(2);
Ksglu = theta(3);
Kilac  = theta(4);
Kiamm = theta(5);
qXd  = theta(6);
qglc  = theta(7);
qglu = theta(8);
qgln1 = theta(9);
qgln2 = theta(10);
qgln3 = theta(11);
qlac1 = theta(12);
qlac2  = theta(13);
qamm = theta(14);
qmAb  = theta(15);


dcdt = zeros(8,1);

mu = mumax;%mumax*(c(3)/(Ksglc+c(3)))*(c(4)/(Ksglu+c(4)))*(Kilac/(Kilac+c(6)))*(Kiamm/(Kiamm+c(7)));

% Control option A:
% if c(1) <= Xv_target
%     Fmedia = CSPR*c(1)*V; % [mL/d]
% else
%     Fmedia = PR*V;
% end
% Fharvest = Fmedia-Fbleed; % This holds because of the assumptionof constsnt bioreactor volume!

% Control option B:
if c(1) <= Xv_target
    Fharvest = CSPR*c(1)*V; % [mL/d], CSPR nL/cell/day
else
    Fharvest = PR*V; % vvd, i.e., vessel volume per day
end
Fmedia = Fharvest+Fbleed; % This holds because of the assumptionof constsnt bioreactor volume!

% mumax = Viable cells specific growth rate, 1/time
% Ksglc = ;
% Ksglu = ;
% Kmlac  = ;
% Kmamm = ;
% qXd = Viable cells specific death rate, 1/time
% qglc = Glucose consumption per viable cell, glc/Xv
% qglu = Glutamate consumption per viable cell, glu/Xv
% qgln1 = Glutamine production per glutamate, gln/glu
% qgln2 = Glutamine consumption per viable cell, gln/Xv
% qgln3 = Specific glutamine dissociation rate, 1/time
% qglc1 = Lactate production per glucose, lac/glc
% qglc2 = Lactate consumption per viable cell, lac/Xv
% qamm = Ammonia production per glutamine, amm/gln
% qmAb = mAb production per viable cell, mAb/Xv

r = [(mu-qXd)*c(1);qXd*c(1);-qglc*c(1);-qglu*c(1);qgln1*c(4)-qgln2*c(1)-qgln3*c(5);qlac1*c(3)-qlac2*c(1);qamm*c(5);qmAb*c(1)];

dcdt(1) = (Fmedia/V)*cfeed(1)-(Fbleed/V)*c(1)-(Fharvest/V)*c(1)*proportioninHarvest(1)+r(1); % Viable cells
dcdt(2) = (Fmedia/V)*cfeed(2)-(Fbleed/V)*c(2)-(Fharvest/V)*c(2)*proportioninHarvest(2)+r(2); % Dead cells
dcdt(3) = (Fmedia/V)*cfeed(3)-(Fbleed/V)*c(3)-(Fharvest/V)*c(3)*proportioninHarvest(3)+r(3); % Glucose
dcdt(4) = (Fmedia/V)*cfeed(4)-(Fbleed/V)*c(4)-(Fharvest/V)*c(4)*proportioninHarvest(4)+r(4); % Glutamate
dcdt(5) = (Fmedia/V)*cfeed(5)-(Fbleed/V)*c(5)-(Fharvest/V)*c(5)*proportioninHarvest(5)+r(5); % Glutamine
dcdt(6) = (Fmedia/V)*cfeed(6)-(Fbleed/V)*c(6)-(Fharvest/V)*c(6)*proportioninHarvest(6)+r(6); % Lactate
dcdt(7) = (Fmedia/V)*cfeed(7)-(Fbleed/V)*c(7)-(Fharvest/V)*c(7)*proportioninHarvest(7)+r(7); % Ammonia
dcdt(8) = (Fmedia/V)*cfeed(8)-(Fbleed/V)*c(8)-(Fharvest/V)*c(8)*proportioninHarvest(8)+r(8); % mAb product
end

% Back propagation function here:
function muFun = muFunction(thetamu,cmu)
mumaxmu = thetamu(1);
Ksglcmu = thetamu(2);
Ksglumu = thetamu(3);
Kilacmu = thetamu(4);
Kiammmu = thetamu(5);
%muFun = mumaxmu*(cmu(1)/(Ksglcmu+cmu(1)))*(cmu(2)/(Ksglumu+cmu(2)))*(Kilacmu/(Kilacmu+cmu(3)))*(Kiammmu/(Kiammmu+cmu(4)));
%muFun = mumaxmu*(cmu(1)/(Ksglcmu+cmu(1))).*(cmu(2)/(Ksglumu+cmu(2))).*(Kilacmu/(Kilacmu+cmu(3))).*(Kiammmu/(Kiammmu+cmu(4)));
%muFun = mumaxmu*(cmu(:,1)./(Ksglcmu+cmu(:,1))).*(cmu(:,2)./(Ksglumu+cmu(:,2)));%.*(Kilacmu./(Kilacmu+cmu(:,3))).*(Kiammmu./(Kiammmu+cmu(:,4)));
%muFun = mumaxmu*(cmu(:,1)./(Ksglcmu+cmu(:,1))).*(cmu(:,2)./(Ksglumu+cmu(:,2))).*(Kilacmu./(Kilacmu+cmu(:,3))).*(Kiammmu./(Kiammmu+cmu(:,4)));
muFun = mumaxmu*(cmu(:,1)./(Ksglcmu+cmu(:,1))).*(cmu(:,2)./(Ksglumu+cmu(:,2))).*(Kilacmu./(Kilacmu+cmu(:,3))).*(Kiammmu./(Kiammmu+cmu(:,4)));
%muFun = mumaxmu*(cmu(:,1)./(Ksglcmu+cmu(:,1))).*(cmu(:,2)./(Ksglumu+cmu(:,2))).*(Kilacmu./(Kilacmu+cmu(:,3))).*(Kiammmu./(Kiammmu+cmu(:,4))).*(cmu(:,5)./(Ksglnmu+cmu(:,5)));
end