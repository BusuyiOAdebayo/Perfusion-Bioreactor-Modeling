clear all; close all; clc
% Define simulation span
tinit = 0;
tfin = 30;

delta_t = 1;
noOfIntervals = (tfin-tinit)/delta_t;
tVector = linspace(tinit,tfin,noOfIntervals+1)';

% Define the initial conditions
% Viable cells, Dead cells, Lysed cells, Glucose, Glutamate, Glutamine, Other amino acids, Lactate, Ammonia, mAb product
%cinit = [1.75,0.00,0.00,4.94,2.46,0.00,5.00,0.00,0.00,0.00]; % Initial concentrations, for lysed
cinit = [1.75,0.00,4.94,2.46,0.00,0.00,0.00,0.00]; % Initial concentrations

% Define simulation options
optionsSim = []; %odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.1); % []

% theta(1) = 10; % mumax
% theta(2) = 10; % Ksglc
% theta(3) = 10; % Ksglu
% theta(4) = 10; % Kilac
% theta(5) = 10; % Kiamm
% theta(6) = 10; % qXd
% theta(7) = 10; % qglc
% theta(8) = 10; % qglu
% theta(9) = 10; % qgln1
% theta(10) = 10; % qgln2
% theta(11) = 10; % qgln3
% theta(12) = 10; % qlac1
% theta(13) = 10; % qlac2
% theta(14) = 10; % qamm
% theta(15) = 10; % qmAb
%theta0(:,16) = 10; % Ksgln

SimTime = tinit;%zeros(0,1);%[]
SimState = cinit;%zeros(0,3);%[]
%22.4955708947415,1.63149375575746,41.1367681530088,59.0149418412106,3.75430879675822,2.50000000000000
theta = [22.4955708947415,1.63149375575746,41.1367681530088,59.0149418412106,3.75430879675822,0.005526637,0.123553246,0.020039192,83.72545907,0.55704426,87.65361892,-1.251269118,-0.007744341,136.8218673,0.021678423];

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

loopCounter = 0;
for i = 1:noOfIntervals
    tSpan = [tVector(i),tVector(i+1)];
    [t,c] = ode45(@(t,c)perfusionBioreactor(t,c,theta,Fbleed(i),Xv_target,cfeed,V,CSPR,PR,proportioninHarvest),tSpan,cinit,optionsSim);%,theta);
    SimTime = [SimTime;t(end)];
    SimState = [SimState;c(end,:)];
    cinit = c(end,:);
    loopCounter = loopCounter+1;
end

t = SimTime;
Xv = SimState(:,1);
Xd = SimState(:,2);
Glc = SimState(:,3);
Glu = SimState(:,4);
Gln = SimState(:,5);
Lac = SimState(:,6);
Amm = SimState(:,7);
mAb = SimState(:,8);

% Plot the results
figure(1)
subplot(2,4,1)
plot(t,Xv); %xlim([tinit tfin]) %axis([0 t(end) 0 inf])
xlabel('Time (d)');
ylabel('Viable Cell Density (E6 Cells/mL)');
subplot(2,4,2)
plot(t,Xd); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Dead Cell Density (E6 Cells/mL)');
subplot(2,4,3)
plot(t,Glc); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Glucose Concentration (g/L)');
subplot(2,4,4)
plot(t,Glu); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Glutamate Concentration (g/L)');
subplot(2,4,5)
plot(t,Gln); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Glutamine Concentration (g/L)');
subplot(2,4,6)
plot(t,Lac); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Lactate Concentration (g/L)');
subplot(2,4,7)
plot(t,Amm); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('Ammonia Concentration (g/L)');
subplot(2,4,8) 
plot(t,mAb); %xlim([tinit tfin])
xlabel('Time (d)');
ylabel('mAb Concentration (g/L)');

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

mu = mumax*(c(3)/(Ksglc+c(3)))*(c(4)/(Ksglu+c(4)))*(Kilac/(Kilac+c(6)))*(Kiamm/(Kiamm+c(7)));

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