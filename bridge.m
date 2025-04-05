clc; clear; close all;
L = 1200; % bridge length
n = 1200; % subsections of L
P = 400; % total weight of train
x = linspace(0, L, n+1); % x-axis for SFD & BMD

%%
x_train = [52 228 392 568 732 908]; % train load locations (wheels)
P_train = [1 1 1 1 1 1] * -P/6; % load for each point

n_train = 1200; % num of train locations

Wi = zeros(n_train, n+1); %applied load
SFDi = zeros(n_train, n+1); % a SFD for each train location
BMDi = zeros(n_train, n+1); % same as above but BMD

for i = 1:n_train
    traincmloc = 1200/(n_train+1) * i; % place center of train at even spots along the bridge 
    flocs = x_train - (480-traincmloc); % make array of force locations
    temp_P_train = P_train; % copy force loads

    for j = 1:6
        if flocs(j) < 0 || flocs(j) > 1200 % check if any force locations are off bridge
            temp_P_train(j) = 0; % if so make force = 0
            
        end
    end

    By = sum(temp_P_train .* flocs) / -L; % solve for By with sum of moments from left
    Ay = -sum(temp_P_train) - By; % solve for Ay by Fy = 0

    Wi(i, 1) = Ay; % add support
    Wi(i, n+1) = By; % add support

    for j = 1:length(flocs) 
        if flocs(j) <= L && flocs(j) >= 0 % make sure is on bridge
            [~, idx] = min(abs(x - flocs(j)));
            Wi(i, idx) = temp_P_train(j);
        end
    end
    
    SFDi(i, :) = cumsum(Wi(i, :)); % calculate SFD
    BMDi(i, :) = cumsum(SFDi(i, :)); % calculate BMD



   %figure; plot(Wi(i, :)) % plot FBD
   %figure; plot(SFDi(i, :)) % plot SFD
   %figure; plot(BMDi(i, :)) % plot BMD

end

SFD = max(abs(SFDi)); % SFD envelope
SFD_for_graphing = max(SFDi)
BMD = max(BMDi); % BMD envelope

%figure; plot(SFD(1, :))
%figure; plot(BMD(1, :))
%%
% Bridge Parameters

% xc = location of cross-section change
% Btf = top flange width
% Ttf = top flange thickness (height)
% Bgt = glue tab width
% Tgt = glue tab thickness
% Bst = sides width 
% Tst = sides thickness (in original this is 75 - 1.27)
% Bb = bottom width
% Tb = bottom thickness



% = xc, Btf, Ttf, Bgt, Tgt, Bst, Tst, Bb, Tb 

% default: param = [0, 100, 1.27, 5, 1.27, 1.27, 73.73, 80, 1.27]
param = [0, 100, 1.27, 5, 1.27, 1.27, 73.73, 80, 1.27]

number_of_webs = 4
% assuming webs are placed evenly

areas_of_param = [param(2)*param(3), param(4)*param(5), param(4)*param(5), param(6)*param(7), param(6)*param(7), param(8)*param(9)]
y_values_from_bot = [param(7)+param(9)+param(3)/2, param(7)+param(9)-param(5)/2, param(7)+param(9)-param(5)/2, param(9) + param(7)/2, param(9) + param(7)/2, param(9)/2]


% bft = interp1(param(:,1), param(:,2), x)


%%
% y value of centriodal axis
ybar = sum(areas_of_param .* y_values_from_bot) / sum(areas_of_param)

% calculate I with parallel axis thereom (note the + Ad^2 functions could
% be wrong but I'm returning the correct value for design 0 so

I = (param(2)*param(3)^3)/12 + 2 * (param(4)*param(5)^3)/12 + 2*(param(6)*param(7)^3)/12 + (param(8)*param(9)^3)/12 + sum(areas_of_param .* (ybar - y_values_from_bot).^2)

%assuming that the centroid is never below the top of the bottom bar (which
%likely should be the case)

%calculate Q's
Qcent = areas_of_param(6) * (ybar - y_values_from_bot(6)) + 2 * param(6)*(ybar - param(9))*(ybar - param(9))/2

Qglue = areas_of_param(1) * (y_values_from_bot(1) - ybar)

M = max(BMD)

V = max(SFD)

%calculate flexural stresses and shear stresses

Flexural_Stress_Top = M*(param(9)+param(7)+param(3)-ybar)/I

Flexural_Stress_Bot = M*(ybar)/I

Shear_Stress_Max = V*Qcent/(I*2*param(6))

Shear_Stress_Glue = V*Qglue/(I*2*(param(4)+param(6)))

% initialize material variables
E = 4000;
mu = 0.2;
S_tens = 30;
S_comp = 6;
T_max = 4;
T_glue = 2;

% calculate maximum buckling stresses based on geometry

S_buck1 = 4*(pi^2)*E/(12*(1-mu^2))*(param(3)/param(8))^2
S_buck2 = 0.4254*(pi^2)*E/(12*(1-mu^2))*(param(3)/((param(2)-param(8))/2))^2
S_buck3 = 6*(pi^2)*E/(12*(1-mu^2))*(param(6)/(param(7)+param(9)-ybar))^2
T_buck = 5*(pi^2)*E/(12*(1-mu^2))*((param(6)/param(7))^2+(param(6)/(L/(number_of_webs-1)))^2)

% report out FOS

FOS_tension = S_tens / Flexural_Stress_Bot
FOS_compression = S_comp / Flexural_Stress_Top
FOS_buck1 = S_buck1 / Flexural_Stress_Top % top flange mid buckling
FOS_buck2 = S_buck2 / Flexural_Stress_Top % top flange sides buckling
FOS_buck3 = S_buck3 / Flexural_Stress_Top
FOS_shear = T_max / Shear_Stress_Max
FOS_glue = T_glue / Shear_Stress_Glue
FOS_shear_buck = T_buck / Shear_Stress_Max

FOS = [FOS_tension FOS_compression FOS_shear FOS_glue FOS_buck1 FOS_buck2 FOS_buck3 FOS_shear_buck];
FOSindex = ["FOS_tension" "FOS_compression" "FOS_shear" "FOS_glue" "FOS_buck1" "FOS_buck2" "FOS_buck3" "FOS_shear_buck"];

% find lowest FOS and failure and print it out

[minFOS, index] = min(FOS);
minFOS
Pf = minFOS*max(SFD)
reason = FOSindex(index)

% find failure loads for all types for graphing

Mf_tens = FOS_tension * M
Mf_comp = FOS_compression * M
Mf_buck1 = FOS_buck1 * M
Mf_buck2 = FOS_buck2 * M
Mf_buck3 = FOS_buck3 * M

Vf_shear = FOS_shear * V
Vf_glue = FOS_glue * V
Vf_buck = FOS_shear_buck * V

%graph plots

subplot(2,3,1)
hold on; grid on; grid minor;
yline(Vf_shear, "r")
plot(x, SFD, 'k');
plot([0, 1200], [0,0],'k','LineWidth',2)
legend("Matboard Shear Failure", 'Location', 'northoutside', 'Orientation', 'vertical')
xlabel("Distance along bridge (mm)")
ylabel("Shear Force (N)")

subplot(2,3,2)
hold on; grid on;, grid minor;
yline(Vf_glue, "r")
plot(x, SFD, 'k');
plot([0,1200],[0,0],'k','LineWidth',2)
legend("Glue Shear Failure", 'Location', 'northoutside', 'Orientation', 'vertical')
xlabel("Distance along bridge (mm)")
ylabel("Shear Force (N)")

subplot(2,3,3)
hold on; grid on; grid minor;
yline(Vf_buck, "r")
plot(x, SFD, 'k');
plot([0, 1200], [0,0],'k','LineWidth',2)
legend("Matboard Shear Buckling Failure", 'Location', 'northoutside', 'Orientation', 'vertical')
xlabel("Distance along bridge (mm)")
ylabel("Shear Force (N)")

subplot(2,3,4)
hold on; grid on; grid minor;
set(gca, 'Ydir','reverse')
axis([0 1200 (-0.5 * 10^5) (3.5*10^5)])
yline(Mf_tens, "r")
yline(Mf_comp, "b")
plot(x, BMD, 'k');
plot([0, 1200], [0,0],'k','LineWidth',2)
legend("Matboard Tension Failure", "Matboard Compression Failure", 'Location', 'northoutside', 'Orientation', 'vertical')
xlabel("Distance along bridge (mm)")
ylabel("Bending Moment (Nmm)")

subplot(2,3,5)
hold on; grid on; grid minor;
set(gca, 'Ydir','reverse')
axis([0 1200 (-0.5 * 10^5) (3.5*10^5)])
yline(Mf_buck1, "r")
yline(Mf_buck2, "b")
plot(x, BMD, 'k');
plot([0, 1200], [0,0],'k','LineWidth',2)
legend("Matboard Buckling Failure, Top Flange - Mid", "Matboard Buckling Failure, Top Flange - Sides", 'Location', 'northoutside', 'Orientation', 'vertical')
xlabel("Distance along bridge (mm)")
ylabel("Bending Moment (Nmm)")

subplot(2,3,6)
hold on; grid on; grid minor;
set(gca, 'Ydir','reverse')
axis([0 1200 (-0.5 * 10^5) (4*10^5)])
yline(Mf_buck3, "r")
plot(x, BMD, 'k');
plot([0, 1200], [0,0],'k','LineWidth',2)
legend("Matboard Buckling Failure, Webs", 'Location', 'northoutside', 'Orientation', 'vertical')
xlabel("Distance along bridge (mm)")
ylabel("Bending Moment (Nmm)")