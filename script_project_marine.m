NETaddAssembly('DHI.Generic.MikeZero.DFS.dll');
import DHI.Generic.MikeZero.DFS.*;


dfsu3 = DfsFileFactory.DfsuFileOpen('project.m3fm - Result files/Cadmium_conc0205_0_1.dfsu');
%change according to concentration you want, then run this + calculate
%average
data=ReadItemTimeStep(dfsu3,2,14);

% Node coordinates
xn = double(dfsu3.X);
yn = double(dfsu3.Y);
zn = double(dfsu3.Z);

% Create element table in Matlab format

tn3D = mzNetFromElmtArray(dfsu3.ElementTable);
% also calculate element center coordinates
[xe,ye,ze] = mzCalcElmtCenterCoords(tn3D,xn,yn,zn);

% Read some item information
items = {};
for i = 0:dfsu3.ItemInfo.Count-1
   item = dfsu3.ItemInfo.Item(i);
   items{i+1,1} = char(item.Name);
   items{i+1,2} = char(item.Quantity.Unit);
   items{i+1,3} = char(item.Quantity.UnitAbbreviation); 
end

nsteps = dfsu3.NumberOfTimeSteps;

nlayers = dfsu3.NumberOfLayers;

 % Load 3D data

SHM=[];
XHM=[];
XSS=[];
SED=[];
for i=1:nsteps-1
  shm   = double(dfsu3.ReadItemTimeStep(2,i).Data);
  xhm   = double(dfsu3.ReadItemTimeStep(3,i).Data);
  xss   = double(dfsu3.ReadItemTimeStep(4,i).Data);
  sed   = double(dfsu3.ReadItemTimeStep(5,i).Data);
  SHM=[SHM;shm];
  XHM=[XHM; xhm];
  XSS=[XSS;xss];
  SED=[SED;sed];
end
%% 0.001 mg/L
%take average over the 74 locations and 10 depths
average_XHM=mean(XHM,2); %adsorbed HM (mg/L)
average_SHM=mean(SHM,2); %dissolved HM (mg/L) %compare with PNEC (0.19 µg/L)
average_XSS=mean(XSS,2); %suspended solids (mg/L)
average_SED=mean(SED,2);
%% 0.002 mg/L
%take average over the 74 locations and 10 depths
average_XHM_0_002=mean(XHM,2); %adsorbed HM (mg/L)
average_SHM_0_002=mean(SHM,2); %dissolved HM (mg/L) %compare with PNEC (0.19 µg/L)
average_XSS_0_002=mean(XSS,2); %suspended solids (mg/L)
average_SED_0_002=mean(SED,2);
%% 0.005 mg/L
%take average over the 74 locations and 10 depths
average_XHM_0_005=mean(XHM,2); %adsorbed HM (mg/L)
average_SHM_0_005=mean(SHM,2); %dissolved HM (mg/L) %compare with PNEC (0.19 µg/L)
average_XSS_0_005=mean(XSS,2); %suspended solids (mg/L)
average_SED_0_005=mean(SED,2);
%% 0.01 mg/L
average_XHM_0_01=mean(XHM,2); %adsorbed HM (mg/L)
average_SHM_0_01=mean(SHM,2); %dissolved HM (mg/L) %compare with PNEC (0.19 µg/L)
average_XSS_0_01=mean(XSS,2); %suspended solids (mg/L)
average_SED_0_01=mean(SED,2);

%% 0.03 mg/L
average_XHM_0_03=mean(XHM,2); %adsorbed HM (mg/L)
average_SHM_0_03=mean(SHM,2); %dissolved HM (mg/L) %compare with PNEC (0.19 µg/L)
average_XSS_0_03=mean(XSS,2); %suspended solids (mg/L)
average_SED_0_03=mean(SED,2);
%% 0.05 mg/L
average_XHM_0_05=mean(XHM,2); %adsorbed HM (mg/L)
average_SHM_0_05=mean(SHM,2); %dissolved HM (mg/L) %compare with PNEC (0.19 µg/L)
average_XSS_0_05=mean(XSS,2); %suspended solids (mg/L)
average_SED_0_05=mean(SED,2);

%% 0.1 mg/L
average_XHM_0_1=mean(XHM,2); %adsorbed HM (mg/L)
average_SHM_0_1=mean(SHM,2); %dissolved HM (mg/L) %compare with PNEC (0.19 µg/L)
average_XSS_0_1=mean(XSS,2); %suspended solids (mg/L)
average_SED_0_1=mean(SED,2);

%% plot of average Cd concentration in water
map=colormap(hot);
colors= [map(1,:);map(50,:); map(80,:); map(120,:); map(160,:);map(200,:);map(230,:)];

t = datetime(2021,1,1)+hours(1:2160);
y=zeros(length(average_SHM_0_01),1);
y(:)=0.19;
figure(1)
plot(t,y, "LineStyle","--", "Color", "k", "LineWidth", 1.5 )
hold on
plot(t,average_SHM_0_1*10^3,"color",colors(1,:), "LineWidth", 1.5)
hold on
plot(t, average_SHM_0_05*10^3,"color",colors(2,:), "LineWidth", 1.5)
hold on
plot(t, average_SHM_0_03*10^3,"color",colors(3,:), "LineWidth", 1.5)
hold on
plot(t, average_SHM_0_01*10^3,"color",colors(4,:), "LineWidth", 1.5)
hold on
plot(t, average_SHM_0_005*10^3,"color",colors(5,:), "LineWidth", 1.5)
hold on
plot(t, average_SHM_0_002*10^3,"color",colors(6,:), "LineWidth", 1.5)
hold on
plot(t, average_SHM*10^3,"color",colors(7,:), "LineWidth", 1.5)
xlabel("Time", "FontSize", 12)
ylabel("Cd concentration (µg/L)", "FontSize", 12)
sgtitle("Average concentration of dissolved Cd in the water")
legend("PNEC","0.1 mg/L","0.05 mg/L","0.03 mg/L", "0.01 mg/L","0.005 mg/L", "0.002 mg/L", "0.001 mg/L", "location", "northwest")
ylim([0,0.30])

%%
% define parameters
ku = 0.005;%L/g.d
AE = 0.45;                                   
IR = 0.1; %g/g.d
ke = 0.016; %d^-1
g = 0.001665; %d^-1
SHM_0_001 = average_SHM*10^3; %dissolved (µg/L)
Cf_0_001 = (average_XHM*10^3)./(average_XSS*10^-3); %adsorbed (mg/L) --> µg/g
SHM_0_002 = average_SHM_0_002*10^3; %dissolved (µg/L)
Cf_0_002 = (average_XHM_0_002*10^3)./(average_XSS_0_002*10^-3); %adsorbed (mg/L) --> µg/g
SHM_0_005 = average_SHM_0_005*10^3; %dissolved (µg/L)
Cf_0_005 = (average_XHM_0_005*10^3)./(average_XSS_0_005*10^-3); %adsorbed (mg/L) --> µg/g
SHM_0_01 = average_SHM_0_01*10^3; %dissolved (µg/L)
Cf_0_01 = (average_XHM_0_01*10^3)./(average_XSS_0_01*10^-3); %adsorbed (mg/L) --> µg/g
SHM_0_03 = average_SHM_0_03*10^3; %dissolved (µg/L)
Cf_0_03 = (average_XHM_0_03*10^3)./(average_XSS_0_03*10^-3); %adsorbed (mg/L) --> µg/g
SHM_0_05 = average_SHM_0_05*10^3; %dissolved (µg/L)
Cf_0_05 = (average_XHM_0_05*10^3)./(average_XSS_0_05*10^-3); %adsorbed (mg/L) --> µg/g
SHM_0_1 = average_SHM_0_1*10^3; %dissolved (µg/L)
Cf_0_1 = (average_XHM_0_1*10^3)./(average_XSS_0_1*10^-3); %adsorbed (mg/L) --> µg/g
%SHM_1 = average_SHM_1*10^3; %dissolved (µg/L)
%Cf_1 = (average_XHM_1*10^3)./(average_XSS_1*10^-3); %adsorbed (mg/L) --> µg/g

%% 0.001 mg/L
initial =0;
Ct = zeros(length(average_SHM_0_1),1); %µg/g

for i = 1:length(average_SHM_0_1)
    Ct(i) = ku*SHM_0_001(i)+AE*IR*Cf_0_001(i)-(ke+g)*initial;
    initial = Ct(i);
end
Ct_0_001=Ct;

% 0.002 mg/L
initial =0;
Ct = zeros(length(average_SHM_0_1),1);

for i = 1:length(average_SHM_0_1)
    Ct(i) = ku*SHM_0_002(i)+AE*IR*Cf_0_002(i)-(ke+g)*initial;
    initial = Ct(i);
end
Ct_0_002=Ct;

% 0.005 mg/L
initial =0;
Ct = zeros(length(average_SHM_0_1),1);

for i = 1:length(average_SHM_0_1)
    Ct(i) = ku*SHM_0_005(i)+AE*IR*Cf_0_005(i)-(ke+g)*initial;
    initial = Ct(i);
end
Ct_0_005=Ct;

% 0.01 mg/L
initial =0;
Ct = zeros(length(average_SHM_0_1),1);

for i = 1:length(average_SHM_0_1)
    Ct(i) = ku*SHM_0_01(i)+AE*IR*Cf_0_01(i)-(ke+g)*initial;
    initial = Ct(i);
end
Ct_0_01=Ct;

% 0.03 mg/L
initial =0;
Ct = zeros(length(average_SHM_0_1),1);

for i = 1:length(average_SHM_0_1)
    Ct(i) = ku*SHM_0_03(i)+AE*IR*Cf_0_03(i)-(ke+g)*initial;
    initial = Ct(i);
end
Ct_0_03=Ct;

% 0.05 mg/L
initial =0;
Ct = zeros(length(average_SHM_0_1),1);

for i = 1:length(average_SHM_0_1)
    Ct(i) = ku*SHM_0_05(i)+AE*IR*Cf_0_05(i)-(ke+g)*initial;
    initial = Ct(i);
end
Ct_0_05=Ct;

% 0.1 mg/L
initial =0;
Ct = zeros(length(average_SHM_0_1),1);

for i = 1:length(average_SHM_0_1)
    Ct(i) = ku*SHM_0_1(i)+AE*IR*Cf_0_1(i)-(ke+g)*initial;
    initial = Ct(i);
end
Ct_0_1=Ct;

%% plot of concentration in fish

map=colormap(hot);
colors= [map(1,:);map(50,:); map(80,:); map(120,:); map(160,:);map(200,:);map(230,:)];

t2 = datetime(2021,1,1)+hours(1:2160);
y2=zeros(length(average_SHM_0_1),1);
y2(:)=0.05;
figure(2)
plot(t2,y2, "LineStyle","--", "Color", "k", "LineWidth", 1.5 )
hold on
plot(t2,Ct_0_1,"color",colors(1,:), "LineWidth", 1.5)
hold on
plot(t2, Ct_0_05,"color",colors(2,:), "LineWidth", 1.5)
hold on
plot(t2, Ct_0_03,"color",colors(3,:), "LineWidth", 1.5)
hold on
plot(t2, Ct_0_01, "color",colors(4,:),"LineWidth", 1.5)
hold on
plot(t2, Ct_0_005, "color",colors(5,:),"LineWidth", 1.5)
hold on
plot(t2, Ct_0_002, "color",colors(6,:),"LineWidth", 1.5)
hold on
plot(t2, Ct_0_001, "color",colors(7,:),"LineWidth", 1.5)
xlabel("Time", "FontSize", 12)
ylabel("Cd concentration (µg/g)", "FontSize", 12)
sgtitle("Average concentration of Cd in the fish")
legend("limit", "0.1 mg/L","0.05 mg/L","0.03 mg/L", "0.01 mg/L","0.005 mg/L", "0.002 mg/L", "0.001 mg/L")
ylim([0,0.5])


%compare Ct with 0.05 mg/kg wet weight of fish = 0.05 µg/g


