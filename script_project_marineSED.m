NETaddAssembly('DHI.Generic.MikeZero.DFS.dll');
import DHI.Generic.MikeZero.DFS.*;


dfsu3 = DfsFileFactory.DfsuFileOpen('project.m3fm - Result files/Cadmium_conc1305_0_03.dfsu');
%change according to concentration you want, then run this + calculate
%average

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
average_XHM=mean(XHM,2); %adsorbed HM (mg/L)
average_SHM=mean(SHM,2); %dissolved HM (mg/L)
average_XSS=mean(XSS,2); %suspended solids (mg/L)
average_SED=mean(SED,2);
%% 0.002 mg/L
average_XHM_0_002=mean(XHM,2); %adsorbed HM (mg/L)
average_SHM_0_002=mean(SHM,2); %dissolved HM (mg/L) 
average_XSS_0_002=mean(XSS,2); %suspended solids (mg/L)
average_SED_0_002=mean(SED,2);

%% 0.005 mg/L 
average_XHM_0_005=mean(XHM,2); %adsorbed HM (mg/L)
average_SHM_0_005=mean(SHM,2); %dissolved HM (mg/L)
average_XSS_0_005=mean(XSS,2); %suspended solids (mg/L)
average_SED_0_005=mean(SED,2);
%% 0.01 mg/L
average_XHM_0_01=mean(XHM,2); %adsorbed HM (mg/L)
average_SHM_0_01=mean(SHM,2); %dissolved HM (mg/L) 
average_XSS_0_01=mean(XSS,2); %suspended solids (mg/L)
average_SED_0_01=mean(SED,2);


%% plot of average Cd concentration in sediment

map=colormap(hot);
colors= [map(120,:);map(160,:); map(200,:); map(230,:)];

t = datetime(2021,1,1)+hours(1:2160);
y=zeros(length(average_SHM_0_002),1);
y(:)=0.19;
figure(1)
plot(t,y, "LineStyle","--", "Color", "k", "LineWidth", 1.5 )
hold on
plot(t, average_SED_0_01*10^3,"Color", colors(1,:), "LineWidth", 1.5)
hold on
plot(t, average_SED_0_005*10^3,"Color", colors(2,:), "LineWidth", 1.5)
hold on
plot(t, average_SED_0_002*10^3,"Color", colors(3,:), "LineWidth", 1.5)
hold on
plot(t, average_SED*10^3, "Color", colors(4,:),"LineWidth", 1.5)
xlabel("Time", "FontSize", 12)
ylabel("Cd concentration (µg/L)", "FontSize", 12)
sgtitle("Total concentration of Cd in the sediment")
legend("PNEC","0.01 mg/L","0.005 mg/L","0.002 mg/L", "0.001 mg/L")
ylim([0,0.8])

