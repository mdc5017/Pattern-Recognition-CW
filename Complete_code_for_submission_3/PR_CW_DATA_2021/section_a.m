%%%%%%%%% Pattern Recognition Coursework %%%%%%%%%%%

% Section: Section A - Data Preparation
% Start Date: 24/Feb/2021

%% Section A - 1


t = linspace(1,1000,1000);

%% Object: kitchen sponge
figure;

% Trial 1
load('PR_CW_DATA_2021/kitchen_sponge_114_01_HOLD')

subplot(3,4,1); plot(t,F0pdc); title('Trial 1 - Pressure'); set(gca,'Fontsize',18);
subplot(3,4,2); plot(t,F0pac(2,:)); title('Trial 1 - Vibrations'); set(gca,'Fontsize',18);
subplot(3,4,3); plot(t,F0tdc);  title('Trial 1 - Temperature');set(gca,'Fontsize',18);
subplot(3,4,4); plot(t,F0Electrodes);  title('Trial 1 - Electrodes');set(gca,'Fontsize',18);

% Trial 4
load('PR_CW_DATA_2021/kitchen_sponge_114_04_HOLD')

subplot(3,4,5); plot(t,F0pdc); title('Trial 4 - Pressure');set(gca,'Fontsize',18);
subplot(3,4,6); plot(t,F0pac(2,:));title('Trial 4 - Vibrations'); set(gca,'Fontsize',18);
subplot(3,4,7); plot(t,F0tdc); title('Trial 4 - Temperature');set(gca,'Fontsize',18);
subplot(3,4,8); plot(t,F0Electrodes); title('Trial 4 - Electrodes');set(gca,'Fontsize',18);

% Trial 7
load('PR_CW_DATA_2021/kitchen_sponge_114_07_HOLD')

subplot(3,4,9); plot(t,F0pdc); title('Trial 7 - Pressure');set(gca,'Fontsize',18);
subplot(3,4,10); plot(t,F0pac(2,:)); title('Trial 7 - Vibrations');set(gca,'Fontsize',18); 
subplot(3,4,11); plot(t,F0tdc); title('Trial 7 - Temperature');set(gca,'Fontsize',18);
subplot(3,4,12); plot(t,F0Electrodes); title('Trial 7 - Electrodes');set(gca,'Fontsize',18);
sgtitle('\fontsize{22}Kitchen sponge')


%% Object: flour sack
figure;

% Trial 1
load('PR_CW_DATA_2021/flour_sack_410_01_HOLD')

subplot(3,4,1); plot(t,F0pdc); title('Trial 1 - Pressure');set(gca,'Fontsize',18);
subplot(3,4,2); plot(t,F0pac(2,:)); title('Trial 1 - Vibrations'); set(gca,'Fontsize',18);
subplot(3,4,3); plot(t,F0tdc); title('Trial 1 - Temperature');set(gca,'Fontsize',18);
subplot(3,4,4); plot(t,F0Electrodes); title('Trial 1 - Electrodes');set(gca,'Fontsize',18);

% Trial 4
load('PR_CW_DATA_2021/flour_sack_410_04_HOLD')

subplot(3,4,5); plot(t,F0pdc); title('Trial 4 - Pressure');set(gca,'Fontsize',18);
subplot(3,4,6); plot(t,F0pac(2,:)); title('Trial 4 - Vibrations');set(gca,'Fontsize',18); 
subplot(3,4,7); plot(t,F0tdc); title('Trial 4 - Temperature');set(gca,'Fontsize',18);
subplot(3,4,8); plot(t,F0Electrodes); title('Trial 4 - Electrodes');set(gca,'Fontsize',18);

% Trial 7
load('PR_CW_DATA_2021/flour_sack_410_06_HOLD')

subplot(3,4,9); plot(t,F0pdc); title('Trial 6 - Pressure');set(gca,'Fontsize',18);
subplot(3,4,10); plot(t,F0pac(2,:)); title('Trial 6 - Vibrations'); set(gca,'Fontsize',18);
subplot(3,4,11); plot(t,F0tdc); title('Trial 6 - Temperature');set(gca,'Fontsize',18);
subplot(3,4,12); plot(t,F0Electrodes); title('Trial 6 - Electrodes');set(gca,'Fontsize',18);

sgtitle('\fontsize{22}Flour sack')

%% Object: steel vase
figure;

% Trial 1
load('PR_CW_DATA_2021/steel_vase_702_01_HOLD')

subplot(3,4,1); plot(t,F0pdc); title('Trial 1 - Pressure');set(gca,'Fontsize',18);
subplot(3,4,2); plot(t,F0pac(2,:)); title('Trial 1 - Vibrations'); set(gca,'Fontsize',18);
subplot(3,4,3); plot(t,F0tdc); title('Trial 1 - Temperature');set(gca,'Fontsize',18);
subplot(3,4,4); plot(t,F0Electrodes); title('Trial 1 - Electrodes');set(gca,'Fontsize',18);

% Trial 4
load('PR_CW_DATA_2021/steel_vase_702_04_HOLD')

subplot(3,4,5); plot(t,F0pdc); title('Trial 4 - Pressure');set(gca,'Fontsize',18);
subplot(3,4,6); plot(t,F0pac(2,:)); title('Trial 4 - Vibrations'); set(gca,'Fontsize',18);
subplot(3,4,7); plot(t,F0tdc); title('Trial 4 - Temperature');set(gca,'Fontsize',18);
subplot(3,4,8); plot(t,F0Electrodes); title('Trial 4 - Electrodes');set(gca,'Fontsize',18);

% Trial 7
load('PR_CW_DATA_2021/steel_vase_702_07_HOLD')

subplot(3,4,9); plot(t,F0pdc); title('Trial 7 - Pressure');set(gca,'Fontsize',18);
subplot(3,4,10); plot(t,F0pac(2,:)); title('Trial 7 - Vibrations'); set(gca,'Fontsize',18);
subplot(3,4,11); plot(t,F0tdc); title('Trial 7 - Temperature');set(gca,'Fontsize',18);
subplot(3,4,12); plot(t,F0Electrodes); title('Trial 7 - Electrodes');set(gca,'Fontsize',18);


sgtitle('\fontsize{22}Steel vase')


%% Quantitative Analysis

% average across trials for all the objects

% acrylic
sum_press=zeros(10,1000); sum_vib=zeros(10,1000); sum_temp=zeros(10,1000); 
for i=1:9
    name_file = sprintf('PR_CW_DATA_2021/acrylic_211_0%d_HOLD',i);
    load(name_file);
    sum_press(i,:)=F0pdc;
    sum_vib(i,:) = F0pac(2,:);
    sum_temp(i,:)= F0tdc;
end

load('PR_CW_DATA_2021/acrylic_211_10_HOLD');
sum_press(10,:)=F0pdc;
sum_vib(10,:) = F0pac(2,:);
sum_temp(10,:)= F0tdc;

avg_press = mean(sum_press);
avg_vib = mean(sum_vib);
avg_temp = mean(sum_temp);

acrylic_avg = [avg_press; avg_vib; avg_temp];

% black foam
sum_press=zeros(10,1000); sum_vib=zeros(10,1000); sum_temp=zeros(10,1000); 
for i=1:9
    name_file = sprintf('PR_CW_DATA_2021/black_foam_110_0%d_HOLD',i);
    load(name_file);
    sum_press(i,:)=F0pdc(1:1000);
    sum_vib(i,:) = F0pac(2,1:1000);
    sum_temp(i,:)= F0tdc(1:1000);
end

load('PR_CW_DATA_2021/black_foam_110_10_HOLD');
sum_press(10,:)=F0pdc;
sum_vib(10,:) = F0pac(2,:);
sum_temp(10,:)= F0tdc;

avg_press = mean(sum_press);
avg_vib = mean(sum_vib);
avg_temp = mean(sum_temp);

black_avg = [avg_press; avg_vib; avg_temp];


% car sponge
sum_press=zeros(10,1000); sum_vib=zeros(10,1000); sum_temp=zeros(10,1000); 
for i=1:9
    name_file = sprintf('PR_CW_DATA_2021/car_sponge_101_0%d_HOLD',i);
    load(name_file);
    sum_press(i,:)=F0pdc;
    sum_vib(i,:) = F0pac(2,:);
    sum_temp(i,:)= F0tdc;
end

load('PR_CW_DATA_2021/car_sponge_101_10_HOLD');
sum_press(10,:)=F0pdc;
sum_vib(10,:) = F0pac(2,:);
sum_temp(10,:)= F0tdc;

avg_press = mean(sum_press);
avg_vib = mean(sum_vib);
avg_temp = mean(sum_temp);

car_avg = [avg_press; avg_vib; avg_temp];



% flour sack
sum_press=zeros(10,1000); sum_vib=zeros(10,1000); sum_temp=zeros(10,1000); 
for i=1:9
    name_file = sprintf('PR_CW_DATA_2021/flour_sack_410_0%d_HOLD',i);
    load(name_file);
    sum_press(i,:)=F0pdc;
    sum_vib(i,:) = F0pac(2,:);
    sum_temp(i,:)= F0tdc;
end

load('PR_CW_DATA_2021/flour_sack_410_10_HOLD');
sum_press(10,:)=F0pdc;
sum_vib(10,:) = F0pac(2,:);
sum_temp(10,:)= F0tdc;

avg_press = mean(sum_press);
avg_vib = mean(sum_vib);
avg_temp = mean(sum_temp);

flour_avg = [avg_press; avg_vib; avg_temp];


% kitchen sponge
sum_press=zeros(10,1000); sum_vib=zeros(10,1000); sum_temp=zeros(10,1000); 
for i=1:9
    name_file = sprintf('PR_CW_DATA_2021/kitchen_sponge_114_0%d_HOLD',i);
    load(name_file);
    sum_press(i,:)=F0pdc;
    sum_vib(i,:) = F0pac(2,:);
    sum_temp(i,:)= F0tdc;
    %sum_elec(i,:)= F0Electrodes;
end

load('PR_CW_DATA_2021/kitchen_sponge_114_10_HOLD');
sum_press(10,:)=F0pdc;
sum_vib(10,:) = F0pac(2,:);
sum_temp(10,:)= F0tdc;

avg_press = mean(sum_press);
avg_vib = mean(sum_vib);
avg_temp = mean(sum_temp);

kitchen_avg = [avg_press; avg_vib; avg_temp];

% steel vase
sum_press=zeros(10,1000); sum_vib=zeros(10,1000); sum_temp=zeros(10,1000); 
for i=1:9
    name_file = sprintf('PR_CW_DATA_2021/steel_vase_702_0%d_HOLD',i);
    load(name_file);
    sum_press(i,:)=F0pdc;
    sum_vib(i,:) = F0pac(2,:);
    sum_temp(i,:)= F0tdc;
end

load('PR_CW_DATA_2021/steel_vase_702_10_HOLD');
sum_press(10,:)=F0pdc;
sum_vib(10,:) = F0pac(2,:);
sum_temp(10,:)= F0tdc;

avg_press = mean(sum_press);
avg_vib = mean(sum_vib);
avg_temp = mean(sum_temp);

steel_avg = [avg_press; avg_vib; avg_temp];

%% Finding variance across objects for different parameters

all_press = [acrylic_avg(1,:); black_avg(1,:); car_avg(1,:); flour_avg(1,:); kitchen_avg(1,:); steel_avg(1,:)];

var_press = var(all_press,0,1);

all_vib = [acrylic_avg(2,:); black_avg(2,:); car_avg(2,:); flour_avg(2,:); kitchen_avg(2,:); steel_avg(2,:)];

var_vib = var(all_vib,0,1);

all_temp = [acrylic_avg(3,:); black_avg(3,:); car_avg(3,:); flour_avg(3,:); kitchen_avg(3,:); steel_avg(3,:)];

var_temp = var(all_temp,0,1);

% find time instance with highest variance 
[~, t_press] = max(var_press);
[~, t_vib] = max(var_vib);
[~, t_temp] = max(var_temp);
 
%% Section A - 2
%For finger F0, this is the time instance at which the greatest variance
%was observed in data points between different objects. i.e., it is easier
%to distinguish objects from each other at this time instance. 
t = 50;

%this contains the starting of the names of the different objects
nameObj = ["acrylic_211_" "black_foam_110_" "car_sponge_101_" "flour_sack_410_" "kitchen_sponge_114_" "steel_vase_702_"];
%Below are the array of zeros that will later be filled with sampled
%values. The first dimension is the object type in the order in which it is
%called in nameObj and the second dimension is the trial number. For the
%electrodes data, the third dimension is the 19 electrode values at time
%instance 50
sampledData = zeros(6,3,10);
sampledDataPressure = zeros(6, 10);
sampledDataVibration = zeros(6, 10);
sampledDataTemperature = zeros(6, 10);
sampledDataElectrodes = zeros(6, 19, 10); 

%Another way of sampling, this time by object for PVT. The first dimension
%is each objects in the order of nameObj. The second dimension is the
%values of P, V and T respectively. The third dimension is the number of
%trials.
sampledDataObjects_PVT = zeros(6, 3, 10);

%Another way of sampling, this time by object, for electrodes. 
%The first dimension is each objects in the order of nameObj.
%The first dimension is each of the 19 electrodes respectively. 
%The second dimension is the value at t = 999 for each of the 19 trials.
sampledDataObjects_E = zeros(6, 19, 10);

%The first dimension is each of the 19 electrodes respectively. 
%The second dimension is the value at t = 999 for each of the 19 trials.

%The outer for loop loops through nameObj to loop through the different
%objects
for i = 1:6
    %The inner for loop loops through each trial for each object
    for j = 1:10
        %This if statement ensures that if the trial number is a single
        %digit, the string version will contain a '0' in front, just like
        %how the objects are named.
        if j<10
            num = strcat(num2str(0), num2str(j));
        else
            num = num2str(j);
        end
       %This forms the name for each trial of each object and loads it into
       %the workspace
       nameObjtemp = strcat(nameObj(1,i), num, "_HOLD.mat");
       load(nameObjtemp)
       
       %Sampling the data at time instance 999, at which the maximum
       %differentation can be found between objects, as measured by
       %variance. Data sampled for pressure, vibration and temperature
       sampledDataPressure(i, j) = F0pdc(1, t);
       sampledDataVibration(i, j) = F0pac(2, t);
       sampledDataTemperature(i, j) = F0tdc(1, t);
       sampledDataElectrodes(i, :, j) = F0Electrodes(:, t);
       
       %Sampling based on objects
       sampledDataObjects_E(i, :, j) = F0Electrodes(:, t);
       sampledDataObjects_PVT(i, 1, j) = F0pdc(1, t);
       sampledDataObjects_PVT(i, 2, j) = F0pac(1, t);
       sampledDataObjects_PVT(i, 3, j) = F0tdc(1, t);

    end
end

sampledDataAcrylic_PVT = sampledDataObjects_PVT(1, :, :);
sampledDataBF_PVT = sampledDataObjects_PVT(2, :, :);
sampledDataCS_PVT = sampledDataObjects_PVT(3, :, :);
sampledDataFS_PVT = sampledDataObjects_PVT(4, :, :);
sampledDataKS_PVT = sampledDataObjects_PVT(5, :, :);
sampledDataSV_PVT = sampledDataObjects_PVT(6, :, :);

sampledDataAcrylic_E = sampledDataObjects_E(1, :, :);
sampledDataBF_E = sampledDataObjects_E(2, :, :);
sampledDataCS_E = sampledDataObjects_E(3, :, :);
sampledDataFS_E = sampledDataObjects_E(4, :, :);
sampledDataKS_E = sampledDataObjects_E(5, :, :);
sampledDataSV_E = sampledDataObjects_E(6, :, :);

sampledDataPVT(:, 1, :) = sampledDataPressure; 
sampledDataPVT(:, 2, :) = sampledDataVibration; 
sampledDataPVT(:, 3, :) = sampledDataTemperature;

x = reshape(sampledDataPressure', [1, 60]);
y = reshape(sampledDataVibration', [1, 60]);
z = reshape(sampledDataTemperature', [1, 60]);
sampledDataPVTCon = zeros(4,60);
%The First row gives details on which object is being referred to 
sampledDataPVTCon(1, 1:10) = 1;
sampledDataPVTCon(1, 11:20) = 2;
sampledDataPVTCon(1, 21:30) = 3;
sampledDataPVTCon(1, 31:40) = 4;
sampledDataPVTCon(1, 41:50) = 5;
sampledDataPVTCon(1, 51:60) = 6;
sampledDataPVTCon(2, :) = x;
sampledDataPVTCon(3, :) = y;
sampledDataPVTCon(4, :) = z;

%Do the same for electrodes
sampledDataECon = zeros(20, 60);

sampledDataECon(1, 1:10) = 1;
sampledDataECon(1, 11:20) = 2;
sampledDataECon(1, 21:30) = 3;
sampledDataECon(1, 31:40) = 4;
sampledDataECon(1, 41:50) = 5;
sampledDataECon(1, 51:60) = 6;

for i = 2:20
    sampledDataECon(i, :) = reshape(squeeze(sampledDataElectrodes(:, i-1, :))', [1, 60]);
end

%save F0_PVT_V6_t50.mat sampledDataPVTCon sampledDataECon

%% Section A - 3
load("F0_PVT_V5_t50.mat")

%X axis will be pressure, y axis will be vibration and z axis will be
%temperature
sampledDataPressure = squeeze(sampledDataPVT(:, 1, :)); 
sampledDataVibration = squeeze(sampledDataPVT(:, 2, :)); 
sampledDataTemperature = squeeze(sampledDataPVT(:, 3, :));


p = reshape(sampledDataPressure', [1, 60]);
v = reshape(sampledDataVibration', [1, 60]);
t = reshape(sampledDataTemperature', [1, 60]);

colors = ['r', 'g', 'b', 'm', 'k', 'c'];
n=0;

figure;
for i=1:10:60
    n=n+1;
    color = colors(n);
    scatter3(p(1, i:(i+9)),v(1,i:(i+9)),t(1, i:(i+9)),15,'filled',color); grid on; hold on;
end

xlabel('Pressure'); ylabel('Vibration'); zlabel('Temperature');
title('3D PVT scatter plot - Perspective 1')
legend('acrylic','black foam','car sponge', 'flour sack', 'kitchen sponge','steel vase')