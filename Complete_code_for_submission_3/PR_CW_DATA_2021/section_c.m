%%%%%%%%% Pattern Recognition Coursework %%%%%%%%%%%

% Section: Section C - Linear Disctiminant Analysis
% Start Date: 1/Mar/2021
clc
clear
close all
%%

%Note: black foam is 2 and car sponge is 3 in the first row of the .mat
%file loaded below
load('F0_PVT_V6_t50.mat')

%This helps find the indices of values of PVT for both black foam and car
%sponge
objects_indices = find(sampledDataPVTCon(1,:)==2|sampledDataPVTCon(1,:)==3);

%This extracts the pressure, vibration and temperature values for both
%objects
objects_pressure = sampledDataPVTCon(2,objects_indices);
objects_vibration = sampledDataPVTCon(3,objects_indices);
objects_temperature = sampledDataPVTCon(4,objects_indices);

%Linear discriminant analysis

%Step 1 - standardize the data to 0 mean and standard deviation 1
mean_p = mean(objects_pressure);
mean_v = mean(objects_vibration);
mean_t = mean(objects_temperature);

std_p = std(objects_pressure);
std_v = std(objects_vibration);
std_t = std(objects_temperature);

pressure_stand = (objects_pressure-mean_p)/std_p;
vibration_stand = (objects_vibration-mean_v)/std_v;
temperature_stand = (objects_temperature-mean_t)/std_t;

%Step 2 - Calculate means (each class and overall)

%these are the means for pressure
meanP_bf = mean(pressure_stand(1:10));
meanP_cs = mean(pressure_stand(11:20));
meanP_overall = mean(pressure_stand);

%these are the means for vibration
meanV_bf = mean(vibration_stand(1:10));
meanV_cs = mean(vibration_stand(11:20));
meanV_overall = mean(vibration_stand);

%these are the means for temperature
meanT_bf = mean(temperature_stand(1:10));
meanT_cs = mean(temperature_stand(11:20));
meanT_overall = mean(temperature_stand);

%Step 3 - Compute within-class scatter matrix

%This is the raw data concatenated depending on which graph we are looking
%at
X_pv = [vibration_stand; pressure_stand];
X_pt = [temperature_stand; pressure_stand];
X_tv = [vibration_stand; temperature_stand];

%These are the means concatenated depending on which graph we are looking
%at
m_bf_pv = [meanV_bf; meanP_bf];
m_cs_pv = [meanV_cs; meanP_cs];

m_bf_pt = [meanT_bf;meanP_bf];
m_cs_pt = [meanT_cs;meanP_cs];

m_bf_tv = [meanV_bf;meanT_bf];
m_cs_tv = [meanV_cs;meanT_cs];

%Within-class scatter matrix for Pressure vs vibration
tempPVBF = (X_pv(:, 1:10)-m_bf_pv)*(X_pv(:, 1:10)-m_bf_pv)';
tempPVCS = (X_pv(:, 11:20)-m_cs_pv)*(X_pv(:, 11:20)-m_cs_pv)';
Sw_PV = tempPVBF + tempPVCS;

%Within-class scatter matrix for Pressure vs temperature
tempPTBF = (X_pt(:, 1:10)-m_bf_pt)*(X_pt(:, 1:10)-m_bf_pt)';
tempPTCS = (X_pt(:, 11:20)-m_cs_pt)*(X_pt(:, 11:20)-m_cs_pt)';
Sw_PT = tempPTBF + tempPTCS;

%Within-class scatter matrix for temperature vs vibration
tempTVBF = (X_tv(:, 1:10)-m_bf_tv)*(X_tv(:, 1:10)-m_bf_tv)';
tempTVCS = (X_tv(:, 11:20)-m_cs_tv)*(X_tv(:, 11:20)-m_cs_tv)';
Sw_TV = tempTVBF + tempTVCS;

%Step 4 - compute between-class scatter matrix
Sb_pv = (m_bf_pv - m_cs_pv)*(m_bf_pv - m_cs_pv)'; %Pressure vs Vibration
Sb_pt = (m_bf_pt - m_cs_pt)*(m_bf_pt - m_cs_pt)'; %Pressure vs Temperature
Sb_tv = (m_bf_tv - m_cs_tv)*(m_bf_tv - m_cs_tv)'; %Temperature vs Vibration

%Step 5 - combine the scatter matrices
S_pv = inv(Sw_PV)*Sb_pv; % P vs V
S_pt = inv(Sw_PT)*Sb_pt; % P vs T
S_tv = inv(Sw_TV)*Sb_tv; % T vs V

%Step 6 - Compute eigenvectors and eigenvalues
[V_pv, D_pv] = eig(S_pv); % P vs V
[V_pt, D_pt] = eig(S_pt); % P vs T
[V_tv, D_tv] = eig(S_tv); % T vs V

%Step 7 - sort the eigenvalues - going to do so by finding the column where
%the highest eigenvalue is in the diagonal matrix
[row_eig_pv col_eig_pv] = find(D_pv == max(D_pv, [], 'all')); % P vs V
[row_eig_pt col_eig_pt] = find(D_pt == max(D_pt, [], 'all')); % P vs T
[row_eig_tv col_eig_tv] = find(D_tv == max(D_tv, [], 'all')); % T vs V

%Step 8 - Create a feature vector
F_pv = V_pv(:, col_eig_pv); % P vs V
F_pt = V_pt(:, col_eig_pt); % P vs T
F_tv = V_tv(:, col_eig_tv); % T vs V

% Step 9 - Project the data 
X_pv_lda = (F_pv'*X_pv)'; % P vs V
X_pt_lda = (F_pt'*X_pt)'; % P vs T
X_tv_lda = (F_tv'*X_tv)'; % T vs V

%Plot the data for Pressure vs Vibration
colors = ['g', 'b'];
n=0;

figure;
subplot(3,2,1)
for i=1:10:20
    n=n+1;
    color = colors(n);
    scatter(vibration_stand(1,i:(i+9)),pressure_stand(1,i:(i+9)),15,'filled',color); grid on; hold on;
end

xlabel('Vibration'); ylabel('Pressure'); 

for i=1:1
    quiver(F_pv(1, i), F_pv(2, 1), 2,'Linewidth',1.5); hold on;
end
legend('black foam','car sponge', 'LD 1')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
title("LDA for Pressure vs Vibration graph")


subplot(3,2,2)
n = 0;
for i=1:10:20
    n=n+1;
    color1 = colors(n);
    scatter(X_pv_lda(i:(i+9),1),X_pv_lda(i:(i+9),1)*0,15,'filled',color1); grid on; hold on;
end

xlabel('LD 1'); 
title("Numberline showing LD 1")
legend('black foam','car sponge')

%Plot the data for Pressure vs Temperature
n=0;

subplot(3,2,3)
for i=1:10:20
    n=n+1;
    color = colors(n);
    scatter(temperature_stand(1,i:(i+9)),pressure_stand(1,i:(i+9)),15,'filled',color); grid on; hold on;
end

xlabel('Temperature'); ylabel('Pressure'); 

for i=1:1
    quiver(F_pt(1, i), F_pt(2, 1), 2,'Linewidth',1.5); hold on;
end
legend('black foam','car sponge', 'LD 1')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
title("LDA for Pressure vs Temperature graph")


subplot(3,2,4)
n = 0
for i=1:10:20
    n=n+1;
    color1 = colors(n);
    scatter(X_pt_lda(i:(i+9),1),X_pt_lda(i:(i+9),1)*0,15,'filled',color1); grid on; hold on;
end

xlabel('LD 1'); 
title("Numberline showing LD 1")
legend('black foam','car sponge')
%Plot the data for Temperature vs Vibration

n=0;

subplot(3,2,5)
for i=1:10:20
    n=n+1;
    color = colors(n);
    scatter(vibration_stand(1,i:(i+9)),temperature_stand(1,i:(i+9)),15,'filled',color); grid on; hold on;
end

xlabel('Vibration'); ylabel('Temperature'); 

for i=1:1
    quiver(F_tv(1, i), F_tv(2, 1), 2,'Linewidth',1.5); hold on;
end
legend('black foam','car sponge', 'LD 1')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
title("LDA for Temperature vs Vibration graph")


subplot(3,2,6)
n = 0;
for i=1:10:20
    n=n+1;
    color1 = colors(n);
    scatter(X_tv_lda(i:(i+9),1),X_tv_lda(i:(i+9),1)*0,15,'filled',color1); grid on; hold on;
end

xlabel('LD 1'); 
title("Numberline showing LD 1")
legend('black foam','car sponge')

%% C - part b

%Note: black foam is 2 and car sponge is 3 in the first row of the .mat
%file loaded below
load('F0_PVT_V6_t50.mat')

%This helps find the indices of values of PVT for both black foam and car
%sponge
objects_indices = find(sampledDataPVTCon(1,:)==2|sampledDataPVTCon(1,:)==3);

%This extracts the pressure, vibration and temperature values for both
%objects
objects_pressure = sampledDataPVTCon(2,objects_indices);
objects_vibration = sampledDataPVTCon(3,objects_indices);
objects_temperature = sampledDataPVTCon(4,objects_indices);

% 3D Linear Discriminant Analysis

%Step 1 - standardize the data to 0 mean and standard deviation 1
mean_p = mean(objects_pressure);
mean_v = mean(objects_vibration);
mean_t = mean(objects_temperature);

std_p = std(objects_pressure);
std_v = std(objects_vibration);
std_t = std(objects_temperature);

pressure_stand = (objects_pressure-mean_p)/std_p;
vibration_stand = (objects_vibration-mean_v)/std_v;
temperature_stand = (objects_temperature-mean_t)/std_t;

%Step 2 - Calculate means (each class and overall)

%these are the means for pressure
meanP_bf = mean(pressure_stand(1:10));
meanP_cs = mean(pressure_stand(11:20));
meanP_overall = mean(pressure_stand);

%these are the means for vibration
meanV_bf = mean(vibration_stand(1:10));
meanV_cs = mean(vibration_stand(11:20));
meanV_overall = mean(vibration_stand);

%these are the means for temperature
meanT_bf = mean(temperature_stand(1:10));
meanT_cs = mean(temperature_stand(11:20));
meanT_overall = mean(temperature_stand);

%Step 3 - Compute within-class scatter matrix

%This is the raw data concatenated together
X_pvt = [pressure_stand; vibration_stand; temperature_stand];

%These are the means concatenated together
m_bf_pvt = [meanP_bf; meanV_bf; meanT_bf];
m_cs_pvt = [meanP_cs; meanV_cs; meanT_cs];

%Within-class scatter matrix for Pressure vs vibration vs temperature
tempPVTBF = (X_pvt(:, 1:10)-m_bf_pvt)*(X_pvt(:, 1:10)-m_bf_pvt)';
tempPVTCS = (X_pvt(:, 11:20)-m_cs_pvt)*(X_pvt(:, 11:20)-m_cs_pvt)';
Sw_PVT = tempPVTBF + tempPVTCS;

%Step 4 - compute between-class scatter matrix for Pressure vs vibration vs temperature
Sb_pvt = (m_bf_pvt - m_cs_pvt)*(m_bf_pvt - m_cs_pvt)';

%Step 5 - combine the scatter matrices
S_pvt = inv(Sw_PVT)*Sb_pvt; % P vs V

%Step 6 - Compute eigenvectors and eigenvalues
[V_pvt, D_pvt] = eig(S_pvt); % P vs V vs T

%Step 7 and 8 - sort the eigenvalues and make feature vector - going to do so by finding the column where
%the highest eigenvalue is in the diagonal matrix. Then copy and remove
%that column from the matrix and run algorithm again to find LDA 2
[row_eig_pvt_1 col_eig_pvt_1] = find(D_pvt == max(D_pvt, [], 'all')); %This finds coordinates of the highest eigenvalue

%store the eigenvector corresponding to the largest eigenvalue in the
%feature vector - this is LDA 1
F_pvt = [V_pvt(:, col_eig_pvt_1)]; 

%Remove the columns corresponding to the highest eigenvalue in both the
%diagonal and eigenvector matrices
V_pvt(:, col_eig_pvt_1) = []; 
D_pvt(:, col_eig_pvt_1) = []; 

%Find the maximum of this new diagonal matrix that doesn't have the highest
%eigenvalue to find the second highest eigenvalue for LDA 2
[row_eig_pvt_2 col_eig_pvt_2] = find(D_pvt == max(D_pvt, [], 'all')); %This finds coordinates of the highest eigenvalue

%store the eigenvector corresponding to the largest eigenvalue in the
%feature vector - this is LDA 2
F_pvt = [F_pvt V_pvt(:, col_eig_pvt_2)];

% Step 9 - Project the data onto the LDA planes
X_pvt_lda = (F_pvt'*X_pvt)'; %This variable contains projections on both LDAs

%Plotting the Pressure vs Temperature graph
%Let's plot pressure vs vibration as a test to see if things work
colors = ['g', 'b'];

figure
subplot(3,1,1)
n = 0;
for i=1:10:20
    n=n+1;
    color1 = colors(n);
    scatter(X_pvt_lda(i:(i+9),1),X_pvt_lda(i:(i+9),1)*0,15,'filled',color1); grid on; hold on;
end

xlabel('LD 1'); 
title("Numberline showing LD 1")
legend('black foam','car sponge')

subplot(3,1,2)
n = 0;
for i=1:10:20
    n=n+1;
    color1 = colors(n);
    scatter(X_pvt_lda(i:(i+9),1),X_pvt_lda(i:(i+9),2),15,'filled',color1); grid on; hold on;
end

xlabel('LD 1'); ylabel('LD 2') 
title("LD 2 vs LD1")
legend('black foam','car sponge')

n=0;
subplot(3,1,3)
for i=1:10:20
    n=n+1;
    color = colors(n);
    scatter3(pressure_stand(1,i:(i+9)), vibration_stand(1,i:(i+9)), temperature_stand(1,i:(i+9)),15,'filled',color); grid on; hold on;
end

xlabel('Pressure'); ylabel('Vibration'); zlabel('Temperature');

for i=1:2
    quiver3(F_pvt(1, i), F_pvt(2, i), F_pvt(3, i), 2,'Linewidth',1.5); hold on;
end
legend('black foam','car sponge', 'LD 1', 'LD 2')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
title("LDA for Pressure vs Vibration vs Temperature graph - Perspective 1")
%% C - part d
%Note: kitchen sponge is 5 and steel vase is 6 in the first row of the .mat
%file loaded below
load('F0_PVT_V6_t50.mat')

%This helps find the indices of values of PVT for both kitchen sponge and
%steel vase
objects_indices = find(sampledDataPVTCon(1,:)==5|sampledDataPVTCon(1,:)==6);

%This extracts the pressure, vibration and temperature values for both
%objects
objects_pressure = sampledDataPVTCon(2,objects_indices);
objects_vibration = sampledDataPVTCon(3,objects_indices);
objects_temperature = sampledDataPVTCon(4,objects_indices);

%Linear discriminant analysis

%Step 1 - standardize the data to 0 mean and standard deviation 1
mean_p = mean(objects_pressure);
mean_v = mean(objects_vibration);
mean_t = mean(objects_temperature);

std_p = std(objects_pressure);
std_v = std(objects_vibration);
std_t = std(objects_temperature);

pressure_stand = (objects_pressure-mean_p)/std_p;
vibration_stand = (objects_vibration-mean_v)/std_v;
temperature_stand = (objects_temperature-mean_t)/std_t;

%Step 2 - Calculate means (each class and overall)

%these are the means for pressure
meanP_bf = mean(pressure_stand(1:10));
meanP_cs = mean(pressure_stand(11:20));
meanP_overall = mean(pressure_stand);

%these are the means for vibration
meanV_bf = mean(vibration_stand(1:10));
meanV_cs = mean(vibration_stand(11:20));
meanV_overall = mean(vibration_stand);

%these are the means for temperature
meanT_bf = mean(temperature_stand(1:10));
meanT_cs = mean(temperature_stand(11:20));
meanT_overall = mean(temperature_stand);

%Step 3 - Compute within-class scatter matrix

%This is the raw data concatenated depending on which graph we are looking
%at
X_pv = [vibration_stand; pressure_stand];
X_pt = [temperature_stand; pressure_stand];
X_tv = [vibration_stand; temperature_stand];

%These are the means concatenated depending on which graph we are looking
%at
m_bf_pv = [meanV_bf; meanP_bf];
m_cs_pv = [meanV_cs; meanP_cs];

m_bf_pt = [meanT_bf;meanP_bf];
m_cs_pt = [meanT_cs;meanP_cs];

m_bf_tv = [meanV_bf;meanT_bf];
m_cs_tv = [meanV_cs;meanT_cs];

%Within-class scatter matrix for Pressure vs vibration
tempPVBF = (X_pv(:, 1:10)-m_bf_pv)*(X_pv(:, 1:10)-m_bf_pv)';
tempPVCS = (X_pv(:, 11:20)-m_cs_pv)*(X_pv(:, 11:20)-m_cs_pv)';
Sw_PV = tempPVBF + tempPVCS;

%Within-class scatter matrix for Pressure vs temperature
tempPTBF = (X_pt(:, 1:10)-m_bf_pt)*(X_pt(:, 1:10)-m_bf_pt)';
tempPTCS = (X_pt(:, 11:20)-m_cs_pt)*(X_pt(:, 11:20)-m_cs_pt)';
Sw_PT = tempPTBF + tempPTCS;

%Within-class scatter matrix for temperature vs vibration
tempTVBF = (X_tv(:, 1:10)-m_bf_tv)*(X_tv(:, 1:10)-m_bf_tv)';
tempTVCS = (X_tv(:, 11:20)-m_cs_tv)*(X_tv(:, 11:20)-m_cs_tv)';
Sw_TV = tempTVBF + tempTVCS;

%Step 4 - compute between-class scatter matrix
Sb_pv = (m_bf_pv - m_cs_pv)*(m_bf_pv - m_cs_pv)'; %Pressure vs Vibration
Sb_pt = (m_bf_pt - m_cs_pt)*(m_bf_pt - m_cs_pt)'; %Pressure vs Temperature
Sb_tv = (m_bf_tv - m_cs_tv)*(m_bf_tv - m_cs_tv)'; %Temperature vs Vibration

%Step 5 - combine the scatter matrices
S_pv = inv(Sw_PV)*Sb_pv; % P vs V
S_pt = inv(Sw_PT)*Sb_pt; % P vs T
S_tv = inv(Sw_TV)*Sb_tv; % T vs V

%Step 6 - Compute eigenvectors and eigenvalues
[V_pv, D_pv] = eig(S_pv); % P vs V
[V_pt, D_pt] = eig(S_pt); % P vs T
[V_tv, D_tv] = eig(S_tv); % T vs V

%Step 7 - sort the eigenvalues - going to do so by finding the column where
%the highest eigenvalue is in the diagonal matrix
[row_eig_pv col_eig_pv] = find(D_pv == max(D_pv, [], 'all')); % P vs V
[row_eig_pt col_eig_pt] = find(D_pt == max(D_pt, [], 'all')); % P vs T
[row_eig_tv col_eig_tv] = find(D_tv == max(D_tv, [], 'all')); % T vs V

%Step 8 - Create a feature vector
F_pv = V_pv(:, col_eig_pv); % P vs V
F_pt = V_pt(:, col_eig_pt); % P vs T
F_tv = V_tv(:, col_eig_tv); % T vs V

% Step 9 - Project the data 
X_pv_lda = (F_pv'*X_pv)'; % P vs V
X_pt_lda = (F_pt'*X_pt)'; % P vs T
X_tv_lda = (F_tv'*X_tv)'; % T vs V

%Plot the data for Pressure vs Vibration
colors = ['k', 'c'];
n=0;

figure;
subplot(3,2,1)
for i=1:10:20
    n=n+1;
    color = colors(n);
    scatter(vibration_stand(1,i:(i+9)),pressure_stand(1,i:(i+9)),15,'filled',color); grid on; hold on;
end

xlabel('Vibration'); ylabel('Pressure'); 

for i=1:1
    quiver(F_pv(1, i), F_pv(2, 1), 2,'Linewidth',1.5); hold on;
end
legend('kitchen sponge','steel vase', 'LD 1')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
title("LDA for Pressure vs Vibration graph")


subplot(3,2,2)
n = 0;
for i=1:10:20
    n=n+1;
    color1 = colors(n);
    scatter(X_pv_lda(i:(i+9),1),X_pv_lda(i:(i+9),1)*0,15,'filled',color1); grid on; hold on;
end

xlabel('LD 1'); 
title("Numberline showing LD 1")
legend('kitchen sponge','steel vase')

%Plot the data for Pressure vs Temperature
n=0;

subplot(3,2,3)
for i=1:10:20
    n=n+1;
    color = colors(n);
    scatter(temperature_stand(1,i:(i+9)),pressure_stand(1,i:(i+9)),15,'filled',color); grid on; hold on;
end

xlabel('Temperature'); ylabel('Pressure'); 

for i=1:1
    quiver(F_pt(1, i), F_pt(2, 1), 2,'Linewidth',1.5); hold on;
end
legend('kitchen sponge','steel vase', 'LD 1')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
title("LDA for Pressure vs Temperature graph")


subplot(3,2,4)
n = 0;
for i=1:10:20
    n=n+1;
    color1 = colors(n);
    scatter(X_pt_lda(i:(i+9),1),X_pt_lda(i:(i+9),1)*0,15,'filled',color1); grid on; hold on;
end

xlabel('LD 1'); 
title("Numberline showing LD 1")
legend('kitchen sponge','steel vase')
%Plot the data for Temperature vs Vibration

n=0;

subplot(3,2,5)
for i=1:10:20
    n=n+1;
    color = colors(n);
    scatter(vibration_stand(1,i:(i+9)),temperature_stand(1,i:(i+9)),15,'filled',color); grid on; hold on;
end

xlabel('Vibration'); ylabel('Temperature'); 

for i=1:1
    quiver(F_tv(1, i), F_tv(2, 1), 2,'Linewidth',1.5); hold on;
end
legend('kitchen sponge','steel vase', 'LD 1')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
title("LDA for Temperature vs Vibration graph")


subplot(3,2,6)
n = 0;
for i=1:10:20
    n=n+1;
    color1 = colors(n);
    scatter(X_tv_lda(i:(i+9),1),X_tv_lda(i:(i+9),1)*0,15,'filled',color1); grid on; hold on;
end

xlabel('LD 1'); 
title("Numberline showing LD 1")
legend('kitchen sponge','steel vase')

%Note: kitchen sponge is 5 and steel vase is 6 in the first row of the .mat
%file loaded below
load('F0_PVT_V6_t50.mat')

%This helps find the indices of values of PVT for both kitchen sponge and
%steel vase
objects_indices = find(sampledDataPVTCon(1,:)==5|sampledDataPVTCon(1,:)==6);

%This extracts the pressure, vibration and temperature values for both
%objects
objects_pressure = sampledDataPVTCon(2,objects_indices);
objects_vibration = sampledDataPVTCon(3,objects_indices);
objects_temperature = sampledDataPVTCon(4,objects_indices);

% 3D Linear Discriminant Analysis

%Step 1 - standardize the data to 0 mean and standard deviation 1
mean_p = mean(objects_pressure);
mean_v = mean(objects_vibration);
mean_t = mean(objects_temperature);

std_p = std(objects_pressure);
std_v = std(objects_vibration);
std_t = std(objects_temperature);

pressure_stand = (objects_pressure-mean_p)/std_p;
vibration_stand = (objects_vibration-mean_v)/std_v;
temperature_stand = (objects_temperature-mean_t)/std_t;

%Step 2 - Calculate means (each class and overall)

%these are the means for pressure
meanP_bf = mean(pressure_stand(1:10));
meanP_cs = mean(pressure_stand(11:20));
meanP_overall = mean(pressure_stand);

%these are the means for vibration
meanV_bf = mean(vibration_stand(1:10));
meanV_cs = mean(vibration_stand(11:20));
meanV_overall = mean(vibration_stand);

%these are the means for temperature
meanT_bf = mean(temperature_stand(1:10));
meanT_cs = mean(temperature_stand(11:20));
meanT_overall = mean(temperature_stand);

%Step 3 - Compute within-class scatter matrix

%This is the raw data concatenated together
X_pvt = [pressure_stand; vibration_stand; temperature_stand];

%These are the means concatenated together
m_bf_pvt = [meanP_bf; meanV_bf; meanT_bf];
m_cs_pvt = [meanP_cs; meanV_cs; meanT_cs];

%Within-class scatter matrix for Pressure vs vibration vs temperature
tempPVTBF = (X_pvt(:, 1:10)-m_bf_pvt)*(X_pvt(:, 1:10)-m_bf_pvt)';
tempPVTCS = (X_pvt(:, 11:20)-m_cs_pvt)*(X_pvt(:, 11:20)-m_cs_pvt)';
Sw_PVT = tempPVTBF + tempPVTCS;

%Step 4 - compute between-class scatter matrix for Pressure vs vibration vs temperature
Sb_pvt = (m_bf_pvt - m_cs_pvt)*(m_bf_pvt - m_cs_pvt)';

%Step 5 - combine the scatter matrices
S_pvt = inv(Sw_PVT)*Sb_pvt; % P vs V

%Step 6 - Compute eigenvectors and eigenvalues
[V_pvt, D_pvt] = eig(S_pvt); % P vs V vs T

%Step 7 and 8 - sort the eigenvalues and make feature vector - going to do so by finding the column where
%the highest eigenvalue is in the diagonal matrix. Then copy and remove
%that column from the matrix and run algorithm again to find LDA 2
[row_eig_pvt_1 col_eig_pvt_1] = find(D_pvt == max(D_pvt, [], 'all')); %This finds coordinates of the highest eigenvalue

%store the eigenvector corresponding to the largest eigenvalue in the
%feature vector - this is LDA 1
F_pvt = [V_pvt(:, col_eig_pvt_1)]; 

%Remove the columns corresponding to the highest eigenvalue in both the
%diagonal and eigenvector matrices
V_pvt(:, col_eig_pvt_1) = []; 
D_pvt(:, col_eig_pvt_1) = []; 

%Find the maximum of this new diagonal matrix that doesn't have the highest
%eigenvalue to find the second highest eigenvalue for LDA 2
[row_eig_pvt_2 col_eig_pvt_2] = find(D_pvt == max(D_pvt, [], 'all')); %This finds coordinates of the highest eigenvalue

%store the eigenvector corresponding to the largest eigenvalue in the
%feature vector - this is LDA 2
F_pvt = [F_pvt V_pvt(:, col_eig_pvt_2)];

% Step 9 - Project the data onto the LDA planes
X_pvt_lda = (F_pvt'*X_pvt)'; %This variable contains projections on both LDAs

%Plotting the Pressure vs Temperature graph
%Let's plot pressure vs vibration as a test to see if things work
colors = ['k', 'c'];

figure
subplot(3,1,1)
n = 0;
for i=1:10:20
    n=n+1;
    color1 = colors(n);
    scatter(X_pvt_lda(i:(i+9),1),X_pvt_lda(i:(i+9),1)*0,15,'filled',color1); grid on; hold on;
end

xlabel('LD 1'); 
title("Numberline showing LD 1")
legend('kitchen sponge','steel vase')

subplot(3,1,2)
n = 0;
for i=1:10:20
    n=n+1;
    color1 = colors(n);
    scatter(X_pvt_lda(i:(i+9),1),X_pvt_lda(i:(i+9),2),15,'filled',color1); grid on; hold on;
end

xlabel('LD 1'); ylabel('LD 2') 
title("LD 2 vs LD1")
legend('kitchen sponge','steel vase')

n=0;
subplot(3,1,3)
for i=1:10:20
    n=n+1;
    color = colors(n);
    scatter3(pressure_stand(1,i:(i+9)), vibration_stand(1,i:(i+9)), temperature_stand(1,i:(i+9)),15,'filled',color); grid on; hold on;
end

xlabel('Pressure'); ylabel('Vibration'); zlabel('Temperature');

for i=1:2
    quiver3(F_pvt(1, i), F_pvt(2, i), F_pvt(3, i), 2,'Linewidth',1.5); hold on;
end
legend('kitchen sponge','steel vase', 'LD 1', 'LD 2')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
title("LDA for Pressure vs Vibration vs Temperature graph - Perspective 1")