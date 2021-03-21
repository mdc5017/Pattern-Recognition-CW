%%%%%%%%% Pattern Recognition Coursework %%%%%%%%%%%

% Section: Section B - Principal Component Analysis
% Start Date: 1/Mar/2021
clc
clear
close all
load('F0_PVT_V6_t50.mat')


%% Covariance Matrix, Eigenvalues and Eigenvectors
PVT = sampledDataPVTCon(2:4,:)';

% 1. Stardardise the Data
N = zscore(PVT);

% 2. Compute the Covariance Matrix
C = cov(N);

% 3. Calculate Eigenvectors V and Eigenvalues D

[V, D] = eig(C);

% 4. Select a Feature Vector
Feature_vector = fliplr(V);


%% Replot standardised data with principal components

[coeff,score] = pca(N);
colors = ['r', 'g', 'b', 'm', 'k', 'c'];
n=0;
figure;

% Data Points
for i=1:10:60
    n=n+1;
    color = colors(n);
    scatter3(N(i:(i+9),1),N(i:(i+9),2),N(i:(i+9),3),30,'filled',color); grid on; hold on;
end
xlabel('Pressure'); ylabel('Vibration'); zlabel('Temperature');

% PC vectors
for i=1:3
    quiver3(0, Feature_vector(i,1), Feature_vector(i,2), Feature_vector(i,3),2,'Linewidth',3); hold on;
end

legend('acrylic','black foam','car sponge', 'flour sack', 'kitchen sponge','steel vase', 'PC 1','PC 2','PC 3')
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'Fontsize',18)
title('Standardised data with PC')


%% Reduce data to 2D and replot

vbls = ["Pressure","Vibration","Temperature"]; 

figure;
h = biplot(coeff(:,1:2),'Scores',score(:,1:2),'VarLabels',vbls); 
i=0;
set(gca,'Fontsize',18)

for k = 10:69 
    if rem(k,10) == 0 
        i=i+1;
    end
    color = colors(i);
    h(k).MarkerEdgeColor = color;  % Specify color for the observations    
    h(k).MarkerSize=30; 
end
title('PVT data reduced to 2D')

%% Show how data is distributed across all principal components

figure;
for i=1:3
    Y = N*Feature_vector(:,i);
    subplot(1,3,i);
    n=0;
    for j=1:10:60
        n=n+1;
        color = colors(n);
        scatter(Y(j:(j+9)),zeros(10,1),30,'filled',color); grid on; hold on;
    end
    legend('acrylic','black foam','car sponge', 'flour sack', 'kitchen sponge','steel vase')
    title(strcat('Principal Component #',num2str(i)))
    set(gca,'Fontsize',18);
end


%% PCA of the electrodes data
electrodes = sampledDataECon(2:20,:)';
N = zscore(electrodes);
[coeff,score,latent,~,explained] = pca(N);

figure;g=plot(explained,'.-');
set(g,'LineWidth',2,'MarkerSize',25)
set(gca,'XTick',1:19)
title('Explained variance by principal components')
xlabel('Principal component number')
ylabel('Fraction of variation explained [%]')
set(gca,'Fontsize',18);

%% Electrode Data on 3 Principal Components

lab = ["","","","","","","","","","","","","","","","","","",""];

figure;
h1 = biplot(coeff(:,1:3),'Scores',score(:,1:3),'VarLabels',lab);

i=0;
n=0;
set(gca,'Fontsize',18);
for k = 58:117 
    if rem(n,10) == 0 
        i=i+1;
    end
    n=n+1;
    color = colors(i);
    h1(k).MarkerEdgeColor = color;   % Specify color for the observations    
    h1(k).MarkerSize=20; 
    
end

