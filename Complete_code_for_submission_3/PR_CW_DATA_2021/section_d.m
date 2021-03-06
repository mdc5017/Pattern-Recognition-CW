%%%%%%%%% Pattern Recognition Coursework %%%%%%%%%%%

% Section: Section D - Clustering
% Start Date: 8/Mar/2021
clc
clear
close all
load('F0_PVT_V6_t50.mat')


%%
data = zscore(sampledDataPVTCon');
X = zscore(sampledDataPVTCon(2:4,:)');

[idx, C] = kmeans(X,6); % add 'Distance','correlation' for different distance metric

figure;
subplot(2,1,1)
colors = ['r', 'g', 'b', 'm', 'k', 'c'];

% Data Points
for i = 1:6
  color = colors(i);
  scatter3(X(idx==i,1),X(idx==i,2),X(idx==i,3),30,'filled',color); 
  grid on; hold on; 
  
end

% Cluster Centroids
plot3(C(:,1),C(:,2),C(:,3),'o','Color','b','MarkerSize',10,'MarkerFaceColor','y'); hold on;
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6','Centroids',...
       'Location','NW')
title ('Cluster Assignments and Centroids')
xlabel('Pressure'); ylabel('Vibration'); zlabel('Temperature');
set(gca,'Fontsize',18)
hold off

%% Plot combining shapes and colours for different objects

k_data = [X idx];
subplot(2,1,2)
shapes = ['o', '+', '*', 'd', 'x', 's'];
n=0;
j=0;

% Data Points
for i=1:60
    
    if rem(n,10) == 0 
        j=j+1;
    end
    color = colors(j);
    shape = shapes(k_data(i,4));
    plot3(X(i,1),X(i,2),X(i,3),shape,'Color',color,'MarkerSize',10,'MarkerFaceColor',color); hold on;grid on;
    n=n+1;
end

% Cluster Centroids
plot3(C(:,1),C(:,2),C(:,3),'o','Color','b','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','b'); hold on;
xlabel('Pressure'); ylabel('Vibration'); zlabel('Temperature');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'Fontsize',18)
title('K means on PVT Data')

%% D - 2

load('PCA_Electrodes.mat')
rng(3)
%Split up the dataset into training/test data sets using 60/40 split for
%each electrode, ensuring that 60% of the trials for each object go to
%training and 40% of trials for each object go to test. The first 3 PCA axes
%are being used
training_data = zeros(36, 3);
training_output = zeros(36, 1);
test_data = zeros(24, 3);
test_output = zeros(24, 1);

for e = 1:3
    for i = 1:6
        training_data((6*i)-5:6*i, e) = score((10*i)-9:(10*i)-4, e);
        training_output((6*i)-5:6*i) = i;
        test_data((4*i)-3:4*i, e) = score((10*i)-3:10*i, e);
        test_output((4*i)-3:4*i) = i;
    end
end

%% Bagging

%Do bagging
B = TreeBagger(150, training_data, training_output, 'OOBPrediction', 'On');
%%
%Part a - Plot the OOB error to decide the number of trees that will be
%used
figure;
subplot(3,1,1)
oobErrorBaggedEnsemble = oobError(B);
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';
title("OOB Classification error over the number of grown trees for 3 PCA dimensions")

%From the graph above, I think I will work with 25 trees as the OOB error
%seems to plateau at around that point before dipping even further. I think
%if we choose a tree number that is higher than this, then the model risks
%overfitting, as it starts to memorize the training data

%% Part b - Visualize 2 of your generated decision trees
%Do bagging
rng(3)
B1 = TreeBagger(30, training_data, training_output, 'OOBPrediction', 'On');

%B = TreeBagger(30, training_data, training_output, 'OOBPrediction', 'On');
view(B1.Trees{1}, 'Mode', 'graph')
view(B1.Trees{2}, 'Mode', 'graph')

%% Part c - running the trained model with test data

YFIT = predict(B1, test_data);

YFIT = cell2mat(YFIT);
YFIT_d = zeros(24,1);
for i = 1:24
   YFIT_d(i,1) = str2double(YFIT(i));
end
%This is the confusion matrix
C = confusionmat(test_output, YFIT_d);

subplot(3,1,2)
confusionchart(C)
title("Confusion Matrix Chart")
%% Doing bagging with all 19 PCA dimensions
rng(3)
%Split up the dataset into training/test data sets using 60/40 split for
%each electrode, ensuring that 60% of the trials for each object go to
%training and 40% of trials for each object go to test. The first 3 PCA axes
%are being used
training_data = zeros(36, 19);
training_output = zeros(36, 1);
test_data = zeros(24, 19);
test_output = zeros(24, 1);

for e = 1:19
    for i = 1:6
        training_data((6*i)-5:6*i, e) = score((10*i)-9:(10*i)-4, e);
        training_output((6*i)-5:6*i) = i;
        test_data((4*i)-3:4*i, e) = score((10*i)-3:10*i, e);
        test_output((4*i)-3:4*i) = i;
    end
end

%% Bagging

%Do bagging
B2 = TreeBagger(150, training_data, training_output, 'OOBPrediction', 'On');rng(3)
%%
%Part a - Plot the OOB error to decide the number of trees that will be
%used
subplot(3,1,3)
oobErrorBaggedEnsemble = oobError(B2);
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';
title("OOB Classification error over the number of grown trees for 19 PCA Dimensions")