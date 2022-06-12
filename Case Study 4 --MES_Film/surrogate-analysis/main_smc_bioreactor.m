%% Main Sequential Monte Carlo 
close all
clearvars
%% Initialize
uqlab
% Fix the random seed to get repeatable results
rng(100)
%Have precalculated realizations
load('saved_reals.mat')
tind=[7,14,20,27,34,42,49];

%For creating state-dependent noise
std_factor = [0.05,0.05,5e-2];
std_floor = [0.00, 0.00, 0.00];
% std_factor = [2.5e-2];
% std_floor = [ 0.00];

nparam = 6;
nq=3; %Number of measured variables
nq0 = 3; %all vars
nT =length(tind);  %Number of time measurements
tind=tind';
npoints = length(Tlist);
N = 2e4 ;%Particles

load('train_data_os.mat')
Xdata = X;
Ydata = Y;


tic
fact_truth = [1.075,0.925,0.925,1.065,0.975,0.915,1.075,1.065,0.935];
for j=1:nparam
rng(1000)
xs_real(j) =((max(Xdata(:,j+nq))+min(Xdata(:,j+nq)))/2).*fact_truth(j); %eg %0.309   
end
%[y1_true,y2_true,y3_true,y4_true,xs_rel_true,pctime]=forward_evals(1,Tlist,npoints);
[y1_true, y2_true,y3_true]=evalforward(xs_real,1,Tlist,npoints);
Y_true = [y1_true(tind),y2_true(tind),y3_true(tind)];
% Y_true = [y4_true(tind)];


E = zeros(nT,nq);
for i = 1:nT
    for j = 1:nq
        %E(i,j) = (std_factor(j)*Y_true(i,j)+std_floor(j))*rand(1,1);
        E(i,j) = (std_factor(j)*Y_true(i,j)+std_floor(j))*randn;
        %((0.15-0.05).*rand(1,1)+0.05)*Y_true(i,j);%
    end
end
Y = Y_true + E;

close all
for j=1:nq
    figure(j)
    plot(1:length(tind), Y(:,j),'o','MarkerFaceColor','r');
    hold on
    plot(1:length(tind), Y_true(:,j),'-','MarkerFaceColor','b');

end


close all


%% Pre-Compute all Trajectories
Y_predicted = zeros(N,nq,nT);

%% SMC Algorithm
%Define particles
particlesexist = 1;

if particlesexist==0
X=zeros(N,nparam);
load('train_data_os.mat')
Xdata = X;
Ydata = Y;
for nn=1:N
fprintf('Computing trajectory')
nn
tic

for j=1:nparam
rng('shuffle')
xs(j) =(max(Xdata(:,j+nq0))-min(Xdata(:,j+nq0))).*rand(1,1) + min(Xdata(:,j+nq0)); %eg %0.309   
end
X(nn,:) = xs;

[y1_pred, y2_pred,y3_pred,y4_pred]=evalforward(xs,1,Tlist,npoints);
Y_predicted(nn,1,:) = y1_pred(tind);
Y_predicted(nn,2,:) = y2_pred(tind);
Y_predicted(nn,3,:) = y3_pred(tind);
% Y_predicted(nn,1,:) = y1_pred(tind);
% Y_predicted(nn,2,:) = y4_pred(tind);

% Y_predicted(nn,1,:)

prediction_time = toc;
% % figure(10)
% % plot(reshape(Y_predicted(nn,1,:),[1,nT]));
% % hold all
% % figure(20)
% % plot(reshape(Y_predicted(nn,2,:),[1,nT]));
% % hold all

end
else
load('saved_reals.mat')
%tind=tindfull;
X = xs_rel;
q1r = reshape(y1,[size(y1,1),size(y1,3)])';
q2r = reshape(y2,[size(y2,1),size(y2,3)])';
q3r = reshape(y3,[size(y3,1),size(y3,3)])';

q1red = q1r(1:N,tind);
q2red = q2r(1:N,tind);
q3red = q3r(1:N,tind);

Y_predicted(:,1,:) = q1red;
Y_predicted(:,2,:) = q2red;
Y_predicted(:,3,:) = q3red;
% Y_predicted(:,1,:) = q4red;


end


W = ones(N,1)/N;
% Initialize data storage
X_data = zeros(N,nparam,nT+1);
W_data = zeros(N,1,nT+1);
X_data(:,:,1) = X(1:N,:);
W_data(:,:,1) = W;

%Store original simulation matrices
X_g = X;
Y_g = Y_predicted;  



% Print statement
fprintf('Starting Sequential Monte Carlo \n')
% Loop over time


for k = 1:nT
    % Print starting statement
    fprintf('starting sample %g / %g...',k,nT)
    tic
   % Calculate likelihood for samples
   %%%%%%%%%%%%%%%%%%%%%%
   tic
   sim_ind = zeros(N,1);
   for j=1:N
     for i=1:N
     if all(X_g(i,:) == X(j,:))
        sim_ind(j)=i;
     end
     end
   end
   %matchingtime=toc
   %mean = evalmetaModel(x, myMM, k, T); %matrix (N by nq), each column corresponds to a specific
   mean_val = Y_predicted(sim_ind,:,k);
   y=Y_true(k,:);
   for i = 1:nq
   std(:,i) = std_factor(i)*mean_val(:,i)+std_floor(i);
   end
   L = prod(normpdf(repmat(y,[N,1]),mean_val,std),2);
   %%%%%%%%%%%%%%%%%%%%%%
    % Reweighting -- recursively update weights and normalize
    W = W.*L;
    W = W./sum(W);
       
    % Resample if effective number of particles is less than 85% of original
    Neff = 1/sum(W.^2);
    if Neff < 0.85*N
        fprintf('need to resample...')
        [X,W,~] = systematicResample(X,W);
    end
    
    % Store data
    X_data(:,:,k+1) = X(1:N,:);
    W_data(:,:,k+1) = W;
    
    % Print time
    fprintf('took %g seconds\n',toc)
       
end

fprintf('\n')
%}
%% Mean of Posterior 
% Get smc estimates as mean of posterior 
for jj=1:nparam
X_opt(jj) = sum(X(:,jj))/N;
end

custom_cols = [0.7294    0.1412    0.0510;...
     0.0824    0.5216    0.0157;...
     0.1882    0.0902    0.8196;...
     0.9020    0.5294    0.1098;...
     0.6078    0.1882    0.8902]; 
deepred    = custom_cols(1,:);
deepgreen  = custom_cols(2,:);
deepblue   = custom_cols(3,:);
deeporange = custom_cols(4,:);
deeppurple = custom_cols(5,:);


close all
for j=1:nparam
ub(j) =max(Xdata(:,j+nq));
lb(j) =min(Xdata(:,j+nq));%0.309   
end
X_prior = X_data(:,:,1);
X_post = X_data(:,:,end);
X_true = xs_real;
%% Generate plots
% Plot histograms of parameter posterior on top of prior
X_lb = lb;
X_ub = ub;
figure (15) 
for i = 1:nparam
    subplot(nparam, nparam, i+nparam*(i-1))
    histogram(X_prior(:,i), 'BinLimits', [X_lb(i), X_ub(i)], 'BinWidth', (X_ub(i)-X_lb(i))/10, 'FaceColor', deepgreen, 'EdgeColor', 'none', 'Normalization', 'probability')
    hold on;
    numbin=histogram(X_post(:,i), 'BinLimits', [X_lb(i), X_ub(i)], 'BinWidth', (X_ub(i)-X_lb(i))/10, 'FaceColor', deepblue, 'EdgeColor', 'none', 'Normalization', 'probability');
   
   % axis([0.95*X_lb(i), 1.05*X_ub(i), 0, 1.5*max(numbin.Values)])
    line([X_true(i), X_true(i)], [0,1], 'LineWidth', 3, 'Color', deepred);
    line([X_opt(i), X_opt(i)], [0,1], 'LineWidth', 3, 'Color', 'k');
    set(gca,'TickLabelInterpreter', 'latex')
   
    %ax = gca;
    %ax.XAxis.LineWidth = 2;
    %ay = gca;
    %ay.YAxis.LineWidth = 2;
    if i == 1
        set(gca,'xtick',[])    
        set(gca,'ytick',[]) 
        set(gca,'Linewidth',2)
            set(gca,'TickLabelInterpreter', 'latex')

    end
    if i > 1 && i < nparam
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'Linewidth',2)
            set(gca,'TickLabelInterpreter', 'latex')

   
    end
    if i == nparam
        set(gca,'ytick',[])
        set(gca,'Linewidth',2)
            set(gca,'TickLabelInterpreter', 'latex')

    end
    for j = 1:i-1
        subplot(nparam, nparam, j+nparam*(i-1))
        hold on;
        scatter(X_prior(:,j),X_prior(:,i),2,deepgreen)
        %scatter(X_post(:,j),X_post(:,i),'b.')
        scatter(X_post(:,j),X_post(:,i),2,deepblue)
        %box on 
        axis([0.95*X_lb(j), 1.05*X_ub(j), 0.95*X_lb(i), 1.05*X_ub(i)])
            set(gca,'TickLabelInterpreter', 'latex')

        
        if j == 1 && i < nparam
            set(gca,'xtick',[])
            set(gca,'Linewidth',2)

        end
        if j > 1 && i < nparam
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            set(gca,'Linewidth',2)

        end
        if j > 1 && i == nparam
            set(gca,'ytick',[])
            set(gca,'Linewidth',2)
        else
         set(gca,'Linewidth',2)
                set(gca,'TickLabelInterpreter', 'latex')
        end
    end
end
set(gcf,'color','w')

save('inference.mat','X_prior','X_post', 'X_true', 'X_opt')
return