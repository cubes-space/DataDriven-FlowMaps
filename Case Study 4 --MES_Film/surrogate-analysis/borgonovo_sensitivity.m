%% Borgonovo Sensitivity

clearvars
clc
uqlab

load('saved_reals.mat')

nstates = 3; 
nparams = 6;
tind=[49];
%tind=[2:65]';
npoints = length(tind);
nrels = size(y1,3);

sens1 = zeros(nparams,npoints);
sens2 = zeros(nparams,npoints);
sens3 = zeros(nparams,npoints);


for j=1:npoints
i=tind(j);
% i
BorgonovoOpts.Type = 'Sensitivity';
BorgonovoOpts.Method = 'Borgonovo';
BorgonovoOpts.Borgonovo.Method = 'HistBased';
%BorgonovoOpts.Borgonovo.NClasses = 20;
%BorgonovoOpts.Borgonovo.Overlap = 0.15;

values_rel = reshape(y1(i,1,:),[nrels,1]) ;
BorgonovoOpts.Borgonovo.Sample.X = xs_rel;
BorgonovoOpts.Borgonovo.Sample.Y = values_rel;
BorgonovoAnalysis{j,1}= uq_createAnalysis(BorgonovoOpts);
sens1(:,j) = BorgonovoAnalysis{j,1}.Results.Delta;


values_rel = reshape(y2(i,1,:),[nrels,1]) ;
BorgonovoOpts.Borgonovo.Sample.X = xs_rel;
BorgonovoOpts.Borgonovo.Sample.Y = values_rel;
BorgonovoAnalysis{j,2}= uq_createAnalysis(BorgonovoOpts);
sens2(:,j) = BorgonovoAnalysis{j,2}.Results.Delta;


values_rel = reshape(y3(i,1,:),[nrels,1]) ;
BorgonovoOpts.Borgonovo.Sample.X = xs_rel;
BorgonovoOpts.Borgonovo.Sample.Y = values_rel;
BorgonovoAnalysis{j,3}= uq_createAnalysis(BorgonovoOpts);
sens3(:,j) = BorgonovoAnalysis{j,3}.Results.Delta;

% values_rel = reshape(y4(i,1,:),[nrels,1]) ;
% BorgonovoOpts.Borgonovo.Sample.X = xs_rel;
% BorgonovoOpts.Borgonovo.Sample.Y = values_rel;
% BorgonovoAnalysis{j,4}= uq_createAnalysis(BorgonovoOpts);
% sens4(:,j) = BorgonovoAnalysis{j,4}.Results.Delta;
end


ccol(1,:) = [0.7;0.08;0.08]; %wine
ccol(2,:) = [0.11;0.11;0.61]; %dark blue
ccol(3,:) = [0.05;0.39;0.09]; %olive
%ccol(4,:) = [0.9;0.01;0.01]; %red

%}
close all
lwi=3.0
figure(1000)
for j=1:nparams
semilogy(linspace(1,length(tind),length(tind)),sens1(j,1:end),'-.','LineWidth',lwi)
hold on 
end

figure(2000)
for j=1:nparams
semilogy(linspace(1,length(tind),length(tind)),sens2(j,1:end),'-.','LineWidth',lwi)
hold on 
end

figure(3000)
for j=1:nparams
semilogy(linspace(2,length(tind),length(tind)-1),sens3(j,2:end),'-','LineWidth',lwi)
hold on 
end


save('sensitivities_results.mat', 'sens1','sens2','sens3','Tlist')
