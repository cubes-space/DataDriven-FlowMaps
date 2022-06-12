results_exist = 1;
trajectories_used =  [100,150,200,250,300,350];

if results_exist == 0
clc
close all
clearvars
Nrepeats=5;
trajectories_used = [100,150,200,250,300,350];
averaged_error=zeros(Nrepeats,length(trajectories_used))
close all
for ii=1:Nrepeats   
mean_error_ps = [];
for k=1:length(trajectories_used)
    Ntrajectories = trajectories_used(k);
    [pck, y1t_list,y2t_list,y3t_list,y1p_list,...
        y2p_list,y3p_list,time_list,mte] = main_pce(Ntrajectories);
mean_error_ps =[mean_error_ps, mte];
end
averaged_error(ii,:) = mean(mean_error_ps,1);
figure(101001)
plot(trajectories_used, averaged_error(ii,:))
hold on
end

end

for jj=1:Nrepeats
    plot(trajectories_used, averaged_error(jj,:))
    hold on
end


mave = zeros(1,length(trajectories_used))
vave = zeros(1,length(trajectories_used))

for i=1:length(trajectories_used)
    mave(i)=mean(averaged_error(:,i))
    vave(i)=sqrt(var(averaged_error(:,i)))
end



close all
figure(1)
plot(trajectories_used,mave,'-')
hold on
plot(trajectories_used,mave+vave,'--')
hold on
plot(trajectories_used,mave-vave,'--')

save('average_error_results.mat', 'averaged_error', 'mave', 'vave')

