function [DCM] = WSCT_model_inversion(mdp, OPTIONS, Exp_MDP, field_params)
% MODEL INVERSION
%==============================================================
DCM.mdp = mdp;
DCM.OPTIONS = OPTIONS;
DCM.U = {Exp_MDP.o};
DCM.Y = {Exp_MDP.u};
DCM.field = field_params;

disp('Starting...')
DCM = WSCT_Estimate_parameters(DCM); % Run the parameter estimation function
disp('Finished!')

disp('Saving DCM...')
filename = ['WSCT_DCM_' num2str(participant) '_' isBMR '.mat'];
save(filename, 'DCM');
disp('DCM saved!')

subplot(2,2,3)
xticklabels(DCM.field),xlabel('Parameter')
subplot(2,2,4)
xticklabels(DCM.field),xlabel('Parameter')

field = fieldnames(DCM.M.pE);
for i = 1:length(field)
    if strcmp(field{i},'eta')
        prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
        posterior(i) = 1/(1+exp(-DCM.Ep.(field{i}))); 
    elseif strcmp(field{i},'omega')
        prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
        posterior(i) = 1/(1+exp(-DCM.Ep.(field{i})));
    else
        prior(i) = exp(DCM.M.pE.(field{i}));
        posterior(i) = exp(DCM.Ep.(field{i}));
    end
end

figure, set(gcf,'color','white')
subplot(2,1,1),hold on
title('Means')
bar(prior,'FaceColor',[.5,.5,.5]),bar(posterior,0.5,'k')
xlim([0,length(prior)+1]),set(gca, 'XTick', 1:length(prior)),set(gca, 'XTickLabel', DCM.field)
legend({'Prior','Posterior'})
hold off
subplot(2,1,2)
imagesc(DCM.Cp),caxis([0 1]),colorbar
title('(Co-)variance')
set(gca, 'XTick', 1:length(prior)),set(gca, 'XTickLabel', DCM.field)
set(gca, 'YTick', 1:length(prior)),set(gca, 'YTickLabel', DCM.field)