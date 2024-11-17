function WSCT_plot(MDP)

spm_figure('GetWin','WSCT'); clf    % display behavior

if iscell(MDP(1).X)
    Nf = numel(MDP(1).B);                 % number of hidden state factors
    Ng = numel(MDP(1).A);                 % number of outcome factors
else
    Nf = 1;
    Ng = 1;
end


Nt    = length(MDP);               % number of trials
Ne    = size(MDP(1).V,1) + 1;      % number of epochs per trial
Np    = size(MDP(1).V,2) + 1;      % number of policies
f_act = 10;                        % state factor for which to display sampled actions
timestep_of_trial_to_plot = 1;

for i = 1:Nt
    % assemble expectations of hidden states and outcomes
    %----------------------------------------------------------------------
    for j = 1:Ne
        for k = 1:Ne
            for f = 1:Nf
                try
                    x{f}{i,1}{k,j} = gradient(MDP(i).xn{f}(:,:,j,k)')';
                catch
                    x{f}{i,1}{k,j} = gradient(MDP(i).xn(:,:,j,k)')';
                end
            end
        end
    end
    act_prob(:,i) = MDP(i).P(:,:,:,:,:,:,:,:,:,:,timestep_of_trial_to_plot);
    act(:,i) = MDP(i).u(f_act,timestep_of_trial_to_plot);
    w(:,i) = mean(MDP(i).dn,2);
end

% Initial states and expected policies (habit in red)
%--------------------------------------------------------------------------
col   = {'r.','g.','b.','c.','m.','k.'};
subplot(5,1,1)
if Nt < 64
    MarkerSize = 24;
else
    MarkerSize = 16;
end

image(64*(1 - act_prob)),  hold on

plot(act,col{3},'MarkerSize',MarkerSize)

try
    plot(Np*(1 - act_prob(Np,:)),'r')
end
try
    E = spm_softmax(spm_cat({MDP.e}));
    plot(Np*(1 - E(end,:)),'r:')
end
title('Action selection and action probabilities')
xlabel('Trial'),ylabel('Action'), hold off
yticks([1, 2])
yticklabels({'Card 1', 'Card 2'})

%============================================

% % posterior beliefs about hidden states
% figure;
% for f = 1:Nf
%     subplot(Nf, 1, f)
%     image(64*(1 - X{gf(f)})), hold on
%     if size(X{gf(f)},1) > 128
%         spm_spy(X{gf(f)},16,1)
%     end
%     plot(MDP.s(gf(f),:),'.c','MarkerSize',16), hold off
%     if f < 2
%         title(sprintf('Hidden states - %s',MDP.label.factor{gf(f)}));
%     else
%         title(MDP.label.factor{gf(f)});
%     end
% end
% 
% % outcome and preferences
% figure;
% for g = 1:Ng
%     subplot(Ng, 1, g)
%     image(64*(1 - C{gg(g)})), hold on
%     if size(C{gg(g)},1) > 128
%         spm_spy(C{gg(g)},16,1)
%     end
%     plot(MDP.o(gg(g),:),'.c','MarkerSize',16), hold off
%     if f < 2
%         title(sprintf('Hidden states - %s',MDP.label.modality{gg(g)}));
%     else
%         title(MDP.label.modality{gg(g)});
%     end
% end
% 
% 
% % posterior beliefs about control states
% %--------------------------------------------------------------------------
% figure;
% Nu     = find(Nu);
% Np     = length(Nu);
% for f  = 1:Np
%     subplot(Np,1,f)
%     P  = MDP.P;
%     if Nf > 1
%         ind     = 1:Nf;
%         for dim = 1:Nf
%             if dim ~= ind(Nu(f))
%                 P = sum(P,dim);
%             end
%         end
%         P = squeeze(P);
%     end
% 
%     % display
%     %----------------------------------------------------------------------
%     image(64*(1 - P)), hold on
%     plot(MDP.u(Nu(f),:),'.c','MarkerSize',16), hold off
%     if f < 2
%         title(sprintf('Action - %s',MDP.label.factor{Nu(f)}));
%     else
%         title(MDP.label.factor{Nu(f)});
%     end
%     set(gca,'XTickLabel',{});
%     set(gca,'XTick',1:size(X{1},2));
%     set(gca,'YTick',1:numel(MDP.label.action{Nu(f)}));
%     set(gca,'YTickLabel',MDP.label.action{Nu(f)});
% 
% 
% end