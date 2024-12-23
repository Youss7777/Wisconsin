function [MDP] = WSCT_performance(MDP, OPTIONS, N, N_agents, f_state)

hit   = @(MDP) any(MDP.o(5,3) == 2) & ~any(MDP.o(5,3) == 1);


for m = 1:N_agents
    
    % create structure array and solve
    %---------------------------------------------------------------------
    rng(m)
    [MDP_noBMR, ~]  = WSCT_X_tutorial(MDP);
    
    % free energy and confidence
    %----------------------------------------------------------------------
    [F_noBMR,Fu_noBMR,~,~,~,~,p_best_policy_noBMR,~]   = WSCT_spm_MDP_F(MDP_noBMR, f_state);
    Fm_noBMR(:,m)  = F_noBMR(:);
    Fum_noBMR(:,m) = p_best_policy_noBMR(:);
    
    % find run of correct responses
    %----------------------------------------------------------------------
    for i = 1:N
        if hit(MDP_noBMR(i)); hit_noBMR(i,m) = 1; else, hit_noBMR(i,m) = 0;  end
    end

    % repeat with BMR
    %----------------------------------------------------------------------
    rng(m)
    [MDP_BMR, OPTIONS] = WSCT_X_tutorial(MDP,OPTIONS);
    if OPTIONS.BMR0.applied == false
        OPTIONS.BMR0.trials = [OPTIONS.BMR0.trials 0];
    end
    [F_BMR,Fu_BMR,~,~,~,~,p_best_policy_BMR,~]   = WSCT_spm_MDP_F(MDP_BMR, f_state);
    Fm_BMR(:,m)  = F_BMR(:);
    Fum_BMR(:,m) = p_best_policy_BMR(:);

    % find analytical solving trial    

    
    % look for instances of BMR
    %----------------------------------------------------------------------
    % vA    = spm_vec(A{1});
    % for i = 1:N
    %    c(i,1) = corr(vA,(spm_vec(RDP(i).a{1}) > 0));
    % end
    % bmr(:,m) = diff(c);
    clear bmr
    for i=1:m
        bmr(:, i) = OPTIONS.BMR0.trials(i);
    end
    % preferred locations
    %----------------------------------------------------------------------
    for i = 1:N
        if hit(MDP_BMR(i)); hit_BMR(i,m) = 1; else, hit_BMR(i,m) = 0;  end
    end
    
    % find run of correct responses
    %----------------------------------------------------------------------
    for i = 1:N
        if hit(MDP_BMR(i)); hit_BMR(i,m) = 1; else, hit_BMR(i,m) = 0;  end
    end
    
    % find non-learners
    %----------------------------------------------------------------------
    Hit_BMR = spm_conv(hit_BMR,2,0);
    Hit_noBMR = spm_conv(hit_noBMR,2,0);
    
    % show results
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 11');clf
    subplot(4,1,1)
    b = bar(mean(Hit_BMR(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.0),hold on    % black bars for BMR
    b = bar(mean(Hit_noBMR(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.8),hold off % grey bars for no BMR
    xlabel('trial'), ylabel('probability of correct'), axis([1/2 (N + 1/2) 1/4 1]);
    title('Average performance','Fontsize',16)
    
    subplot(4,1,2)
    spm_plot_ci(mean(Fm_BMR'),var(Fm_BMR')), hold on        % bar plot with bar as mean and error bars as variance
    plot(mean(Fm_noBMR'),'r'),              hold off        %
    xlabel('trial'), ylabel('free energy'); set(gca,'XLim',[1 N])   % label compresses the bar in one trial
    title('Average free energy','Fontsize',16)
    
    subplot(4,1,3)
    spm_plot_ci(mean(Fum_BMR'),var(Fum_BMR')), hold on
    plot(mean(Fum_noBMR'),'r'),              hold off
    xlabel('trial'), ylabel('confidence');  set(gca,'XLim',[1 N])
    title('Average confidence','Fontsize',16)
    
    % show individual performance
    %----------------------------------------------------------------------
    subplot(4,1,4), image(32*(hit_BMR'))    % grey for correct trials, rest (incorrect) in black
    xlabel('trial'), ylabel('subject')
    title('Aha moments','Fontsize',16)
    
    % plot model updates
    %----------------------------------------------------------------------
    hold on
    for i = 1:m
        if bmr(:,i) > 0
            plot(bmr(:,i), i, '.r','MarkerSize',16)
        else

    
        end
        %try, plot(j(1),i,'.m','MarkerSize',32), end
        %try, plot(j(2),i,'.r','MarkerSize',32), end
        %j = find(bmr(:,i) < 0);
        %plot(j, (j - j + i),'.b','MarkerSize',32)
    end
    hold off, drawnow
    
    % reset BMR structure
    OPTIONS.BMR0.jmin = [];
    OPTIONS.BMR0.rFmin = [];
    OPTIONS.BMR0.rD = [];
    OPTIONS.BMR0.sD = [];
end

return
end
