function [MDP] = WSCT_performance(mdp, OPTIONS, N, N_agents)

N     = 32;
N_agents = 64;
for m = 1:N_agents
    
    % create structure array and solve
    %----------------------------------------------------------------------
    for i = 1:N
        MDP(i) = mdp;
    end
    rng(m)
    MDP  = WSCT_X_tutorial(MDP);
    
    % free energy and confidence
    %----------------------------------------------------------------------
    [F,Fu]   = spm_MDP_F(MDP);
    Fm(:,m)  = F(:);
    Fum(:,m) = Fu(:);
    
    % find run of correct responses
    %----------------------------------------------------------------------
    for i = 1:N
        if hit(MDP(i)); h(i,m) = 1; else, h(i,m) = 0;  end
    end

    % repeat with BMR
    %----------------------------------------------------------------------
    RDP = MDP;
    RDP = rmfield(RDP,{'u','o'});
    for i = 1:N
        RDP(i).a  = mda.a;
        RDP(i).a0 = mda.a0;
        RDP(i).s  = RDP(i).s(:,1);
    end
    rng(m)
    RDP      = spm_MDP_VB_X(RDP,OPT);
    [F,Fu]   = spm_MDP_F(RDP);
    Rm(:,m)  = F(:);
    Rum(:,m) = Fu(:);
    
    % look for instances of BMR
    %----------------------------------------------------------------------
    vA    = spm_vec(A{1});
    for i = 1:N
       c(i,1) = corr(vA,(spm_vec(RDP(i).a{1}) > 0));
    end
    bmr(:,m) = diff(c);
    
    % preferred locations
    %----------------------------------------------------------------------
    for i = 1:N
        if hit(RDP(i)); r(i,m) = 1; else, r(i,m) = 0;  end
    end
    
    % find run of correct responses
    %----------------------------------------------------------------------
    for i = 1:N
        if hit(RDP(i)); r(i,m) = 1; else, r(i,m) = 0;  end
    end
    
    % find non-learners
    %----------------------------------------------------------------------
    R = spm_conv(r,2,0);
    H = spm_conv(h,2,0);
    
    % show results
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 11');clf
    subplot(4,1,1)
    b = bar(mean(R(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.0),hold on
    b = bar(mean(H(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.8),hold off
    xlabel('trial'), ylabel('probability of correct'), axis([1/2 (N + 1/2) 1/3 1]);
    title('Average performance','Fontsize',16)
    
    subplot(4,1,2)
    spm_plot_ci(mean(Rm'),var(Rm')), hold on
    plot(mean(Fm'),'r'),              hold off
    xlabel('trial'), ylabel('free energy'); set(gca,'XLim',[1 N])
    title('Average free energy','Fontsize',16)
    
    subplot(4,1,3)
    spm_plot_ci(mean(Rum'),var(Rum')), hold on
    plot(mean(Fum'),'r'),              hold off
    xlabel('trial'), ylabel('confidence');  set(gca,'XLim',[1 N])
    title('Average confidence','Fontsize',16)
    
    % show individual performance
    %----------------------------------------------------------------------
    subplot(4,1,4), image(32*(r'))
    xlabel('trial'), ylabel('subject')
    title('Aha moments','Fontsize',16)
    
    % plot model updates
    %----------------------------------------------------------------------
    hold on
    for i = 1:m
        j = find(bmr(:,i) > 0) + 1;
        try, plot(j(1),i,'.m','MarkerSize',32), end
        try, plot(j(2),i,'.r','MarkerSize',32), end
        j = find(bmr(:,i) < 0);
        plot(j, (j - j + i),'.b','MarkerSize',32)
    end
    hold off, drawnow
    
    save paper
    
end
