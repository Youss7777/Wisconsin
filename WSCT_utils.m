

s = 0;
for i=1:N
ss = sum(Exp_MDP(i).u(6,2)==MDP(i).u(6,2));
s = s + ss;
end
res = s/N

for i=1:N
exp_act(i) = Exp_MDP(i).u(6,2);
sim_act(i) = MDP(i).u(6,2);
end

figure
subplot(2, 1, 1), bar(1:N,exp_act),  xlabel('trial'), spm_axis tight, title('exp_act')
subplot(2, 1, 2), bar(1:N,sim_act),  xlabel('trial'), spm_axis tight, title('sim_act')


for i=1:N
corr(i)=MDP(i).o(5,3);
end
sum(corr == 2)/N;


for i=1:N
    pg  = 1;                                    % prior precision
    qg  = MDP(i).w;                                % posterior precision
    pu  = spm_softmax(MDP(i).G*diag(qg));          % prior policies
    qu  = spm_softmax(MDP(i).F + MDP(i).G*diag(qg));  % posterior policies
    fu =  sum(qu.*log(qu));                    % confidence
    fs = -sum(qu.*MDP(i).F);                      % free energy of states
    fq = -sum(qu.*log(pu));                    % free energy of policies
    fg = qg/pg - log(qg);                      % free energy of precision
    ftot = fu + fs + fq + fg;
    ffs{i} = fs;
    ffq{i} = fq;
    ffg{i} = fg;
    Ff{i} = MDP(i).F;
    Gg{i} = MDP(i).G;
    Qu{i} = qu;
    Ppu{i} = pu;
    ffu{i} = fu;
    Ftot(i) = sum(ftot);
    Fu(i) = sum(fu);
    Fs(i) = sum(fs);
    Fq(i) = sum(fq);
    Fg(i) = sum(fg);
    F(i) = sum(sum(MDP(i).F));
    G(i) = sum(sum(MDP(i).G));

end

for i=10:15
disp(['trial' num2str(i)])
disp('qu(i)'), Qu{i}
disp('sum(Qu{i}.*log(Qu{i}))'), sum(Qu{i}.*log(Qu{i}))
disp('Confidence'), sum(sum(Qu{i}.*log(Qu{i})))
end

for i=21:26
disp(['trial' num2str(i)])
disp('EFE(i)'), MDP(i).G
disp('FE(i)'), MDP(i).F
end


figure
subplot(6, 1, 1), plot(1:N,Fu),  xlabel('trial'), spm_axis tight, title('confidence')
subplot(6, 1, 2), plot(1:N,Fs),  xlabel('trial'), spm_axis tight, title('free energy of states')
subplot(6, 1, 3), plot(1:N,Fq),  xlabel('trial'), spm_axis tight, title('free energy of policies')
subplot(6, 1, 4), plot(1:N,Fg),  xlabel('trial'), spm_axis tight, title('free energy of precision')
subplot(6, 1, 5), plot(1:N,F),  xlabel('trial'), spm_axis tight, title('Summed F')
subplot(6, 1, 6), plot(1:N,G),  xlabel('trial'), spm_axis tight, title('Summed G')


figure
plot(1:N, movmean(Fu, 2)), xlabel('trial'), spm_axis tight, title('smooth confidence')