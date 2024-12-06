function [count] = WSCT_nb_correct_actions(MDP, Exp_MDP)
ss = 0;
for i=1:size(MDP, 2)
    s = sum(Exp_MDP(i).u(6,2)==MDP(i).u(6,2));
    ss = ss + s;
end
count = ss/size(MDP, 2);