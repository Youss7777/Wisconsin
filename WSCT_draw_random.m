function [MDP] = WSCT_draw_random(MDP, Ns)
for i=1:size(MDP, 2)
    for feature=1:3
        MDP(i).D{feature} = zeros(Ns(feature), 1);
        rand_idx = randi([1, Ns(feature)]);
        MDP(i).D{feature}(rand_idx) = 1;
    end
end
end