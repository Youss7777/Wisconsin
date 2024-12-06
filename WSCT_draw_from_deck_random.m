function [MDP] = WSCT_draw_from_deck_random(MDP, deck)
for i=1:size(MDP, 2)
    rand_idx = randi([1, size(deck, 2)]);
    for feature=1:3
        MDP(i).D{feature} = deck{rand_idx}{feature};
    end
end
end