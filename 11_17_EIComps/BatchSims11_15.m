similV2 = zeros(21,21);
for ind = 1:21
    for ind2 = 1:21
        similV2(ind,ind2) = findSim(ALLTS2{ind,ind2});
    end
end