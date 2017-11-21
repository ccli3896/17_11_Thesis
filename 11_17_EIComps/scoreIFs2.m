% Script to automatically run the integrate-and-fire network and score
% reactivation based on standard deviation (using procNets.m)

function [ALLTS,meansC,stC,stMeans,stMeds] = scoreIFs2
% for the large world network. 
% Outputs:
% ALL_TIMESTAMPS is a cell(21,21) of all the timestamps recorded, in the same layout
% as all variables below.
% meansC is a cell from (0:.1:2,0:.1:2) (E,I) injected currents. Holds all the
% mean vectors from procNets, which are numECells long.
% stC is a cell, same as above, but holds the standard deviation vectors.
% stMeans is the average standard deviation of all cells from each
% simulation. Also (0:.1:2,0:.1:2).
% stMeds is the median standard deviation; like stMeans.

    ALLTS = cell(21);
    meansC = cell(21);
    stC = cell(21);
    stMeans = zeros(21);
    stMeds = zeros(21);
    
    ind2inj = @(ind) (ind/10)-.1;
    for injE_ind = 1:21
        for injI_ind = 1:21
            injE = ind2inj(injE_ind);
            injI = ind2inj(injI_ind);
            ALLTS{injE_ind,injI_ind} = runIF2('tempScoreLW.mat',injE,injI);
            [meansC{injE_ind,injI_ind},stC{injE_ind,injI_ind}]=procNets2('tempScoreLW.mat');
            stMeans(injE_ind,injI_ind) = mean(stC{injE_ind,injI_ind});
            stMeds(injE_ind,injI_ind) = median(stC{injE_ind,injI_ind});
        end
    end
    
end

