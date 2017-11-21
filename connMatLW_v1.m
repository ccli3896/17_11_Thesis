function connMatLW_v1
    % Makes a matrix of cell# by cell# (600) as a reference for runIF.m,
    % telling which cells are dependent on which. 
    % Columns depend on rows. (3,2) means 3 feeds 2
    %
    % v1 uses the same probability at all levels (see journal 10.30)
    
    E_cells = 500;
    I_cells = 100;
    I_n = 1;

    cMat = zeros(E_cells+I_cells);
    
    TOTALPOS = 249500; % Total connections possible given the clusters, see journal 10.10
    TOTALLINKS = 9413; % Total links from SW network
    
    perc = TOTALLINKS/TOTALPOS;
    
    % Clusters at first level
    count = round(perc*25*20*19);
    rowr = randi(E_cells,[count,1]);
    ind = 1;
    while (ind<=count)
        colr = (floor(rowr(ind)/20)*20)+randi(20,[1,1]);
        if (cMat(rowr(ind),colr)==0 && rowr(ind)~=colr)
            cMat(rowr(ind),colr) = 1;
        else
            ind = ind-1;
        end
        ind = ind+1;
    end
    
    % Clusters at second level
    count = round(perc*5*100*80);
    rowr = randi(E_cells,[count,1]);
    clstr = @(n) floor(mod(n,100)/20);
    ind = 1;
    while (ind<=count)
        colr = (floor(rowr(ind)/100)*100)+randi(100,[1,1]);
        if (cMat(rowr(ind),colr)==0 && clstr(rowr(ind))~=clstr(colr) && rowr(ind)~=colr)
            cMat(rowr(ind),colr) = 1;
        else
            ind = ind-1;
        end
        ind = ind+1;
    end
    
    % Clusters at third level
    count = round(perc*500*400);
    rowr = randi(E_cells,[count,1]);
    clstr = @(n) floor(n/100);
    ind = 1;
    while (ind<=count)
        colr = randi(500,[1,1]);
        if(clstr(rowr(ind))~=clstr(colr) && cMat(rowr(ind),colr)==0 && rowr(ind)~=colr)
            cMat(rowr(ind),colr) = 1;
        else
            ind = ind-1;
        end
        ind = ind+1;
    end
   
    % Connecting all neighboring inhibitory cells 
    for ind = E_cells+1:E_cells+I_cells
        % If very low, close to the edge
        if (ind<=E_cells+I_n)
            cMat(ind,E_cells+1:ind) = 1;
            cMat(ind,I_cells-I_n+ind:E_cells+I_cells) = 1;
        else
            cMat(ind,ind-I_n:ind) = 1;
        end
        
        % If very high, close to the edge
        if (ind>E_cells+I_cells-I_n)
            cMat(ind,ind:E_cells+I_cells) = 1;
            cMat(ind,E_cells+1:E_cells+(E_cells+I_cells-ind)+1) = 1;
        else
            cMat(ind,ind:ind+I_n) = 1;
        end
    end
    
    % Connecting each inhibitory cell to five neighboring excitatory ones
    for ind = 1:I_cells
        cMat(((ind-1)*5)+1:(ind*5),ind+E_cells) = 1;
    end
    
    % Connecting each excitatory cell to ten random inhibitory ones 
    for ind = 1:E_cells
        IMat = E_cells+1:E_cells+I_cells;
        ITen = IMat(randperm(length(IMat),10));
        cMat(ITen,ind) = 1;
    end
    
    % Addition of inhibitory random connections
    % For a ratio of 1:1 local:long-range
    for i_ind = 1:I_cells
        for ind = 1:2
            tmpN = randi([E_cells+1 E_cells+I_cells],1,1);
            if (cMat(i_ind+E_cells,tmpN)==1)
                ind = ind-1;
                continue;
            end
            cMat(i_ind+E_cells,tmpN) = 1;
        end
    end
    
    % Erasing all diagonals
    for c=1:E_cells+I_cells
        cMat(c,c) = 0;
    end
    
    dlmwrite('connMatLW1.dat',cMat,'delimiter','\t');
end
    