function testSolver_v2
% Runs a three-cell IF simulation.
% v2 uses the analytical solution to the differential eqn for voltage.
% There are no current vectors assigned to each cell; each one is recreated
% as needed. Instead, this program stores firing times.

    tiFr = 100; % Timeframe in ms
    dt = .1; % in ms
    %spSp = .001; % Spontaneous spiking probability at any timestep
    refrTi = 3; % in ms NOT SPECIFIED IN PAPER
    injE = .8; % injected current to alter excitability
    injI = .8;
    capac = 1; % capacitance. Compare to NEURON simulations, maybe.
    numC = 3;
    numE = 1;
    
    alphj = 1.15;
    
    AVs = cell(numC,1); % holds v
    APts = cell(numC,1); % holds timestamps for action potentials.
    cMat = [0 1 1; 1 0 1; 1 1 0];
    
    % Making all the zero vectors so space doesn't have to be dynamically
    % allocated
    for ind1 = 1:numC
        AVs{ind1} = zeros(1,(tiFr/dt)+1);
    end
    
    refr = zeros(1,numC); 
    
    % FIRING TIME MATRIX: FOR CHECKING. Random numbers, basically.
    % See journal 10.3
    fireTimes = [5 25 27 28 29 80; 10 50 81 -1 -1 -1; 11 55 -1 -1 -1 -1];
    
    for t = 2:tiFr/dt
        fprintf('%f\n',t*dt);
    % Next section is to update the voltage vectors.
        for c = 1:numC
            
        % If AP is in progress.
            if (refr(c) > 0)
                refr(c) = refr(c)-dt;
                AVs{c,1}(t) = 0;
            end
        % If AP should be fired.
            if (AVs{c,1}(t-1) >= 1)
                refr(c) = refrTi;
                AVs{c,1}(t) = 0;
                apTime = t*dt;
                APts{c} = [APts{c} apTime];
            end
        % For firing at times specified by fireTimes matrix
            if (refr(c)<=0)
                if (sum(fireTimes(c,:)==dt*t))
                    refr(c) = refrTi;
                    AVs{c,1}(t) = 0;
                    apTime = t*dt;
                    APts{c} = [APts{c} apTime]; 
                end
            end
                
        % Normal summing over all cells.
            if (refr(c)<=0)            
                % Adjusting for I vs E cell
                if (c<=numE)
                    AVs{c,1}(t) = cVolt(AVs{c,1}(t-1),t,injE,alphj,cMat,dt,capac,numC,numE,APts,c);
                else
                    AVs{c,1}(t) = cVolt(AVs{c,1}(t-1),t,injI,alphj,cMat,dt,capac,numC,numE,APts,c);
                end
            end
            
        end
    end
    
    %Plot all voltages from all cells
    hold on;
    for c_ind = 1:numC
        plot(AVs{c_ind,1})
    end
end

function newV = cVolt(prevV,t,inj,alphj,cMat,dt,capac,numC,numE,APts,c)
% Calculates the voltage for the next time step. A single number.
    % Inputs:
    % prevV is the t-1 voltage for that cell.
    % t is time, to figure out tk vectors.
    % inj is the injected current for i/e.
    % alphj is pretty much the conductance. Units are 1/R, at least.
    % cMat is the connectivity matrix.
    % dt and capac are defined and explained in the main function.
    % numC and numE are for assigning the correct weights.
    % APts is the cell of AP timestamps. 
    % c is the cell number
    
% This version uses the analytical solution to the differential eqn for
% voltage. 

TAUF = .3; % From Zochowski paper
TAUS = 3;

% First makes the tk vector for each cell based on what timestep the
% simulation is on. Saves two AP times, to fix current jumps.
    tk = zeros(numC,2)-1000;
    for c_i = 1:numC
        if (~isempty(APts{c_i})) % if there's an APts vector for this cell
            if (size(APts{c_i})>1)
                tk(c_i,1) = APts{c_i}(end-1) - (dt*(t-1));
            else
                tk(c_i,1) = APts{c_i}(end) - (dt*(t-1));
            end
            tk(c_i,2) = APts{c_i}(end) - (dt*(t-1));
        end
    end

% Calculates the summation term (see journal 10.16)
    synSum = 0;
    for c2 = 1:numC
        % to skip the same cell
        if (c == c2)
            continue;
        end
        % Check for the projecting cells and set wjk
        if (cMat(c2,c)==1)
            if (c2<=numE && c<=numE) % both E
                wjk = 2;
            elseif (c2<=numE && c>numE) % E to I
                wjk = 4;
            elseif (c2>numE && c<=numE) % I to E
                wjk = -2;
            else % I to I;
                wjk = -10; 
            end
            sumTerm = @(tkTest,TAUS,alphj,capac,dt,TAUF) ((TAUS/(alphj*TAUS-capac))*(exp(tkTest/TAUS))*...
                (exp(-dt/TAUS)-exp(-alphj*dt/capac)) - ...
                (TAUF/(alphj*TAUF-capac))*(exp(tkTest/TAUF))*...
                (exp(-dt/TAUF)-exp(-alphj*dt/capac)));
            synSum = synSum + wjk*max([sumTerm(tk(c2,1),TAUS,alphj,capac,dt,TAUF); sumTerm(tk(c2,2),TAUS,alphj,capac,dt,TAUF)]);
        end
    end

    newV = prevV*exp(-alphj*dt/capac) + inj/alphj - (inj/alphj)*exp(-alphj*dt/capac) + synSum;
end