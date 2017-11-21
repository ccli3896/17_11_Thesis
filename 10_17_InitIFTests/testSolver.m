function testSolver
% Runs a three-cell IF simulation.

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
    
    AVs = cell(numC,2); % 1 holds v, 2 holds syn current
    APts = cell(numC,1); % holds timestamps for action potentials.
    cMat = [0 1 1; 1 0 1; 1 1 0];
    
    % Making all the zero vectors so space doesn't have to be dynamically
    % allocated
    for ind1 = 1:numC
        for ind2 = 1:2
            AVs{ind1,ind2} = zeros(1,(tiFr/dt)+1);
        end
    end
    
    refr = zeros(1,numC); 
    apFlag = 0;
    
    % FIRING TIME MATRIX: FOR CHECKING. Random numbers, basically.
    % See journal 10.3
    fireTimes = [5 25 27 28 29 80;10 50 81 -1 -1 -1;11 55 -1 -1 -1 -1];
    
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
                apFlag = 1;
            end
        % For firing at times specified by fireTimes matrix
            if (refr(c)<=0)
                if (sum(fireTimes(c,:)==dt*t))
                    refr(c) = refrTi;
                    AVs{c,1}(t) = 0;
                    apTime = t*dt;
                    APts{c} = [APts{c} apTime]; 
                    apFlag = 1;
                end
            end
            
        % Makes current if an AP fired
            if (apFlag)
                % 10.23 changed from AVs{c,2}+makeCurr(apTime,tiFr,0,dt);
                AVs{c,2} = max([AVs{c,2};makeCurr(apTime,tiFr,0,dt)]);
                apFlag = 0;
            end
                
        % Normal summing over all cells.
            if (refr(c)<=0)
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
                    synSum = synSum + wjk*AVs{c2,2}(t);
                end
            end
            
            % Adjusting for I vs E cell
            if (c<=numE)
                AVs{c,1}(t) = cVolt(AVs{c,1}(t-1),injE,alphj,synSum,dt,capac);
            else
                AVs{c,1}(t) = cVolt(AVs{c,1}(t-1),injI,alphj,synSum,dt,capac);
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

function newV = cVolt(prevV,inj,alphj,synSum,dt,capac)
% Calculates the voltage for the next time step. A single number.
    % Inputs:
    % prevV is the t-1 voltage for that cell.
    % inj is the injected current for i/e.
    % alphj is pretty much the conductance. Units are 1/R, at least.
    % synSum is the already-calculated summation term.
    % dt and capac are defined and explained in the main function.
% This version uses an ODE solver to check against the old version.
%
    tspan = [0,2*dt];
    v0 = prevV;
    [t,v] = ode45(@(t,v) (1/capac)*(-alphj*v+inj+synSum), tspan, v0);
    newV = interp1(t,v,dt);
%{
   coef = prevV-(inj/alphj)-((1/alphj)*synSum);
   newV = (inj/alphj)+(1/alphj)*synSum+coef*exp((-1*dt)/(capac/alphj));
%}
end

function synVec = makeCurr(APtime,tlength,tstart,dt)
    % APtime is the time a spike fires
    % tlength is the timeframe of the whole simulation
    % tstart is 0: the beginning of the entire current vector
    % dt is dt
    %
    % synVec is a vector with the current trace from Jablonski et al. 2007
    
    tvec = linspace(tstart,tlength,((tlength-tstart)/dt)+1);
    synVec = zeros(1,length(tvec));
    
    synVec = exp((-1*(tvec-APtime))/3) - exp((-1*(tvec-APtime))/.3);
    synVec(tvec<APtime) = 0;
end