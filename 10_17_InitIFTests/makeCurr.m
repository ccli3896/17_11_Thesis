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