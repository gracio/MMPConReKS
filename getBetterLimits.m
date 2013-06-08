function[betterMin, betterMax, betterInc] = getBetterLimits(trueMin,trueMax,numTicks,absRange)

trueInc = (trueMax-trueMin)/(numTicks - 2);
betterInc = round(trueInc/5*(10^ceil(-1*log10(trueInc/5))))/(10^ceil(-1*log10(trueInc/5)))*5;
betterMin = max(absRange(1),floor(trueMin/5*(10^ceil(-1*log10(trueInc/5))))/(10^ceil(-1*log10(trueInc/5)))*5 - betterInc);
betterMax = min(absRange(2),betterMin + numTicks*betterInc);

if trueMax >  betterMax
    betterMax = betterMax + betterInc;
    betterMin = betterMin + betterInc; 
end