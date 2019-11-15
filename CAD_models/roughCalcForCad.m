clear
clc
format compact

%% rough calc for cad
chord = 0.5;
AR = 8;
TR = 0.8;

TE_sweep = (180/pi)*atan(((1-TR)*chord)/(AR*chord/2));

LE_sweep = TE_sweep;

% fit a line to the chord as a function of position aling the span
fitChord = [chord TR*chord];
spanLoc = [0 AR*chord/2];

p = polyfit(spanLoc,fitChord,1);

% spanwise location of aerodynamic center
aeroSpanWise = TR*AR*chord/4;

% chordWise location of aerodynamic center
aeroChordWise = aeroSpanWise*tand(LE_sweep) + 0.25*polyval(p,aeroSpanWise);

% displaced volume
vol = 9500*(1e-3)^3;
rhoWater = 1000;
g = 9.81;

Fb = rhoWater*vol*g;

% weight force
volNew = 6824*(1e-3)^3;
rhoResin = 1000*1.3;

Fg = volNew*rhoResin*g;
percBuoy = Fb/Fg

% polar plots
load('ayazmat.mat')
CL_tot = clLW + clRW + clHS;
CL_wing = clLW + clRW;

wingRatio = CL_wing./CL_tot;
hsRatio = clHS./CL_tot;

% x
x = 20.4*CL_wing - 67.3*clHS;









