function [F_tension,magTension]  = fcn(R1_o,V1_o,Rn_o,Vn_o,Ri_o,Vi_o,thrLength,F_tension,magTension,...
    numNode, massBdy, zetaTether, youngsTether, diaTether)

if numNode == 2
    Ri_o = [R1_o Rn_o];
    Vi_o = [V1_o Vn_o];
else
    Ri_o = [R1_o Ri_o Rn_o];
    Vi_o = [V1_o Vi_o Vn_o];
end

% cross sectional area
CS = (pi/4)*diaTether^2;

% total spring stiffness
k_tot = youngsTether*CS/thrLength;

% total damping coefficient
c_tot = zetaTether*2*sqrt(k_tot*massBdy);

% inidividual spring stiffness and damping
kElem = k_tot*(numNode-1);
cElem = c_tot*(numNode-1);
LiElem = thrLength/(numNode-1);

% determine distance between nodes and the derivatives of each element
di = zeros(numNode-1,1);
direction = zeros(3,numNode-1);
di_dot = zeros(numNode-1,1);
Tk = zeros(numNode-1,1);
Td = zeros(numNode-1,1);
F_t_i = zeros(3,numNode-1);


for ii = 1:numNode-1
    di(ii) = norm(Ri_o(:,ii+1) - Ri_o(:,ii));
    direction(:,ii) = (Ri_o(:,ii+1) - Ri_o(:,ii))/di(ii);
    di_dot(ii) = dot(Ri_o(:,ii+1) - Ri_o(:,ii),Vi_o(:,ii+1) - Vi_o(:,ii))/di(ii);
    
end

% spring and damper forces
for ii = 1:numNode-1
    if di(ii) > LiElem
        Tk(ii) = kElem*(di(ii) - LiElem);
        Td(ii) = cElem*di_dot(ii);
        F_t_i(:,ii) = (Tk(ii) + Td(ii))*direction(:,ii);
        magTension(ii) = norm(F_t_i(:,ii));
    else
        Tk(ii) = 0;
        Td(ii) = 0;
        F_t_i(:,ii) = (Tk(ii) + Td(ii))*direction(:,ii);
        magTension(ii) = norm(F_t_i(:,ii));
        
    end
end

for ii = 1:numNode
    if ii == 1
        F_tension(:,ii) = F_t_i(:,ii);
    elseif ii > 1 && ii < numNode
        F_tension(:,ii) = F_t_i(:,ii) - F_t_i(:,ii-1);
    elseif ii == numNode
        F_tension(:,ii) = -F_t_i(:,ii-1);
    end
end


end







