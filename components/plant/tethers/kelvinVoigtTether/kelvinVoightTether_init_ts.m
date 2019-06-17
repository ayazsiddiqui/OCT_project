x = 1;

for ii = 1:length(tN)   
    temp_cls(ii).initNodePos = [...
        linspace(iniR1_o(1,ii),iniRn_o(1,ii),tN(ii));...
        linspace(iniR1_o(2,ii),iniRn_o(2,ii),tN(ii));...
        linspace(iniR1_o(3,ii),iniRn_o(3,ii),tN(ii))];
    
    temp_cls(ii).initNodePos = temp_cls(ii).initNodePos(:,2:end-1);
    temp_cls(ii).initNodeVel = zeros(size(temp_cls(ii).initNodePos));
    
end

x = 2;

