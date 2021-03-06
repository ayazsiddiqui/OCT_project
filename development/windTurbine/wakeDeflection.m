function [effectiveVel,du]  = wakeDeflection(nD,nTheta,diameter,upstreamTurbFlowVel,...
    upstreamTurbPos,downstreamTurbPos,upstreamTurbYaw,ky,kz,CT,effectiveVel)

yaw = upstreamTurbYaw;
xturb = downstreamTurbPos(1);
yturb = downstreamTurbPos(2);
zturb = downstreamTurbPos(3);

x0 = upstreamTurbPos(1);
zh = upstreamTurbPos(3);

d = diameter;
ths = linspace(0,2*pi,nTheta);
radii = linspace(0,d/2,nD);

yT = NaN(nD,nTheta);
zT = NaN(nD,nTheta);

% create grid
for ii = 1:nTheta
    for jj = 1:nD
        yT(jj,ii) = yturb + radii(jj)*cos(ths(ii));
        zT(jj,ii) = zturb + radii(jj)*sin(ths(ii));
    end
end

x = xturb;
y = yT(:);
z = zT(:);

% wake growth rates
TI = 1;
alpStar = 2.32;
betaStar = 0.154;

% wake widths
sigmaYbyD = ky*(x-x0)/d + (cos(yaw)/sqrt(8));
sigmaZbyD = kz*(x-x0)/d + (1/sqrt(8));

% normalized length of potential core
x0byD = (cos(yaw)*(1 + sqrt(1-CT)))/(sqrt(2)*(alpStar*TI + betaStar*(1-sqrt(1-CT))));

% wake skew angle
theta = ((0.3*yaw)/(cos(yaw)))*(1 - sqrt(1 - CT*cos(yaw)));

% normalized wake deflection\
t1 = theta*x0byD;
t2 = (theta/14.7)*sqrt(cos(yaw)/(ky*kz*CT));
t3 = 2.9 + 1.3*sqrt(1 - CT) - CT;
t4 = (1.6 + sqrt(CT))*(1.6*sqrt((8*sigmaYbyD.*sigmaZbyD)/cos(yaw)) - sqrt(CT));
t5 = (1.6 - sqrt(CT))*(1.6*sqrt((8*sigmaYbyD.*sigmaZbyD)/cos(yaw)) + sqrt(CT));

deltabyD = t1 + t2*t3*log(t4./t5);

% % % far wake velocity
t1n = NaN(length(y),length(z));
t2n = NaN(length(y),length(z));
t3n = NaN(length(y),length(z));
du = NaN(length(y),length(z));

for jj = 1:length(z)
    for ii = 1:length(y)
        
        t1n(ii,jj) = 1 - sqrt(1 - ((CT*cos(yaw))./(8*(sigmaYbyD.*sigmaZbyD))));
        t2n(ii,jj) = exp(-0.5*(((y(ii)/d)-deltabyD)./sigmaYbyD).^2);
        t3n(ii,jj) = exp(-0.5*(((z(jj) - zh)/d)./sigmaZbyD).^2);
        
        du(ii,jj) = t1n(ii,jj).*t2n(ii,jj).*t3n(ii,jj);
    end
    
end

effectiveVel = mean(upstreamTurbFlowVel*(1-du(:)));

end

