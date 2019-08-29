clear
clc
format compact

%% attempt to recreate analytical wake redirection model from Bastankhah paper
d = 1;
yaw = 20*pi/180;
u = 1.5;

% locations
x0 = 0;
xs = 2;
xf = 16;
x = linspace(xs,xf,1);

% make grid about point
xturb = 10;
yturb = 0;
zturb = 0;

nTheta = 20;
nD = 20;

ths = linspace(0,2*pi,nTheta);
radii = linspace(0,d/2,nD);

yT = NaN(nD,nTheta);
zT = NaN(nD,nTheta);

for ii = 1:nTheta
    for jj = 1:nD
        yT(jj,ii) = yturb + radii(jj)*cos(ths(ii));
        zT(jj,ii) = zturb + radii(jj)*sin(ths(ii));
    end
end

y = yT(:)./d;
z = zT(:)./d;
zh = 0;

[X,Y,Z] = meshgrid(x,y,z);

upstreamTurbPos = [0;0;0];
downstreamTurbPos = [xturb;yturb;zturb];

[effectiveVel,du] = wakeDeflection(nD,nTheta,d,u,upstreamTurbPos,downstreamTurbPos,yaw);

%%
fg = figure(1);
for ii = 1:numel(x)
    subplot(2,4,ii)
    [C,h] = contourf(squeeze(Z(:,ii,:)),squeeze(Y(:,ii,:)),u*(1-squeeze(du(:,ii,:))));
    h.LineStyle = 'none';
    caxis([min(u*(1-squeeze(du(:,:,:))),[],'all') max(u*(1-squeeze(du(:,:,:))),[],'all')])
    axis equal
    hold on
    grid on
    colorbar
    colormap jet
    xlabel('y/d')
    ylabel('z/d')
    title(sprintf('x/d = %0.1f, pitch = %0.1f deg',x(ii),yaw*180/pi))
    
    
end

% saveas(fg,sprintf('pitchTurb%0.0f.png',yaw*180/pi));

