function vectorPlotter(time,data,field_names,yAxisLabel)

lwd = 1;

colors = 1/255*[228,26,28
    55,126,184
    77,175,74
    152,78,163
    255,127,0
    255,255,51];

sdata = squeeze(data);
sz = size(sdata);

for ii = 1:sz(1)
    plot(time,sdata(ii,:),'-','linewidth',lwd,'color',colors(ii,:))
    hold on
    if ii == 1
        grid on
        hold on
    end
end
xlabel('Time (s)');
ylabel(yAxisLabel);
legend(field_names)
end

