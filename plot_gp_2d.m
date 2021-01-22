function [fig]=plot_gp_2d(X, z, ptitle, saveplot)

if(nargin >3)
    savep = saveplot;
else
    savep =false;
end
x= X(1, :);
y = X(2, :);
xv = linspace(min(x), max(x), 3000);
yv = linspace(min(y), max(y), 3000);
[xx,yy] = meshgrid(xv, yv);
zz = griddata(x, y, z, xx,yy);
fig= figure(); clf
surf(xx, yy, zz);
xlabel('Latitude', 'FontSize', 18)
ylabel('Longitude', 'FontSize', 18)
zlabel('Temprature', 'FontSize', 18)
grid on
set(gca, 'ZLim',[0 50])
shading interp
view(2)
colorbar;

title(ptitle, 'FontSize', 20)  

if(savep)
fig.PaperPosition = [ 0 0 12 10];
set(fig,'PaperSize',[12 10])
print(fig, ptitle,'-dpdf')
end

end