screen_size = get(0, 'ScreenSize');

buff=zeros(screen_size(3),screen_size(4));

for i=1:50
    buff(screen_size(3)/2-i,screen_size(4)/2+i)=100;

end
f1 = image(buff)
colormap(gray)

set(gcf,'windowstyle','modal');
set(gcf,'OuterPosition', screen_size); 
set(gcf,'position',screen_size); 

set(gcf,'Units','normal', 'outerposition',[0 0 1 1])
set(gca,'Visible', 'Off', 'Position',[0 0 1 1]) 

