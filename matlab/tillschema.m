function tillschema
% TILLSCHEMA  Produce a figure showing how the till water amount W_tl evolves.

set(0,'defaultlinelinewidth',2.0)
da = 0.1;
rht = 1.6;

y = 0.55;
plot([-0.2 rht], [y y], 'k'), hold on,  arrow(rht,y)
tick(0.0,y,'0'),  tick(0.6,y,'\tau W'),  tick(1.0,y,'W_{til}^{max}')
plot([0.2 0.5], [y+da y+da], 'k'),  sharpright(0.5,y+da)
text(0.05,y+1.7*da,'W_{til} increases','fontsize',16)
plot([0.77 0.93], [y+da y+da], 'k'),  sharpleft(0.77,y+da)
text(0.7,y+1.7*da,'W_{til} decreases','fontsize',16)

y = 0.0;
plot([-0.2 rht], [y y], 'k'),  arrow(rht,y)
tick(0.0,y,'0'),  tick(1.0,y,'W_{til}^{max}'),  tick(1.4,y,'\tau W')
plot([0.25 0.75], [y+da y+da], 'k'),  sharpright(0.75,y+da)
text(0.2,y+1.7*da,'W_{til} increases','fontsize',16)

hold off
axis([-0.3 (rht+0.2) -0.2 1.3]), axis off
end

  function tick(x,y,s)
    plot([x x], [y-0.03 y+0.03],'k')
    text(x-0.04,y-0.08,s,'fontsize',16)
  end

  function arrow(x,y)
    plot([x-0.02 x], [y+0.02 y],'k')
    plot([x-0.02 x], [y-0.02 y],'k')
  end
 
  function sharpleft(x,y)
    plot([x+0.03 x], [y+0.015 y],'k')
    plot([x+0.03 x], [y-0.015 y],'k')
  end

  function sharpright(x,y)
    plot([x-0.03 x], [y+0.015 y],'k')
    plot([x-0.03 x], [y-0.015 y],'k')
  end
