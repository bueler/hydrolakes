function makestencil
% MAKESTENCIL   Plot a stencil.  Works under Matlab better than Octave, especially
% because saving the figures sucks under Octave.

msize=16;

basic(1) % stencil for non-constant diffusivity
hold on
plot([-0.5 0.5 0 0], [0 0 -0.5 0.5], '^k', 'Markersize',msize)
plot([-1 0 1 0 0], [0 0 0 -1 1], 'dk', 'Markersize',msize)
hold off
% to create .eps in Octave:
%print -depsc2 ../pdffigs/diffstencil.eps
%print -dpdf diffstencil.pdf

   function basic(i)
   figure(i)
   lwidth =4.0;

   % make the grid
   x=-1.3:0.1:1.3;
   plot(x,[-ones(size(x)); zeros(size(x)); ones(size(x))]','k',...
       'LineWidth',lwidth)
   axis([-1.7 1.7 -1.5 1.5])
   axis off
   hold on
   plot([-ones(size(x)); zeros(size(x)); ones(size(x))]',x,'k',...
       'LineWidth',lwidth)

   % make the cell
   xcell=-0.5:0.1:0.5;
   clwidth=2.0;
   plot(xcell,[-0.5*ones(size(xcell)); 0.5*ones(size(xcell))]','k--',...
        'LineWidth',clwidth)
   plot([-0.5*ones(size(xcell)); 0.5*ones(size(xcell))],xcell','k--',...
        'LineWidth',clwidth)

   % label the regular grid positions
   fsize = 16;
   text(-1,-1.5,'i-1','FontSize',fsize,'FontAngle','italic',...
      'HorizontalAlignment','center')
   text(0,-1.5,'i','FontSize',fsize,'FontAngle','italic',...
      'HorizontalAlignment','center')
   text(1,-1.5,'i+1','FontSize',fsize,'FontAngle','italic',...
      'HorizontalAlignment','center')
   text(-1.6,-1,'j-1','FontSize',fsize,'FontAngle','italic')
   text(-1.5,0,'j','FontSize',fsize,'FontAngle','italic')
   text(-1.62,1,'j+1','FontSize',fsize,'FontAngle','italic')
   hold off

