clf 
xu=[0 .1 .11 .4 .79 .8 .88 .97 0.99 1];
 yu=[1.3 2.2 2.2  1.4 2.2  2.2 1.5 1 1 1];
 xv=[0 .2 .4 .6  .7 .8  .9 0.98 1];
 yv=[1  1  1  1  .9  .1  0  0   0];
 pu0=spline(xu,yu);
 pv0=spline(xv,yv);
 x=0:.01:1;
 %plot(x,ppval(pu0,x),'-',xu,yu,'o');
  plot(x,ppval(pv0,x),'-',xv,yv,'o');