1C:
Files still are named hw1b so a new cmake is not needed. Just in case the list is still included.

Demo Screencapture Link: 
https://raw.githubusercontent.com/NotoriousENG/ogl_tutorial/master/ogl_doubleViewInteraction.gif?token=AHWX43FMUQNUYMS47C6MEW26LR5CQ

*****************************************************************************
1B:
Demo Screencapture Link: 
https://raw.githubusercontent.com/NotoriousENG/ogl_tutorial/master/ogl_interactableCurves.gif?token=AHWX43DB5XSZZH2Y27DZSVC6J5L4Y

Modulo operator % is used as we want 
to get only possible control points

We will go to point 0 when exceeding the bounds (positively)
We will go to the final point when exceeding the bounds (negatively)

These points are computed as the average of the interior points 
of their neighbors

c[i][0] = (c[i][1] + c[(i-1)%n][2] ) /2
c[i][3] = (c[i][2] + c[(i+1)%n][1] ) /2