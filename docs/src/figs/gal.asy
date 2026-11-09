// `gal.asy` Asymptote Figure Source https://asymptote.sourceforge.io/
// 2025 Tom C. Gorordo - `StandardHaloModelsPP.jl` https://github.com/tgorordo/StandardHaloModelsPP.jl
// adapted from https://asymptote.sourceforge.io/gallery/3Dgraphs/near_earth.asy
// Galactocentric Spherical Coordinates, the LSR tangent and solar peculiar velocity.

import three;
import math;
texpreamble("\usepackage{bm}");
texpreamble("\usepackage{starfont}");
size(300,0);

pen thickp=linewidth(0.5mm);
real d0=0.8, lambda=(90 - 10), aux=60;

currentprojection=perspective(4,1,1);

// Draw and label the Galactic Center
draw(scale3(0.036) * unitsphere, surfacepen = emissive(black));
label("Sgr. A$^\star$",(0,0,0),4*W);


// Draw and Label the Location of the Sun
triple rO=d0*dir(lambda,aux); // Sun Position
draw(shift(rO) * scale3(0.02) * unitsphere , surfacepen = emissive(orange));
label("\Sun",rO, NE);

triple rP=d0*dir(90,aux); // Proj. Sun Position in Plane
draw(Label("$R_0$", 0.55), O--rP, dashed, align = S);
label("$\theta=0$",d0*dir(90,aux), align=SE + E);

draw(Label("$\eta_0$", MidPoint), arc(O,0.3rO,0.3*d0*dir(90,aux)), E); // altitude of the Sun

// Draw Spherical octant
real rr=sqrt(rO.x^2+rO.y^2);
draw(arc((0,0,rO.z),(rr,0,rO.z),(0,rr,rO.z)),dashed);
draw(arc(O,d0*Z,d0*dir(90,aux)),dashed);
draw(arc(O,d0*Z,d0*X),thickp);
draw(arc(O,d0*Z,d0*Y),thickp);
draw(arc(O,d0*X,d0*Y),thickp);
draw(Label("N",1),d0*Z + 0.02*Z,red);

triple i=dir(90+lambda,aux);
triple k=unit(rO);
triple j=cross(k,i);

triple phihat = 0.2*i;
triple thetahat = 0.32*j;
triple rhat = 0.26*k;

triple xhat = -0.2*unit(rP);
triple yhat = -thetahat;
triple zhat = cross(xhat, yhat);

pen p = rgb(0.0, 0.7, 0.0);
triple rB = rO - rO.z * Z;
draw(Label("$\hat{x}$",1),rB--rB+xhat,2S,p,Arrow3);
draw(Label("$\hat{y}$",1),rB--rB+yhat,2N,p,Arrow3);
draw(Label("$\hat{z}$",1),rB--rB+zhat,2E,p,Arrow3);

draw(Label("$\hat{\phi}$",1),rO--rO+phihat,2W,red,Arrow3);
draw(Label("$\hat{\theta}$",1),rO--rO+thetahat,N,red,Arrow3);
draw(Label("$\hat{r}$",1),rO--rO+rhat,N,red,Arrow3);

triple vLSR = -0.2*j;
// (U, V, W) = (11.1, 12.2, 7.3) km/s along (x̂, ŷ, ẑ): toward Sgr. A*, with the rotation, north.
triple vPEC = 0.1*unit(11.1*unit(xhat) + 12.2*unit(yhat) + 7.3*unit(zhat));

draw(Label("$\vec{v}_{LSR}$",1),rO--rO+vLSR,blue+linewidth(0.6pt),Arrow3, align=NW);
draw(Label("$\vec{v}_{PEC}$",1),rO--rO+vPEC, blue+linewidth(0.6pt),Arrow3, align=NW);

//draw(Label("$\vec{v}_\odot=\vec{v}_{LSR} +\vec{v}_{PEC}$",1),rO--rO+vLSR+vPEC, blue+linewidth(0.6pt),Arrow3, align=2*NW + 2*N);

draw("$\bm{d}_0$",O--rO,PenMargin3, align=N);
draw("$\bm{z}_0$",rO--rO - rO.z * Z,PenMargin3, align=W);
