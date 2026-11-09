// `eph.asy` Asymptote Figure Source https://asymptote.sourceforge.io/
// 2025 Tom C. Gorordo - `StandardHaloModelsPP.jl` https://github.com/tgorordo/StandardHaloModelsPP.jl
// The Earth's orbit about the Sun, in Galactic-Rectangular coordinates.
//
// Oriented with the package's ecliptic basis (ϵ1, ϵ2 at the epoch of teq; McCabe 2014,
// Lee, Lisanti & Safdi 2013). Their normal ϵ1 × ϵ2 ≈ (-0.10, 0.86, 0.50) points to galactic
// (l, b) ≈ (96°, 30°), so the ecliptic is inclined ≈ 60° to the galactic plane. The Earth is
// drawn at the March equinox teq and at T0 (≈ 72° of orbit later), where its velocity is most
// closely aligned with the Sun's. Velocities share one scale, |v_ES| / |v_⊙| ≈ 30 / 251,
// independent of the 1 au orbit.

import three;
size(320, 0);
currentprojection = perspective(6, -3.5, 3, up = Z);

triple e1 = unit((0.9941, 0.1086, 0.0047));   // Earth's velocity direction at teq
triple e2 = unit((-0.0500, 0.4947, -0.8677)); // Earth -> Sun direction at teq
triple n = unit(cross(e1, e2));               // ecliptic north

real ph0 = 2pi * 73.0 / 365.25;               // orbital phase of T0 after teq
triple rEq = -e2;                             // Sun -> Earth at teq
triple rT0 = -cos(ph0) * e2 + sin(ph0) * e1;  // Sun -> Earth at T0
triple vdir(real ph) { return cos(ph) * e1 + sin(ph) * e2; }

real vscale = 1.4;                            // |v_⊙| as drawn
triple vsun = vscale * unit((11.1, 250.2, 7.3));
real ves = vscale * 29.8 / 250.6;

pen gp = rgb(0.0, 0.6, 0.0);
pen vp = blue + linewidth(0.6pt);

// Galactic plane (reference) and the ecliptic, both at 1 au
draw(circle(O, 1, Z), gray(0.65) + dashed);
label("galactic plane", (0.0, -1.0, 0), gray(0.45), align = SW);
draw(circle(O, 1, n), black);
label("ecliptic (1 au)", e2, align = 3E + S);

// Direction of the Earth's motion along the orbit
draw(arc(O, -cos(ph0 + 2.2) * e2 + sin(ph0 + 2.2) * e1,
         -cos(ph0 + 2.6) * e2 + sin(ph0 + 2.6) * e1, n), black, Arrow3);

// Galactic-Rectangular axes at the Sun
draw(Label("$\hat{x}$", 1, align = NE), O--0.55X, gp, Arrow3);
draw(Label("$\hat{y}$", 1), O--0.55Y, gp, Arrow3);
draw(Label("$\hat{z}$", 1), O--0.55Z, gp, Arrow3);
label("to Sgr. A$^\star$", 0.35X, gp, align = 3S + 0.6E);

// The Sun and its galactic velocity
draw(scale3(0.045) * unitsphere, surfacepen = emissive(orange));
label("$\odot$", O, align = 2NW);
draw(Label("$\vec{v}_\odot$", 1), O--vsun, vp, Arrow3);

// The Earth at the March equinox, with the ecliptic basis
draw(shift(rEq) * scale3(0.03) * unitsphere, surfacepen = emissive(blue));
label("$\oplus\ (t_\mathrm{eq})$", rEq, align = 2N);
draw(Label("$\hat\epsilon_1$", 1), rEq--rEq + 0.4 * e1, red, Arrow3);
draw(Label("$\hat\epsilon_2$", 1), rEq--rEq + 0.4 * e2, red, Arrow3);

// The Earth at T0, where its velocity is most closely aligned with the Sun's
draw(shift(rT0) * scale3(0.03) * unitsphere, surfacepen = emissive(blue));
label("$\oplus\ (T_0)$", rT0, align = 2E);
draw(Label("$\vec{v}_{ES}$", 1, align = E), rT0--rT0 + ves * vdir(ph0), vp, Arrow3);
