# ODEPACK

The double precision ODEPACK solvers, unmodified, from
<https://www.netlib.org/odepack/> (`opkdmain.f`, `opkda1.f`, `opkda2.f`).
The reference engine calls `DLSODA`, which switches between Adams and BDF
methods as the problem turns stiff. It is the LSODA of A. C. Hindmarsh and
L. R. Petzold (Lawrence Livermore National Laboratory), the solver ADAPT
and S-ADAPT use.

ODEPACK is in the public domain.
