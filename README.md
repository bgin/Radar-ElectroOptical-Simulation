# Radar-IR-EOS-Simulation

This project attempts to model and simulate the innerworking of Radar system, Electro-Optical active and passive sources of IR radiation.
The main purpose is to achieve a realistic system modeling and simulation as much as possible, hence the main sources of knowledge
is a engineering and technical literature reaching design level.

The second firm foundation which this project stand upon is being wholly optimized at basic level of massive manual vectorization by 
leveraging Intel Intrinsic programming i.e. usage of SSE/AVX/AVX2/AVX512 code path for almost every algorithm which is vectorizable.
Compiler-level autovectorization is of secondary importance and is being inserted mainly to vectorize descriptive statistics routines
and profiling metrics calculations.

The second code path beside the SIMD  is the GPGPU Cuda implementation counting so far close to 15000 lines of code of computational
and helper routines and kernels.

I envision five main components:
1) Radar system modeling and simulation.
2) Radio altimeter modeling and simulation.
3) Propagation of laser and IR radiation through the turbulent atmospheric channels.
4) Optical signals processing (background noise extraction).
5) Electro-optical sensor modeling and simulation.
   
The main structure of the projects is a collection of free standing 'modules' programmatically describing
various modelled components.
It is a software library of framework and may be used as computational backend of larger program of be
connected to GUI front-end.
Currently only a few hundreds (circa 1000) kernels belonging to SSE/AVX/AVX512 double and single precision executing path
were implemented.
All of these kernels compute analytical Radar Cross Section of simple and to lesser extent complex objects.
Beside aforestated, the Antennae modeling computational kernels were developed which describe various radiation
characteristics of different antennae.






