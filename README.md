# Radar-IR-EOS-Simulation

This project attempts to model and simulate the innerworking of Radar system, Electro-Optical active and passive sources of IR radiation
and missile guidance system.
The main purpose is to achieve a realistic system modeling and simulation as much as possible, hence the main sources of knowledge
are russian-language (soviet-era) engineering and technical literature reaching design level (candidate of science and experienced system
engineers and designers).

The second firm foundation which this project stand upon is being wholly optimized at basic level of massive manual vectorization by 
leveraging Intel Intrinsic programming i.e. usage of AVX/AVX2/AVX512 code path for almost every algorithm which is vectorizable.
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
Currently only hundreds (circa 400) kernels belonging to AVX512 double and single precision executing path
were implemented.
All of these kernels compute analytical Radar Cross Section of simple and to lesser extent complex objects.





List of references (incomplete):

Detection, Estimation, and Modulation Theory Part III: Radar-Sonar Signal Processing and Gaussian Signals in Noise Harry L. van Trees ISBN-10: 047110793X ISBN-13: 978-0471107934

Detection Estimation and Modulation Theory, Part I: Detection, Estimation, and Filtering Theory Harry L. van Trees ISBN-10: 9780470542965 ISBN-13: 978-0470542965 ASIN: 0470542969

Automatic Control of Aircraft and Missiles John H. Blakelock ASIN: B01FJ0JOU2

Principles of High-Resolution Radar (Artech House Radar Library), August W. Rihaczek ISBN-10: 089006900X
ISBN-13: 978-0890069004

Леонов А.И. - Моделирование в радиолокации (1979, Сов. радио ) 

Кремер И.Я. - Модулирующие помехи и прием радиосигналов (1972)

Abramowitz M., Stegun I.A. (eds.) - Handbook of mathematical functions (1972, NBS) 

Шифрин Я.С. - Вопросы статистической теории антенн 

Горяинов В.Т., Журавлев А.Г., Тихонов В.И - Статистическая радиотехника. Примеры и задачи (1980, Советское радио)

George T. Ruck, Donald E. Barrick , William D. Stuart , - Radar Cross Section Handbook (1970, Kluwer Academic Plenum Publishers)

Тихонов В. И. Статистическая радиотехника. «Сов. радио», 1966

Вайнштейн Л. А., Зубаков В. Д. Выделение сигналов на фоне случайных помех. «Сов. радио», 1960.

Зубковнч С. Г. Статистические характеристики радиосигналов, отраженных от земной поверхности. «Сов. радио», 1968

