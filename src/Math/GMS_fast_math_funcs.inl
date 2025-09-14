
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T>  &MatFunctions<T>::fastcot(const T arg)
{

	 constexpr  T PI = static_cast<T>(3.14159265358979323846264);
	 constexpr  T NEG_PI = -PI;
	 T& sum = 0.0;
	 if (arg > PI)
	 {
		 return std::numeric_limits<T>::quiet_NaN();
	 }
	 else if (-arg < NEG_PI)
	 {
		 return std::numeric_limits<T>::quiet_NaN();
	 }
	 else
	 {
		 constexpr T coeff1 = static_cast<T>(-0.333333333333333333333333);
		 constexpr T coeff2 = static_cast<T>(-0.022222222222222222222222);
		 constexpr T coeff3 = static_cast<T>(-0.002116402116402116402116);
		 constexpr T coeff4 = static_cast<T>(-0.000211640211640211640211);
		 constexpr T coeff5 = static_cast<T>(-0.000021377799155576933354);
		 constexpr T coeff6 = static_cast<T>(-2.16440428080639720851361e-6);
		 constexpr T coeff7 = static_cast<T>(-2.19259478518737777997037e-7);
		 constexpr T coeff8 = static_cast<T>(-2.22146087899796790760584e-8);
		 constexpr T coeff9 = static_cast<T>(-2.25078465168089928542084E-9);
		 constexpr T coeff10 = static_cast<T>(-2.28051512045921828658638e-10);
		 constexpr T coeff11 = static_cast<T>(-2.31064325990026240965325e-11);
		 constexpr T coeff12 = static_cast<T>(-2.34117068198248839592094e-12);
		 T invx = 1.0 / arg;
		 T sqrx = arg * arg;
		 sum = invx + arg * (coeff1 + sqrx * (coeff2 + sqrx * (coeff3 + sqrx * (coef4 + sqrx * (coeff5 + sqrx * (coeff6 + sqrx * (coeff7 + sqrx * (coeff8 + sqrx * (coeff9 + sqrx * (coeff10 + sqrx * (coeff11 + sqrx * (coeff12 + sqrx))))))))))));
		 return sum;
	 }
}

//------------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------------//

template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T> &MatFunctions<T>::fastsin(const T arg)
{

	constexpr   T PI          = static_cast<T>(3.14159265358979323846264);
	constexpr   T NEG_PI      = -PI;
	constexpr   T HALF_PI     = 0.5 * PI;
	constexpr   T NEG_HALF_PI = -HALF_PI;

	if (arg > HALF_PI)
	{
		return ::sin(arg);
	}
	else if (-arg < NEG_HALF_PI)
	{
		return ::sin(-arg);
	}
	else
	{

		constexpr T coeff1 = static_cast<T>(-0.16666666666666666);
		constexpr T coeff2 = static_cast<T>(0.00833333333333333);
		constexpr T coeff3 = static_cast<T>(-1.9841269841269841);
		constexpr T coeff4 = static_cast<T>(2.75573192239858906);
		constexpr T coeff5 = static_cast<T>(-2.5052108385441718);
		constexpr T coeff6 = static_cast<T>(1.60590438368216145);
		constexpr T coeff7 = static_cast<T>(-7.6471637318198164);
		constexpr T coeff8 = static_cast<T>(2.81145725434552076);
		constexpr T coeff9 = static_cast<T>(-8.2206352466243297);
		constexpr T coeff10 =static_cast<T>( 1.95729410633912612);
		constexpr T coeff11 = static_cast<T>(-3.8681701706306840);
		T arg_squared = arg * arg;
		T& sine = arg + arg*arg_squared*(coeff1 + arg_squared*(coeff2 + arg_squared*(coeff3 + arg_squared*(coeff4 + arg_squared
			*(coeff5 + arg_squared*(coeff6 + arg_squared*(coeff7 + arg_squared*(coeff8 + arg_squared*(coeff9 + arg_squared
			*(coeff10 + arg_squared*(coeff11 + arg_squared)))))))))));
		return sine;
	}
}
//------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------//
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T> &MatFunctions<T>::fastcos(const T arg)
{

	constexpr   T PI = static_cast<T>(3.14159265358979323846264);
	constexpr   T NEG_PI = -PI;
	constexpr   T HALF_PI = 0.5 * PI;
	constexpr   T NEG_HALF_PI = -HALF_PI;
	if (std::fabs(arg ) > HALF_PI)
	{
		return ::cos(arg);
	}
	else if (std::fabs(arg) < 0.0)
	{
		return ::cos(arg);
	}
	else
	{
		constexpr T coeff1 = static_cast<T>(-0.50000000000000000);
		constexpr T coeff2 = static_cast<T>(0.041666666666666666);
		constexpr T coeff3 = static_cast<T>(-0.00138888888888888);
		constexpr T coeff4 = static_cast<T>(0.000024801587301587);
		constexpr T coeff5 = static_cast<T>(-2.75573192239858906e-7);
		constexpr T coeff6 = static_cast<T>(2.087675698786809897e-9);
		constexpr T coeff7 = static_cast<T>(-1.14707455977297247e-11);
		constexpr T coeff8 = static_cast<T>(4.77947733238738529e-14);
		constexpr T coeff9 = static_cast<T>(-1.56192069685862264e-16);
		constexpr T coeff10 =static_cast<T>( 4.11031762331216485e-19);
		constexpr T coeff11 =static_cast<T>( -8.8967913924505732e-22);
		constexpr T one = 1.0e+00;
		T squared_arg = arg * arg;
		T& sum = one + squared_arg*(coeff1 + squared_arg*(coeff2 + squared_arg*(coeff3 + squared_arg*(coeff4 + squared_arg
			*(coeff5 + squred_arg*(coeff6 + squared_arg*(coeff6 + squred_arg*(coeff7 + squred_arg
			*(coeff8 + squred_arg*(coeff9 + squred_arg*(coeff9 + squred_arg*(coeff10 + squred_arg*(coeff11 + squred_arg)))))))))))));
		return sum;
	}
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T>  &MatFunctions<T>::fasttan(const T arg)
{
	constexpr   T PI      = static_cast<T>(3.14159265358979323846264);
	constexpr   T HALF_PI = 0.5 * PI;
	if (fabs(arg) > HALF_PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}else 
	if (arg < 0.0)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else
	{
		constexpr T coeff1 = static_cast<T>(0.33333333333333333);
		constexpr T coeff2 = static_cast<T>(0.13333333333333333);
		constexpr T coeff3 = static_cast<T>(0.05396825396825396);
		constexpr T coeff4 = static_cast<T>(0.02186948853615520);
		constexpr T coeff5 = static_cast<T>(0.00886323552990219);
		constexpr T coeff6 = static_cast<T>(0.00359212803657248);
		constexpr T coeff7 = static_cast<T>(0.00145583438705131);
		constexpr T coeff8 = static_cast<T>(0.00059002744094558);
		constexpr T coeff9 = static_cast<T>(0.00023912911424355);
		constexpr T coeff10 = static_cast<T>(0.0000969153795692);
		constexpr T coeff11 = static_cast<T>(0.00003927832388331);
		T squared_arg = arg * arg;
		T& sum = arg + arg*squared_arg*(coeff1 + squred_arg*(coeff2 + squred_arg*(coeff3 + squred_arg*(coeff4 + squared_arg*(coeff5 + squared_arg
			*(coeff6 + squared_arg*(coeff7 + squared_arg*(coeff8 + squared_arg*(coeff9 + squared_arg*(coeff10 + squared_arg*(coeff11 + squared_arg)))))))))));
		return sum;
	}
	
	

}
//-----------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------//
template<typename T> inline static  std::enable_if<std::is_floating_point<T>::value, T> &MatFunctions<T>::fastexp(const T arg)
{
	if (fabs(arg) < 3.0)
		return std::exp(arg);

	
	constexpr T coeff1 = static_cast<T>(1.0e+00);
	constexpr T coeff2 = static_cast<T>(0.5);
	constexpr T coeff3 = static_cast<T>(0.166666666666666666);
	constexpr T coeff4 = static_cast<T>(0.041666666666666666);
	constexpr T coeff5 = static_cast<T>(0.008333333333333333);
	constexpr T coeff6 = static_cast<T>(0.001388888888888888);
	constexpr T coeff7 = static_cast<T>(0.001388888888888888);
	constexpr T coeff8 = static_cast<T>(0.000198412698412698);
	constexpr T coeff9 = static_cast<T>(0.000024801587301587);
	constexpr T coeff10 = static_cast<T>(2.755731922398589065e-6);
	constexpr T coeff11 = static_cast<T>(2.755731922398589065e-7);
	return  1.0 + arg*coeff1 + arg*(coeff2 + arg*(coeff3 + arg*(coeff4 + arg*(coeff5 + arg*(coeff6 +
		arg*(coeff7 + arg*(coeff8 + arg*(coeff9 + arg*(coeff10 + arg*(coeff11 + arg))))))))));
	
}
//-----------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------//
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T> &MatFunctions<T>::fastcsc(const T arg)
{
	constexpr  T PI = static_cast<T>(3.14159265358979323846264);
	constexpr  T NEG_PI = -PI;
	if (arg > PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else if (-arg < NEG_PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else
	{
		constexpr T coeff1 = static_cast<T>(0.166666666666666);
		constexpr T coeff2 = static_cast<T>(0.019444444444444);
		constexpr T coeff3 = static_cast<T>(0.002050264550264);
		constexpr T coeff4 = static_cast<T>(0.000209986772486);
		constexpr T coeff5 = static_cast<T>(0.000021336045641);
		constexpr T coeff6 = static_cast<T>(2.163347442778659e-6);
		constexpr T coeff7 = static_cast<T>(2.192327134456764e-7);
		constexpr T coeff8 = static_cast<T>(2.2213930853920414e-8);
		constexpr T coeff9 = static_cast<T>(2.2507674795567867e-9);
		constexpr T coeff10 = static_cast<T>(2.280510770721821e-10);
		constexpr T coeff11 = static_cast<T>(2.310642158099696e-11);
		constexpr T coeff12 = static_cast<T>(2.341170402893194e-12);
		T& inv_arg = 1.0 / arg;
		T arg_squared = arg * arg;
		return inv_arg + arg*(coeff1 + arg_squared*(coeff2 + arg_squared*(coeff3 + arg_squared*(coeff4 + arg_squared*(coeff5 + arg_squared
			*(coeff6 + arg_squared*(coeff7 + arg_squared*(coeff8 + arg_squared*(coeff9 + arg_squared*(coeff10 + arg_squared
			*(coeff11 + arg_squared*(coeff12 + arg_squared))))))))))));
	}
		
}
//--------------------------------------------------------------------------//

//--------------------------------------------------------------------------//
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T> &MatFunctions<T>::fastsec(const T arg)
{
	constexpr   T PI = static_cast<T>(3.14159265358979323846264);
	constexpr   T NEG_PI = -PI;
	constexpr   T HALF_PI = 0.5 * PI;
	constexpr   T NEG_HALF_PI = -HALF_PI;
	if (arg > HALF_PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else if (arg < NEG_HALF_PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else
	{
		constexpr T coeff1 = static_cast<T>(0.500000000000000);
		constexpr T coeff2 = static_cast<T>(0.208333333333333);
		constexpr T coeff3 = static_cast<T>(0.084722222222222);
		constexpr T coeff4 = static_cast<T>(0.034350198412698);
		constexpr T coeff5 = static_cast<T>(0.013922233245149);
		constexpr T coeff6 = static_cast<T>(0.005642496810031);
		constexpr T coeff7 = static_cast<T>(0.002286819095164);
		constexpr T coeff8 = static_cast<T>(0.000926812927377);
		constexpr T coeff9 = static_cast<T>(0.000375623133852);
		constexpr T coeff10 = static_cast<T>(0.00015223432221);
		constexpr T coeff11 = static_cast<T>(0.00006169824687);
		constexpr T coeff12 = static_cast<T>(0.00002500535760);
		T& squared_arg = arg * arg;
		return 1.0 + squared_arg*(coeff1 + squared_arg*(coeff2 + squared_arg*(coeff3 + squared_arg*(coeff4 + squared_arg
			*(coeff5 + squared_arg*(coeff6 + squared_arg*(coeff7 + squared_arg*(coeff8 + squared_arg*(coeff9 + squared_arg
			*(coeff10 + squared_arg*(coeff11 + squared_arg*(coeff12 + squared_arg))))))))))));

	}
}
//-------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------//
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T> &MatFunctions<T>::fastarcsin(const T arg)
{

	if (arg > 1.0)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else if (-arg < -1.0)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else
	{
		constexpr T coeff1 = static_cast<T>(0.16666666666666666);
		constexpr T coeff2 = static_cast<T>(0.07500000000000000);
		constexpr T coeff3 = static_cast<T>(0.04464285714285714);
		constexpr T coeff4 = static_cast<T>(0.03038194444444444);
		constexpr T coeff5 = static_cast<T>(0.02237215909090909);
		constexpr T coeff6 = static_cast<T>(0.01735276442307692);
		constexpr T coeff7 = static_cast<T>(0.01396484375000000);
		constexpr T coeff8 = static_cast<T>(0.01155180089613970);
		constexpr T coeff9 = static_cast<T>(0.00976160952919407);
		constexpr T coeff10 =static_cast<T>(0.0083903358096168);
		constexpr T coeff11 =static_cast<T>(0.0073125258735988);
		constexpr T coeff12 =static_cast<T>(0.0064472103118896);
		T& squared_arg;
		return arg + arg*squared_arg*(coeff1 + squared_arg*(coeff2 + squared_arg*(coeff3 + squared_arg*(coeff4 + squared_arg*(coeff5 + squared_arg
			*(coeff6 + squared_arg*(coeff7 + squared_arg*(coeff8 + squared_arg*(coeff9 + squared_arg*(coeff10 + squared_arg*(coeff11 + squared_arg
			*(coeff12 + squared_arg))))))))))));
	}
}
//-----------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------//
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T > &MatFunctions<T>::fastarctan(const T arg)
{
	constexpr T one = static_cast<T>(1.0e+00);
	constexpr T neg_one = (-one);
	if (arg > one)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else if (arg < neg_one)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else
	{
		constexpr T coeff1 = static_cast<T>(-0.33333333333333);
		constexpr T coeff2 = static_cast<T>(0.200000000000000);
		constexpr T coeff3 = static_cast<T>(-0.14285714285714);
		constexpr T coeff4 = static_cast<T>(0.111111111111111);
		constexpr T coeff5 = static_cast<T>(-0.09090909090909);
		constexpr T coeff6 = static_cast<T>(0.07692307692307);
		constexpr T coeff7 = static_cast<T>(-0.06666666666666);
		constexpr T coeff8 = static_cast<T>(0.058823529411764);
		constexpr T coeff9 = static_cast<T>(-0.05263157894736);
		constexpr T coeff10 =static_cast<T>(0.04761904761904);
		constexpr T coeff11 =static_cast<T>(-0.0434782608695);
		constexpr T coeff12 =static_cast<T>(0.04000000000000);
		T& squared_arg = arg * arg;
		return arg + arg*squared_arg*(coeff1 + squared_arg*(coeff2 + squared_arg*(coeff3 + squared_arg*(coeff4 + squared_arg*(coeff5 + squared_arg
			*(coeff6 + squared_arg*(coeff7 + squared_arg*(coeff8 + squared_arg*(coeff9 + squared_arg*(coeff10 + squared_arg*(coeff11 + squared_arg
			*(coeff12 + squared_arg))))))))))));
	}
}
//----------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------//
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T> &MatFunctions<T>::fastsinh(const T arg)
{
	constexpr T coeff1 = static_cast<T>(0.166666666666666);
	constexpr T coeff2 = static_cast<T>(0.008333333333333);
	constexpr T coeff3 = static_cast<T>(0.000198412698412);
	constexpr T coeff4 = static_cast<T>(2.755731922398589E-6);
	constexpr T coeff5 = static_cast<T>(2.505210838544171E-8);
	constexpr T coeff6 = static_cast<T>(1.605904383682161E-10);
	constexpr T coeff7 = static_cast<T>(7.647163731819816E-13);
	constexpr T coeff8 = static_cast<T>(2.811457254345520E-15);
	constexpr T coeff9 = static_cast<T>(8.220635246624329E-18);
	constexpr T coeff10 =static_cast<T>( 1.95729410633912E-20);
	constexpr T coeff11 =static_cast<T>( 3.86817017063068E-23);
	constexpr T coeff12 =static_cast<T>( 6.44695028438447E-26);
	T& squared_arg = arg * arg;
	return arg + arg*squared_arg*(coeff1 + squared_arg*(coeff2 + squared_arg*(coeff3 + squared_arg*(coeff4 + squared_arg*(coeff5 + squared_arg
		*(coeff6 + squared_arg*(coeff7 + squared_arg*(coeff8 + squared_arg*(coeff9 + squared_arg*(coeff10 + squared_arg*(coeff11 + squared_arg
		*(coeff12 + squared_arg))))))))))));

}
//--------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------//
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T> &MatFunctions<T>::fastcosh(const T arg)
{
	constexpr T coeff1 = static_cast<T>(0.500000000000000);
	constexpr T coeff2 = static_cast<T>(0.041666666666666);
	constexpr T coeff3 = static_cast<T>(0.001388888888888);
	constexpr T coeff4 = static_cast<T>(0.000024801587301);
	constexpr T coeff5 = static_cast<T>(2.755731922398589E-7);
	constexpr T coeff6 = static_cast<T>(2.087675698786809E-9);
	constexpr T coeff7 = static_cast<T>(1.147074559772972E-11);
	constexpr T coeff8 = static_cast<T>(4.779477332387385E-14);
	constexpr T coeff9 = static_cast<T>(1.561920696858622E-16);
	constexpr T coeff10 =static_cast<T>(4.11031762331216E-19);
	constexpr T coeff11 =static_cast<T>(8.89679139245057E-22);
	constexpr T coeff12 =static_cast<T>(1.61173757109611E-24);
	T& squared_arg = arg * arg;
	return 1.0 + squared_arg*(coeff1 + squared_arg*(coeff2 + squared_arg*(coeff3 + squared_arg*(coeff4 + squared_arg*(coeff5 + squared_arg
		*(coeff6 + squared_arg*(coeff7 + squared_arg*(coeff8 + squared_arg*(coeff9 + squared_arg*(coeff10 + squared_arg*(coeff11 + squared_arg
		*(coeff12 + squared_arg))))))))))));
}
//------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------//
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T> &MatFunctions<T>::fasttanh(const T arg)
{
	constexpr   T PI = static_cast<T>(3.14159265358979323846264);
	constexpr   T NEG_PI = -PI;
	constexpr   T HALF_PI = 0.5 * PI;
	constexpr   T NEG_HALF_PI = -HALF_PI;
	if (arg > HALF_PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else if (-x < NEG_HALF_PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else
	{
		constexpr T coeff1 = static_cast<T>(-0.333333333333333);
		constexpr T coeff2 = static_cast<T>(0.133333333333333);
		constexpr T coeff3 = static_cast<T>(-0.05396825396825);
		constexpr T coeff4 = static_cast<T>(0.021869488536155);
		constexpr T coeff5 = static_cast<T>(-0.00886323552990);
		constexpr T coeff6 = static_cast<T>(0.003592128036572);
		constexpr T coeff7 = static_cast<T>(-0.00145583438705);
		constexpr T coeff8 = static_cast<T>(0.000590027440945);
		constexpr T coeff9 = static_cast<T>(-0.00023912911424);
		constexpr T coeff10 =static_cast<T>( 0.00009691537956);
		constexpr T coeff11 =static_cast<T>(-0.0000392783238);
		constexpr T coeff12 =static_cast<T>( 0.00001591890506);
		T& squared_arg = arg * arg;
		return arg + arg*squared_arg*(coeff1 + squared_arg*(coeff2 + squared_arg*(coeff3 + squared_arg*(coeff4 + squared_arg*(coeff5 + squared_arg
			*(coeff6 + squared_arg*(coeff7 + squared_arg*(coeff8 + squared_arg*(coeff9 + squared_arg*(coeff10 + squared_arg*(coeff11 + squared_arg
			*(coeff12 + squared_arg))))))))))));
	}
}
//-----------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------//
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T> &MatFunctions<T>::fastcsch(const T arg)
{
	constexpr   T PI = static_cast<T>(3.14159265358979323846264);
	constexpr   T NEG_PI = -PI;
	constexpr   T HALF_PI = 0.5 * PI;
	constexpr   T NEG_HALF_PI = -HALF_PI;
	if (arg > HALF_PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else if (-arg < NEG_HALF_PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else
	{

		constexpr T coeff1 = static_cast<T>(-0.1666666666666);
		constexpr T coeff2 = static_cast<T>(0.01944444444444);
		constexpr T coeff3 = static_cast<T>(-0.0020502645502);
		constexpr T coeff4 = static_cast<T>(0.0002099867724);
		constexpr T coeff5 = static_cast<T>(-0.0000213360456);
		constexpr T coeff6 = static_cast<T>(2.16334744277865E-6);
		constexpr T coeff7 = static_cast<T>(-2.192327134456E-7);
		constexpr T coeff8 = static_cast<T>(2.221393085392E-8);
		constexpr T coeff9 = static_cast<T>(-2.2507674795567E-9);
		constexpr T coeff10 =static_cast<T>( 2.2805107707218E-10);
		constexpr T coeff11 =static_cast<T>( -2.31064215809E-11);
		constexpr T coeff12 =static_cast<T>( 2.341170402893E-12);
		T inv = 1.0 / arg;
		T& squared_arg = arg * arg;
		return inv + arg*(coeff1 + squared_arg*(coeff2 + squared_arg*(coeff3 + squared_arg*(coeff4 + squared_arg*(coeff5 + squared_arg
			*(coeff6 + squared_arg*(coeff7 + squared_arg*(coeff8 + squared_arg*(coeff9 + squared_arg*(coeff10 + squared_arg
			*(coeff11 + squared_arg*(coeff12 + squared_arg))))))))))));
	}
}
//-------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------//
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T> &MatFunctions<T>::fastsech(const T arg)
{
	constexpr   T PI = static_cast<T>(3.14159265358979323846264);
	constexpr   T NEG_PI = -PI;
	constexpr   T HALF_PI = 0.5 * PI;
	constexpr   T NEG_HALF_PI = -HALF_PI;

	if (arg > HALF_PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else  if (-arg < NEG_HALF_PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else
	{
		constexpr T coeff1 = static_cast<T>(-0.5000000000000);
		constexpr T coeff2 = static_cast<T>(0.2083333333333);
		constexpr T coeff3 = static_cast<T>(-0.084722222222);
		constexpr T coeff4 = static_cast<T>(0.0343501984126);
		constexpr T coeff5 = static_cast<T>(-0.013922233245);
		constexpr T coeff6 = static_cast<T>(0.0056424968100);
		constexpr T coeff7 = static_cast<T>(-0.002286819095);
		constexpr T coeff8 = static_cast<T>(0.0009268129273);
		constexpr T coeff9 = static_cast<T>(-0.000375623133);
		constexpr T coeff10 =static_cast<T>(0.000152234322);
		constexpr T coeff11 =static_cast<T>(-0.00006169824);
		constexpr T coeff12 =static_cast<T>(0.000025005357);
		T& squared_arg = arg * arg;
		return 1.0 + squared_arg*(coeff1 + squared_arg*(coeff2 + squared_arg*(coeff3 + squard_arg*(coeff4 + squared_arg*(coeff5 + squared_arg
			*(coeff6 + squared_arg*(coeff7 + squared_arg*(coeff8 + squared_arg*(coeff9 + squared_arg*(coeff10 + vsquared_arg*(coeff11 + squared_arg
			*(coeff12 + squared_arg))))))))))));
	}
}
//------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------//
template<typename T> inline static std::enable_if<std::is_floating_point<T>::value, T> &MatFunctions<T>::fastcoth(const T arg)
{
	constexpr   T PI = static_cast<T>(3.14159265358979323846264);
	constexpr   T NEG_PI = -PI;
	if (x > PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else if (-x < NEG_PI)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	else
	{
		constexpr T coeff1 = static_cast<T>(0.33333333333333);
		constexpr T coeff2 = static_cast<T>(-0.0222222222222);
		constexpr T coeff3 = static_cast<T>(0.0021164021164);
		constexpr T coeff4 = static_cast<T>(-0.000211640211);
		constexpr T coeff5 = static_cast<T>(0.0000213777991);
		constexpr T coeff6 = static_cast<T>(-2.1644042808063E-6);
		constexpr T coeff7 = static_cast<T>(2.1925947851873E-7);
		constexpr T coeff8 = static_cast<T>(-2.2214608789979E-8);
		constexpr T coeff9 = static_cast<T>(2.2507846516808E-9);
		constexpr T coeff10 =static_cast<T>(-2.28051512045E-10);
		T inv = 1.0 / arg;
		T& squared_arg = arg * arg;
		return inv + arg*(coeff1 + squared_arg*(coeff2 + squared_arg*(coeff3 + squared_arg*(coeff4 + squared_arg*(coeff5 + squared_arg
			*(coeff6 + squared_arg*(coeff7 + squared_arg*(coeff8 + squared_arg*(coeff9 + squared_arg*(coeff10 + squared_arg))))))))));
	}
}
