

template<typename T> inline std::enable_if<ad::is_double_precision<T>::value ||
	ad::is_single_precision<T>::value, T>::type atmosphere::MetricViscosity<T>::operator()(_In_ const T theta) {

		constexpr static  T TZERO(288.15);
		constexpr static  T BETAVISC(1.458E-6);
		constexpr static  T SUTH(110.4);
		T t = theta * TZERO;
		return (BETAVISC*std::sqrt(t*t*t) / (t + SUTH));
	}