
template<class Coords>  inline gms::math::Cartesian3D<Coords>::Cartesian3D(_In_ const Coords &X_Axis, _In_ const Coords &Y_Axis,
	_In_ const Coords &Z_Axis, _In_ const std::string &Type, const bool isFixed) :
	mX_Axis{ X_Axis },
	mY_Axis{ Y_Axis },
	mZ_Axis{ Z_Axis },
	m_Type{ Type },
	m_isFixed{ isFixed },
	mX_Origin{ mX_Axis.getPosition()[0] },
	mY_Origin{ mY_Axis.getPosition()[1] },
	mZ_Origin{ mZ_Axis.getPosition()[2] }
{
	this->m_OrientationQuaternions = std::array<gms::math::Quaternion, 3U>{};
	this->m_OrientationQuaternions.operator[](0) = mX_Axis.getOrientation();
	this->m_OrientationQuaternions.operator[](1) = mY_Axis.getOrientation();
	this->m_OrientationQuaternions.operator[](2) = mZ_Axis.getOrientation();
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//

template<class Coords> template<class OtherCoords> inline gms::math::Cartesian3D<Coords>::Cartesian3D(_In_ const OtherCoords &X,
	_In_ const OtherCoords &Y, _In_ const OtherCoords &Z, _In_ const std::string &Type, const bool isFixed) :
	mX_Axis{ X },
	mY_Axis{ Y },
	mZ_Axis{ Z },
	m_Type{ Type },
	m_isFixed{ isFixed },
	mX_Origin{ X.getPosition()[0] },
	mY_Origin{ Y.getPosition()[1] },
	mZ_Origin{ Z.getPosition()[2] }

{
	this->m_OrientationQuaternions = std::array<gms::math::Quaternion, 3U>{};
	this->m_OrientationQuaternions.operator[](0) = mX_Axis.getOrientation();
	this->m_OrientationQuaternions.operator[](1) = mY_Axis.getOrientation();
	this->m_OrientationQuaternions.operator[](2) = mZ_Axis.getOrientation();
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords>  inline gms::math::Cartesian3D<Coords>::Cartesian3D(_In_ const Cartesian3D &rhs) :
mX_Axis{ rhs.mX_Axis },
mY_Axis{ rhs.mY_Axis },
mZ_Axis{ rhs.mZ_Axis },
m_Type{ rhs.m_Type },
m_isFixed{ rhs.m_isFixed },
mX_Origin{ rhs.mX_Origin },
mY_Origin{ rhs.mY_Origin },
mZ_Origin{rhs.mZ_Origin}
{

	this->m_OrientationQuaternions.operator[](0).operator=( rhs.m_OrientationQuaternions[0]);
	this->m_OrientationQuaternions.operator[](1).operator=(rhs.m_OrientationQuaternions[1]);
	this->m_OrientationQuaternions.operator[](2).operator=(rhs.m_OrientationQuaternions[2]);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords>  inline gms::math::Cartesian3D<Coords>::Cartesian3D(_In_ Cartesian3D &&rhs) :
mX_Axis{ std::move(rhs.mX_Axis) },
mY_Axis{ std::move(rhs.mY_Axis) },
mZ_Axis{ std::move(rhs.mZ_Axis) },
m_Type{ std::move(rhs.m_Type) },
m_isFixed{ std::move(rhs.m_isFixed) },
mX_Origin{ std::move(rhs.mX_Origin) },
mY_Origin{ std::move(rhs.mY_Origin) },
mZ_Origin{ std::move(rhs.mZ_Origin) }
{

	this->m_OrientationQuaternions.operator[](0).operator=(std::move(rhs.m_OrientationQuaternions[0]));
	this->m_OrientationQuaternions.operator[](1).operator=(std::move(rhs.m_OrientationQuaternions[1]));
	this->m_OrientationQuaternions.operator[](2).operator=(std::move(rhs.m_OrientationQuaternions[2]));
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::operator=(_In_ const Cartesian3D &rhs)->gms::math::Cartesian3D<Coords> & {

	if (this == &rhs) return (*this);

	this->mX_Axis.operator=(rhs.mX_Axis);
	this->mY_Axis.operator=(rhs.mY_Axis);
	this->mZ_Axis.operator=(rhs.mZ_Axis);
	this->m_Type.operator=(rhs.m_Type);
	this->m_isFixed = rhs.m_isFixed;
	this->mX_Origin = rhs.mX_Origin;
	this->mY_Origin = rhs.mY_Origin;
	this->mZ_Origin = rhs.mZ_Origin;
	this->m_OrientationQuaternions[0].operator=( rhs.m_OrientationQuaternions[0]);
	this->m_OrientationQuaternions[1].operator=( rhs.m_OrientationQuaternions[1]);
	this->m_OrientationQuaternions[2].operator=(rhs.m_OrientationQuaternions[2]);
	return (*this);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
 template<class Coords> template<class OtherCoords>  inline auto gms::math::Cartesian3D<Coords>
	::operator=(_In_ const OtherCoords &rhs)->gms::math::Cartesian3D<Coords> & {
		if (this == &rhs) return (*this);
		/*
		rhs must have X(),Y() and Z() accessors.
		rhs must have X().getPosition().operator[](), Y().getPosition().operator[](),Z().getPosition().operator[]() accessors.
		*/
		this->mX_Axis = rhs.X();
		this->mY_Axis = rhs.Y();
		this->mZ_Axis = rhs.Z();
		this->m_Type.operator=(rhs.m_Type);
		this->m_isFixed = rhs.isFixed;
		this->mX_Origin = rhs.X().getPosition()[0];
		this->mY_Origin = rhs.Y().getPosition()[1];
		this->mZ_Origin = rhs.Z().getPosition()[2];
		this->m_OrientationQuaternions[0].operator=(rhs.m_OrientationQuaternions[0]);
		this->m_OrientationQuaternions[1].operator=(rhs.m_OrientationQuaternions[1]);
		this->m_OrientationQuaternions[2].operator=(rhs.m_OrientationQuaternions[2]);
		return (*this);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::operator=(_In_ Cartesian3D &&rhs)->gms::math::Cartesian3D<Coords> & {

	if (this == &rhs) return (*this);
	
	this->mX_Axis.operator=(std::move(rhs.mX_Axis));
	this->mY_Axis.operator=(std::move(rhs.mY_Axis));
	this->mZ_Axis.operator=(std::move(rhs.mZ_Axis));
	this->m_Type.operator=(std::move(rhs.m_Type));
	this->m_isFixed = std::move(rhs.m_isFixed);
	this->mX_Origin = std::move(rhs.mX_Origin);
	this->mY_Origin = std::move(rhs.mY_Origin);
	this->mZ_Origin = std::move(rhs.mZ_Origin);
	this->m_OrientationQuaternions[0].operator=(std::move(rhs.m_OrientationQuaternions[0]));
	this->m_OrientationQuaternions[1].operator=(std::move(rhs.m_OrientationQuaternions[1]));
	this->m_OrientationQuaternions[2].operator=(std::move(rhs.m_OrientationQuaternions[2]));
	return (*this);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//

template<class Coords> inline auto gms::math::Cartesian3D<Coords>::operator==(_In_ Cartesian3D &rhs)->std::array<__m256d, 6U> {
	constexpr static const int N{ 6 };
	std::array<__m256d, N> ret_val = std::array<__m256d, N>();
	for (int i{ 0 }; i != N / 2; ++i) ret_val[i] = this->m_OrientationQuaternions[i].operator==(rhs.getOrientQuaternions()[i]);

	ret_val[3] = this->mX_Axis.getPosition().operator==(rhs.getX_Axis().getPosition());
	ret_val[4] = this->mY_Axis.getPosition().operator==(rhs.getY_Axis().getPosition());
	ret_val[5] = this->mZ_Axis.getPosition().operator==(rhs.getZ_Axis().getPosition());
	return (ret_val);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::operator!=(_In_ Cartesian3D &rhs)->std::array<__m256d, 6U> {
	constexpr static const int N{ 6 };
	std::array<__m256d, N> ret_val = std::array<__m256d, N>();
	for (int i{ 0 }; i != N / 2; ++i) ret_val[i] = this->m_OrientationQuaternions[i].operator!=(rhs.getOrientQuaternions()[i]);

	ret_val[3] = this->mX_Axis.getPosition().operator!=(rhs.getX_Axis().getPosition());
	ret_val[4] = this->mY_Axis.getPosition().operator!=(rhs.getY_Axis().getPosition());
	ret_val[5] = this->mZ_Axis.getPosition().operator!=(rhs.getZ_Axis().getPosition());
	return (ret_val);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::getX_Axis() const->Coords {

	return (this->mX_Axis);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::getY_Axis() const->Coords {

	return (this->mY_Axis);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::getZ_Axis() const->Coords {

	return (this->mZ_Axis);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::getType() const->std::string {

	return (this->m_Type);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::getOrientQuaternions() const->std::array<gms::math::Quaternion, 3U> {

	return (this->m_OrientationQuaternions);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::isFixed() const->bool {

	return (this->m_isFixed);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::getXOrigin() const->double {

	return (this->mX_Origin);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::getYOrigin() const->double {

	return (this->mY_Origin);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::getZOrigin() const->double {

	return (this->mZ_Origin);
}
//-----------------------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------------------------------//
template<class Coords> inline auto gms::math::Cartesian3D<Coords>::displayCartesian3D() const->void {

	std::cout << std::setprecision(16) << std::fixed << std::showpoint
		<< "Cartesian3D object context \n\n"
		<< "Position vectors and Direction quaternions: \n";
	this->mX_Axis.displayAxis();
	this->mY_Axis.displayAxis();
	this->mZ_Axis.displayAxis();
	std::cout << "Name:" << this->m_Type.c_str() << std::endl;
    std::cout << "Is fixed:" << this->m_isFixed << std::boolalpha;
}

template<class Coords> inline auto gms::math::Cartesian3D<Coords>::cross_prod()->gms::math::Cartesian3D<Coords> & {

	//Coords axis = Coords{ this->getZ_Axis().getPosition().cross(this->getX_Axis().getPosition(), this->getY_Axis().getPosition()), this->getX_Axis().getOrientation() };
	//axis.displayAxis();
	this->getZ_Axis().getPosition() = this->getZ_Axis().getPosition().cross(this->getX_Axis().getPosition(), this->getY_Axis().getPosition());
	return (*this);
}
		
		







