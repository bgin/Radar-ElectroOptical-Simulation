
#include <iostream>
#include "GMS_arma_model.h"
#include "GMS_timsac_iface.h"
#include "GMS_fpexcept.h"
#if defined _WIN64
    #include "../GMS_common.h"
#elif defined __linux
    #include "GMS_common.h"
#endif

#if !defined (GMS_ARMA_MODEL_FP_EXCEPT_PROLOG)
#define GMS_ARMA_MODEL_FP_EXCEPT_PROLOG(status)             \
	(status) = clear_fpexcepts();                           \
    if (0 != (status))  {                                   \
	  std::cerr << "clear_fpexcpets failed to clear FP environment!!!  Detection of floating-point exceptions will not be possible!!!\n"; \
	}
#endif

#if !defined (GMS_ARMA_MODEL_FP_EXCEPT_EPILOG)
#define GMS_ARMA_MODEL_FP_EXCEPT_EPILOG(status,fpexcepts)  \
	do {		                                                    \
									                                \
		if (0 != (status)) {                                        \
		  (fpexcepts[0]) = test_feinvalid(FE_INVALID);		        \
		  (fpexcepts[1]) = test_feinexact(FE_INEXACT);              \
		  (fpexcepts[2]) = test_fedivbyzero(FE_DIVBYZERO);          \
		  (fpexcepts[3]) = test_feunormal(FE_DENORMAL);             \
		  (fpexcepts[4]) = test_feoverflow(FE_OVERFLOW);            \
		  (fpexcepts[5]) = test_feunderflow(FE_UNDERFLOW);          \
		}											                \
  } while (0); 
#endif

#if !defined (ARMA_MODEL_COPY_OP_BODY)
#define ARMA_MODEL_COPY_OP_BODY      \
	delete m_iq;				     \
	delete m_b2;				     \
	delete m_ip;				     \
	delete m_a2;				     \
	delete m_std;				     \
	delete m_cxx2;				     \
	delete m_g;					     \
	delete m_saic;	                 \
	m_newn = x.m_newn;			     \
	m_nmax = x.m_nmax;			     \
	m_mmax = x.m_mmax;			     \
	m_kq = x.m_kq;				     \
        m_kp = x.m_kp;				     \
        m_aicm = x.m_aicm;			     \
        m_iq = new VAi32(*(x.m_iq));     \
        m_b2 = new VAf64(*(x.m_b2));     \
        m_ip = new VAi32(*(x.m_ip));     \
        m_a2 = new VAf64(*(x.m_a2));     \
        m_std = new VAf64(*(x.m_std));   \
        m_cxx2 = new VAf64(*(x.m_cxx2)); \
        m_g = new VAf64(*(x.m_g));       \
        m_saic = new VAf64(*(m_saic));
#endif

#if !defined (ARMA_MODEL_MOVE_CTOR_BODY)
#define ARMA_MODEL_MOVE_CTOR_BODY        \
	m_nmax = x.m_nmax;					\
        m_mmax = x.m_mmax;					\
	m_iq = x.m_iq;						\
	m_b2 = x.m_b2;						\
	m_ip = x.m_ip;						\
	m_a2 = x.m_a2;						\
	m_std = x.m_std;				    \
	m_cxx2 = x.m_cxx2;				    \
	m_g = x.m_g;						\
	m_saic = x.m_saic;					\
	x.m_nmax = 0;						\
	x.m_mmax = 0;						\
	x.m_iq   = nullptr;				    \
	x.m_b2   = nullptr;					\
	x.m_ip  = nullptr;				    \
	x.m_a2  = nullptr;					\
	x.m_std = nullptr;					\
	x.m_cxx2 = nullptr;					\
	x.m_g    = nullptr;                 \
	x.m_saic = nullptr;
#endif

#if !defined (ARMA_MODEL_MOVE_OP_BODY)
#define ARMA_MODEL_MOVE_OP_BODY         \
	delete m_iq;				        \
	delete m_b2;				        \
	delete m_ip;				        \
	delete m_a2;				        \
	delete m_std;				        \
	delete m_cxx2;				        \
	delete m_g;					        \
	delete m_saic;	                    \
	m_newn = x.m_newn;			        \
	m_nmax = x.m_nmax;			        \
	m_mmax = x.m_mmax;			        \
	m_kq = x.m_kq;				        \
	m_kp = x.m_kp;				        \
	m_aicm = x.m_aicm;                  \
	m_iq = x.m_iq;						\
	m_b2 = x.m_b2;						\
	m_ip = x.m_ip;						\
	m_a2 = x.m_a2;						\
	m_std = x.m_std;				    \
	m_cxx2 = x.m_cxx2;				    \
	m_g = x.m_g;						\
	m_saic = x.m_saic;					\
	x.m_nmax = 0;						\
	x.m_mmax = 0;						\
	x.m_iq = nullptr;				    \
	x.m_b2 = nullptr;					\
	x.m_ip = nullptr;				    \
	x.m_a2 = nullptr;					\
	x.m_std = nullptr;					\
	x.m_cxx2 = nullptr;					\
	x.m_g = nullptr;                    \
	x.m_saic = nullptr;
#endif


gms::math::stat
::ArmaModel
::ArmaModel()
:
m_newn{},
m_nmax{},
m_mmax{},
m_kq{},
m_kp{},
m_aicm{},
m_iq(  nullptr),
m_b2(  nullptr),
m_ip(  nullptr),
m_a2(  nullptr),
m_std( nullptr),
m_cxx2(nullptr),
m_g(   nullptr),
m_saic(nullptr)
{}



gms::math::stat
::ArmaModel
::ArmaModel(const int32_t nmax,
	    const int32_t mmax,
	    const int32_t ival,
	    const double  dval)
:
m_newn{},
m_nmax{ nmax },
m_mmax{ mmax },
m_kq{},
m_kp{},
m_aicm{},
m_iq(    new VAi32(ival,static_cast<size_t>(m_nmax))),
m_b2(    new VAf64(dval,static_cast<size_t>(m_mmax*m_nmax))),
m_ip(    new VAi32(ival,static_cast<size_t>(m_nmax))),
m_a2(    new VAf64(dval,static_cast<size_t>(m_mmax*m_nmax))),
m_std(   new VAf64(dval,static_cast<size_t>(m_mmax*m_nmax))),
m_cxx2(  new VAf64(dval,static_cast<size_t>(m_nmax))),
m_g(     new VAf64(dval,static_cast<size_t>(m_mmax*m_nmax))),
m_saic(  new VAf64(dval,static_cast<size_t>(m_nmax))) {}

gms::math::stat
::ArmaModel
::ArmaModel(const ArmaModel &x)
:
m_newn{ x.m_newn },
m_nmax{ x.m_nmax },
m_mmax{ x.m_mmax },
m_kq{   x.m_kq },
m_kp{   x.m_kp },
m_aicm{ x.m_aicm},
m_iq(   new VAi32(*(x.m_iq))),
m_b2(   new VAf64(*(x.m_b2))),
m_ip(   new VAi32(*(x.m_ip))),
m_a2(   new VAf64(*(x.m_a2))),
m_std(  new VAf64(*(x.m_std))),
m_cxx2( new VAf64(*(x.m_cxx2))),
m_g(    new VAf64(*(x.m_g))),
m_saic( new VAf64(*(x.m_saic)))  {}



gms::math::stat
::ArmaModel
::ArmaModel(ArmaModel &&x)
:
m_newn{ x.m_newn },
m_nmax{ 0 },
m_mmax{ 0 },
m_kq{   x.m_kq },
m_kp{   x.m_kp },
m_aicm{ x.m_aicm},
m_iq(   nullptr),
m_b2(   nullptr),
m_ip(   nullptr),
m_a2(   nullptr),
m_std(  nullptr),
m_cxx2( nullptr),
m_g(    nullptr),
m_saic( nullptr)  {
	
	ARMA_MODEL_MOVE_OP_BODY
}

gms::math::stat
::ArmaModel
::~ArmaModel() {
	if (nullptr != m_iq)   delete m_iq;     m_iq   = nullptr;
	if (nullptr != m_b2)   delete m_b2;     m_b2   = nullptr;
	if (nullptr != m_ip)   delete m_ip;     m_ip   = nullptr;
	if (nullptr != m_a2)   delete m_a2;     m_a2   = nullptr;
	if (nullptr != m_std)  delete m_std;    m_std  = nullptr;
	if (nullptr != m_cxx2) delete m_cxx2;   m_cxx2 = nullptr;
	if (nullptr != m_g)    delete m_g;      m_g = nullptr;
	if (nullptr != m_saic) delete m_saic;   m_saic = nullptr;
	
 }

gms::math::stat::ArmaModel &
gms::math::stat::ArmaModel
::operator=(const ArmaModel &x) {

	if (this == &x) return (*this);
	ARMA_MODEL_COPY_OP_BODY
	  
	 return (*this);
}	

	


gms::math::stat::ArmaModel &
gms::math::stat::ArmaModel
::operator=(ArmaModel &&x) {
	
	if (this == &x) return (*this);

	m_newn = x.m_newn, m_nmax = x.m_nmax,
	m_mmax = x.m_mmax, m_kq = x.m_kq,
	m_kp = x.m_kp, m_aicm = x.m_aicm;

	

	return (*this);
}

void
gms::math::stat
::ArmaModel::compute_arma_model(int32_t n,
				int32_t lagh01,
				VAf64 & cyy1,
				int32_t newl1,
				VAi32 & iqi1,
				VAf64 & b1,
				VAi32 & ipi1,
				VAf64 & a1,
				int32_t lmax,
				int32_t fpexcepts[6]) {

	int32_t status{};

	GMS_ARMA_MODEL_FP_EXCEPT_PROLOG(status)

	MOD_TIMSAC_mp_AUTARMF(&n,
						  &lagh01,
						  &cyy1[0],
						  &newl1,
						  &iqi1[0],
						  &b1[0],
						  &ipi1[0],
						  &a1[0],
						  &m_newn,
						  &m_iq->operator[](0),
						  &m_b2->operator[](0),
						  &m_ip->operator[](0),
						  &m_a2->operator[](0),
						  &m_std->operator[](0),
						  &m_cxx2->operator[](0),
						  &m_g->operator[](0),
						  &m_saic->operator[](0),
						  &m_aicm,
						  &m_kq,
						  &m_kp,
						  &lmax,
						  &m_mmax,
						  &m_nmax );
						  
	GMS_ARMA_MODEL_FP_EXCEPT_EPILOG(status,fpexcepts)
		
}









