/*
 * Qmc Single Header
 * Commit: a0afeded6daf591e19aba26cc7a48c97dd596c2a
 * Generated: 16-10-2021 19:29:24
 *
 * ----------------------------------------------------------
 * This file has been merged from multiple headers.
 * Please don't edit it directly
 * ----------------------------------------------------------
 */
#ifndef QMC_H
#define QMC_H

#include <mutex>
#include <random> // mt19937_64, uniform_real_distribution
#include <vector>
#include <map>
#include <set>

#include <gsl/gsl_multifit_nlinear.h>

namespace integrators
{
    using U = unsigned long long int;
}

// Custom Types
#ifndef QMC_LOGGER_H
#define QMC_LOGGER_H

#include <functional> // reference_wrapper
#include <ostream> // ostream
#include <chrono> // chrono
#include <string> // to_string

namespace integrators
{
    struct Logger : public std::reference_wrapper<std::ostream>
    {
        //        using std::reference_wrapper<std::ostream>::reference_wrapper;

        bool display_timing;
        std::chrono::steady_clock::time_point reference_time; // time logger was initialised

        template<typename T> std::ostream& operator<<(T arg) const
        {
            std::string time_string = "";
            if( display_timing )
            {
                std::chrono::steady_clock::time_point now_time = std::chrono::steady_clock::now();
                time_string = "[" + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(now_time - reference_time).count()) + " ms] ";
            }
            return this->get() << time_string << arg;
        }
        std::ostream& operator<<(std::ostream& (*arg)(std::ostream&)) const { return this->get() << arg; }

        Logger(std::ostream& stream) : reference_wrapper(std::ref(stream)), display_timing(true), reference_time(std::chrono::steady_clock::now()) {}

    };
};

#endif
#ifndef QMC_RESULT_H
#define QMC_RESULT_H

namespace integrators
{
    template <typename T>
    struct result
    {
        T integral;
        T error;
        U n;
        U m;
        U iterations;
        U evaluations;
    };
};

#endif
#ifndef QMC_SAMPLES_H
#define QMC_SAMPLES_H

#include <vector> // vector
#ifndef QMC_MATH_MUL_MOD_H
#define QMC_MATH_MUL_MOD_H

#include <type_traits> // make_signed

namespace integrators
{
    namespace math
    {
        template <typename R, typename D>
#ifdef __CUDACC__
        __host__ __device__
#endif
        R mul_mod(U a, U b, U k)
        {
            // Computes: (a*b % k) correctly even when a*b overflows std::numeric_limits<typename std::make_signed<U>>
            // Assumes:
            // 1) std::numeric_limits<U>::is_modulo
            // 2) a < k
            // 3) b < k
            // 4) k < std::numeric_limits<typename std::make_signed<U>::type>::max()
            // 5) k < std::pow(std::numeric_limits<D>::radix,std::numeric_limits<D>::digits-1)
            using S = typename std::make_signed<U>::type;
            D x = static_cast<D>(a);
            U c = static_cast<U>( (x*b) / k );
            S r = static_cast<S>( (a*b) - (c*k) ) % static_cast<S>(k);
            return r < 0 ? static_cast<R>(static_cast<U>(r)+k) : static_cast<R>(r);
        };
    };
};

#endif

namespace integrators
{
    template <typename T, typename D>
    struct samples
    {
        std::vector<U> z;
        std::vector<D> d;
        std::vector<T> r;
        U n;

        D get_x(const U sample_index, const U integration_variable_index) const
        {
            using std::modf;

            D mynull;
            return modf( integrators::math::mul_mod<D,D>(sample_index,z.at(integration_variable_index),n)/(static_cast<D>(n)) + d.at(integration_variable_index), &mynull);
        }
    };
};

#endif
#ifndef QMC_ERRORMODE_H
#define QMC_ERRORMODE_H

namespace integrators
{
    enum ErrorMode : int
    {
        all = 1,
        largest = 2
    };
};

#endif

// Fit Functions
#ifndef QMC_FITFUNCTIONS_NONE_H
#define QMC_FITFUNCTIONS_NONE_H

#include <cstddef> // nullptr_t
#include <stdexcept> // logic_error

#ifndef QMC_HASBATCHING_H
#define QMC_HASBATCHING_H

#include <utility> // declval
#include <type_traits> // true_type, false_type

namespace integrators
{
    namespace core
    {
        template <typename I, typename T, typename D, typename U, typename = void>
        struct has_batching_impl : std::false_type {};
        template <typename I, typename T, typename D, typename U>
        struct has_batching_impl<I,T,D,U,std::void_t<decltype(std::declval<I>().operator()(std::declval<D*>(),std::declval<T*>(),std::declval<U>()))>> : std::true_type {};

        // Helper function for detecting if the user's functor has a operator(D* x, T* r, U batchsize) used for computing batches of points on CPU
        template <typename I, typename T, typename D, typename U> inline constexpr bool has_batching = has_batching_impl<I,T,D,U>::value;

    };
};

#endif

namespace integrators
{
    namespace fitfunctions
    {

        template <typename D>
        struct NoneFunction
        {
            static const int num_parameters = 0;
            const std::vector<std::vector<D>> initial_parameters = {};

            D operator()(const D x, const double* p) const
            {
                throw std::logic_error("fit_function called");
            }
        };
        template<typename I, typename D, U M>
        struct NoneTransform
        {
            static const U num_parameters = 0;

            I f; // original function
            const U number_of_integration_variables;
            D p[M][1]; // fit_parameters

            NoneTransform(const I& f) : f(f), number_of_integration_variables(f.number_of_integration_variables) {}

#ifdef __CUDACC__
            __host__ __device__
#endif
            auto operator()(D* x) -> decltype(f(x))
            {
                return f(x);
            }
            void operator()(D* x, decltype(f(x))* res, U count)
            {
                if constexpr (integrators::core::has_batching<I, decltype(f(x)), D, U>) {
                    f(x, res, count);
                } else {
                    for (U i = U(); i != count; ++i) {
                        res[i] = operator()(x + i * f.number_of_integration_variables);
                    }
                }
            }
        };

        template<typename I, typename D, U M>
        struct NoneImpl
        {
            using function_t = NoneFunction<D>;
            using jacobian_t = std::nullptr_t;
            using hessian_t = std::nullptr_t;
            using transform_t = NoneTransform<I,D,M>;
        };

        struct None
        {
            template<typename I, typename D, U M> using type = NoneImpl<I, D, M>;
        };
    };
};

#endif
#ifndef QMC_FITFUNCTIONS_POLYSINGULAR_H
#define QMC_FITFUNCTIONS_POLYSINGULAR_H

#include <stdexcept> // std::domain_error
#include <cmath> // abs
#include <vector>

// (Included Above): #include "../core/has_batching.hpp"

namespace integrators
{
    namespace fitfunctions
    {
        template <typename D>
        struct PolySingularFunction
        {
            static const int num_parameters = 6;
            const std::vector<std::vector<D>> initial_parameters = { {1.1,-0.1, 0.1,0.1, 0.9,-0.1} };

            D operator()(const D x, const double* p) const
            {
                using std::abs;
                // constraint: no singularity and singular terms have positive coefficients
                if (p[0]<=static_cast<D>(1.001) or p[1]>=static_cast<D>(-0.001))
                    return D(10.); // std::numeric_limits<D>::max() will sometimes result in fit parameters being NaN
                
                D p2 = abs(p[2]);
                D p3 = abs(p[3]);
                if(p2<1e-4) p2=0.;
                if(p3<1e-4) p3=0.;
                D y = p2*(x*(p[0]-D(1)))/(p[0]-x) + p3*(x*(p[1]-D(1)))/(p[1]-x)  + x*(p[4]+x*(p[5]+x*(D(1)-p2-p3-p[4]-p[5])));

                // constraint: transformed variable within unit hypercube
                if ( y<static_cast<D>(0) || y>static_cast<D>(1) )
                    return std::numeric_limits<D>::max();

                return y;
            }
        };

        template <typename D>
        struct PolySingularJacobian
        {
            static const int num_parameters = 6;

            D operator()(const D x, const double* p, const size_t parameter) const
            {
                using std::abs;
                if (parameter == 0) {
                    if(abs(p[2])<1e-4) return D(0);
                    return abs(p[2])*((D(1) - x)*x)/(x - p[0])/(x - p[0]);
                } else if (parameter == 1) {
                    if(abs(p[3])<1e-4) return D(0);
                    return abs(p[3])*((D(1) - x)*x)/(x - p[1])/(x - p[1]);
                } else if (parameter == 2) {
                    if(abs(p[2])<1e-4) return D(0);
                    return ((x*(p[0]-D(1)))/(p[0]-x) -x*x*x) * ((p[2] < 0) ? D(-1) : D(1));
                } else if (parameter == 3) {
                    if(abs(p[3])<1e-4) return D(0);
                    return ((x*(p[1]-D(1)))/(p[1]-x) -x*x*x) * ((p[3] < 0) ? D(-1) : D(1));
                } else if (parameter == 4) {
                    return  x*(D(1)-x*x);
                } else if (parameter == 5) {
                    return  x*x*(D(1)-x);
                } else {
                    throw std::domain_error("fit_function_jacobian called with invalid parameter: " + std::to_string(parameter));
                }
            }
        };
        template <typename D>
        struct PolySingularHessian
        {
            D operator()(const D x, const double* v, const double* p) const
            {
                using std::abs;
                D xmp0 = x-p[0];
                D xmp1 = x-p[1];
                return x*(D(1)-x)*D(2)*(v[0]*(abs(p[2])*v[0]+(x - p[0])*v[2])/xmp0/xmp0/xmp0 + v[1]*(abs(p[3])*v[1]+(x - p[1])*v[3])/xmp1/xmp1/xmp1);
          }
        };

        template<typename I, typename D, U M>
        struct PolySingularTransform
        {
            static const U num_parameters = 6;

            I f; // original function
            const U number_of_integration_variables;
            D p[M][num_parameters]; // fit_parameters

            PolySingularTransform(const I& f) : f(f), number_of_integration_variables(f.number_of_integration_variables) {}

#ifdef __CUDACC__
            __host__ __device__
#endif
            auto operator()(D* x) -> decltype(f(x)) 
            {
                using std::abs;
                D wgt = 1;
                for (U d = 0; d < number_of_integration_variables ; ++d)
                {
                    D p2 = abs(p[d][2]);
                    D p3 = abs(p[d][3]);
                    wgt *= p2*p[d][0]*(p[d][0]-D(1))/(p[d][0]-x[d])/(p[d][0]-x[d]) + p3*p[d][1]*(p[d][1]-D(1))/(p[d][1]-x[d])/(p[d][1]-x[d]) + p[d][4] + x[d]*(D(2)*p[d][5]+x[d]*D(3)*(D(1)-p2-p3-p[d][4]-p[d][5]));
                    x[d] = p2*(x[d]*(p[d][0]-D(1)))/(p[d][0]-x[d]) + p3*(x[d]*(p[d][1]-D(1)))/(p[d][1]-x[d])  + x[d]*(p[d][4]+x[d]*(p[d][5]+x[d]*(D(1)-p2-p3-p[d][4]-p[d][5])));
                    if ( x[d] > D(1) || x[d] < D(0) ) return D(0);
                }
                return wgt * f(x);
            }
            void operator()(D* x, decltype(f(x))* res, U count)
            {
                if constexpr (integrators::core::has_batching<I, decltype(f(x)), D, U>) {
                    auto xx = x;
                    D* wgts = new D[count];
                    for (U i = 0; i!= count; ++i, xx+=number_of_integration_variables) {
                        wgts[i] = 1;
                        for (U d = 0; d < number_of_integration_variables ; ++d)
                        {
                            D p2 = abs(p[d][2]);
                            D p3 = abs(p[d][3]);
                            wgts[i] *= p2*p[d][0]*(p[d][0]-D(1))/(p[d][0]-xx[d])/(p[d][0]-xx[d]) + p3*p[d][1]*(p[d][1]-D(1))/(p[d][1]-xx[d])/(p[d][1]-xx[d]) + p[d][4] + xx[d]*(D(2)*p[d][5]+xx[d]*D(3)*(D(1)-p2-p3-p[d][4]-p[d][5]));
                            xx[d] = p2*(xx[d]*(p[d][0]-D(1)))/(p[d][0]-xx[d]) + p3*(xx[d]*(p[d][1]-D(1)))/(p[d][1]-xx[d])  + xx[d]*(p[d][4]+xx[d]*(p[d][5]+xx[d]*(D(1)-p2-p3-p[d][4]-p[d][5])));
                            if ( xx[d] > D(1) || xx[d] < D(0) ) wgts[i] = D(0);
                        }
                    }    
                    f(x, res, count);
                    for (U i = 0; i!= count; ++i, xx+=number_of_integration_variables) {
                        res[i] = wgts[i] * res[i];
                    }
                    delete[] wgts;
                } else {
                    for (U i = U(); i != count; ++i) {
                        res[i] = operator()(x + i * f.number_of_integration_variables);
                    }
                }
            }
        };

        template<typename I, typename D, U M>
        struct PolySingularImpl
        {
            using function_t = PolySingularFunction<D>;
            using jacobian_t = PolySingularJacobian<D>; // set to std::nullptr_t to compute numerically
            using hessian_t = PolySingularHessian<D>; // set to std::nullptr_t to compute numerically (also set fitparametersgsl.trs = gsl_multifit_nlinear_trs_lm);
            using transform_t = PolySingularTransform<I,D,M>;
        };
        
        struct PolySingular
        {
            template<typename I, typename D, U M> using type = PolySingularImpl<I, D, M>;
        };
        
    };
};

#endif

// Periodizing Transforms
#ifndef QMC_TRANSFORMS_NONE_H
#define QMC_TRANSFORMS_NONE_H

#include <cstddef> // nullptr_t

// (Included Above): #include "../core/has_batching.hpp"

namespace integrators
{
    namespace transforms
    {
        template<typename I, typename D>
        struct NoneImpl
        {
            I f; // original function
            const U number_of_integration_variables;

            NoneImpl(I f) : f(f), number_of_integration_variables(f.number_of_integration_variables) {}

#ifdef __CUDACC__
            __host__ __device__
#endif
            auto operator()(D* x) -> decltype(f(x)) const
            {
                return f(x);
            }
            void operator()(D* x, decltype(f(x))* res, U count)
            {
                if constexpr (integrators::core::has_batching<I, decltype(f(x)), D, U>) {
                    f(x, res, count);
                } else {
                    for (U i = U(); i != count; ++i) {
                        res[i] = operator()(x + i * f.number_of_integration_variables);
                    }
                }
            }
        };
        struct None
        {
            template<typename I, typename D, U M> using type = NoneImpl<I, D>;
        };
    };
};

#endif
#ifndef QMC_TRANSFORMS_BAKER_H
#define QMC_TRANSFORMS_BAKER_H

// (Included Above): #include "../core/has_batching.hpp"

namespace integrators
{
    namespace transforms
    {
        template<typename I, typename D>
        struct BakerImpl
        {
            I f; // original function
            const U number_of_integration_variables;

            BakerImpl(I f) : f(f), number_of_integration_variables(f.number_of_integration_variables) {}

#ifdef __CUDACC__
            __host__ __device__
#endif
            auto operator()(D* x) -> decltype(f(x)) const
            {
                D wgt = 1;
                for (U s = 0; s < number_of_integration_variables; s++)
                {
                    x[s] = D(1) - fabs(D(2)*x[s]-D(1)) ;
                    // loss of precision can cause x < 0 or x > 1 must keep in x \elem [0,1]
                    if (x[s] > D(1)) x[s] = D(1);
                    if (x[s] < D(0)) x[s] = D(0);
                }
                return wgt * f(x);
            }
            void operator()(D* x, decltype(f(x))* res, U count)
            {
                if constexpr (integrators::core::has_batching<I, decltype(f(x)), D, U>) {
                    auto xx = x;
                    D wgt = 1;
                    for (U i = 0; i!= count; ++i, xx+=number_of_integration_variables) {
                        for(U s = 0; s<number_of_integration_variables; s++)
                        {
                            xx[s] = D(1) - fabs(D(2)*xx[s]-D(1)) ;
                            // loss of precision can cause x < 0 or x > 1 must keep in x \elem [0,1]
                            if (xx[s] > D(1)) xx[s] = D(1);
                            if (xx[s] < D(0)) xx[s] = D(0);
                        }
                    }
                    f(x, res, count);
                    for (U i = 0; i!= count; ++i, xx+=number_of_integration_variables) {
                        res[i] = wgt * res[i];
                    }
                } else {
                    for (U i = U(); i != count; ++i) {
                        res[i] = operator()(x + i * f.number_of_integration_variables);
                    }
                }
            }
        };
        struct Baker
        {
            template<typename I, typename D, U M> using type = BakerImpl<I, D>;
        };
    };
};

#endif
#ifndef QMC_TRANSFORMS_KOROBOV_H
#define QMC_TRANSFORMS_KOROBOV_H

#include <type_traits> // integral_constant

#ifndef QMC_TRANSFORMS_DETAIL_IPOW_H
#define QMC_TRANSFORMS_DETAIL_IPOW_H

#include <type_traits> // enable_if

namespace integrators
{
    namespace transforms
    {
        namespace detail
        {
            /*
             * Power function: IPow<D,i>(d) raises the D d to the U power i
             */
            template<typename D, U n, typename = void>
            struct IPow // n%2 == 0 && n != 0
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                static D value(D base)
                {
                    D power = IPow<D,n/2>::value(base);
                    return power * power;
                }
            };
            template<typename D, U n>
            struct IPow<D, n, typename std::enable_if< n%2 != 0 && n != 0>::type>
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                static D value(D base)
                {
                    D power = IPow<D,(n-1)/2>::value(base);
                    return base * power * power;
                }
            };
            template<typename D, U n>
            struct IPow<D, n, typename std::enable_if< n == 0>::type>
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                static D value(D base)
                {
                    return D(1);
                }
            };
        };
    };
};

#endif
#ifndef QMC_TRANSFORMS_DETAIL_BINOMIAL_H
#define QMC_TRANSFORMS_DETAIL_BINOMIAL_H

#include <type_traits> // enable_if

namespace integrators
{
    namespace transforms
    {
        namespace detail
        {
            /*
             * Binomial Coefficients: Binomial<n,k>::value gives the type U binomial coefficient (n k)
             */
            template<U n, U k, typename = void>
            struct Binomial
            {
                constexpr static U value = (Binomial<n-1,k-1>::value + Binomial<n-1,k>::value);
            };

            // Note: potential for optimisation
            // k > n-k ? Binomial<n,n-k> : Binomial<n,k>

            template<U n, U k>
            struct Binomial<n, k, typename std::enable_if<n < k>::type>
            {
                constexpr static U value = 0;
            };

            template<U n, U k>
            struct Binomial<n, k, typename std::enable_if<k == 0>::type>
            {
                constexpr static U value = 1;
            };

            template<U n, U k>
            struct Binomial<n, k, typename std::enable_if<n == k && k != 0>::type>
            {
                constexpr static U value = 1;
            };

            // require declaration
            template<U n, U k, typename T>
            constexpr U Binomial<n,k,T>::value;
        };
    };
};

#endif
#ifndef QMC_TRANSFORMS_DETAIL_KOROBOV_COEFFICIENT_H
#define QMC_TRANSFORMS_DETAIL_KOROBOV_COEFFICIENT_H

#include <type_traits> // enable_if

// (Included Above): #include "binomial.hpp"

namespace integrators
{
    namespace transforms
    {
        namespace detail
        {
            /*
             * Korobov Coefficients
             */
            template<typename D, U k, U a, U b, typename = void>
            struct KorobovCoefficient
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                const static D value()
                {
                    return (D(-1)*(D(b)-D(k)+D(1))*D(KorobovCoefficient<D,k-1,a,b>::value())*(D(a)+D(k)))/(D(k)*(D(a)+D(k)+D(1)));
                }
            };

            template<typename D, U k, U a, U b>
            struct KorobovCoefficient<D, k, a, b, typename std::enable_if<k == 0>::type>
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                const static D value()
                {
                    return ((D(a)+D(b)+D(1))*D(Binomial<a+b,b>::value))/(D(a)+D(1));
                }
            };
        };
    };
};

#endif
#ifndef QMC_TRANSFORMS_DETAIL_KOROBOV_TERM_H
#define QMC_TRANSFORMS_DETAIL_KOROBOV_TERM_H

#include <type_traits> // enable_if

// (Included Above): #include "binomial.hpp"
// (Included Above): #include "korobov_coefficient.hpp"

namespace integrators
{
    namespace transforms
    {
        namespace detail
        {
            /*
             * Korobov Transform Terms
             */
            template<typename D, U k, U a, U b, typename = void>
            struct KorobovTerm
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                const static D value(const D& x)
                {
                    return KorobovTerm<D,k-1,a,b>::value(x)*x+KorobovCoefficient<D,b-k,a,b>::value();
                }
            };
            template<typename D, U k, U a, U b>
            struct KorobovTerm<D,k, a, b, typename std::enable_if<k == 0>::type>
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                const static D value(const D& x)
                {
                    return KorobovCoefficient<D,b,a,b>::value();
                }
            };
        };
    };
};

#endif
// (Included Above): #include "../core/has_batching.hpp"

namespace integrators
{
    namespace transforms
    {
        /*
         * Korobov Transform: Korobov<D,r0,r1>(func) takes the weight r0,r1 Korobov transform of func
         */
        template<typename I, typename D, U r0, U r1 = r0>
        struct KorobovImpl
        {
            I f; // original function
            const U number_of_integration_variables;

            KorobovImpl(I f) : f(f), number_of_integration_variables(f.number_of_integration_variables) {}

#ifdef __CUDACC__
            __host__ __device__
#endif
            auto operator()(D* x) -> decltype(f(x)) const
            {
                D wgt = 1;
                const D prefactor = (D(r0)+D(r1)+D(1))*detail::Binomial<r0+r1,r0>::value;
                for(U s = 0; s<number_of_integration_variables; s++)
                {
                    wgt *= prefactor*detail::IPow<D,r0>::value(x[s])*detail::IPow<D,r1>::value(D(1)-x[s]);
                    x[s] = detail::IPow<D,r0+1>::value(x[s])*detail::KorobovTerm<D,r1,r0,r1>::value(x[s]);
                    // loss of precision can cause x < 0 or x > 1 must keep in x \elem [0,1]
                    if (x[s] > D(1)) x[s] = D(1);
                    if (x[s] < D(0)) x[s] = D(0);
                }
                return wgt * f(x);
            }
            void operator()(D* x, decltype(f(x))* res, U count)
            {
                if constexpr (integrators::core::has_batching<I, decltype(f(x)), D, U>) {
                    auto xx = x;
                    D* wgts = new D[count];
                    for (U i = 0; i!= count; ++i, xx+=number_of_integration_variables) {
                        wgts[i] = 1;
                        const D prefactor = (D(r0)+D(r1)+D(1))*detail::Binomial<r0+r1,r0>::value;
                        for(U s = 0; s<number_of_integration_variables; s++)
                        {
                            wgts[i] *= prefactor*detail::IPow<D,r0>::value(xx[s])*detail::IPow<D,r1>::value(D(1)-xx[s]);
                            xx[s] = detail::IPow<D,r0+1>::value(xx[s])*detail::KorobovTerm<D,r1,r0,r1>::value(xx[s]);
                            // loss of precision can cause x < 0 or x > 1 must keep in x \elem [0,1]
                            if (xx[s] > D(1)) xx[s] = D(1);
                            if (xx[s] < D(0)) xx[s] = D(0);
                        }
                    }
                    f(x, res, count);
                    for (U i = 0; i!= count; ++i, xx+=number_of_integration_variables) {
                        res[i] = wgts[i] * res[i];
                    }
                    delete[] wgts;
                } else {
                    for (U i = U(); i != count; ++i) {
                        res[i] = operator()(x + i * f.number_of_integration_variables);
                    }
                }
            }

        };
        template<U r0, U r1 = r0>
        struct Korobov
        {
            template<typename I, typename D, U M> using type = KorobovImpl<I, D, r0, r1>;
        };
    };
};

#endif
#ifndef QMC_TRANSFORMS_SIDI_H
#define QMC_TRANSFORMS_SIDI_H

#include <type_traits> // enable_if
#include <cmath> // sin, acos

// (Included Above): #include "detail/binomial.hpp"
#ifndef QMC_TRANSFORMS_DETAIL_FACTORIAL_H
#define QMC_TRANSFORMS_DETAIL_FACTORIAL_H

#include <type_traits> // enable_if

namespace integrators
{
    namespace transforms
    {
        namespace detail
        {
            /*
             * Factorial: Factorial<U,n>::value gives the type U factorial of n
             */
            template<U n, typename = void>
            struct Factorial
            {
                constexpr static U value = n*Factorial<n-1>::value;
            };

            template<U n>
            struct Factorial<n, typename std::enable_if<n == 0>::type>
            {
                constexpr static U value = U(1);
            };

            template<U n, typename T>
            constexpr U Factorial<n,T>::value;
        };
    };
};

#endif
// (Included Above): #include "detail/ipow.hpp"
#ifndef QMC_TRANSFORMS_DETAIL_SIDI_COEFFICIENT_H
#define QMC_TRANSFORMS_DETAIL_SIDI_COEFFICIENT_H

#include <type_traits> // enable_if

// (Included Above): #include "binomial.hpp"
// (Included Above): #include "ipow.hpp"

namespace integrators
{
    namespace transforms
    {
        namespace detail
        {
            /*
             * Sidi Coefficients
             */
            template<typename D, U k, U r, typename = void>
            struct SidiCoefficient{};

            // Odd r
            template<typename D, U k, U r>
            struct SidiCoefficient<D, k, r, typename std::enable_if<(r % 2) != 0>::type>
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                const static D value()
                {
                    return IPow<D,(r-U(1))/U(2)-k>::value(D(-1))*D(Binomial<r,k>::value)/(D(2)*k-r);
                }
            };

            // Even r
            template<typename D, U k, U r>
            struct SidiCoefficient<D, k, r, typename std::enable_if<(r % 2) == 0>::type>
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                const static D value()
                {
                    return IPow<D,r/U(2)-k>::value(D(-1))*D(Binomial<r,k>::value)/(D(2)*k-r);
                }
            };
        };
    };
};

#endif
#ifndef QMC_TRANSFORMS_DETAIL_SIDI_TERM_H
#define QMC_TRANSFORMS_DETAIL_SIDI_TERM_H

#include <type_traits> // enable_if
#include <cmath> // cos, sin

// (Included Above): #include "binomial.hpp"
// (Included Above): #include "sidi_coefficient.hpp"

namespace integrators
{
    namespace transforms
    {
        namespace detail
        {
            /*
             * Sidi Transform Terms
             */
            template<typename D, U k, U r, typename = void>
            struct SidiTerm{};

            // Odd r
            template<typename D, U k, U r>
            struct SidiTerm<D, k, r, typename std::enable_if<(r % 2) != 0 && (k != 0)>::type>
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                const static D value(const D& x, const D pi)
                {
                    using std::cos;
                    return SidiCoefficient<D,k,r>::value()*(cos((D(2)*D(k)-D(r))*pi*x) - D(1)) + SidiTerm<D,k-U(1),r>::value(x,pi);
                }
            };
            template<typename D, U k, U r>
            struct SidiTerm<D, k, r, typename std::enable_if<((r % 2) != 0) && (k == 0)>::type>
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                const static D value(const D& x, const D pi)
                {
                    using std::cos;
                    return SidiCoefficient<D,0,r>::value()*(cos(-D(r)*pi*x) - D(1));
                }
            };

            // Even r
            template<typename D, U k, U r>
            struct SidiTerm<D, k, r, typename std::enable_if<(r % 2) == 0 && (k != 0)>::type>
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                const static D value(const D& x, const D pi)
                {
                    using std::sin;
                    return SidiCoefficient<D,k,r>::value()*sin((D(2)*D(k)-D(r))*pi*x) + SidiTerm<D,k-U(1),r>::value(x,pi);
                }
            };
            template<typename D, U k, U r>
            struct SidiTerm<D, k, r, typename std::enable_if<((r % 2) == 0) && (k == 0)>::type>
            {
#ifdef __CUDACC__
                __host__ __device__
#endif
                const static D value(const D& x, const D pi)
                {
                    using std::sin;
                    return SidiCoefficient<D,0,r>::value()*sin(-D(r)*pi*x) + D(1)/D(2)*pi*Binomial<r,r/2>::value*x;
                }
            };

        };
    };
};

#endif
// (Included Above): #include "../core/has_batching.hpp"

namespace integrators
{
    namespace transforms
    {
        /*
         * Sidi Transform: Sidi<D,r>(func) takes the weight r Sidi transform of func
         */
        template<typename I, typename D, U r, typename = void>
        struct SidiImpl{};

        // Odd r
        template<typename I, typename D, U r>
        struct SidiImpl<I, D, r, typename std::enable_if<(r % 2) != 0 && (r != 0)>::type>
        {
            I f; // original function
            const U number_of_integration_variables;
            const D pi = acos( D(-1) );

            SidiImpl(I f) : f(f), number_of_integration_variables(f.number_of_integration_variables) {}

#ifdef __CUDACC__
            __host__ __device__
#endif
            auto operator()(D* x) -> decltype(f(x)) const
            {
                using std::sin;

                D wgt = 1;

                const D fac1 = detail::Factorial<r>::value;
                const D fac2 = detail::Factorial<(r-U(1))/U(2)>::value;

                const D wgt_prefactor = pi/detail::IPow<D,r>::value(D(2))*fac1/fac2/fac2;
                const D transform_prefactor = D(1)/detail::IPow<D,U(2)*r-U(1)>::value(D(2))*fac1/fac2/fac2;
                for(U s = 0; s<number_of_integration_variables; s++)
                {
                    wgt *= wgt_prefactor*detail::IPow<D,r>::value(sin(pi*x[s]));
                    x[s] = transform_prefactor*detail::SidiTerm<D,(r-U(1))/U(2),r>::value(x[s],pi);
                    // loss of precision can cause x < 0 or x > 1 must keep in x \elem [0,1]
                    if (x[s] > D(1)) x[s] = D(1);
                    if (x[s] < D(0)) x[s] = D(0);
                }
                return wgt * f(x);
            }
            void operator()(D* x, decltype(f(x))* res, U count)
            {
                if constexpr (integrators::core::has_batching<I, decltype(f(x)), D, U>) {
                    auto xx = x;
                    D* wgts = new D[count];
                    for (U i = 0; i!= count; ++i, xx+=number_of_integration_variables) {
                        wgts[i] = 1;
                        const D fac1 = detail::Factorial<r>::value;
                        const D fac2 = detail::Factorial<(r-U(1))/U(2)>::value;

                        const D wgt_prefactor = pi/detail::IPow<D,r>::value(D(2))*fac1/fac2/fac2;
                        const D transform_prefactor = D(1)/detail::IPow<D,U(2)*r-U(1)>::value(D(2))*fac1/fac2/fac2;
                        for(U s = 0; s<number_of_integration_variables; s++)
                        {
                            wgts[i] *= wgt_prefactor*detail::IPow<D,r>::value(sin(pi*xx[s]));
                            xx[s] = transform_prefactor*detail::SidiTerm<D,(r-U(1))/U(2),r>::value(xx[s],pi);
                            // loss of precision can cause xx < 0 or xx > 1 must keep in xx \elem [0,1]
                            if (xx[s] > D(1)) xx[s] = D(1);
                            if (xx[s] < D(0)) xx[s] = D(0);
                        }
                    }
                    f(x, res, count);
                    for (U i = 0; i!= count; ++i, xx+=number_of_integration_variables) {
                        res[i] = wgts[i] * res[i];
                    }
                    delete[] wgts;
                } else {
                    for (U i = U(); i != count; ++i) {
                        res[i] = operator()(x + i * f.number_of_integration_variables);
                    }
                }
            }
        };

        // Even r
        template<typename I, typename D, U r>
        struct SidiImpl<I, D, r, typename std::enable_if<(r % 2) == 0 && (r != 0)>::type>
        {
            I f; // original function
            const U number_of_integration_variables;
            const D pi = acos( D(-1) );

            SidiImpl(I f) : f(f), number_of_integration_variables(f.number_of_integration_variables) {}

#ifdef __CUDACC__
            __host__ __device__
#endif
            auto operator()(D* x) -> decltype(f(x)) const
            {
                using std::sin;

                D wgt = 1;

                const D fac1 = detail::Factorial<r/U(2)-U(1)>::value;
                const D fac2 = detail::Factorial<r-U(1)>::value;

                const D wgt_prefactor = detail::IPow<D,r-U(2)>::value(D(2))*D(r)*fac1*fac1/fac2;
                const D transform_prefactor = D(r)/D(2)/pi*fac1*fac1/fac2;
                for(U s = 0; s<number_of_integration_variables; s++)
                {
                    wgt *= wgt_prefactor*detail::IPow<D,r>::value(sin(pi*x[s]));
                    x[s] = transform_prefactor*detail::SidiTerm<D,r/U(2)-U(1),r>::value(x[s],pi);
                    // loss of precision can cause x < 0 or x > 1 must keep in x \elem [0,1]
                    if (x[s] > D(1)) x[s] = D(1);
                    if (x[s] < D(0)) x[s] = D(0);
                }
                return wgt * f(x);
            }
            void operator()(D* x, decltype(f(x))* res, U count)
            {
                if constexpr (integrators::core::has_batching<I, decltype(f(x)), D, U>) {
                    auto xx = x;
                    D* wgts = new D[count];
                    for (U i = 0; i!= count; ++i, xx+=number_of_integration_variables) {
                        wgts[i] = 1;
                        const D fac1 = detail::Factorial<r/U(2)-U(1)>::value;
                        const D fac2 = detail::Factorial<r-U(1)>::value;

                        const D wgt_prefactor = detail::IPow<D,r-U(2)>::value(D(2))*D(r)*fac1*fac1/fac2;
                        const D transform_prefactor = D(r)/D(2)/pi*fac1*fac1/fac2;
                        for(U s = 0; s<number_of_integration_variables; s++)
                        {
                            wgts *= wgt_prefactor*detail::IPow<D,r>::value(sin(pi*xx[s]));
                            xx[s] = transform_prefactor*detail::SidiTerm<D,r/U(2)-U(1),r>::value(xx[s],pi);
                            // loss of precision can cause xx < 0 or xx > 1 must keep in xx \elem [0,1]
                            if (xx[s] > D(1)) xx[s] = D(1);
                            if (xx[s] < D(0)) xx[s] = D(0);
                        }
                    }
                    f(x, res, count);
                    for (U i = 0; i!= count; ++i, xx+=number_of_integration_variables) {
                        res[i] = wgts[i] * res[i];
                    }
                    delete[] wgts;
                } else {
                    for (U i = U(); i != count; ++i) {
                        res[i] = operator()(x + i * f.number_of_integration_variables);
                    }
                }
            }
        };

        // r == 0
        template<typename I, typename D, U r>
        struct SidiImpl<I, D, r, typename std::enable_if<r == 0>::type>
        {
            I f; // original function
            const U number_of_integration_variables;

            SidiImpl(I f) : f(f), number_of_integration_variables(f.number_of_integration_variables) {}

#ifdef __CUDACC__
            __host__ __device__
#endif
            auto operator()(D* x) -> decltype(f(x)) const
            {
                return f(x);
            }
            void operator()(D* x, decltype(f(x))* res, U count)
            {
                f(x, res, count);
            }
        };

        template<U r0>
        struct Sidi
        {
            template<typename I, typename D, U M> using type = SidiImpl<I, D, r0>;
        };

    };
};

#endif

namespace integrators
{
    template <
                 typename T, typename D, U M,
                 template<typename,typename,U> class P = transforms::Korobov<3>::template type,
                 template<typename,typename,U> class F = fitfunctions::None::template type,
                 typename G = std::mt19937_64, typename H = std::uniform_real_distribution<D>
             >
    class Qmc
    {

    private:

        H uniform_distribution{0,1};

        void init_z(std::vector<U>& z, const U n, const U number_of_integration_variables) const;
        void init_d(std::vector<D>& d, const U m, const U number_of_integration_variables);
        void init_r(std::vector<T>& r, const U m, const U r_size_over_m) const;

        template <typename I> void sample_worker(const U thread_id,U& work_queue, std::mutex& work_queue_mutex, const std::vector<U>& z, const std::vector<D>& d, std::vector<T>& r, const U total_work_packages, const U n, const U m,  I& func, const int device, D& time_in_ns, U& points_computed) const;
        template <typename I> void evaluate_worker(const U thread_id,U& work_queue, std::mutex& work_queue_mutex, const std::vector<U>& z, const std::vector<D>& d, std::vector<T>& r, const U n, I& func, const int device, D& time_in_ns, U& points_computed) const;
        template <typename I> result<T> sample(I& func, const U n, const U m, std::vector<result<T>> & previous_iterations);
        void update(const result<T>& res, U& n, U& m) const;
        template <typename I> result<T> integrate_no_fit_no_transform(I& func);

    public:

        Logger logger;
        G randomgenerator;
        U minn;
        U minm;
        D epsrel;
        D epsabs;
        U maxeval;
        U maxnperpackage;
        U maxmperpackage;
        ErrorMode errormode;
        U cputhreads;
        U cudablocks;
        U cudathreadsperblock;
        std::set<int> devices;
        std::map<U,std::vector<U>> generatingvectors;
        U verbosity;

        bool batching;

        U evaluateminn;

        size_t fitstepsize;
        size_t fitmaxiter;
        double fitxtol;
        double fitgtol;
        double fitftol;
        gsl_multifit_nlinear_parameters fitparametersgsl;

        U get_next_n(U preferred_n) const;

        template <typename I> result<T> integrate(I& func);
        template <typename I> samples<T,D> evaluate(I& func);
        template <typename I> typename F<I,D,M>::transform_t fit(I& func);
        Qmc();
        virtual ~Qmc() {}
    };
};

// Implementation
// (Included Above): #include "math/mul_mod.hpp"
#ifndef QMC_ARGSORT_H
#define QMC_ARGSORT_H

#include <algorithm> // std::sort
#include <numeric> // std::iota
#include <vector> // std::vector

namespace integrators
{
    namespace math
    {
        // argsort as suggested in https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
        template <typename T>
        std::vector<size_t> argsort(const std::vector<T> &v) {

          // initialize original index locations
          std::vector<size_t> idx(v.size());
          std::iota(idx.begin(), idx.end(), 0);

          // sort indexes based on comparing values in v
          std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

          // return vector of indices
          return idx;
        };
    };
};

#endif
#ifndef QMC_OVERLOADS_REAL_H
#define QMC_OVERLOADS_REAL_H

#include <cmath> // abs, sqrt

namespace integrators
{
    namespace overloads
    {
        template <typename T>
        T compute_variance(const T& mean, const T& variance, const T& sum, const T& delta )
        {
            return variance + delta*(sum - mean);
        };

        template <typename T>
        T compute_error(const T& variance)
        {
            using std::sqrt;
            using std::abs;
            return T(sqrt(abs(variance)));
        };

        template <typename T>
        T compute_variance_from_error(const T& error)
        {
            return T(error*error);
        };

        template <typename T, typename D>
        D compute_error_ratio(const result<T>& res, const D& epsrel, const D&epsabs, const ErrorMode errormode)
        {
            using std::abs;

            #define QMC_ABS_CALL abs(res.error/(res.integral*epsrel))

            static_assert(std::is_same<decltype(QMC_ABS_CALL),D>::value, "Downcast detected in integrators::overloads::compute_error_ratio. Please implement \"D abs(D)\".");
            return std::min(res.error/epsabs, QMC_ABS_CALL);

            #undef QMC_ABS_CALL
        };
    };
};

#endif
#ifndef QMC_OVERLOADS_COMPLEX_H
#define QMC_OVERLOADS_COMPLEX_H

#include <complex>
#include <cmath> // abs, sqrt
#include <stdexcept> // invalid_argument
#include <string> // to_string

#ifdef __CUDACC__
#include <thrust/complex.h>
#endif

namespace integrators
{
    namespace overloads
    {
        // Implementation
        template <typename T>
        T compute_variance_complex(const T& mean, const T& variance, const T& sum, const T& delta )
        {
            return variance + T(delta.real()*(sum.real() - mean.real()), delta.imag()*(sum.imag() - mean.imag()));
        }

        template <typename T>
        T compute_error_complex(const T& svariance)
        {
            using std::sqrt;
            using std::abs;
            return T(sqrt(abs(svariance.real())), sqrt(abs(svariance.imag())));
        };

        template <typename T>
        T compute_variance_from_error_complex(const T& error)
        {
            return T(error.real()*error.real(), error.imag()*error.imag());
        }

        template <typename T, typename D>
        D compute_error_ratio_complex(const result<T>& res, const D& epsrel, const D& epsabs, const ErrorMode errormode)
        {
            using std::abs;
            if( errormode == all )
            {
                return std::max(
                                std::min(res.error.real()/epsabs, res.error.real()/abs(res.integral.real()*epsrel)),
                                std::min(res.error.imag()/epsabs, res.error.imag()/abs(res.integral.imag()*epsrel))
                                );
            }
            else if ( errormode == largest )
            {
                return std::min(
                                std::max(res.error.real(),res.error.imag())/epsabs,
                                std::max(res.error.real(),res.error.imag())/(std::max(abs(res.integral.real()),abs(res.integral.imag()))*epsrel)
                                );
            }
            else
            {
                throw std::invalid_argument("Invalid errormode = " + std::to_string(errormode) + " passed to compute_error_ratio.");
            }
        };

        // Overloads (std::complex)
        template <typename T> std::complex<T> compute_variance(const std::complex<T>& mean, const std::complex<T>& variance, const std::complex<T>& sum, const std::complex<T>& delta ) { return compute_variance_complex(mean,variance,sum,delta); };
        template <typename T> std::complex<T> compute_error(const std::complex<T>& svariance) { return compute_error_complex(svariance); };
        template <typename T> std::complex<T> compute_variance_from_error(const std::complex<T>& error) { return compute_variance_from_error_complex(error); };
        template <typename T, typename D> D compute_error_ratio(const result<std::complex<T>>& res, const D& epsrel, const D& epsabs, const ErrorMode errormode) { return compute_error_ratio_complex(res, epsrel, epsabs, errormode); };

#ifdef __CUDACC__
        // Overloads (thrust::complex)
        template <typename T> thrust::complex<T> compute_variance(const thrust::complex<T>& mean, const thrust::complex<T>& variance, const thrust::complex<T>& sum, const thrust::complex<T>& delta ) { return compute_variance_complex(mean,variance,sum,delta); };
        template <typename T> thrust::complex<T> compute_error(const thrust::complex<T>& svariance) { return compute_error_complex(svariance); };
        template <typename T> thrust::complex<T> compute_variance_from_error(const thrust::complex<T>& error) { return compute_variance_from_error_complex(error); };
        template <typename T, typename D> D compute_error_ratio(const result<thrust::complex<T>>& res, const D& epsrel, const D& epsabs, const ErrorMode errormode) { return compute_error_ratio_complex(res, epsrel, epsabs, errormode); };
#endif
    };
};

#endif
#ifndef QMC_GENERATINGVECTORS_CBCPT_DN1_100_H
#define QMC_GENERATINGVECTORS_CBCPT_DN1_100_H

#include <vector>
#include <map>

namespace integrators
{
    namespace generatingvectors
    {
        inline std::map<U,std::vector<U>> cbcpt_dn1_100()
        {

            // Vectors generated using Dirk Nuyens' fastrank1pt.m tool https://people.cs.kuleuven.be/~dirk.nuyens/fast-cbc
            // Settings:
            // s = 100
            // omega=inline('2*pi^2*(x.^2-x+1/6)')
            // gamma=1/s
            // beta=1

            std::map<U,std::vector<U>> generatingvectors;
            generatingvectors[1021]={1,374,421,220,482,449,309,72,382,152,328,247,212,414,119,503,196,438,418,452,350,208,92,84,316,59,192,405,338,342,478,160,272,427,147,248,429,268,224,314,86,498,127,54,183,493,174,78,263,454,372,470,329,479,281,83,139,11,188,45,240,256,466,284,159,10,17,305,409,229,245,31,430,91,443,378,111,317,463,6,410,471,462,336,486,131,75,195,507,19,199,364,487,100,23,180,151,474,122,324};
            generatingvectors[1123]={1,438,413,324,169,121,429,395,127,510,531,110,90,185,131,236,462,553,538,547,336,520,322,35,342,433,424,160,119,546,379,294,19,44,227,26,300,537,408,92,495,312,89,74,41,299,376,456,522,385,214,145,142,283,203,218,549,331,197,465,20,478,113,128,63,240,28,278,486,68,471,305,391,8,258,392,493,85,504,542,207,67,153,390,225,344,282,422,200,382,158,36,375,209,473,554,430,387,224,380};
            generatingvectors[1237]={1,457,545,284,422,471,375,514,194,200,491,214,577,349,307,258,503,105,506,72,547,329,234,550,102,370,430,455,499,66,278,372,449,180,603,561,80,335,485,406,469,308,273,208,354,519,443,61,126,111,362,161,470,476,347,140,293,174,32,321,326,357,279,287,464,300,356,17,12,38,10,213,584,265,202,441,262,433,415,314,271,318,500,164,436,27,130,46,612,145,342,592,55,589,133,344,496,52,340,511};
            generatingvectors[1361]={1,380,405,533,214,172,611,280,481,251,549,649,131,96,447,635,99,562,436,54,147,646,426,622,625,200,671,229,58,302,158,469,312,354,105,573,422,30,256,308,34,186,349,643,389,638,395,148,513,620,123,596,7,261,381,582,545,416,22,84,412,609,534,277,152,134,225,72,248,260,235,471,264,71,391,224,498,339,598,362,306,160,119,326,434,341,173,209,259,184,161,559,404,548,623,585,658,16,233,227};
            generatingvectors[1499]={1,629,551,586,614,452,114,574,312,660,489,648,714,581,91,126,255,656,400,508,207,409,236,22,149,98,493,192,180,536,423,524,560,547,188,666,471,413,506,680,706,231,635,116,494,215,567,55,65,103,175,484,285,556,210,276,249,475,692,41,28,197,710,53,161,266,303,637,366,529,246,107,219,234,568,102,652,436,448,370,152,227,176,565,502,144,569,213,184,89,602,396,444,15,674,125,171,490,135,289};
            generatingvectors[1657]={1,612,494,270,770,722,168,203,563,314,656,129,784,679,521,475,341,100,797,327,396,80,699,402,52,763,746,714,384,819,472,122,410,176,618,369,776,93,278,459,171,413,211,748,113,147,583,431,766,404,12,217,74,377,65,281,162,757,806,737,257,106,155,411,85,136,438,378,446,133,733,561,334,693,660,464,706,515,169,625,182,79,585,805,375,397,473,72,356,502,790,287,348,558,569,107,23,313,617,354};
            generatingvectors[1811]={1,669,509,700,536,373,524,770,564,663,879,96,785,690,407,304,515,551,801,749,329,415,493,437,298,673,695,843,543,387,491,251,757,231,622,721,435,150,168,424,572,410,601,583,806,463,882,66,35,353,226,540,146,244,430,716,778,19,200,45,727,652,240,377,117,102,246,93,457,662,854,194,767,12,273,404,796,222,380,197,561,14,308,301,489,483,170,23,699,518,409,167,901,852,88,726,442,481,323,439};
            generatingvectors[1993]={1,835,758,447,542,466,794,867,685,816,631,848,601,782,345,644,611,528,986,72,56,510,39,370,491,475,620,959,789,823,309,768,921,876,577,991,107,961,897,678,703,318,125,285,739,141,402,377,364,294,565,20,312,343,650,903,590,97,680,407,737,162,924,483,64,319,38,502,23,180,676,918,265,235,473,378,260,552,626,177,709,399,629,743,713,381,829,521,444,663,426,496,174,244,571,730,186,383,592,942};
            generatingvectors[2203]={1,669,970,462,996,644,299,777,1059,369,587,687,954,434,874,489,921,819,447,853,833,129,621,471,985,407,656,1004,501,581,674,821,887,794,1064,553,513,606,228,201,61,292,509,639,1039,528,940,666,158,806,715,524,381,112,1067,766,919,139,848,405,951,418,136,52,327,31,47,40,216,84,236,192,830,70,1000,1044,824,976,742,747,459,399,1050,1053,1085,22,354,779,636,193,1099,624,733,771,927,453,7,751,536,957};
            generatingvectors[2411]={1,923,681,985,1094,704,388,1077,592,559,433,765,847,532,628,663,236,670,146,188,214,787,1121,957,160,646,567,31,395,380,499,840,634,52,366,833,180,702,903,928,1161,966,550,1194,408,342,243,363,873,1136,999,1040,601,675,688,98,879,919,229,105,259,1101,1010,252,400,762,1197,1168,10,982,832,553,462,888,617,165,810,285,412,854,1073,72,489,167,750,492,39,861,264,238,1115,852,1045,379,276,803,613,22,599,65};
            generatingvectors[2647]={1,972,1160,679,1257,383,718,462,634,1195,83,322,408,1100,1235,739,569,782,1309,341,1209,1175,813,316,72,1014,313,1118,280,234,192,141,933,920,1038,265,839,772,890,963,1134,528,1277,677,714,273,1032,1315,50,196,835,909,723,218,511,823,27,1177,744,1047,539,202,153,487,164,37,678,281,863,544,520,953,1076,1102,504,183,663,858,1141,779,22,1239,1181,211,783,905,287,129,1133,292,1215,194,683,1053,825,439,45,131,615,345};
            generatingvectors[2917]={1,1265,909,754,1375,1025,279,598,468,783,336,1216,1280,646,554,106,266,659,1388,223,1199,1168,543,1345,767,1416,546,1261,204,403,381,145,1446,1145,295,422,516,1323,305,171,392,214,1336,822,1231,916,906,327,233,201,444,1210,1021,226,369,992,1343,1082,348,192,586,18,1034,857,620,1194,706,1173,1431,949,368,130,417,741,1385,1291,300,352,460,429,188,1379,66,794,249,651,1390,1236,14,908,1047,23,722,871,689,607,228,331,20,512};
            generatingvectors[3203]={1,1183,1316,991,432,709,1536,911,604,740,1033,1581,979,1431,963,857,161,504,1504,1132,1396,1471,891,660,588,103,1290,694,1013,1464,1507,413,671,646,815,1230,1097,125,895,371,1576,1487,873,1027,1243,1276,92,1424,1206,1159,856,1429,989,926,1089,290,200,467,476,455,668,1254,1178,431,1492,973,77,733,172,1498,30,127,226,79,1407,874,792,653,1124,1218,880,937,143,986,1301,1452,1131,457,1295,893,85,109,824,1139,1148,420,587,744,74,1448};
            generatingvectors[3527]={1,1368,952,1110,1564,1459,1654,1476,840,1424,1297,229,418,385,488,338,1418,829,1547,301,1729,427,1737,1022,1515,1605,1748,667,381,145,1373,102,1757,597,202,1408,551,749,670,131,927,945,1574,649,1327,693,722,767,499,645,245,60,774,858,973,617,424,1557,54,502,1347,288,581,160,470,911,861,1479,537,1393,444,274,208,426,674,558,1265,742,968,642,47,813,38,239,65,412,151,1388,126,1084,1183,788,904,265,1188,1703,1168,383,1375,1170};
            generatingvectors[3877]={1,1469,850,1018,1189,523,1828,702,1509,1378,1272,1230,475,1485,872,1052,370,1321,817,1584,1098,1216,781,420,539,801,1457,438,1608,1336,500,1717,1784,1517,1663,1614,1435,1174,976,1732,207,362,1070,1358,353,822,1867,1307,1626,38,694,63,112,348,641,1067,1360,1533,727,770,928,1806,1265,92,1280,599,829,225,473,364,749,1461,444,1914,1477,980,1127,993,692,395,1048,1339,744,869,36,1131,1882,685,1787,1376,854,13,1387,1699,1191,1248,1637,999,584,1463};
            generatingvectors[4261]={1,1648,1902,1757,1533,2032,646,1043,366,719,432,285,469,1401,1783,960,1722,230,1319,168,1745,585,1571,1537,1108,1156,306,866,1272,1941,670,75,1024,618,43,844,880,1228,607,2064,71,1968,444,1748,1453,1547,1521,233,2026,2115,1377,684,847,1234,1847,1161,1736,200,1359,829,950,706,104,763,1905,1959,1424,1123,737,748,1711,396,1281,735,411,357,1998,1531,652,1667,752,2012,2127,1143,135,339,1915,458,1056,1994,1839,1392,1327,1865,496,315,1818,1230,2004,59};
            generatingvectors[4691]={1,1372,2097,1122,709,1034,659,1752,1462,835,1001,886,2286,1952,1912,509,1098,1847,346,2305,727,1356,1081,1540,351,1650,966,455,240,1184,1257,174,1882,1936,2138,328,624,875,2060,988,2333,133,1520,2094,1160,2024,542,365,1706,1375,1326,207,1087,677,1986,217,1558,1363,950,1648,1071,791,101,2130,370,1915,2027,1171,816,1176,2032,1022,1118,1113,531,807,2187,2047,1892,229,1266,822,560,788,1918,898,1107,2173,754,982,216,1445,1581,1838,1284,1846,641,1452,1276,1359};
            generatingvectors[5167]={1,1428,2393,982,521,1681,1839,2483,663,229,2135,1894,2437,2379,1057,1607,931,81,2331,2401,1267,2414,1659,1919,1353,502,267,2105,684,889,2382,966,2443,1566,1541,597,218,290,1982,2167,1575,1301,1037,694,214,727,2276,1763,317,2262,1927,1675,2531,1790,438,1605,160,853,2339,1694,591,1626,785,1021,2070,2099,1650,1257,354,47,993,56,1367,845,608,1482,1543,2492,620,1437,2205,26,679,92,1708,2169,452,263,1150,2015,532,2456,447,1883,2348,1087,2545,1536,1032,1439};
            generatingvectors[5683]={1,2169,2042,1538,1280,1113,2535,1385,834,1879,395,226,981,2237,342,873,2467,81,2449,650,125,1313,1270,413,1955,2719,2197,1845,1397,1341,1646,2587,165,1685,253,1259,137,1090,2013,1747,1804,758,2003,359,1692,229,60,1853,1170,1602,881,722,2790,680,281,693,1870,1076,1158,786,191,626,2514,2053,1775,2095,517,1410,2786,689,1464,926,1689,1698,1730,1495,1338,2576,453,1901,1989,1182,1935,2325,219,293,2835,1431,1298,134,2320,16,535,964,2241,950,19,571,1056,2838};
            generatingvectors[6247]={1,2423,2709,1622,2906,674,720,2221,910,2285,2593,2862,1459,897,1406,1588,804,2250,1867,2459,3034,939,1176,2110,949,789,1258,920,2751,1077,2385,1152,156,1210,2116,2534,1424,2445,1783,2898,2738,3069,1030,699,1583,1933,2990,1391,499,3110,277,2469,1278,2673,2048,3095,586,652,205,2800,1489,778,3006,952,1083,1300,420,982,2511,2425,2327,2273,614,1467,2462,1352,866,877,2951,2928,736,775,2869,1370,2766,1205,2506,1295,2211,220,2728,1573,1170,713,1444,1745,597,2440,404,336};
            generatingvectors[6863]={1,2626,2990,1543,2081,2892,641,1514,789,1269,348,1333,2226,2484,229,2027,1638,2268,608,867,929,3095,629,1025,1646,3001,2373,3366,2014,1300,523,2853,2586,1607,3403,1051,2916,1709,334,425,863,387,566,2926,2884,2846,2240,1255,69,1586,1879,1129,1204,1428,880,1522,1678,2534,1849,1167,447,197,2354,669,3353,2075,1219,2094,3102,3423,984,711,2821,1436,1984,1185,192,1232,901,465,1635,2751,1107,2261,2231,1148,743,2116,643,637,1400,3031,944,1375,564,762,1689,1641,1308,1563};
            generatingvectors[7549]={1,3121,1635,2730,2280,732,3488,880,3209,1489,2485,3324,590,2555,3260,649,748,1344,2100,526,3223,675,256,2108,1977,2885,2041,2931,1842,1650,1884,2897,1075,972,1126,2970,857,927,762,2211,803,1336,356,2337,2677,3601,3041,3461,1084,665,202,1415,2061,3552,2275,1289,638,516,106,1194,3622,2010,2182,2311,407,1145,2728,2466,2049,843,3563,3522,193,1460,1594,3591,2502,2527,2805,251,2365,1019,2627,2505,1569,2235,143,3088,483,1861,1557,2191,1620,18,3656,3347,2941,452,2514,1729};
            generatingvectors[8311]={1,3068,1811,1128,1964,516,4056,3313,3787,3725,2917,539,3189,3548,888,1243,2715,3678,410,2972,1120,595,3226,2395,1324,637,3735,3947,1850,3434,3172,1094,1652,2688,1288,823,3490,2923,1515,393,3921,3239,3247,3890,1491,1755,3449,750,981,63,781,3411,2781,216,1670,3013,1388,913,3217,1431,1696,2836,1052,2140,2299,46,326,678,720,2164,2654,3445,1821,1686,2415,3715,801,3397,1721,1193,1116,3504,1458,903,3039,1186,943,1103,680,1773,3437,4078,1677,1912,3593,152,2672,922,3414,1719};
            generatingvectors[9137]={1,2515,3818,2072,3705,3561,3756,976,1200,2629,2965,2893,3433,147,3783,1791,4216,1023,1616,4134,4486,2339,1045,391,1454,1777,2195,2474,4085,1745,1712,805,3595,1113,3507,4150,3082,1645,260,116,2431,2098,1586,4004,2442,3228,3612,3098,3052,681,1299,4311,1902,4399,86,2331,237,3707,382,3537,1443,1089,108,169,2986,3148,4533,831,2819,1499,708,894,1689,4168,2719,178,1073,77,4102,998,2190,1651,1875,552,60,2570,422,1757,1272,3205,1458,1915,181,1249,4465,2318,1864,1461,2300,465};
            generatingvectors[10061]={1,3850,2758,1065,676,4224,1699,4881,3515,2413,2652,1876,4935,4408,2057,644,4328,4629,2145,2634,1280,1140,1825,2906,122,1499,516,1313,575,2815,2478,1103,1856,3335,2087,4728,376,1087,1157,1628,790,4703,309,4447,1867,3835,1421,3179,4649,3648,2745,3893,2673,3266,3878,3975,1895,1585,3544,3944,3080,1973,1492,314,3193,648,2797,4945,2044,222,694,2915,4608,2831,3522,45,982,208,1175,548,4561,4997,923,2060,3610,3775,2883,1194,4780,4597,4978,1672,4542,1848,3796,2690,3995,4720,3364,4846};
            generatingvectors[11057]={1,4276,3002,4932,2671,1415,1956,520,1858,2433,3868,1723,2535,1169,188,4046,5425,1700,2029,1245,5389,1026,2003,4202,2166,1291,4544,2759,4709,1826,1071,250,4487,1836,894,4303,2525,3612,4989,2705,2810,4855,1947,201,5471,2824,1050,1396,1190,3935,2175,3458,4197,2360,3090,1091,3658,4444,946,592,2322,4950,3749,2865,1060,580,3231,3175,905,2845,1164,1068,1017,1582,4892,3319,1555,2619,1964,1327,2388,429,826,3455,3182,142,4327,2519,5486,4109,2430,209,4223,1315,4037,3145,4089,3991,2400,363};
            generatingvectors[12157]={1,3396,4694,3778,5270,5576,4206,5753,703,2405,3323,4365,4165,5410,3220,5289,3803,5063,1777,5870,4475,5926,1297,3985,3934,4180,2478,1040,1531,3156,4251,3723,1599,319,3967,5475,3578,1659,1651,1769,1721,1221,2545,1320,1162,1807,3593,1631,5638,425,5663,397,1823,3073,3927,5691,4121,3137,251,4215,3684,1264,2121,565,2468,6002,4564,1644,4595,1301,3468,4383,3481,5778,379,2775,5981,1119,2701,3695,1020,5162,3504,4418,5891,1908,4278,3543,2049,4465,1708,870,4537,1441,2794,2682,4316,342,2447,4085};
            generatingvectors[13381]={1,5532,5254,2015,2116,871,3259,1348,6468,2477,2395,2942,763,5593,4739,3982,1164,1032,1893,4115,2553,1816,4307,341,3746,971,5802,5853,3302,800,2349,3771,3357,641,1669,5389,3531,2203,5616,6541,842,560,2970,4630,4316,2429,5114,6182,5642,5247,6488,370,772,1791,748,1108,3566,4927,3033,4359,1148,4249,932,1227,6628,736,2168,6559,3197,721,6236,2470,1982,3972,202,241,5711,3009,3616,5603,6368,2747,4201,2258,3238,5813,850,4847,6406,4239,3019,5367,248,2567,2505,3722,5237,2667,2407,2673};
            generatingvectors[14713]={1,5637,4311,3905,2313,1626,4488,6284,7249,3846,3359,3078,3656,5295,4351,2525,4444,2728,5846,4699,5399,3425,1961,1323,6749,1607,5574,3305,1237,6016,6989,5328,4827,1577,754,2640,4418,452,364,5822,4770,3629,5101,532,5728,391,7237,6702,6394,760,4304,863,4135,4146,1229,6381,4116,1256,1774,242,6832,5953,6535,6714,5939,2724,4582,1170,112,2115,3661,3785,5011,6515,1133,3672,5127,4848,2663,3534,2908,3733,4425,2434,4478,2411,5500,6242,6764,285,1503,256,102,7329,5385,3106,6047,2032,263,6452};
            generatingvectors[16183]={1,3689,6135,4805,4694,3861,3429,7078,3175,2620,2547,3769,7146,4185,7233,2176,586,5572,704,6044,1204,6828,547,6961,5909,2691,1432,4845,1795,4383,4556,3524,1022,2964,3006,5265,1591,5989,6400,6240,639,1313,1903,7699,980,1109,6169,2841,6195,2365,5483,2062,5618,385,5763,287,3412,425,1847,1498,5149,5109,525,1815,4039,4135,4644,5704,8034,5565,7430,4024,3922,3765,7007,7303,7670,758,488,2368,2425,4181,5754,4922,4321,566,2435,7001,7300,7403,7695,6760,2299,4717,2346,7636,3029,5248,3747,1653};
            generatingvectors[17807]={1,6801,7999,5312,2438,2316,4090,538,1119,5801,7358,2074,1861,5201,1330,2784,1432,3462,2891,352,5650,1288,5558,6020,5828,5630,7028,1461,7507,6198,6498,457,4349,6840,5001,1248,1830,2303,2014,656,278,6721,2748,5293,2673,683,3041,6217,2042,6597,7575,6315,8374,1642,6347,4115,6529,5018,3296,3956,2684,5072,3795,6416,2843,8784,5117,4599,6983,7937,5913,2190,5381,520,7041,1627,6559,2019,6063,1394,2650,1467,4778,4422,5974,4706,753,2514,2566,3773,8557,8132,8860,4982,157,1129,8781,7694,6906,5416};
            generatingvectors[19583]={1,5421,4085,8266,5060,7621,2346,747,4782,3849,1991,6283,5315,7991,6480,8415,954,9142,9203,3503,3433,5542,7256,4217,1024,5243,3008,1419,1146,5101,8013,6340,9104,4767,8710,1548,217,6437,3176,5956,1196,4112,1846,7781,2333,9692,9128,2133,6036,5263,984,288,4175,4952,5495,7789,3390,8046,1253,1767,8386,5486,2553,6177,2437,5911,2179,8761,5035,7586,8358,3763,4260,8891,3713,4682,5600,234,5590,2412,5414,7575,4714,2685,1424,472,7166,4294,2928,2665,5433,5295,2242,6816,6088,7906,600,4229,879,6851};
            generatingvectors[21557]={1,5940,4705,6313,10084,1700,6932,1908,5592,3331,793,4576,7405,1007,7025,4102,5827,3408,4039,9483,9705,1633,1953,8974,1328,4959,2268,2423,5015,3094,10483,240,6464,3822,7468,5424,3604,6976,8199,5716,390,7328,4506,3852,10171,1445,5990,2102,2414,1064,10644,8438,7203,10184,2753,1617,4191,5305,4943,2808,2302,7672,1743,4730,4818,6755,6492,10555,8610,9623,1557,6234,8457,915,6197,9122,8778,2054,5078,4064,9145,3385,10146,9129,5998,7286,8704,10283,8550,5481,8106,9226,9831,2059,7856,1093,4351,8374,3647,5377};
            generatingvectors[23719]={1,8796,9824,11301,6623,8429,7176,10057,9901,7010,1908,8956,1462,3099,4373,7248,6379,6592,642,8072,5959,9171,5895,11669,4464,4827,9577,616,2943,4063,3053,9725,7850,6882,5051,8411,9789,4116,7670,11589,11686,2408,8374,1061,8226,11561,4273,9508,800,3164,10516,4854,10994,10353,6006,6816,11653,8878,6442,10255,7349,7148,9065,4925,4426,5512,11628,5561,5772,1665,8507,4636,4356,9892,4547,3565,9403,407,9072,7398,3442,8705,490,593,11393,2297,6969,9033,11817,687,6666,7336,5554,4264,7319,11104,6803,10338,9274,11462};
            generatingvectors[26083]={1,11458,9987,7869,4195,12308,4720,5040,12139,11872,11705,6389,10484,2958,5274,8120,3218,8421,1370,2547,6280,3100,11968,5499,6668,10170,11205,11624,6175,131,2003,3777,5145,9911,9435,11043,3813,2449,8279,11994,921,12028,2845,2210,3844,7711,11946,1595,9092,9006,8413,11515,12235,756,451,7503,916,7442,1756,278,4072,10364,2599,505,960,2913,2368,5488,10234,7151,2311,8083,1539,5899,10578,3572,5402,5994,6825,9699,12775,5755,594,3584,4397,328,1769,282,7562,9849,203,12414,11106,3073,11297,8862,1206,5334,3748,5738};
            generatingvectors[28669]={1,10974,7676,11986,8567,3929,5838,2524,9098,4861,4326,6173,6389,2145,7490,3546,12311,8673,3957,10411,13577,980,8064,13966,7191,2822,5447,6267,8252,3002,8695,4579,13286,13097,7139,6784,6855,910,8853,6691,12022,7031,4493,10116,11420,10182,6908,2083,12336,7588,6584,996,10917,8942,5895,10073,9125,5298,8169,2535,2208,11727,1112,5938,895,10641,5186,2733,234,6562,12200,13049,5235,14224,8447,8644,450,6344,2896,6817,9888,13915,8805,2543,6545,10278,12425,3046,12731,7709,11982,10011,10992,11553,6251,6964,10616,10392,14169,4773};
            generatingvectors[31531]={1,12194,8335,9171,10806,3885,6682,9497,14707,4689,11424,2366,13854,5092,13770,13533,15221,4444,9191,1423,875,10602,9633,13204,10435,5598,11035,12566,2705,12859,10758,7125,13508,14170,7846,14352,2955,1755,13414,7703,9869,4083,9697,4907,14060,9532,10975,12364,13794,8520,10428,3813,2745,13603,14202,5405,8495,9254,5881,3308,14341,13090,2463,991,10708,12474,15701,6419,14275,8587,3539,12101,10795,10139,815,4194,13712,13936,12752,1015,14312,14033,5985,1335,9428,1096,5115,5901,2889,15688,14724,3529,5223,3224,5840,10652,1359,8740,835,5355};
            generatingvectors[34687]={1,14564,10745,12209,7027,16829,2849,1670,11838,9148,7184,16714,6810,4805,9677,7532,13551,2643,14914,4452,13783,12355,12080,2606,5830,1775,6461,11719,12796,2459,1148,15281,8755,11780,12398,13444,10560,2215,3422,4057,7059,14638,1315,5582,406,4510,13668,10120,12551,5043,9345,17294,15225,2987,8021,1407,4735,2749,17149,4223,3681,12854,8284,14526,9375,12132,362,8253,6408,482,9552,3169,16609,15394,8123,1101,1033,186,16984,2018,4727,1626,7512,14417,6595,3462,15768,11429,16846,15241,10291,4482,1382,8511,4312,15791,1939,42,15149,17199};
            generatingvectors[38153]={1,16115,11664,6826,10005,15783,5686,1793,7091,13954,3220,15485,16583,15055,11919,18659,16998,14031,5413,10869,11207,11070,3572,12017,8580,16903,18448,4516,15694,10415,7643,3185,2621,6650,12157,2539,7894,9432,11486,10038,2223,12230,4780,1106,13552,5120,11939,7385,16472,16973,8215,662,9047,10394,9181,5698,12400,7935,1942,16954,15533,8660,2705,6286,7594,9306,13098,7118,2343,10318,619,1711,3488,10745,16106,12031,10461,2137,6660,3239,18067,13685,7808,3341,18049,18124,861,8379,1306,14996,5982,18715,13765,3477,15845,5584,14391,15725,7675,9103};
            generatingvectors[41969]={1,16266,15435,12445,9588,7262,19534,14710,19662,18967,4125,11679,13381,5439,3072,7454,13020,8632,6959,17226,15379,1656,3772,6544,14421,9733,1197,16751,4104,17894,5730,15986,19129,13885,12762,8050,737,10543,2697,11695,19601,1899,15644,12722,7727,1918,5637,2679,4076,16165,16695,9406,9264,18623,10378,9139,10104,9520,12565,19335,6214,10928,7371,19104,19186,7250,13547,10035,11890,6034,10414,8733,9727,14029,15886,10262,15596,17578,18057,12364,19153,11401,18780,8009,9619,14834,14805,5600,9153,14158,10332,13823,401,12026,15616,5390,7540,10299,8602,6988};
            generatingvectors[46171]={1,17070,10462,13109,8323,2352,11130,8456,5666,18333,12943,19864,4103,3681,13517,5872,10383,14993,6872,6662,1497,3027,8891,12391,3095,8046,444,14621,14161,3993,11706,8975,21506,21320,21278,3477,21487,20973,20568,15324,21980,16378,15707,7739,18699,1437,21178,18106,1832,2300,11480,22903,629,14343,18130,22276,18576,9114,12184,18972,19151,7765,5003,14501,10121,4772,18177,14197,11934,21578,5214,18076,11303,19599,4398,14207,8826,9034,20209,17462,5048,18431,16071,14308,21912,16716,1064,2028,22882,17033,9970,20273,12640,20166,22672,21862,9277,22753,3959,5451};
            generatingvectors[50789]={1,19431,22336,13991,4310,11463,10686,6652,18997,1642,24047,5547,17961,24507,3896,8889,2077,6448,21871,9104,20162,10618,19936,15144,20785,22028,21745,12336,6562,5213,9990,4191,12378,13869,20678,6880,15418,20811,8149,7761,7856,9809,8005,23987,24778,17941,5043,17513,24102,2413,23807,17328,16547,21196,15042,2374,18664,7480,10825,19490,21629,18016,22196,5005,3161,10066,11725,19086,17630,11225,4705,24518,18536,19556,8553,1389,9458,15355,23719,2883,14642,13235,19379,16670,21396,16767,23536,8352,6578,8126,23695,9906,25019,19071,15407,1335,5525,2044,3372,15092};
            generatingvectors[55871]={1,20670,15630,21442,16345,11511,12525,19046,6438,26082,14658,16875,24105,10257,2867,9223,18038,6691,23063,26719,19832,17342,728,15112,3690,5188,25404,4081,14054,14084,12129,17155,7222,16130,10523,17115,25622,19842,1817,11215,7637,14231,13469,8450,20639,10613,2207,19593,22600,11235,18730,4581,16972,9015,23705,3037,2676,7017,5279,26518,3903,4479,8640,26458,4420,9500,8996,9259,27104,25021,27525,17326,23275,12190,11025,1030,20330,17041,19922,23618,21707,13462,19522,4543,18390,24669,20442,26629,21315,6000,25884,17128,12118,8745,1056,9808,14768,718,8415,19329};
            generatingvectors[61463]={1,23526,9889,22496,24939,18003,5689,8563,1740,25618,10768,28000,18562,12548,30335,14818,28362,8123,2228,28848,12967,2679,3778,17909,4863,27528,10661,11235,6179,19716,19109,9634,969,26142,14076,21486,17363,7659,29034,15616,27905,16283,26366,10983,19872,22011,19829,13942,2832,443,24143,29918,18293,16766,14702,22319,9207,11724,1009,21025,27131,10629,10011,12055,25222,1266,26041,3755,25814,24043,26718,10348,7110,3097,25183,21481,11937,23745,14197,26629,18453,29859,26021,10020,6611,17211,10930,10494,21823,8429,16079,2727,28422,22577,7478,29894,25872,16198,22761,25565};
            generatingvectors[67601]={1,24821,28748,19803,25712,17700,21990,6942,10654,25539,9992,6363,29102,5706,17283,33486,12505,16046,22212,29524,14188,9775,27755,13823,6619,22780,31533,2769,4366,15671,24486,14437,21863,30042,10870,4200,19451,1329,24299,7982,18006,7340,10088,20269,869,7419,12381,27894,29799,6086,33248,28358,16753,27560,26782,27295,1628,12306,13675,12820,10556,14521,27111,32242,16290,25278,18099,4491,21235,6224,735,30026,33578,22722,26685,8414,25890,5781,1280,32215,22310,8381,12642,154,20511,27330,32819,20718,8314,2032,8678,12672,5276,20651,27666,24568,6817,21365,4721,25323};
            generatingvectors[74353]={1,27461,32398,11960,22089,18900,21142,12839,16356,28614,4586,35428,12731,5068,36316,13365,26653,4190,31774,6973,32995,36525,13037,33603,33392,17996,5011,4763,24005,17289,17726,16610,2720,36920,27570,23952,11530,8545,5151,12081,2281,31399,7997,30773,31808,19542,2412,23532,35178,7124,1853,4070,7794,25320,20030,20360,32909,21780,352,36110,16733,15681,10386,21699,25874,32433,35955,12300,21802,9444,126,4935,6606,35008,24882,32632,17568,19705,37153,34730,17673,2918,6612,21539,25912,6707,25356,14630,25335,10742,20286,667,28650,7863,11789,1083,34834,23772,25648,14816};
            generatingvectors[81799]={1,31315,23929,17708,36771,9473,10311,3965,26857,35596,3185,35140,16868,38615,14026,6539,38811,16930,40350,25153,10071,30791,22450,28184,5751,22934,16280,8831,12911,2446,4699,38874,18132,3810,21353,1423,28275,27465,13468,5711,11599,3574,23466,31291,24063,28394,21686,36141,25727,36678,8087,33514,16076,14922,11780,33626,35830,3952,2623,20001,9776,517,9367,6321,27746,18695,31364,9076,28646,2353,28532,5530,11379,1991,6124,29837,32582,30539,22411,38574,29355,28428,18293,10053,13647,25578,2439,25409,11315,1920,26977,8398,3483,23784,458,27859,16535,25532,33531,1695};
            generatingvectors[89963]={1,34364,24782,19475,41115,10440,2878,42592,35659,24026,20775,10685,29883,2622,42267,34741,33316,14640,4691,22304,32303,32640,9143,30130,16070,16409,29453,11125,12582,29007,27482,37093,35407,18629,34560,37888,43118,37713,30365,22460,24213,33087,28052,24266,19778,29829,4110,5216,40244,8985,4772,1216,3095,23389,3376,40083,29653,16475,10551,15076,11644,41138,13802,42004,6262,19989,10054,26174,23838,17482,24566,9655,28498,22100,11880,14666,38685,14435,30294,5360,17937,38148,9730,31513,14192,1658,28646,36648,40167,34978,22939,3238,16398,21814,25510,14660,40049,6911,29922,6447};
            generatingvectors[98963]={1,27258,23359,21445,20691,9170,36578,19051,46475,39818,16333,5906,30040,13999,21102,31074,41878,5683,6472,39424,25671,4189,1472,37142,41111,44913,24825,46348,7915,8297,43651,38814,16421,2307,21393,12716,38244,14504,32233,19426,49402,46515,26770,6968,46251,31824,33837,7611,20961,6046,43098,876,30625,35282,27330,8442,14034,36216,18296,12694,20219,35412,3734,2209,46083,18929,47666,13651,32624,7283,16555,37277,29857,39334,8606,15758,16266,4913,3715,10718,22322,28505,30598,13037,18612,29117,33526,47207,31526,38531,19316,9609,14169,31572,30658,34145,24962,1626,19008,20616};
            generatingvectors[108863]={1,42235,34421,15792,38785,22724,52662,6712,44274,46020,24694,51652,10403,53143,39258,24175,23642,6904,23339,48108,19473,2363,5000,21618,5415,20187,26152,7273,12188,34155,30651,51603,48490,39806,18673,4270,46786,7435,24613,33142,10759,32396,23033,36633,22498,3578,45737,47310,53734,1282,33493,48158,9012,32139,9791,13124,51276,54337,32551,7760,32623,28984,15637,9594,9198,4352,52422,31657,13317,44225,32335,5686,3403,47739,10060,49967,49437,19267,8833,6547,52563,40432,24102,27664,39931,13152,26977,8219,41082,37098,43163,41506,4855,54244,43191,17380,42270,7889,53257,15691};
            generatingvectors[119747]={1,50752,48825,27409,57182,16469,31715,12959,14656,17949,18621,41641,29371,30714,3739,25020,16881,52884,25743,12456,56389,9080,53374,2652,41567,3365,5537,26128,33535,38125,55325,34924,29257,35734,33802,42174,41058,59747,5776,51684,43142,49476,26232,31794,14047,47814,26529,37714,16784,36939,19810,33910,1096,12756,57576,54483,29847,28807,975,33995,768,10020,30834,54886,23985,39225,5331,37305,37110,23049,52057,59229,52976,311,44511,41982,4433,12855,25910,26942,57331,57024,57395,13139,43254,35407,21142,40426,18970,53695,27339,49863,35217,482,37000,34401,38363,47738,25250,12075};
            generatingvectors[131713]={1,57442,59612,50816,22396,45300,25780,17851,31187,63229,17537,15241,14267,28333,60285,4600,10705,30615,48583,42423,62911,57710,7906,54864,32765,2888,63075,47957,3976,55349,52779,17092,55493,29928,21566,4714,18361,63828,64310,47675,32800,45792,55598,53916,47566,34856,6827,52148,405,16406,3834,63154,15687,41590,44529,8752,48229,2545,16882,29739,27034,23932,62874,37219,7067,11831,35817,16015,23346,50030,41372,52857,45463,37911,6814,41088,9472,8507,53689,63353,1241,34058,59571,18440,22354,61666,65353,14229,25416,48247,9332,33000,50293,41535,35031,34369,32450,18907,30140,61255};
            generatingvectors[144887]={1,56003,51930,45573,19701,32059,21288,68871,29363,13618,50176,4326,10133,31819,16166,23212,31166,38585,9974,45476,3632,22136,59336,35925,41073,35308,56163,59941,15277,10616,58386,50789,26154,9317,65704,26084,62131,39986,16041,39816,30693,9605,18587,59004,40467,54605,12738,403,56789,66301,3861,67800,33531,4544,26238,40901,10708,13967,7455,18985,8484,31348,66494,52730,5838,2112,60281,24687,68844,49545,27427,39320,5730,59921,18755,5471,13031,55463,49662,71315,10806,69255,42358,15921,16450,33925,32137,37621,11878,26293,29168,59964,47424,41491,7247,67704,48072,49924,25223,58823};
            generatingvectors[159389]={1,65858,61603,71847,34917,57110,7043,9682,41237,71046,7633,59287,55648,37656,26391,30743,35443,16797,41571,41124,75757,32365,14909,21483,8208,46850,63355,43699,68432,34733,15712,57956,55831,41853,12175,46019,32319,51040,67461,10873,68154,57172,59686,25336,55216,32851,11679,29225,60148,23822,41337,64171,7445,30902,13464,69888,34038,61197,55554,21833,13781,57386,56904,31288,18367,62337,60673,45661,59925,29517,74262,33914,67255,3875,4311,64591,1960,73770,59193,20221,55795,70216,1432,70684,4580,4232,17829,41077,21910,46002,71401,59088,21406,53244,36437,75452,60388,73206,9381,44887};
            generatingvectors[175327]={1,72749,67074,40138,49165,37045,33902,31742,61180,31315,33358,43363,55882,57680,3572,63976,78458,84966,28550,46640,18675,34842,62463,37998,55557,41919,53428,62220,83443,53569,28658,32325,63854,7202,68764,44939,86017,26225,11097,3362,36348,75579,80690,4287,77395,44700,32490,57452,26314,11996,23370,55941,64919,74045,56533,64128,21207,34544,58597,42065,24534,313,69127,60402,12304,10684,2487,24788,15984,86089,85552,30613,7304,37549,62479,4620,71146,15791,47467,73535,86385,81340,35610,72818,61424,80873,32652,14143,29839,78592,6386,40824,81929,5855,64703,20478,61484,65962,75606,82419};
            generatingvectors[192847]={1,71225,45638,28634,30069,53664,18291,93699,43196,93193,73985,61429,80687,24403,89376,22228,74473,9874,78728,10960,55042,10134,80206,5174,68008,21163,18713,89971,85424,73357,29448,7631,84844,94904,28857,29073,23310,64000,92012,32803,72662,14622,58796,90327,20100,92307,65016,46057,88722,59594,60860,71374,5616,73251,12317,65132,84135,76583,2120,38484,67326,64602,59628,88857,13238,69150,76644,10076,89924,27105,83415,35199,43984,69097,29294,33775,46813,55359,61411,38346,27271,27540,42888,30256,12397,66530,18637,56627,71537,82478,34732,63625,70707,40855,12701,82565,49139,34167,16938,54334};
            generatingvectors[212131]={1,62758,92550,34565,58633,78848,10305,56208,25635,72841,41554,36765,91073,16120,36003,54893,66039,41820,26602,14517,95113,74399,95370,2937,30176,93765,80678,104918,103276,68169,68749,39605,66238,89501,93527,96827,60633,53352,96542,65405,94423,97178,58915,35606,37077,15681,17837,14712,19603,44831,25900,54585,2423,105338,67440,13533,27914,90805,6175,34036,85344,2087,24934,22105,64602,30459,32537,79898,1412,11965,7139,94046,56799,100333,87218,21116,52432,62067,25726,34604,56314,78145,54730,95654,81162,9099,4660,92370,47609,27601,101296,105293,27262,54142,7651,33489,37962,22937,85866,35180};
            generatingvectors[233341]={1,62944,83087,100587,37193,55183,11321,109505,65566,24695,49950,97560,4813,72115,112356,107241,90081,55738,45749,63823,26314,90843,43335,82482,15595,64735,21350,60099,22479,71127,8912,99185,64564,114241,11398,30780,7428,69562,22619,77312,96048,36042,106499,31685,60910,5344,72983,103055,110008,3885,90713,49438,50491,23505,105646,27930,108162,12671,7873,110287,21971,19990,13293,93394,97113,80676,94981,79595,10862,48454,110113,77085,80936,74103,61829,17495,31236,85491,14342,17745,68062,86713,26043,10930,40395,115769,56702,35550,93994,3040,55494,59579,58088,106582,53988,101350,28907,61642,92883,75534};
            generatingvectors[256687]={1,111833,61875,53746,81793,68891,11059,117959,50877,35111,75056,122947,89961,45166,54718,40426,4899,48243,43631,51180,65846,98173,61092,107119,80581,119011,25477,13248,23442,99151,91852,80660,53810,6899,12940,35940,88216,48336,8547,99778,84419,20262,22837,1299,8305,21024,17055,78996,8804,89733,70548,45359,45091,4854,102974,9217,13331,29225,9397,55961,40930,107050,66215,117500,7887,59465,84920,18850,59360,12636,65312,95331,42470,50596,121441,88249,96377,100545,25215,61708,72506,123178,122749,80333,106065,68541,93753,92406,107608,33580,48561,73544,49212,108531,90062,71493,33425,67939,1138,52634};
            generatingvectors[282349]={1,108093,118683,104484,121537,123012,77186,28988,98500,62024,30580,9050,76934,30157,50241,103985,69212,48112,106249,97281,38898,105643,140323,79702,30817,139551,51776,114844,19375,32386,55678,45092,65651,117234,79826,121446,55020,118436,44312,9331,69082,119200,23926,7670,140088,129868,127186,101129,78764,75980,58193,55815,19109,98685,27108,112158,33612,106971,99533,114161,43602,14787,70854,140451,123561,15607,15788,34709,95338,115605,107615,88997,111949,3045,70147,7771,90070,23751,3684,24676,135179,80596,6267,56522,29884,20472,111914,14704,137494,69311,16007,13524,24981,89418,41372,125091,50598,42732,79160,87702};
            generatingvectors[310577]={1,128243,57866,65494,137609,44142,100788,37493,27901,76814,142641,119283,9036,52584,147422,26734,98944,125173,57016,139703,44709,99374,33297,152321,144399,147938,123030,89304,30355,136240,82730,50607,134704,122836,61784,65809,66766,20320,63341,29620,101648,90122,58996,838,82609,74311,113711,32859,112329,121712,123976,76735,2290,7552,39652,94811,29075,17545,103108,53778,54962,5718,24935,66073,132570,91615,122958,110139,31946,9978,121240,124342,47137,133280,132534,40456,31656,42193,141832,64002,113046,72214,74826,138850,7113,58311,113253,18885,108044,95566,88542,61763,58387,114652,128158,73716,14890,47932,101143,105806};
            generatingvectors[341629]={1,132408,126212,91264,96594,53545,105887,130817,164982,69008,93464,7486,73002,167920,41541,57833,36713,48618,149351,138405,39013,95760,49083,125909,51908,156782,9449,72873,25722,153062,52618,101834,41059,52573,31866,63904,144456,131539,71341,95418,157619,129053,136747,59942,58658,77390,124005,66946,87395,45932,66263,152630,125435,140707,50721,82530,2418,76173,103097,122014,123379,80023,68703,94454,150199,35274,5676,134820,56822,42519,150347,51490,142483,128574,37733,136018,168770,82159,156474,160536,67571,150945,131838,90591,86585,24980,27793,150706,123917,115849,24820,95485,88511,161624,2149,153886,166625,132158,164246,40514};
            generatingvectors[375799]={1,155714,81324,69481,144359,76046,168096,31557,128091,30046,174948,33319,104471,73059,116140,22966,152713,28487,47232,22575,26858,90789,85548,65892,154500,90392,136707,153489,15587,146745,177141,79815,106190,34650,86398,13324,81250,184029,162045,52219,127857,2164,27842,115020,100832,111307,53628,113496,41720,178658,174549,119884,106958,39417,186481,166313,91268,142813,49943,136858,5624,157270,151523,70738,54456,62541,151667,33591,81999,105938,185283,157496,69135,90916,141152,62982,144093,15046,38321,7586,8924,104254,75926,115376,178510,42756,80303,121702,5606,17129,85221,107879,76830,149133,19107,37894,10787,17960,33502,135127};
            generatingvectors[413411]={1,151354,120785,71720,94683,170501,80046,169643,153751,87766,112929,122068,86283,110285,205229,113272,186204,152027,84054,79447,184821,24186,4894,89691,194663,78866,157083,131226,144108,180913,182711,91116,73358,49126,165160,24659,159726,71380,7838,52176,133994,19460,130489,73732,29810,76828,201685,57276,9152,37132,196814,81774,76418,70967,41768,85068,107414,49002,205131,139825,110887,128390,124105,62783,23681,124980,61293,167340,174098,29252,4071,125695,30697,163373,109834,136869,3664,134867,190542,3745,113658,118201,3791,20163,74567,175877,89503,50446,14464,182845,98195,70753,194799,112700,13335,39835,37085,134395,45440,185670};
            generatingvectors[454709]={1,166951,108693,53112,203016,163284,141944,210460,56052,189008,160237,92501,147541,119507,143036,97919,31905,178805,79944,20902,61944,69472,65350,29242,45619,205807,104596,126995,66647,146506,172404,78683,144232,189694,48124,103564,6890,83295,24474,95239,107954,108392,202025,139408,202188,87385,106323,127032,74156,57775,54578,163621,171930,56957,68387,146997,201996,75168,129841,140955,148424,166037,127818,190178,38967,216221,142197,98164,8353,156308,227141,8421,124095,133586,87705,125001,136584,155980,136139,77454,120255,207040,105799,215821,173634,47672,202255,218068,54043,70334,94998,104120,25916,140793,126772,18168,196860,53863,26498,100557};
            generatingvectors[500179]={1,207183,190227,226446,116452,108276,172060,19889,7849,184428,69392,221377,167536,62880,35048,49881,158387,176742,176088,88707,41350,203332,237239,93265,2199,128475,237963,184959,143136,149740,119166,192316,96330,224641,97682,31219,161582,137843,5979,27429,97086,10818,213981,68963,7005,138004,216891,139232,185571,112003,111019,81464,130952,47387,171797,40574,202285,242618,110502,97935,233921,216199,80079,109349,197423,165234,7639,184522,186353,149689,91574,29172,162000,33233,186558,58997,119528,182387,56187,234516,212250,74069,56659,173584,79349,219473,10585,203289,25524,123213,94449,43998,61429,133618,133025,222379,25464,8879,19138,189825};
            generatingvectors[550211]={1,232176,248304,191290,218032,213615,130977,126475,85262,41700,25179,243030,86064,47194,190404,240830,184632,193609,46755,19116,66664,19985,127443,249974,14829,260144,55645,263177,17679,177995,197018,256966,106565,3828,140684,122663,8205,95534,30249,216761,31272,174524,41135,253816,189404,100787,32663,85325,25961,229937,190991,32596,13893,199896,52258,36743,214452,210696,158052,236354,210292,139789,274056,48611,171912,14276,162594,218210,97662,126990,104636,90664,128954,139203,270890,201316,271927,48149,48788,116251,251113,74195,117525,180414,165272,94766,6075,3283,247348,55038,14900,106106,272083,84059,136328,170538,36943,141561,161382,153240};
            generatingvectors[605221]={1,231183,157309,107336,261394,89794,98821,198659,147622,232289,160659,103694,145028,83236,148862,173421,227706,14611,8845,70488,185599,9091,206301,216815,176273,115679,47333,76595,246713,232696,68146,86650,44746,21069,106870,105253,63926,265013,179742,123484,91446,218485,184376,55000,108457,196308,38476,285932,26852,264638,189374,126264,62159,140595,185981,195905,68663,149861,158764,85209,234399,25608,154459,259782,88432,174810,85294,68826,14895,91660,3410,42999,257164,122360,79283,200139,189449,1167,111753,164323,273330,143424,172813,70722,50078,39852,269285,252893,221140,81936,116422,56743,35693,67859,159122,299279,136765,269988,297560,140117};
            generatingvectors[665747]={1,195038,205880,238868,105284,311952,226771,128607,88134,152937,31541,175307,320356,258870,151564,132815,121655,309036,269222,177782,160468,289279,125419,114740,130448,140846,170948,36268,294606,25351,147458,232896,213951,207697,270747,139543,324488,78607,5275,310665,212731,158066,301718,144443,290649,205735,235227,99767,60466,280032,145292,284579,302159,104794,331667,203342,100849,233927,242570,116481,146764,317512,94429,267079,86487,101823,264570,296631,202222,239396,316880,326118,27565,189768,276597,200880,73171,51899,306713,87728,191470,48478,65393,135003,44076,179677,207483,53438,110135,125147,314817,236632,56748,111425,321138,264739,36873,300222,21855,116590};
            generatingvectors[732311]={1,309022,136073,177211,261033,85922,126346,277049,295477,321423,287914,65526,334495,346435,200348,236754,74517,166033,199179,344540,16434,162616,63527,35792,16562,272652,314346,105470,173040,124545,175978,47940,274320,121808,132447,41467,217175,237473,345931,309539,69019,297595,168071,280952,39386,49199,216451,94656,353076,161356,81815,135263,249816,338224,360549,254287,276866,360874,190235,120826,238120,6435,227936,146597,252467,30974,78437,56163,161023,226360,208747,21813,113624,20021,194340,261355,90440,245781,352838,286767,17790,75933,95679,360055,37631,111890,285225,252970,186477,203698,5757,279664,308695,28987,26428,321514,361008,68777,159747,60337};
            generatingvectors[805559]={1,297501,363557,108581,125660,165132,352786,172283,26922,70622,277999,216181,106560,397188,93737,31719,248080,365340,210207,206218,181335,38526,148510,400763,34913,323819,113492,43956,287018,221866,340639,270086,290522,246415,238172,179138,306907,247067,160764,121040,75083,72309,138664,78755,223896,252495,60371,2944,134877,36363,197519,289532,240480,269360,199194,236378,186345,79885,165425,332558,184340,234927,191906,311139,140347,291154,160072,400288,305201,117157,185975,391179,42271,27331,105881,274109,281880,186713,53169,350885,114191,41577,139029,122183,251262,116475,275180,331362,217031,114234,1027,345449,364432,350099,115319,183078,76869,207120,220869,161232};
            generatingvectors[886097]={1,327320,261356,199253,410738,40562,336990,273446,146951,69544,362852,303271,326371,299510,216706,378641,92255,195436,158513,89742,33399,305379,222391,314488,113933,22382,421119,414278,431121,322069,29444,1932,133241,143936,172069,199562,102705,339733,17917,117588,14592,33781,102603,167607,95127,142609,306834,341700,214800,20404,301753,419157,120676,140817,9571,438244,79949,56627,164378,44267,16116,438752,416293,176883,69167,259237,409982,227733,215045,208051,217964,309587,108672,105665,228226,191431,357168,219498,29039,191529,354155,211743,404492,94726,373416,423730,172836,55469,287644,266216,22990,404089,22440,317789,3100,86466,53696,433073,299684,211832};
            generatingvectors[974707]={1,403710,339829,273290,413164,141506,251926,177965,390797,373603,270385,67460,179451,51865,331151,456985,380501,126299,54820,38677,134277,263915,248890,110956,136039,363667,374316,5736,195218,40355,118962,27344,358751,185561,308491,96465,313504,156833,138417,142080,304612,123881,310325,482912,167746,285335,402466,252995,14725,219485,166086,143011,271543,220189,449500,121440,382796,10597,65265,323649,2576,323506,187536,390688,92378,402788,167714,58944,294825,183765,355631,217759,258056,400440,359722,250389,260292,353777,24684,153164,235954,232438,391522,285893,422707,345047,336954,370012,106103,402223,399904,263031,415146,136284,234524,7815,376536,374862,232327,411604};
            generatingvectors[1072187]={1,410484,298878,451876,506612,235277,34002,41778,303570,219276,244124,170428,322068,520517,480596,103775,349439,226652,13553,233788,217777,282230,20381,196262,498185,424948,166900,212116,157463,325034,441861,172160,53038,10304,474888,492563,289576,378769,391829,150648,491524,476100,376299,348894,286410,418789,326085,352094,276286,447611,177346,49991,259694,292823,165087,52471,271489,147060,271759,189517,101619,200324,336799,527316,182325,389399,457712,421344,39092,341716,130213,60535,378145,167619,399534,110307,467835,216676,31879,327422,375853,439815,179968,153382,365636,487842,302064,352487,292004,368142,137652,483960,243265,505144,490837,40728,139402,138902,199426,63685};
            generatingvectors[1179403]={1,519571,319093,151286,407304,389165,256337,489745,516880,133775,35311,317209,506277,75986,194265,86336,453380,28320,390907,93763,65882,75468,356708,410522,12722,524135,341833,150132,55363,381586,461201,64933,355509,366493,453746,157825,62650,99678,344673,197491,297819,366985,7755,119643,191010,262781,154733,372589,12957,51202,396006,441753,340621,69385,131888,161533,146076,1287,91354,283089,352139,33525,24129,203695,374488,55292,79854,452309,531751,229428,224790,146247,275168,415263,177371,527719,439464,136144,502512,463290,587133,260934,219449,53310,499775,509976,132300,582819,495629,357616,481808,231698,181539,334056,336293,311213,421022,388144,315069,159046};
            generatingvectors[1297337]={1,543372,403797,609097,626657,186198,275802,452730,246775,267283,36610,338846,594379,137051,411564,424117,420192,476370,243554,111517,93659,84808,35081,249968,641191,95783,244013,386928,177167,236420,4962,210196,395579,192769,312825,229884,603097,307759,349206,369538,419808,426679,552468,297506,548608,56441,31612,461526,643502,78902,41102,325050,224055,477305,123084,520086,462603,200374,6151,42655,44085,88871,472614,143146,313914,89458,633267,106855,415669,485164,455368,217465,390574,34873,124165,96112,622196,634145,608881,101911,160274,469533,296938,329653,155284,544437,165648,457107,99801,424719,596014,390275,174648,140862,16670,289455,237991,346891,12291,240520};
            generatingvectors[1427089]={1,546067,401730,598008,567323,492221,411122,661391,314683,689785,554343,171073,156008,454945,124900,82203,295201,501969,266734,132614,464800,650065,310275,430286,396609,607019,316826,277871,618013,541726,33418,630579,101785,581095,450397,190333,355763,364463,458890,327515,14927,390511,678182,229161,306768,212018,389726,625071,180606,202378,602060,341418,410261,629634,32707,604052,337776,407893,428837,41301,98679,190417,481646,155249,86945,607325,250860,155519,510150,146516,541906,631176,91379,681241,424439,159382,385226,112137,393929,487839,405599,21671,122054,85828,445689,664950,410727,529384,428219,146408,493778,83930,333252,237879,355424,182350,580331,688408,206055,627036};
            generatingvectors[1569781]={1,434434,327582,249381,384243,716016,759441,90863,728366,227793,753838,165162,334940,63704,777391,100566,196025,481589,648211,704331,114893,129446,374283,486349,281274,18626,491866,325181,146904,23448,351601,74060,490922,136022,590531,107163,262316,600556,374797,167292,197965,696261,433207,605049,737063,443156,258961,279237,503538,779464,144641,460363,8803,380379,500431,228816,97713,460725,11774,413943,472532,351105,308283,423318,474500,153188,113491,660405,223550,650490,608450,424744,644596,84317,352590,174450,427359,712960,467154,30042,120864,174302,133366,463898,769655,240993,310667,692743,2327,648614,409260,271431,22390,186509,57942,643544,170775,213215,516685,89020};
            generatingvectors[1726757]={1,669335,734324,593797,222277,155429,400435,492534,714256,147298,80949,119640,586566,17230,426894,322599,829254,126060,452400,605831,461467,758656,361039,624907,217613,39075,517837,95472,124400,206250,659133,294201,672148,608833,858781,626800,515531,23252,62370,274873,295334,92817,859951,678760,526295,832763,494857,390123,857526,725333,377191,264177,838565,686304,384778,321353,742020,839446,565497,738387,535570,570808,199765,81567,862211,470570,228801,456988,701500,588951,468381,500327,153099,194020,177220,404249,149879,785502,707036,697444,523700,10209,567676,531882,494984,87538,822348,466714,228420,613667,490811,709764,42406,734947,861515,445496,845749,698823,626558,824591};
            generatingvectors[1899437]={1,686669,415887,567267,704393,169303,462417,101217,191750,93901,20400,120233,562724,344239,345829,511681,550329,630213,424803,744807,321534,271662,917698,859924,824057,766840,908115,512491,541310,174100,809026,687370,858111,873348,419222,623756,715190,759488,281602,450453,909662,378487,261552,402044,736244,929806,78323,546597,508655,194032,698103,144854,730065,19738,580832,904005,938038,534133,484079,158668,107044,681363,280948,397809,721735,847639,877151,5913,854234,333734,7953,352534,618553,797001,388974,517723,629239,104739,121123,821644,109034,723788,847396,916180,253581,818559,627026,915014,167956,261614,666446,486516,855845,949464,750786,768583,198379,502033,569822,298208};
            generatingvectors[2089379]={1,807731,735419,858357,551129,257803,753211,52764,763034,296585,607688,213252,684574,828071,8355,416560,960537,40843,749989,477988,605527,67933,633626,933878,279991,377526,564603,978706,514277,485620,545201,770448,14909,126877,269735,719282,529331,209564,968634,838967,885568,191213,124096,302510,853996,480388,348409,191029,926216,679165,951437,207043,228680,996415,602190,1030389,544874,552710,189841,459126,209253,133912,166402,838467,1011758,592303,99707,178823,711466,94756,654307,774534,329938,906920,53768,296470,248646,550513,628474,953068,235365,84375,881475,1029426,583921,488500,636891,431165,113467,460670,691050,519702,24894,770630,433645,496778,615037,964939,1020249,359886};
            generatingvectors[2298311]={1,877283,974174,805527,178426,938256,479221,247898,1102963,95243,830469,896761,1125956,768647,451570,633146,1116372,654891,583590,749868,693010,382758,52480,898273,285216,196456,1052514,582382,472514,212627,64608,655146,860806,77901,153415,663115,617777,770013,1036387,10868,1072788,139807,1002181,736520,1067213,807710,588843,370623,134216,676541,973275,976388,868506,461093,196135,348111,911341,2603,258576,714476,801867,444527,781877,441142,312230,1097846,726423,46789,1039794,1040971,438977,508970,481070,152097,1081892,819743,135041,276191,978549,233097,11099,17983,516288,97661,666882,457594,206509,358254,1042792,1041672,895322,459452,563946,381308,558566,63452,525555,1147260,670603,994218};
            generatingvectors[2528147]={1,926051,1057418,1132363,325583,764520,1017199,948670,664141,540765,227158,349232,627785,154484,429050,910355,1021934,1031686,8807,410424,1070642,576554,335603,5631,1104922,962510,1236100,237386,1161410,1181158,88444,327226,489502,987978,16696,367752,794924,1013184,326129,1181555,566623,546378,677959,390361,788712,1076172,178135,694183,1072359,1172307,787799,903332,379475,732017,241392,435141,886713,953515,907261,1244062,513449,985653,950369,852990,326824,545275,328113,249587,579072,1012904,960791,414758,572805,698908,127417,1430,528233,957208,304689,1261180,103077,774223,2018,58694,872363,836749,1135733,550921,778616,718522,139821,1223024,378303,850572,1065979,757175,355386,977606,198833,861902};
            generatingvectors[2780951]={1,1056494,975582,1094842,417955,1199150,322265,1103564,1370795,59522,1380070,1334806,159876,545787,219798,592045,1014098,1134261,911186,277473,193311,549567,672273,750101,999238,1118021,1139880,498666,268800,1217265,541489,1196318,1204622,839883,1051778,360199,527970,1232393,1061578,302255,162880,233950,896961,871939,279613,524336,743658,228771,567922,683150,242307,482373,713992,1194135,1054708,1317060,1303328,1027455,1173281,753426,1260685,1214482,330036,1012829,339575,652982,242034,1060823,1198460,1250650,1122103,1324485,250208,1306296,309547,188174,264628,930532,611316,456416,956528,204365,1214305,101031,404166,99837,1074058,1312623,349903,1189169,345523,967707,1273100,1197068,195159,898146,1031195,1358625,847757,1359507};
            generatingvectors[3059047]={1,845491,1341331,638496,498477,548471,58619,1207277,311370,1180745,366377,949246,224352,1359941,1117981,814123,1427857,1190101,640652,1515073,1068946,826947,1344360,639446,896486,1103991,484572,507131,1511083,37419,589346,1505370,1202791,908784,1395223,1000532,1169341,1096115,1164224,170165,847968,255162,755393,1350376,1352868,120662,583891,1222583,1480944,672779,466785,239516,980684,920235,876914,662680,34141,707178,744945,1394310,827108,1483267,219579,298450,206134,97614,791066,691196,730949,378543,1222887,164628,157523,748163,104725,912022,685228,422447,85302,1119774,1375009,1420640,287520,1173339,1031558,109295,204434,1073895,1475964,540714,398602,157775,506171,262761,156744,1457509,19222,1273762,627124,113310};
            generatingvectors[3364951]={1,992620,1167926,632784,465048,888467,722404,1389613,882265,1425624,1180398,685696,1409313,964656,148305,75378,849070,497060,765889,1355990,1673699,717198,758531,998240,481717,651073,372472,356924,1407952,1333839,1635084,589221,1475691,1390165,213380,1193533,117913,1558099,549637,1604427,919140,956181,515042,397974,1376519,357990,146103,1312764,1669293,599822,1177250,227288,1424809,1093803,920502,694354,57597,1208376,500746,237229,1557278,992153,991,592850,259255,923326,529989,1620565,1338716,1255514,677657,266343,143597,338161,1084614,662041,1591499,100161,1163647,599951,229283,1110040,87211,670320,367192,573564,55665,567317,1142899,18858,79012,1139835,766368,66806,1507581,678196,1226410,510326,1004092,622861};
            generatingvectors[3701471]={1,1095045,1553899,1724445,1374455,786242,449772,200996,1019900,663970,866905,878744,1027316,1014736,264023,1501883,1581261,684152,1006063,776149,1189629,596338,533751,48561,1270557,978491,492328,847286,150697,1326873,1537042,1001352,1453852,228647,1469966,1476053,1438467,1529194,771499,415555,483597,1022609,667867,604617,1561147,1599234,875080,1609723,1152257,1456651,290799,1420035,304506,1569761,1188830,1191645,567197,22299,566201,255117,1418978,718411,1305707,530449,874609,1752200,1207484,346631,143436,373262,1001972,1065353,1632265,839707,1799565,1191150,1058758,313438,840613,498311,641178,239024,793875,1316323,1387341,1368987,332241,1165895,1753674,841351,664294,197938,1703447,11228,1420377,102543,707757,582667,842186,979799};
            generatingvectors[4071589]={1,1506285,723200,1549904,860038,883415,1215466,1299515,1613010,2022076,748627,663621,434573,943366,456324,1802963,1955441,1761198,513239,279152,369161,1434080,402472,231115,1696349,282244,425725,1596430,306327,1685700,90259,1985333,580051,1141155,1248177,1883806,1436578,1063977,815979,958573,803844,1374448,1569313,1411010,1076237,573122,658779,636605,1292978,216545,1936498,1543242,200779,782138,1136125,1099524,1637486,1574239,1796501,1804802,1640801,398108,1862254,1962954,1487837,1510936,253805,1093720,608575,661736,1855358,332158,989651,1886292,1602429,1883989,195111,1177477,673831,1517680,1546017,775593,1150843,1066998,819476,1746136,1921977,965683,1692817,621458,25486,1772637,908919,1094560,167133,299285,576936,1549185,95006,784719};
            generatingvectors[4478777]={1,1639710,1902601,1025044,578986,420676,1555336,1703070,79722,405345,721535,195805,1731163,1456977,279711,2085169,633048,212799,682546,991875,1758392,860550,1161273,1550125,668520,1685682,160229,1571436,841738,417544,1570364,382578,1021593,1735833,750481,377393,2084496,1263674,1864679,1717691,943169,314960,1754799,1974653,952154,2049302,273821,783633,751324,1580552,696691,777062,219545,1197385,822949,591907,2032323,281714,620590,1714427,1727585,1766471,377981,2182962,1204021,1531202,1143822,2097872,1139192,74561,2165760,2142932,1810053,873274,2091672,633534,1121909,993020,67968,444085,1258468,1821636,1821005,470101,1103689,1141083,604723,1612816,477119,1747026,343645,1448085,308162,2071394,981570,2159961,1597816,948457,1908329,2038905};
            generatingvectors[4926629]={1,2022206,1929938,2310795,487927,901663,438140,1242350,2176919,2052324,691022,1692067,2254841,215469,1582783,1749506,189131,1987615,669158,2395780,139186,900240,1002617,1739038,905910,1612788,674731,1535719,75797,2270426,1774677,1067352,446136,584390,1419432,757771,1095219,1119585,1059237,1061613,2143548,2068262,1186893,981780,63204,1890944,907798,181894,1936495,338946,2082766,1895178,223604,2402451,1484552,135534,832641,1151648,1437361,2439699,1848940,2369843,109679,876089,2122744,1013851,1119916,1879267,9252,2186477,559765,174865,919926,807417,1274837,341359,934703,2377291,107897,1307348,2093183,994097,2428252,1597678,317235,1357431,1760823,2285344,841897,1152761,1767557,239794,1009512,2340333,913682,515849,1104833,357632,1247050,1293119};
            generatingvectors[5419291]={1,2094979,2244538,2014676,1540387,757273,2121065,2086602,548870,1265547,2346608,88717,2454878,339194,245344,1723001,2017985,1085082,1308638,192844,703863,2041258,808186,565310,1842503,1485144,864636,830556,2504062,196311,2507231,1191694,1057763,2282479,1867082,1586074,1803230,553442,2600783,755565,1445921,1269065,2571159,473397,2616451,438455,2516009,1349509,592078,941101,568588,1252320,2191204,306512,791199,1287304,1881112,756648,2663032,2605683,606231,2283898,813499,1015909,1365687,2648358,243977,1370100,1938965,2058794,2330060,2270244,2401264,973619,601898,1249961,1007879,602644,1698709,670938,634004,27516,2390709,1033152,2645681,2169485,957825,1423103,2439143,1480039,2149050,1741361,1633121,539889,1712407,1333471,1436208,423685,147702,275224};
            generatingvectors[5961217]={1,2463643,2210464,949289,578486,2761084,379424,1717565,130224,1695736,1355207,170164,1973231,2566607,1793556,1877408,912485,364462,947032,2745186,820679,2671141,224887,93894,1681730,2431249,1402340,1621322,156307,2511798,1714042,1109551,1291119,695189,1102650,2599922,2545156,1874052,880086,2881305,2099792,309187,1464867,1189046,2122993,2396238,1968741,1576188,1462529,405012,162338,311407,191352,525162,306387,1360476,927650,215490,2804151,193106,1971516,751720,2194248,1953062,1772893,2546709,608122,623762,1061869,2828768,1511190,735599,2202972,1150112,1388467,2834521,2252663,2074470,568305,2314092,2680831,2387143,609457,662872,2067675,2574857,1877654,1538802,2372413,2411433,442669,2147756,1032248,267354,2467923,968271,1931958,498527,1714981,1746039};
            generatingvectors[6557333]={1,1828212,1252634,710617,3074263,2347881,2065166,2262378,2045958,2784894,2127821,1296015,390133,963454,2881580,2988023,735017,185510,2270835,68171,1264654,1730554,114280,559151,975199,421887,32643,27631,2875324,1532732,587500,3141919,1399407,2825338,1559169,1083008,440642,622966,2256759,185979,2071311,2077728,2062473,69947,1694870,3269659,3167978,261788,2326041,899050,2178512,2826143,1910389,66364,573935,826507,2786790,2294876,2802843,437649,875974,1018736,3103425,16089,253412,543357,863885,1336269,1145271,415535,81735,2211606,2404421,1814788,2609120,458123,3103975,3150886,171631,1910930,1254400,3241871,1816073,718415,1932085,911633,1672807,2862606,989138,2659462,451850,2727977,3111926,3216708,1529678,560360,2647447,2585383,1417220,1376602};
            generatingvectors[7213069]={1,2640928,2106983,3399952,1643646,1060162,985126,3016135,86050,2179273,2708937,2821310,1722479,97575,2449818,1676575,1953659,1987752,2041662,1735988,936753,170940,3157409,2073471,2587930,1120291,446215,594536,3324496,11995,1813329,2291333,2572222,3160004,2154287,3466323,2072157,1727006,871573,1129211,1177293,149173,3050072,2458062,2199073,1726294,3281220,3286212,1183764,2298517,2083545,1592305,1527831,859161,1083764,2637983,2295762,1004857,3148463,1022231,626709,1369369,1133034,2610702,1564962,1019223,1234738,1269322,1913334,539096,1420832,1963383,164737,1761543,2333871,841925,3574940,3276348,3039187,699700,317299,231309,1429131,1574802,526521,1366341,1820940,885694,1385732,2833795,915926,2697664,26948,1303998,1300769,1648732,1903231,2806620,647832,2647190};
            generatingvectors[7934383]={1,3030685,1643900,2422625,3492292,2152956,2705390,2928149,3219204,3517631,1052564,3412745,318840,2605951,2733211,3935908,3313161,3466769,2438604,723674,1986595,2264533,2794782,3092981,2251466,1780323,3471622,2518403,701324,1603176,2145169,3408902,3427063,1373356,1934951,2358369,2789790,273187,1884397,2524096,766337,3774437,2868116,3751273,921547,3446296,3953257,1875078,3215511,3248600,2367207,3214115,196392,1944420,1765584,1640880,1088018,1566662,2559213,1599355,2953256,384719,1972699,679648,3941213,1042036,195635,3913985,652320,56821,314770,3889998,670092,3508507,420614,2127320,895995,914167,575331,1347753,2671486,3038218,2424173,114117,3664837,3255371,1655468,1029087,1312416,3199432,3663488,3816560,2965536,438828,979949,2336263,888543,1768506,1439858,1177657};
            generatingvectors[8727857]={1,3340892,2551955,2747068,1866142,2512036,523044,1602327,4332958,3687899,2303211,2340268,1937003,2405297,639634,551915,3535135,2116118,1542590,4313421,3550299,1304644,401225,490306,880681,1644089,4335317,1214083,1282179,1274744,3544610,758323,853915,2950410,3685316,2100006,3686851,1852370,3714148,2501486,322119,1119586,505748,685756,3944203,3932926,2673905,2577880,719406,1140162,2169798,489215,1461568,797085,2773439,1288755,3497720,2791018,1346229,2763415,2392378,3465304,3496615,3134500,3527617,200880,4287702,1767491,352525,333245,3565708,1322184,1184943,4211691,1324034,4336226,2653082,745423,2538857,2309146,2380318,2683572,1531676,1618502,2663666,2991340,241977,1411333,1741900,610281,3232359,1943516,4316888,2047495,2387072,238141,2481516,2656280,496419,87033};
            generatingvectors[9600599]={1,2811835,3568964,879235,1222167,2150269,1349231,4081243,1732900,4706891,1620321,1870795,3993785,1595692,1852273,1248766,117572,3022372,2343236,489956,1276692,1813904,3996830,1546847,3638131,4036649,3935124,2659072,2110198,4792295,2752319,2466294,2088488,32435,1408115,3307284,3894923,336033,4678051,581219,3257626,4012953,1566082,1540774,1861361,3276447,2441629,4159973,1016760,2152183,3532974,1821624,1855105,1347493,2152760,1951899,2281340,3492182,3366746,2234227,381328,2658887,608059,1734574,3479744,1455894,2051581,3168168,2652154,2583578,675578,337648,2919864,664840,4149085,259105,4445479,660526,1899965,2719528,119151,2384625,4307695,2351543,433920,662896,3455328,1575073,813140,3093133,4419237,2317357,3161460,4137687,2309362,1882686,1827489,3130898,2395059,791020};
            generatingvectors[10560653]={1,4433656,3997810,975466,2507010,1912833,2448047,4473771,171099,2291232,1975834,2092962,4448240,948620,4574116,1349110,494945,689167,4193220,1944759,469841,2413942,3541398,5020509,1250938,2028792,4532090,1450218,4140305,993025,3007010,2264831,1305047,2593087,2374592,3768785,2814034,1224439,4981633,4056561,1200983,2653730,2519483,4688970,4930021,500890,3512733,1993193,5251284,1135398,2675483,4537962,2849005,2941490,4218237,308727,1713756,4947218,5138723,64248,1354881,2713096,206894,3095535,4041680,2378986,2680167,4585137,1063252,4076169,3287907,3885040,4457061,318966,486872,424912,2302721,3228740,2611697,225896,3819115,4659578,4040803,3818603,176381,4488902,1131429,5152411,3316975,2506860,3077501,4291451,1497168,3340489,757278,442582,2728918,2987852,3315684,2357346};
            generatingvectors[11616721]={1,3397023,5199573,5302266,5585339,5458202,1830574,1077491,1483613,3904482,1545389,369865,2415088,2275691,1897452,3722274,3765540,1629532,5406942,4345052,1743391,5636096,577404,5146084,4424953,3568253,2789476,5506240,2376181,2979515,168958,882065,841932,4118181,163610,5554081,3703993,1060858,2699876,3986861,2348285,272984,1413230,5745491,1571449,3269155,5664145,2271303,1369274,4290463,3452415,4883572,5113230,1066728,1302546,989010,271286,2743665,3906264,1622061,3297661,1601979,5241641,4458292,2602646,66852,869462,4619649,3575826,3349828,2484567,1178780,5629698,5263005,4907334,917593,2685561,4772795,510114,567952,1124649,4042863,840620,1254178,3888076,4016963,826545,636395,435104,3277526,4436201,1705097,4803277,3589087,2485424,2158478,495372,668040,2497034,2802924};
            generatingvectors[12778391]={1,3568997,2894329,5505658,2416170,4468481,3439422,4314535,311019,3124554,1591286,142361,5534431,6062840,669943,1930185,459496,5400627,4206578,6356485,5577699,1898100,1856037,3780469,1667434,3474926,268190,1313021,1396140,1794514,5718333,1160391,1081766,1995181,1023078,2553003,865910,4387306,1462988,691276,4296654,628627,4007687,4235815,1881756,104491,5758369,264561,5992121,4276720,298527,5878246,4524351,301069,827555,1790241,4307254,5338618,5760835,5857747,1327221,1215631,1514952,4693912,1179967,1926342,2522202,2765252,3644243,170772,2567866,2128970,1376775,1035679,188869,1313180,4781841,784484,1778805,779742,3485623,4689033,2389779,271904,6916,1345600,913298,4638719,5245491,5358138,3177595,339652,4756587,427365,2137191,6224898,3061227,5489692,4392063,1546293};
            generatingvectors[14056241]={1,3891084,2669825,4967508,5943388,1834345,6954888,2397255,2131264,3031290,5288397,5954295,5056426,1375992,4065132,6435777,5175554,6356022,6176851,2346060,6818479,4661703,6150488,4137414,3899926,1684879,5802587,6774532,319644,4183904,3193002,2236652,5351493,1881398,4586250,2828955,1019043,3314634,1777643,4664722,938632,3313814,6386971,4324313,2350910,2154484,927651,4298269,3786204,5664199,1576632,4888469,2224117,3081303,1382757,4459151,1345997,3862905,195813,5784416,4444952,6569625,6224707,1370296,2555078,3425186,4722881,26291,6105737,283997,780665,6346105,6429242,4567364,1800224,4424026,3271113,1461625,6104774,393053,5719071,1737697,207073,1330629,79626,6764640,1701767,2984648,4287369,4261359,2564367,4469134,374095,4353711,5782953,4226371,4918275,4743490,2351593,2761376};
            generatingvectors[15461861]={1,5916663,2785206,6671271,7039101,7569342,1551650,5780231,5302787,6457205,4143101,7683249,338005,1225058,7430603,3215498,7065632,6566620,610912,295286,1169493,1961240,1989995,4792972,2464598,3360888,6543705,1746192,215622,2212134,5817906,1579474,6413492,7722499,1657217,5552443,1723840,1887145,4802808,4324007,7563651,657378,2228103,7658571,715509,4287927,2477849,1504615,4059682,1095054,6860155,7195311,1549841,5535904,6034372,1245801,5980219,827362,6233296,2039938,5073759,5379897,2979172,4991711,4736207,132049,6283498,3376865,3157092,108762,1021502,4783552,3057940,1857293,4830368,6647098,5981749,2662325,6615901,608656,5521515,2359986,656932,463314,7273961,4421917,6514412,2367554,444172,6483375,926006,6134428,5502116,4221873,3669001,2809119,2922048,812826,6209756,4210967};
            generatingvectors[17008067]={1,7137894,3692315,4716893,2178460,1388112,2589332,2250607,2881047,6427811,8040791,8465407,8101732,2926680,1579392,7308475,1751104,894551,7586743,5636357,6598814,6467150,6525343,1587799,4154851,1521995,141659,2218766,7663671,6716508,5803603,3301328,1656367,1282667,5755529,2272016,6983949,6346662,581232,6239682,6661766,6125700,2991674,5298363,1014695,3427823,8416785,7143326,2946996,5640103,5729331,5886285,6447343,2303367,4327605,6903247,3688850,1563939,6793001,1899588,3714303,564517,7979873,2189485,3623361,392026,5058955,6885578,7830586,315299,5201561,8331566,2077844,6540752,7103672,6878412,7503300,2785144,85993,5678850,4297196,951392,3590977,1474037,8410572,845028,1758310,6355761,6276579,4706614,5020661,4507899,5410393,3175752,6878733,8286441,6478153,420734,3476614,1204008};
            generatingvectors[18708839]={1,5177807,8085631,6543608,7706728,1960246,6409782,4841532,4256542,5020768,822844,800764,8728475,1779422,106183,1562670,1141955,815341,1935729,1113281,2153696,8150006,1824468,5810488,8013682,7644193,2498782,8351250,1785687,9337255,7878520,9115424,615250,7346119,7871131,6373287,1140374,948405,8496491,8165933,8429775,6061535,5082516,4766622,5864549,3374899,3642868,354756,6403243,5060187,4233091,4992691,9029574,5375163,3369781,3949191,2064905,7403141,5029122,2730846,7818450,3848671,9207452,7672285,6955468,5576320,5322863,6646025,6967428,7020709,1513455,3546528,932235,1324070,4782624,8024993,3021434,9089448,4589610,5665745,863422,6524435,7136821,7894387,2176550,138715,1936560,7737727,4834378,2411645,5088971,651991,3513550,9187614,7771697,5240114,6976757,205186,2795820,7619289};
            generatingvectors[20579719]={1,9066949,4429647,8461836,5318823,6556516,2680484,1505194,2799249,1273010,6290980,393129,1598026,8404095,7366112,4843911,7141210,8609510,1615848,599732,8202441,4639140,9712271,1920621,6892967,4937467,4273723,9644653,9664286,131866,1971874,8391183,4513277,2962850,3463671,5540979,8321182,418080,5634995,1046882,225317,4684058,4348123,3972951,4656148,6650953,3643129,7117908,159513,10103087,9739379,3324731,8799009,552063,4763199,9434443,9777782,5350658,4675596,6123714,3980091,5188559,3110853,5908232,4966497,3853915,1664375,2453545,603218,3588848,9908464,4151566,7761420,3677136,7222228,8713108,2584029,7793638,1176143,3832926,6217511,4473972,4936184,6026760,5236155,1529656,9404172,9233690,9981365,5982350,6966660,7158737,241028,3494151,1686277,7695520,9368031,3208360,5737420,810425};
            generatingvectors[22637707]={1,8641615,6874463,9776537,1905288,2188647,7081400,1674658,8866038,4774985,6837910,6049221,8570124,4130580,7727102,10199043,4675350,2040603,2913559,11260827,9006072,4633203,929676,6925545,4007593,360412,5880933,9184276,5932445,6028330,877536,985010,523100,2645831,9178211,4499942,7175189,734687,3787140,950090,10478008,9583215,5683908,5262111,8805684,7000959,1656358,957626,3334709,6120723,3392105,5459401,11258555,3586231,4248097,3866365,10820202,350380,5127469,2239228,1501280,4999305,6276909,7391463,2121603,5493604,9479931,6993724,3897644,2016705,11043337,2391238,707454,3832301,3949210,1887988,5178760,9443920,791305,4649152,3952889,6352741,5951382,1675737,3013103,2031009,59260,2259075,10436707,10286425,10664544,9598843,6071767,7129228,5887842,10844157,1891625,2191126,6476123,6025832};
            generatingvectors[24901507]={1,6942646,10335048,8022314,4835473,9599986,6512354,6372927,9270851,1767343,9135269,6643169,9325047,1217950,160204,7790023,867048,7773673,4359149,8778884,4128525,40543,5424580,5632415,6361790,8308212,2467399,2531121,1329597,11101328,505324,9821564,10579701,1488468,11980042,2289990,11679908,3059538,2122597,9778465,6871839,5157634,1472654,714561,9291772,2808001,9075622,2094623,7754622,11592532,9921925,8047068,8851145,2841540,8280054,37944,1020778,7851181,2312489,10996689,361332,6504046,9956838,5085261,3408540,8762694,6857982,10044256,27236,7288125,8252445,9590539,6374837,7662674,2119453,2791203,1304783,1833228,3866848,264795,10473436,3720422,4765167,8058559,6565524,6923086,9287347,6943631,7879115,8343118,12404562,10764648,2493451,10477896,3454874,7332201,11272473,11398658,10418436,8214921};
            generatingvectors[27391613]={1,10403242,11223619,11927608,3026389,7348001,7094042,5424820,3695262,8856059,2349148,11826234,13268747,490794,4062996,3401808,2084154,1990755,6722454,11769733,12566322,7207705,9613707,828970,674963,4227985,1240920,6139855,2563010,3237738,7734596,2517947,13642594,8998742,583625,8310660,10080989,7957725,10404708,6502962,3049985,5995386,9307391,5575047,4624377,8422207,12214734,13457876,251624,10496668,8170380,2946058,12645191,12605596,6086822,1009032,1528406,9817570,3659483,4252821,13577371,13637348,12745727,8336339,3050535,12511972,4093875,8690392,9761726,9787961,10637126,3254728,13441056,1252748,501473,1698814,8625108,7780816,13524122,7314570,3893378,7313210,13316874,745157,7351764,11971883,12350841,2546499,2219733,9985596,2450917,6484194,13336602,5913732,407705,11515959,235720,10742751,12314606,5489989};
            generatingvectors[30130781]={1,11677788,12704278,2084242,13084982,7754140,6495065,13144736,6358170,418574,5964183,9030569,8434757,1452660,5381633,11170046,2781008,5870313,14894166,2564993,4960290,1910158,1373893,2947748,11198477,4788992,1594555,10872660,10330798,11599394,6437326,14640099,942556,11679451,1446941,2746936,9676702,5091296,13456509,5701049,14785481,8125575,2193027,14342540,6011514,13485221,3108063,6125649,9955534,171860,56814,12682844,325469,2057326,3394493,4425689,10442230,13324068,2683872,13810576,14744917,605311,5895101,13238300,11451260,4265748,2046180,7527131,8945354,6248863,9889743,4223932,6578343,8513985,10733522,8676002,11053451,11132807,12410083,6958013,6231764,14446957,5768849,12726062,1947289,9266171,69318,8255348,1815052,13918653,6371041,5923967,12751988,5173597,2051111,7416194,13222957,9022118,6537765,3808318};
            generatingvectors[33143849]={1,12799279,13622836,5355010,6552126,13958384,2113057,7040964,8077812,2654877,5016566,5112726,9619589,4987987,3619893,1169759,9909451,8441201,11802697,4250672,269738,15913390,13195328,11130576,4166374,14292945,9700477,7217297,7548906,11563651,3102217,13292040,159849,4531824,9807462,14151838,4106980,11255024,14652731,281618,12030220,14689385,12146612,2667814,4782714,9037254,4733393,5446498,8462506,7970441,11769917,14363482,8074234,9250186,5392259,10681984,1711022,6192366,13127923,3043429,2369484,7634840,13639823,15606019,12562953,2541718,6877861,8701883,14331707,11309858,6084784,7719703,781727,14278886,628117,6048910,3253473,954857,14058540,8692066,10349587,9304482,11071804,6845295,707929,2957434,15114069,12186719,13223076,14188527,8160058,714292,6396427,76543,12391386,4998677,9108529,5428087,2908479,11151630};
            generatingvectors[36458269]={1,13479692,15517872,6313237,13144471,16969280,9438726,12363917,1864592,6964280,1525489,9145210,4002754,16158029,14038024,12347223,1992223,10743655,10026953,9181359,6958934,12834270,17283117,11072823,15431100,2456033,1462300,10833018,14738840,13128070,17039110,9162183,9378894,17396935,4853893,1822870,1538772,11318835,3767620,3128770,5055324,2209725,17241811,5484467,15573677,7755070,14224116,4584731,16088161,4114497,7773060,990415,8714452,10125409,13684743,8632180,11672965,11464335,4332182,10480289,7569472,15166376,11863672,17484964,10331936,7598044,6699999,15143677,1516829,17070982,11184442,14461746,228604,7108818,14635202,8824078,12859667,9838750,12857969,10466596,3524291,15560054,7988942,5257370,1276315,13154864,17818897,443046,6989592,5838448,8454755,6993055,9800271,9405748,9401657,1047852,652083,16419358,11447063,15306749};
            generatingvectors[40104059]={1,16578933,9136081,15113323,6371958,11518076,8709306,17282522,19421394,3826451,19624057,8195276,3448383,16745518,763621,13628973,5444874,10943332,13414491,3383395,13850157,15831125,12431169,16748482,4523211,12758663,18284239,14343072,8806790,12120009,18113591,15408581,3476313,3270880,14799131,5069633,8037516,3331512,9271862,14420308,13250586,2618836,3702795,9235295,10336432,14027918,12771738,11331682,5658179,408989,12211522,17736256,392585,13879279,3028642,696277,15125237,13759710,13871315,2868882,17728556,2961857,4946943,19494936,14457483,10178167,7492204,13041360,14132252,18638470,9904948,12481563,8441105,15165404,11493497,1478529,15319474,18322239,135251,12721492,3152044,2062327,18156108,12781865,8608352,141616,3200652,3517752,6220683,8002816,2131292,5035457,9565380,15193509,1004143,650100,13067209,11661986,6981668,13515293};
            generatingvectors[44114459]={1,16755765,11843835,7799672,12540272,10741518,13516905,1178893,10821502,21064559,10296324,12245749,14239633,6316488,17343844,15854989,15149877,9332317,10945300,8992824,16539663,8892609,179052,19222365,14406523,8045600,21606258,21740627,14677814,1961164,3853384,8824193,4888116,7690402,18429166,11403356,6678093,3719141,20213505,9696293,21809408,10356383,5287705,10035898,17908300,5248148,19275848,4288228,16007198,6994394,9380076,17412425,17989930,20309683,20450728,16223679,14285566,15015537,20537408,10709575,1493343,2137941,16672999,11351863,21326079,10889495,16899431,8924889,811281,14374879,18506745,12469,11368785,2747616,3155425,12901229,7219259,3253071,6105824,7641979,71029,17274071,14499261,16964531,8722466,7072008,18759657,6539962,1932656,1381394,14980329,11710662,10767441,8124228,19150388,16897123,15494620,4415566,11490731,7430349};
            generatingvectors[48525913]={1,17922744,21172742,8589079,5413399,10141394,14404190,8165463,4531787,7490973,19169336,2460174,1225410,10160652,901192,13723315,14283422,11675467,24185551,17868025,13685569,15127533,21678674,16945510,11160049,14021392,23525134,11534739,21574710,18941478,23663225,24054019,15477676,11207687,6445936,13833537,11223181,13964014,12425438,21902119,8243694,7521984,23427758,20278842,8032488,3518268,7272990,23983825,8071147,21893834,14421913,16833332,16908318,11861300,22464150,22217599,7104418,19770089,21975898,13139647,23946666,7722450,112626,15725638,13559078,19863491,22845642,20500887,18382143,16773210,4561655,24228732,17875855,1594862,7976807,991016,2300193,1328088,17724119,22409162,16855153,21982656,22008215,20374995,20277502,19465239,15399516,18541339,18113115,8689829,22601983,18798262,1870740,2289368,11194001,1877447,23567937,8096805,212952,13719543};
            generatingvectors[53378527]={1,20690513,12630266,4949397,5605448,16102386,3243066,9440020,6958426,14168856,3993914,20055148,8139291,21000201,24930417,23878866,8577464,1682886,13996161,6717263,25502680,24389156,9265831,6338259,9522775,5518556,26430337,7179359,6040204,20441785,9355929,5416153,2703859,10066277,23528439,6492776,9667598,11755288,22776132,7880172,4646307,3801861,1762870,24685642,14081188,780678,24680010,24483937,17217644,6611587,18526245,24259554,12561287,16387687,20579614,20774552,19473543,8502078,2340204,1160665,4251237,15843321,23134135,24605772,4226785,24965207,23254969,8901045,17817239,8619535,2022950,9110947,5639027,22219770,19126527,4028430,20504927,23440807,20670265,4214924,23077111,2098915,8912048,7791075,10580349,13169805,23683945,17359021,9065573,12294301,5608915,7570531,26135684,17827956,13492113,59880,19137852,8319895,20838650,14559104};
            generatingvectors[58716341]={1,22399796,21156547,19126263,17908970,11348310,26485881,6651123,29004893,28105313,22520325,8497127,6137112,15570406,17458207,4225304,14755057,15319824,10318628,15043301,8908773,17371662,13403897,23764272,11628540,4892326,7229223,10424482,6590544,1018294,6862150,11147807,16365839,24161059,11887194,11380491,28721859,20716852,6829786,27965527,27070253,17874735,24989542,8153812,6386629,1987668,5891209,9033450,12656163,17303772,23805533,27852706,11281660,18736612,23884786,2077321,22986660,26637114,29012376,7629544,4997229,7683438,5373394,18243674,24314617,24055634,14347024,4052084,3100861,18397923,136462,7422653,14021491,15479612,28809212,23090083,6749487,16378457,22524132,27556778,950914,12765669,27674179,27728176,14082819,20764561,14159554,14390335,23779039,19570254,20598764,27330797,3290642,3628714,21329340,27461448,3641753,14911280,10671225,11523676};
            generatingvectors[64588003]={1,24650655,27949566,7901483,22363671,4754358,14574987,3370778,1736647,15859802,13829543,6108170,14107025,1478435,16755534,19622271,8856971,29684649,21431565,24235692,10983665,22620856,10155579,29641344,18376495,12175246,14387779,22966208,11937663,25934731,24242117,15521344,19213561,16326157,27976663,6269462,32080264,155678,22562722,4160910,18983988,13321165,7010396,26299678,2784426,26023116,2364392,26666101,11150395,31700409,16445009,160866,7951910,21331634,22346900,11768116,21226132,17094645,21612569,27254989,7215218,15535102,1407084,4075737,6334023,11117323,7445934,28497179,13647428,27700406,25596754,147780,24554140,1964365,23971458,16214834,26752877,17463721,7474981,7706988,12770920,30484562,18900224,573913,12119314,31256558,8154464,670973,2606532,8469901,9101495,15997269,27915960,12320141,23363288,31199084,15509174,12249682,14648471,7269866};
            generatingvectors[71046779]={1,29426512,8245344,32030003,21609904,4883419,23374394,25413292,35111400,33705548,12925014,748412,11005002,19711989,9186359,24990289,18971066,8815569,21070753,16811516,22657188,34419262,3811195,19804574,6346249,26473602,24040388,8401008,31032341,18432612,21742198,26861627,1375926,2315331,231177,8173975,6756761,20503307,2131077,24213789,6378935,19355677,34007911,28969261,13489669,1743822,32936440,6903682,20685431,34085796,6229727,4343479,31685702,13145956,20102675,29193697,12029896,24784262,19326261,30553194,2299562,27635937,27623657,27524955,22009506,23854314,5802895,17774124,26421125,19357773,18192999,8203354,32144031,28327139,32499328,6957851,27861149,6546469,17415380,35053718,32353883,13974023,15948082,28623484,26192793,23290848,5591252,4503431,3444523,32150532,12593893,33960854,8054034,26902518,11109490,12719667,21655413,8022321,5239449,32460549};
            generatingvectors[78151457]={1,23108336,32788258,21357460,37999541,26665269,9452590,5937391,33368699,6021764,20921772,6980992,1021467,5372133,37423180,35033943,21166666,34227055,19093559,11359530,1227023,3112139,6092612,31201939,12557882,14070320,37994345,3324527,12279752,21344729,21270023,38954626,2738248,26173121,13511050,7042242,33400146,21022772,7690174,12312503,37669398,38725008,7076607,37408115,8520785,13439775,15026894,3225527,27871478,29639369,32818034,35226597,36781257,15240414,19460367,37644389,37085449,27535210,15524118,8760890,5336799,27415898,2004484,22337342,33045583,15350671,21725127,38809628,2314999,7926022,32045018,31927352,2999132,10081391,29321458,32684604,27012211,7154530,26739564,21733599,8618935,3185137,34902101,13662432,11838995,27601655,21752747,22209718,16231077,35810257,25466087,15134350,28321942,14120887,26997353,15136106,4847487,18639678,23785713,23519238};
            generatingvectors[85966597]={1,36091127,23830247,20618233,14078812,11431985,19390573,32365389,1537028,10166204,23331900,9527822,38743374,10671758,38179743,17828963,25822767,7964104,30114051,36278641,16839710,40458376,22442502,28450071,10400594,30018624,29009448,19793460,25596669,33781477,17333960,40546441,27872176,22060821,28217591,12942726,13574304,14025139,26935049,38475352,2815550,10012611,12639866,9117174,9263110,41263712,17518671,35417539,5928532,19268923,31383301,39654236,35665713,9979085,27118978,31616433,5909544,26502849,19979450,24041987,8227889,4799955,42533026,21618746,33180583,4577770,12905446,9727903,28774968,23376075,24555894,28517509,903634,19156316,37008115,15705889,20387834,31354118,31045494,13933014,36566148,28441621,18670504,37232070,33925007,1307342,18003494,2303963,34255802,20456651,36071661,16885465,11080003,20468491,31576318,21777832,10747023,7801981,16996146,16961906};
            generatingvectors[94563257]={1,39904777,42349334,28264604,25903384,38903155,6906953,7930359,2934498,9591184,14876740,2893702,5814246,34082053,42688142,18791697,11306501,44337723,8793666,26193319,34519987,11235018,16651032,22218534,15632814,13634143,35452370,18193236,15546273,36306162,8945552,4573999,6282514,13337458,29746517,17298508,16764538,23515118,24170648,9642633,24067257,19608813,22940534,10241718,12870191,23841967,9139322,21870343,7853675,43243383,42462062,27520439,27248907,42847567,14842339,20644133,37472202,24114224,33829855,1221579,35344925,37571932,19216195,24650813,26673411,31262377,27090559,34376402,8482155,39650543,11312993,12145610,15344523,22375271,26116112,23688895,46333400,44420130,9537560,39522902,10964059,9579552,21659808,31603907,27283252,22168169,23393157,23103525,19382288,13140423,24777610,13248076,35279517,9590198,6245440,27864723,2384365,7068809,2565496,45971195};
            generatingvectors[104019611]={1,39508188,30677879,33133318,30959245,40637279,35228206,8307376,36579559,38488549,21631708,45641706,50148830,34380998,37681209,45288612,3775892,47577856,1237635,38014355,8005054,51356630,23067796,15614964,35206497,51574666,1558253,30182659,5271738,21034019,49775821,11525364,16496860,33457496,45944552,20301464,35434138,33863006,2665843,26613849,34486725,41801091,27688813,30556279,28910581,29156722,20741880,20217621,13657322,8920890,14637195,19826461,20130882,23616555,33530525,29800087,35086725,49193493,3249779,43125291,42383305,30048838,49511882,911546,32676362,23952362,10415284,43640898,15078067,5018394,13310060,42432979,45423629,1999413,46106159,31053108,48728175,11465676,14376083,17784023,9228875,10713460,23979266,20271585,28421336,41437820,20609176,8852907,48766851,48564200,50985909,42327638,32807435,49664321,48851729,20433886,18183959,10055595,36941686,22327146};
            generatingvectors[114421537]={1,44260893,26596153,46838239,48874554,7259483,35767990,44852313,34251722,48252555,23579486,33611998,23022808,50071524,49003702,54774662,40321884,34405356,57066038,35584530,6758529,35408020,48184000,39560691,32260886,7051093,35065994,45274117,38477902,2751299,18337163,49660460,32112782,55507488,1816181,10737940,10837287,3444708,22613266,41069852,51176477,32843212,38843898,4009514,4928597,41041640,18384909,43573907,15927652,96705,47950086,41781076,28506027,8467964,10484247,54442578,19935485,32466839,39441999,7164031,15500305,24696301,45719859,15144432,49755832,4733104,29371321,45461,2863263,56371010,54775599,7950420,41426023,15334154,16685728,14938611,46382030,1943098,9973209,17749672,55543344,24703740,20685760,45969984,29467899,30335772,52404141,22088241,21172320,54330037,24130522,41851917,1587214,41108211,35540942,26522445,39911074,7973046,36996668,28932539};
            generatingvectors[125863729]={1,52152094,56539146,15228441,37313646,50845266,57068015,61417219,4859725,25802191,19820857,6385285,27545916,45799236,10742811,11622145,47880118,54768693,52000748,3203940,41217073,47632155,6608554,35805817,62322478,150751,61866750,58244216,38287712,32332996,60640487,11722914,5676549,54969596,4008377,59569991,34740812,14027995,12894840,31460298,37406002,53973410,51062870,11142439,61320152,47191413,14558829,707122,61029416,50887167,52600271,24275730,13667082,23389794,55688939,62501909,17241076,46915180,41988520,62811353,7684258,9054610,49733834,14783729,23984717,39442215,14613681,34896225,32860823,8603741,48654245,54471248,11690038,1433555,51040406,53152465,48922313,2961228,57637739,1682248,20030471,34837915,22500690,30333555,18756586,58329083,56022979,28529030,48468952,52144812,26110694,15336390,62070811,45823644,28839256,7015892,5323591,16955711,21558755,48289859};
            generatingvectors[138450077]={1,40549133,50919072,35544732,51299374,53501111,66501369,59798896,59072861,49211176,4827458,55274730,50118524,7264488,24654385,63350942,62514239,39040507,51353924,52077889,34967431,48247267,24051214,4320211,7742713,9751542,48962214,18962855,49187737,65404270,8086936,63621261,16458084,48474378,17107804,26149386,47727257,66826288,52527568,2515957,58759146,11740924,6293307,43890134,9494747,25884698,11625150,29055381,29467996,67930312,1866509,52264529,15453919,40397049,67881984,11429564,12672111,46883153,31314038,33635233,39704259,50124361,6630542,13384392,42183027,19118180,65148557,60051688,494476,29693131,57525480,22594575,16526463,27133194,15513745,27565217,57176154,1932701,28112070,57294411,40816178,33877150,55093435,13191560,34443939,38816053,58945849,34244083,63048017,61611939,44743765,36790967,34184016,48127466,3429378,34692716,43653434,32705557,54905583,2837138};
            generatingvectors[152295083]={1,62936497,55913687,65911936,71176893,2776851,24407447,72407375,67911672,9619484,11832787,56395280,42098523,1467095,18912192,1212638,4290935,3295229,20731922,55287403,51546423,51510781,43561384,17125449,72177117,51964378,39392979,69984204,10143907,25827531,65923225,15509538,65427452,59963659,43964232,61067597,56808116,19450789,32006316,58133716,48914799,33639372,22271105,65740894,70219408,13234103,34587515,57984625,33949798,69719368,56195572,71998395,67102275,19588072,26625981,72853055,61319902,12002698,75167207,50512720,30673458,22982705,26473186,31326272,25280389,14564852,68519544,48740053,68631107,70243005,54432986,35224352,55422020,31844458,8770990,61453128,20433811,45538668,38550867,11373173,61200170,29661649,54723362,20506258,3433238,12083630,1038206,1412412,40597457,47337678,68437490,2897170,43907316,30836515,69713590,18034305,25077649,23491583,61695432,54483225};
            generatingvectors[167524601]={1,61947680,38502197,39812844,59595373,53240402,78977342,46619166,68042566,52958220,19794178,72234915,35831884,42000229,24466048,5077018,50233083,40830504,31697735,36076574,4243144,18331393,33745470,30259805,42501096,425302,25995969,61710146,53649555,60926316,71052827,44900430,68404611,61050964,43544054,3732705,32475003,40364693,45534658,11721783,71552847,67260436,611640,11751338,52508962,42534776,45542674,41584679,44588203,44085741,44666564,42282736,71556154,37429727,28879491,17108581,65744298,9282075,55191467,58582235,30867218,39021155,49403134,68642565,64303054,67413918,80587196,30683811,71422495,69804039,25757439,25887576,17899851,36983272,75421253,21744988,83193730,65688978,54029025,44828892,41119858,10001109,46667435,66708596,9409685,14339664,71909622,1737031,28890700,42597287,16343501,17183146,65916162,39066851,47598310,19604756,43723623,15014636,27542670,4538605};
            generatingvectors[184277053]={1,70500799,77437348,54917980,22011523,50582326,64805099,10035848,86521380,62841509,43917104,46387017,11139899,91829161,85384689,8996961,79841116,12761940,62961243,1976286,9860264,41126490,3447304,11090359,71422424,34427490,65956055,30591248,48089575,52486766,34335138,42340432,66358102,36982787,19649996,26947332,90148180,18751149,55834829,21015562,17906662,23808463,11616622,87786317,30137267,20418331,82398803,32842907,61854080,61004694,76458235,44742275,83397707,43934467,44976559,55133124,61452301,17877005,61906641,71065289,71621863,14008538,7172558,52171491,49378283,68194538,69680381,46165535,46499974,8485473,83033085,2287502,1504856,66966839,91569672,71425976,56597340,5992876,64074758,63805438,48698649,9632997,32796266,13205216,24218801,32861104,40675488,26117740,80180493,31707136,45023857,11084610,77208964,39221289,29454235,65115614,25411206,35822098,53710793,2960441};
            generatingvectors[202704781]={1,83766937,26100604,58456766,91182822,48375604,29330081,17432448,68609601,6280099,83519555,14181575,42252176,41039249,49998921,84838677,64700595,2431355,9668919,76767441,91136679,18060585,14661474,53021802,28053584,62845694,62504732,2506861,62016217,2728423,67479390,29662895,41107232,36309883,46258549,74162877,98141174,2920918,6741539,73690542,39177775,98537505,95963094,27845971,68798809,7169696,36349257,86151867,20644445,12062355,59038023,68981888,59490836,52389946,84784378,32407123,88027836,17940504,69072013,18820001,41576617,98722620,86043282,10940792,38342647,33695755,96802763,15700623,40556726,28384910,26514627,84302957,42645116,82683378,23232463,48115131,96709355,53944298,44195910,27713859,24222384,58880162,29515791,19576345,55876587,81899662,6240580,83793824,82883679,86306765,3975540,1965642,96136191,33114101,35021115,97413552,66930096,68619331,41650229,81466869};
            generatingvectors[222975199]={1,82345163,28732093,17716118,52282666,101921477,68012400,71966705,17196879,100745937,58414361,87793701,4101222,11179058,79854696,27212570,47414309,57765133,106339973,27638088,58668091,123196,33898081,49502102,49133603,110874073,108166564,85238784,70583022,88514055,16084622,110466256,69746380,21500317,32773304,105629923,5782493,70086376,34239457,103789917,3745235,29891168,8851552,73031433,105845227,102838202,32327109,38485608,47480704,104868647,50149591,102214030,64037767,95304992,92264115,71046561,58012727,78376744,25037641,14594449,111155273,72031135,17725601,103131391,85721470,61873351,110683930,18519075,26888489,48105350,67952016,6979839,66171489,72868220,99480377,60745997,68934599,48954062,110074167,10028155,105051438,44312805,50948578,90966101,64899871,37075954,38902413,23239672,14992062,83011378,16357295,107234143,102121258,99636147,106878143,82557835,94995718,41398491,32377574,91484649};
            generatingvectors[245272721]={1,93684110,102761406,65581484,116518380,90711223,69143571,56649403,117028489,115588861,67554297,80812307,107882213,26252134,98789473,60159838,976372,21873315,29028382,4748339,16202506,115702956,858633,43256340,107479249,35987811,15208207,112614628,28336197,42158452,100641559,33227122,82132920,31799015,14151097,75529631,38782764,106901471,77536378,113520523,68469382,113132076,68132640,73767769,65070454,27627891,38673463,7932639,20278249,24280657,108563756,5027537,43928345,60184504,17038719,75928574,108412847,6936764,88349400,119860871,61041371,54142112,20424288,4503073,93217394,15122794,98715577,94305040,120178499,64736946,97660759,78487288,104932565,104787798,94088976,12498822,22102409,45010972,102938107,90052592,50906434,41686973,37773601,100057449,71765185,77424971,113637364,49801683,26946096,97934102,51766751,18860668,79569643,39968647,76960164,30178023,35704793,906629,27633304,99329486};
            generatingvectors[269799989]={1,111799261,72220899,127878270,85969912,94974205,113200563,8219641,92700415,115309270,90370976,97288027,59374082,74521842,55836535,29993212,82088334,29579038,61832751,109657264,3425576,87125288,17356130,67059478,89320087,59391365,70350160,8763279,13163888,35280111,62374472,128722774,48322550,64945146,19877811,7805784,10564143,60901994,92600704,52895559,110620836,20834581,134516287,72584448,108614385,113340605,72094923,127508098,43371237,14224580,109671373,39160875,83229024,73395229,18879631,99254895,35664921,37367467,44518757,3288470,121459975,134623297,108558583,46006119,53883309,75533587,34234994,70313669,55666718,123830807,110790506,58068298,15065734,71138860,101166091,97075132,119835030,124056614,32130359,33401297,39990205,33885859,132028532,133653669,103636674,108231415,53452171,63188103,18629805,123954262,125813018,53074456,79736177,126216144,94497771,90772841,58059329,88680913,118947243,32458514};
            generatingvectors[296780009]={1,110122803,80101225,67894282,52177070,129304725,62637592,108411702,135575588,92382671,81912819,88557378,1326746,84687204,58134090,98848773,18237697,35245055,72565967,67262880,77834332,41712522,5093519,4309373,77157226,73383875,142339584,115313189,127551753,127379930,33752523,109879518,23952998,826299,108202104,2356236,74029967,36898429,1935626,62913522,64910241,21734495,15469499,121307747,3311425,92019234,142917300,105806367,114058205,111434064,59947714,31968775,121236459,8134245,108185305,117093536,72777398,44374665,76378150,91648892,18950422,103948161,22154171,56934942,11051788,122914563,14941051,24267681,13380888,104419872,142639252,12842392,32441249,7299022,11600626,41038078,102009167,23784786,29018262,134948764,108301708,104550608,32744964,81748438,112201723,84101564,99339372,118952387,29567659,98835874,51884310,52954798,120710160,91994029,135804636,85770353,74242307,84589909,2776501,98259299};
            generatingvectors[326457997]={1,137972599,79470984,99083423,152558666,15337332,101652179,97548121,154691812,115660218,22143813,139175535,31320175,67349019,156258080,99500518,159896489,160971640,101766957,45598350,15905425,80115342,41184727,136111748,132334434,7587326,71590635,45871755,118463107,137177913,150914889,1835554,76080346,45155323,159985361,33366646,131312446,104692196,153459953,124758897,151503017,88219651,101907041,53732324,57835702,122128223,52109822,126032378,56961558,125679023,20656811,138680841,70938313,141311631,17402910,92793788,127775605,11930410,98653065,116532499,127179590,32676293,73699358,30255884,83294533,43334115,78810607,90347209,55323475,67871752,37126988,85474625,22774391,89988285,81634570,144871736,27591607,6953420,88215955,34721228,43744968,18636994,41663928,80119051,38870998,116414598,158838123,135372742,142473787,29439859,47016389,119668139,44191093,143757304,114423641,43225752,104383988,70410771,63193197,54253926};
            generatingvectors[359103791]={1,150758595,160864047,95400372,52940266,65019235,141123687,55072142,102569203,135312360,65718811,133515277,105362929,150149894,66829307,79516670,62456471,24048249,134172717,109987408,65429208,139455496,178344421,148575496,97712370,164281964,28326225,160473145,76040926,141651261,28026135,103726689,150975884,96852322,109715443,176426314,69341750,84119227,21950578,168612196,2215023,70067409,129706170,109293072,86833265,112410951,4902537,50513827,152998773,28550572,38089959,139198033,108070151,76317260,97618440,60056129,155120177,14228041,20764834,89417891,79593237,77302719,91698216,106164150,165841886,173198312,102461495,93039083,59830136,40956149,167822796,60665589,96024075,89134549,175731761,35866940,82411851,89761372,50089614,144192897,42981067,67031413,57628487,75016689,141459738,102380505,64444159,85163959,16749558,177547072,27371836,118438592,86780773,173603197,89927400,65741379,2879199,71444996,78235763,101944181};
            generatingvectors[395014229]={1,165835390,155350671,40623321,25021470,37337465,57508209,14886332,120782143,168986285,145381607,115842618,54054301,24777180,176743742,178279483,173789851,189835111,179263020,69591428,147872376,117033882,93013066,45274405,96834794,161625230,180966454,111325836,51889650,189507182,176335781,61444217,87843486,77542616,108840283,105124609,48616302,16407887,19466903,180723980,67321866,189324673,145753567,76807003,96923800,3617674,42359965,192342441,99213142,51444102,55825556,124492034,14695926,24792561,194795046,163053695,65319015,134596933,163725667,56250485,148710164,39508432,176611481,97281353,19633007,12640495,65917919,121658529,48890018,193810125,134861608,75233070,16288555,127831481,45829677,87074264,97688621,115031197,35518300,64801213,150474000,8217808,185971129,18122247,41777424,72142858,35787447,181666925,10009445,55149667,88313975,164833518,162686366,129127223,183838958,155488748,30729380,142659055,10192706,161049054};
            generatingvectors[434515607]={1,166351266,151614025,44805218,63708583,60202802,74093393,161100295,42929266,15804059,135941639,145172136,87027571,40504348,94340940,184695275,26030596,42709085,118766629,55965461,155395660,216669659,125465133,145604031,216644257,159615112,87563233,140740725,58985190,117610781,18053293,5618563,156992489,186332921,118493087,21016668,82043964,120831553,189528488,60533158,159877254,56822309,148393713,215835898,171966159,81747962,26083874,163550754,121550294,46997259,180238095,95364214,206432867,87679839,18913476,136184525,98680598,104611035,11077626,193442737,109000283,53814871,9929071,97248876,61983810,168617512,121993822,80119377,175198964,206120301,61171880,43987861,170728149,82162944,171450637,197989711,206471621,29471037,136698399,65519453,169087133,153175998,49128655,95726400,93415308,183719241,115566232,174012482,105561104,38232162,148606679,151690430,191201420,177588713,51447637,128209088,82133825,50435230,16978168,210636076};
            generatingvectors[477967141]={1,200157411,171466394,92207201,149906013,196722201,141952899,197573823,41069031,133526464,86040753,184016500,176666108,225392736,49859439,212836170,113753110,45723344,101981249,188165239,16404144,173420111,39278578,164541524,118596618,114618549,7451592,24160523,120728683,82534055,211041657,12400891,48681299,14428282,10383537,44189447,129814895,98865232,88768475,236234333,37651107,66994822,218441701,237783264,35108119,195681176,757954,179611928,119530633,4344788,139716007,95241027,11316082,208903914,106324451,136626106,191886589,98800621,22436221,204236792,115218282,108374711,35443336,234550860,40653536,207511151,74276718,95129279,164739399,667119,59929890,70577238,194833679,130422901,35975405,4684799,206119123,189311705,64427131,167468106,183542821,205843245,208831107,31565055,68947614,13077056,202380658,41565052,234018614,213161523,34511850,227435356,108491820,211542727,226020787,226985739,221073320,180850564,74467107,186723121};
            generatingvectors[525763853]={1,220725468,97085888,204067241,191608665,162075429,60109632,46911150,48677981,98402411,113228841,144542804,168307581,225920720,219551888,142260280,157387281,183568851,53335968,261542088,173346842,147433864,212178126,103373235,9227723,157894269,212894228,88671041,226107153,174916063,116098042,135205423,137901397,64660619,260774838,122452750,120230755,64838525,123604262,261937036,114458311,261490444,178191154,172587766,84013798,209660016,238027544,80382871,23767638,84436491,140709636,101943345,137865624,123866427,185628572,39602666,145261540,84539172,116322658,27048392,254639422,59304199,56300265,261801679,76780689,42603458,208618539,28203038,220681316,160521008,198214055,81194238,122558838,151161366,68856939,138788275,52141217,159340094,184328842,65741741,170219038,90631182,10004134,46707337,70405553,245724347,40027264,84199164,246501320,209621101,193976004,215555602,743103,129714059,126434578,122383574,146455177,95291167,54864280,243626377};
            generatingvectors[578340271]={1,223580344,255221567,61135149,237910307,215696204,24108990,72707703,226747657,287734524,136987827,105309942,185575712,37668057,25305112,51387834,20805516,161804701,218616316,198076701,104626653,271869080,194261251,62299883,123443335,242514931,179036597,26891333,205281505,1665388,232073755,222085870,179554529,177418531,191200421,29598755,62992424,116766940,168948900,249103367,182242294,213582850,265383895,192973212,168480816,91186995,4212353,162859454,2406449,133445817,206677087,189688611,100744337,154361244,248778799,43341449,279998018,45858171,62749018,215366294,97425435,284580729,19661897,254684612,224400488,132323638,120943482,263399561,119075906,87832587,163756809,262569003,159935881,181619204,93385358,196385271,194259684,286772202,289107855,48289838,274346344,196108555,173312264,281955969,74238431,21805481,149548532,285658387,130917277,10377439,2347653,280770315,231711591,151494569,60925907,176126773,27709046,284520201,276139633,134667151};
            generatingvectors[636174277]={1,233583642,264147005,171723434,136987674,227895839,33163918,62436115,196495378,255116372,282537586,13324266,147548699,66100394,93804788,24104168,251101283,298304799,250690044,62813806,255715170,288298483,255904373,263273365,195534968,8183365,103644541,306744483,138548317,108676414,105711171,85784098,154210390,244243070,206948944,40612228,87690643,315101662,311678909,38259601,106388199,122090374,278664935,67706114,278918402,144483203,168284924,130908946,39038337,281108170,96323117,50471424,118580917,147730407,92824573,53180039,305877712,312379496,262052172,111346088,148134627,309755537,142144438,36248881,315228454,229648665,220681711,173404252,261809820,38008220,117018154,313189576,16029802,111383592,94240336,40696271,242847540,58164721,224466955,119943390,260625900,168792115,100616384,203189184,300714627,38694617,12534861,176833666,162608057,65671861,129601409,83997353,254383761,143740761,223456646,103739607,273501126,232829282,190987668,302146270};
            generatingvectors[699791699]={1,289182067,149558242,209004870,225148822,90235026,189596948,119200095,252194363,93636917,39028021,321152999,241440208,8392949,201981504,115552195,45954029,334130969,313309767,78348929,114433898,134067042,151426667,312735150,11634415,29795926,226119318,218838475,338081291,169119524,116855513,3773431,120657694,343183385,178822991,171576199,286048836,140796453,125155401,84171746,311508434,153503451,98701069,68092700,83805822,118447744,340634990,56708087,33469928,35666770,89624779,41951962,104463865,191533874,24835987,187416652,342784815,65254272,287485837,153450530,37794182,80965195,117955610,173691527,121940999,62117481,182203738,30188191,248518207,209246157,215347498,88577405,291097544,288318196,228000461,106949023,100769971,37975579,184610492,70867370,265194610,126801139,290870526,329666188,122169042,26351753,33234396,317627249,152073140,266622206,60727612,66045284,87142639,148606742,134966614,71123889,232379154,233196439,346726224,24590142};
            generatingvectors[769770899]={1,213089604,314339689,99629012,252287069,253980247,240257823,149985253,6717673,257650494,287432443,98708807,354991005,30325101,71848601,79206898,94210615,166744887,275629296,228604224,363300724,174299313,241879002,245085660,179786765,214091373,2019751,209539848,179553980,208360705,30639296,88061725,312911521,138220544,336187925,129241702,77857601,133366564,201529109,345204378,327883143,31497991,281546823,206334016,191508308,164918564,362385552,140561394,169184083,341692873,123545561,139893708,340360281,25134134,308596548,236711134,191149325,271597725,113945291,187292875,323538050,69506582,371481325,72821341,218800033,50815912,203490480,32381435,181379847,241407656,343791513,3809652,317088702,104627273,288010741,44674037,179122035,382876448,162329644,291765788,381544361,48734648,202532957,38375015,215172467,78399652,210373275,97619403,238267331,250021108,338528344,80989749,348264831,10729372,368239350,161268596,195040739,137852786,159123242,326991864};
            generatingvectors[846747947]={1,350870958,220567109,309110445,412258224,207854662,49751542,75718830,384156095,242156755,28780682,396519422,21671003,259497204,409390478,77959402,148406614,170503946,193708562,235608769,283620492,353833100,54089321,29956417,370313225,170122761,6975804,44322450,36141931,165546179,403965808,128108418,23850300,171058471,383275956,127338694,140648174,288948157,405740922,97685719,116355301,410900004,216343883,175539252,147111446,166156314,173072466,419831795,113461135,370066129,371806506,318721397,326081976,405721753,404258775,297908161,152180228,15634647,273337486,6735982,131168077,403472788,34011948,365001906,419250516,136730422,46960750,128186883,402137310,410749322,22856688,164415384,339733963,374381461,285456638,112572635,115353156,149672167,410709538,392965956,272465934,17836250,90297169,60922074,323003954,100155569,406997228,67869260,172146142,179247024,148490721,302441799,125304481,260306297,371772966,222261850,84449727,131036755,379505296,192285863};
            generatingvectors[931422787]={1,353827333,284062076,435895381,149983806,245424487,21742684,317212988,451491767,85142648,448501108,319629934,413148323,152471882,27255078,255860609,389083995,333750357,257451574,247192110,429703892,259593395,2537753,164999808,266744147,180919638,348171039,104463975,178083092,123767621,328507941,413699393,184508290,341872204,211701050,311422677,417412269,451287765,1651002,114821851,193859882,377723350,286950330,370117346,253825368,292604616,404758375,453340817,424245046,340008818,293407735,464934119,40039920,134603873,36295952,189429510,70558533,441021664,350727526,138271474,161094713,49526973,117057925,271915174,370458223,267244106,196722307,80429500,394205775,196182360,435339152,308058557,440056112,401877298,30674510,267600526,132313810,328773640,101580404,348886769,293251113,268095629,446909334,379303934,404034228,58422927,96901022,291436651,85950259,252533711,17278554,254704187,174086397,219205409,448441056,347700013,71999362,176762655,205489652,277881104};
            generatingvectors[1024565021]={1,397147529,375038078,139165758,245032651,95458928,424241884,294719718,133423427,296722256,285990824,37594659,270120280,487944603,218354727,70167966,317284685,94701448,337356302,240138903,157415573,183419100,36292955,356862098,192639550,510580387,226018514,50418443,128044737,259773369,146668195,39841400,82910666,97483391,389008142,260814056,258834742,110523831,67931621,44131386,21618695,262616887,391490998,494402522,20943852,297023132,424984942,256620536,203690075,153178014,90605562,162866959,116971619,134608197,17195762,105360710,4130142,511042998,200230654,329180222,219267827,388993056,440516573,384018929,315708474,220024489,239630052,400202869,77991659,386745459,427768485,325341530,231862463,70873152,487503388,334402182,87921022,231742084,99960831,402698385,9404234,449616207,262695722,390026096,376114875,499438713,49546832,273665427,386169766,359636313,12894182,148121116,17557730,34030561,435178964,481989487,508283684,11603784,336990021,109226416};
            generatingvectors[1127021513]={1,431407286,201486647,234008540,471119775,480710433,320666353,120943383,254399571,314307910,138186816,352614454,174139792,291757978,177595608,408483158,538539158,251969451,395136606,207311509,119901567,499176643,50966039,517466417,285406440,358694676,123941213,112941719,176298884,293020885,268362557,557851423,124487352,32939233,498182530,310078280,182114392,103252654,323772609,281594468,170384892,396297814,241686872,18645191,306156073,245100395,391139307,549203345,164176630,510763573,309842725,475834424,289967214,98663679,343891840,340296341,309414450,345255695,393608239,379289571,343136297,261257070,46460420,444394124,190197128,362839724,270012672,508162023,454520427,340890461,258257618,110463138,259773940,549336640,296263857,95563480,330498685,342360489,127908040,420871947,261094422,319244355,97075509,131839095,252801684,10999707,215966259,514001001,369382853,281575075,153161996,503739329,96206709,195678754,2277466,464508400,312580557,446818136,551796850,533904603};
            generatingvectors[1239723671]={1,453165313,344765187,236657535,590043539,479772617,494033696,315938335,341317120,147968482,239329805,301951473,263399566,597369955,149929082,399278371,245137476,366312393,130992986,182850050,435104108,39583064,316520267,93133742,304343850,259733156,584763042,511030959,550336174,416113547,195090035,94186090,573686938,344517521,486306318,362402600,153336105,230444086,430251621,595420120,118577115,350726857,435590980,38867063,478921409,194576566,10792671,600460211,614318545,346934338,609465277,551244473,434066826,571709246,341363653,471505529,417953675,318130051,278635137,138876415,95712007,31516301,170030560,60075780,130627614,245133311,608434590,376831270,556889012,527401851,77786213,206048342,348097406,475331131,17035964,388388477,213946100,364033453,486295235,258094314,598662638,257222745,320471731,383255053,196800123,266247287,55840513,430852116,198546414,611484999,492063895,508041971,351339377,151262849,417218225,485170452,494218333,509164744,369809898,75574868};
            generatingvectors[1363696049]={1,532546482,205326264,623977315,353018806,283642084,473343852,214741848,301082673,191425505,324000265,497152484,639291789,467388015,197078724,342262347,348897601,239335398,44591292,449386960,522476650,541676998,379027144,462153709,367783905,421376367,499375029,569856140,299189876,318043788,116405354,122516591,352279767,67677419,526344623,143216870,465127576,214228168,503188620,385744818,387027476,445088837,389824821,462570752,325931332,251888954,646386571,3006673,194496724,666747911,39494956,446005763,607976218,36923825,516307396,546617751,411131835,551432369,165817225,380327419,568569381,195029423,59722657,119421183,448609162,5164702,408579242,554576512,603564936,294906226,231882913,87684521,610101941,63295258,37756263,233304143,605580924,477760427,263118323,624743315,27172322,322412088,291675744,514687174,568957877,173923641,397550756,526709972,420549959,427135468,189577883,341195574,521729153,426671263,601137524,19861764,416072677,38755337,156092803,270309965};
            generatingvectors[1500065653]={1,445152050,311193132,275684529,722769685,230384853,482279756,604790637,718498725,216840847,638911676,325528253,28822404,672373440,686446111,30339087,51177992,745281171,174349131,664160052,81927029,83615313,234403722,416036217,209500337,568006775,83551757,582863623,277823102,504173603,23266483,358907868,115727035,670986297,691495797,380687306,688368817,325452214,351385122,626551787,659510645,120141821,531203884,699789227,114900138,472848053,25533981,409375577,589141650,211477977,125223575,366832757,283375767,430777691,334828592,571982228,686969363,496243275,643164568,728055814,145930477,367910416,482543277,165944981,166071818,177599250,190400298,690590359,402431083,252066148,303742554,447842698,96371097,521066125,543733516,64236468,187876233,205785696,390280644,650169239,380087749,12995838,229018544,647710354,547645019,743232125,536881336,591481697,114707063,365148345,65824418,413557282,61242245,249815581,187709083,469207384,125669462,157267564,612996995,646478704};
            generatingvectors[1650072199]={1,627357855,730064373,484284357,656815925,691073993,497891490,486117798,788485553,426990475,230832880,317159873,182166088,422903305,497297307,773926834,709295337,285815216,768307860,499916774,138857159,283191183,220187164,215832320,254018855,641134087,28671731,243848522,247434944,300559455,428654010,213454768,728073937,455391973,586130998,774320938,201722485,770684676,363471102,165654003,794439222,193128023,419578805,1982649,723991295,447349211,92923860,324259388,659992849,650668668,480325532,316957438,76481252,414829915,467410042,274494288,92351030,63363309,538016980,487819595,758526884,227665621,126239047,568410545,469052400,126254327,475481388,41501023,704956384,7502854,364409654,13095077,343856776,785488951,376348630,32666002,822169125,363593058,191344411,556874691,653809407,462572350,285727369,426248450,706112658,313718149,70821941,54985795,639492547,529957065,458117251,161487791,47543546,663671958,310343602,151072638,467366285,551330564,776846432,319554389};
            generatingvectors[1815079421]={1,692426799,380140186,842893478,653499579,803769147,153416438,640182760,339988640,862245859,139803454,255741738,725607449,332442437,859693709,136246910,234064213,769607618,406736973,212943763,62409788,117176832,702635281,341414654,497133474,736666040,162077250,626088072,101937192,875292861,643231932,474858156,635992539,748448181,423428347,868736245,466694266,366767649,65039816,652667126,744819776,743376798,577276706,516841577,43410049,694554692,801934014,338111253,95111209,481513615,741435010,282883446,865502391,715239522,442037302,383661621,327100397,238566755,455770648,755747208,464501050,859660215,222998774,685922662,556642128,630315335,334473923,94702317,868291024,188370821,878761992,325873205,510987664,799464060,747895998,714475359,40839372,718238115,523396053,570745997,120620478,383177185,837745897,731766917,439161423,579594262,612157674,29210531,797063218,716235526,331848078,666181694,324992528,53850912,879084460,590859237,611264409,761157593,39864155,783540475};
            generatingvectors[1996587361]={1,738296812,589511649,106668437,729374048,948522830,471683709,451896876,826201231,978618770,269686254,61396354,154270332,529322286,287834673,519245380,959525685,93351830,190196121,615928705,383872519,867261353,921245530,392575456,47174967,716583629,989934616,65814333,845169993,378360913,158818377,498550281,706807089,147779999,79462943,487491804,935153549,945179220,981351995,859441592,185281746,952831765,105167667,735124439,892675712,727729928,29881118,942103364,81415104,383929389,227450244,332620347,352823927,416312549,766388464,644875416,717970696,270224749,455183632,784260870,749091128,740784796,916134153,407103435,356594700,90136595,351454286,685081457,842763617,288867145,834877188,933518818,763436417,413501482,129947223,732329923,529309927,399429682,841834737,641424819,419391311,593602825,850682478,179756529,279525819,257932720,990908387,120066569,224422790,713942922,832741985,30656547,868541979,163800221,460440031,753409442,859834032,152244385,361543234,275057607};
            generatingvectors[2147483647]={1,367499618,943314825,398462436,111463890,251546619,416692448,474655545,974979524,44547761,305844442,714785972,1039898094,662132929,280201912,497391858,178010392,537431569,556169982,246219650,449294340,548185113,817699382,310001399,614710621,399511030,882993569,1056870146,71368987,333222182,1020980125,764556999,58967594,753972837,556745317,506343439,140520151,308901573,948653264,457392340,2082238,558005132,1052473394,399198201,471837708,792442062,479579819,452840212,25344030,1054177295,65350767,16366279,655758979,101840147,751287227,644149708,913889773,749115983,927845151,700507179,405455845,60972350,341628273,844353996,933046539,923820109,783709236,482979181,603839889,49858221,100486019,769455829,769504197,722954331,974399700,800297415,746147283,984345387,918849071,1052383354,750599997,148589686,420091285,823797948,339997452,509374923,301812454,696799242,134896690,225747989,744320441,696764539,53183285,857228966,945665079,540559905,470168551,321475774,164004175,778243523}; // INT_MAX for int32

            return generatingvectors;

        }
    };
};
#endif
#ifndef QMC_GENERATINGVECTORS_CBCPT_DN2_6_H
#define QMC_GENERATINGVECTORS_CBCPT_DN2_6_H

#include <vector>
#include <map>

namespace integrators
{
    namespace generatingvectors
    {
        inline std::map<U,std::vector<U>> cbcpt_dn2_6()
        {

            // Vectors generated using Dirk Nuyens' fastrank1pt.m tool https://people.cs.kuleuven.be/~dirk.nuyens/fast-cbc
            // Settings:
            // s = 100
            // omega=inline('2*pi^2*(x.^2-x+1/6)')
            // gamma=1/s
            // beta=1

            // Used for integration in arXiv:1608.04798, arXiv:1604.06447, arXiv:1802.00349

            std::map<U,std::vector<U>> generatingvectors;

            generatingvectors[65521] = {1,18303,27193,16899,31463,13841};
            generatingvectors[131071] = {1,49763,21432,15971,52704,48065};
            generatingvectors[196597] = {1,72610,13914,40202,16516,29544};
            generatingvectors[262139] = {1,76811,28708,119567,126364,5581};
            generatingvectors[327673] = {1,125075,70759,81229,99364,145331};
            generatingvectors[393209] = {1,150061,176857,160143,163763,27779};
            generatingvectors[458747] = {1,169705,198529,128346,134850,173318};
            generatingvectors[524287] = {1,153309,134071,36812,159642,245846};
            generatingvectors[655357] = {1,253462,69526,294762,304980,238532};
            generatingvectors[786431] = {1,300187,232015,63830,343869,39791};
            generatingvectors[982981] = {1,360960,73766,194632,51680,293702};
            generatingvectors[1245169] = {1,368213,239319,593224,147860,546740};
            generatingvectors[1572853] = {1,459925,736430,70288,373919,634109};
            generatingvectors[1966079] = {1,826127,686058,926897,417836,183049};
            generatingvectors[2359267] = {1,696379,1060519,640757,812754,262923};
            generatingvectors[2949119] = {1,1090495,479029,595914,64689,895947};
            generatingvectors[3670013] = {1,1357095,1026979,857015,644825,1129417};
            generatingvectors[4587503] = {1,1742417,1399874,672080,1827715,1488353};
            generatingvectors[5767129] = {1,2210064,514295,1675989,137965,1611055};
            generatingvectors[7208951] = {1,3144587,1709091,872489,489266,2288306};
            generatingvectors[8933471] = {1,3453775,2110983,261972,2555740,1086124};
            generatingvectors[12506869] = {1,3663001,2298621,853317,2983823,5576578};
            generatingvectors[17509627] = {1,6426637,1486159,2528828,866597,1015123};
            generatingvectors[24513479] = {1,9363157,10935868,7904120,7202893,10833044};
            generatingvectors[34318871] = {1,14408021,2474791,14056163,13619371,8142161};
            generatingvectors[48046423] = {1,17766606,14390535,18752150,15489536,22204578};
            generatingvectors[67264993] = {1,26061396,18907982,30760850,28273663,360289};
            generatingvectors[94170997] = {1,34493943,45822183,33604771,17761662,27235450};
            generatingvectors[131839397] = {1,50467100,12217927,32766578,62069641,43610269};
            generatingvectors[184575163] = {1,70104887,26696463,66178896,33835785,44887749};
            generatingvectors[258405233] = {1,93229880,121240317,81359405,13132851,3566987};
            generatingvectors[361767331] = {1,151870194,126971921,92157910,95131599,23325957};
            generatingvectors[506474291] = {1,187358266,98361527,30873241,97613632,101655550};
            generatingvectors[658416589] = {1,241782268,273651328,292254504,308395171,128636084};
            generatingvectors[855941587] = {1,330894152,231467230,413890922,328644676,277687765};
            generatingvectors[1112724071] = {1,469999985,500083067,198639601,547255654,47506555};
            generatingvectors[1446541331] = {1,597774846,391859496,244464942,540060630,431672305};
            generatingvectors[1735849631] = {1,728493464,221397948,296306937,119145578,528265440};
            generatingvectors[2083019591] = {1,863176824,162064764,556032311,529803849,822974245};
            generatingvectors[2499623531] = {1,946563710,446148663,365195622,977705092,1214380109};

            return generatingvectors;

        }
    };
};
#endif
#ifndef QMC_GENERATINGVECTORS_CBCPT_CFFTW1_6_H
#define QMC_GENERATINGVECTORS_CBCPT_CFFTW1_6_H

#include <vector>
#include <map>

namespace integrators
{
    namespace generatingvectors
    {
        inline std::map<U,std::vector<U>> cbcpt_cfftw1_6()
        {

            // Vectors generated using custom CBC tool based on FFTW
            // Settings:
            // s = 100
            // omega=inline('2*pi^2*(x.^2-x+1/6)')
            // gamma=1/s
            // beta=1

            std::map<U,std::vector<U>> generatingvectors;

            generatingvectors[2500000001]={1,1056092002,604902782,1140518443,1168484358,678540231};
            generatingvectors[3010560001]={1,1265039176,710224583,570392900,246175051,776375237};
            generatingvectors[3527193601]={1,1477280710,1631679535,500900763,1337951012,990240443};
            generatingvectors[4046192641]={1,1545350222,1821675932,852726071,257150351,1540501786};
            generatingvectors[4515840001]={1,1895743314,1978110099,1051107732,1249084094,95135867};
            generatingvectors[5040947521]={1,2084035980,585653597,448523180,856444223,2389197079};
            generatingvectors[5505024001]={1,2280282288,503769990,2547746687,2668753240,2100976149};
            generatingvectors[6165626881]={1,2360163115,1727923807,3043833953,2316665784,2702804871};
            generatingvectors[6561000001]={1,1812543072,1410669934,1037177071,1156985284,3184493703};
            generatingvectors[7112448001]={1,2716396872,3010889702,2894020956,1485836748,1959799747};
            generatingvectors[7501410001]={1,2757177450,879378109,2460097563,3069036981,3181359993};
            generatingvectors[10088401751]={1,4180383535,3234338499,4352977744,1736039557,4095101376};
            generatingvectors[15173222401]={1,5795573206,4481927636,1112677318,2916664894,2804891062};

            return generatingvectors;

        }
    };
};
#endif
#ifndef QMC_GENERATINGVECTORS_CBCPT_CFFTW2_10_H
#define QMC_GENERATINGVECTORS_CBCPT_CFFTW2_10_H

#include <vector>
#include <map>

namespace integrators
{
    namespace generatingvectors
    {
        inline std::map<U,std::vector<U>> cbcpt_cfftw2_10()
        {

            // Vectors generated using custom CBC tool based on FFTW
            // Settings:
            // s = 10
            // omega=inline('2*pi^2*(x.^2-x+1/6)')
            // gamma=1/s
            // beta=1

            std::map<U,std::vector<U>> generatingvectors;
            generatingvectors[2147483659]={1,594332840,987962412,92850199,232950545,732223484,216871252,167983424,134396467,967335352};
            generatingvectors[3037000507]={1,1154488078,1364044675,403532140,1418532776,1211730398,1291040262,1460469075,8245344,1211134081};
            generatingvectors[4294967311]={1,1782332935,2038179138,312057978,558140820,1823542953,1209402799,2050709695,826107915,1042894327};
            generatingvectors[6074001001]={1,2515792464,2214346390,2112169118,1128189490,1484488734,1224278437,2351616022,1611384151,1114862459};
            generatingvectors[8589934609]={1,3321942526,1019713590,1739656802,1346968868,647002714,2524234491,3572784402,3780296281,1932783779};
            generatingvectors[12148002047]={1,4460034585,4192279828,956270804,4233805881,531494255,4314610224,5815452895,2216195588,991235505};
            generatingvectors[17179869209]={1,6509448386,4956350883,7081563299,3187015543,3057004275,1370423947,2853256536,6529660412,6947593721};
            generatingvectors[24296004011]={1,10063714232,11598913330,2292935963,1047062119,1472104940,758619429,10725375394,1245101568,3691451466};
            generatingvectors[34359738421]={1,13284241015,4747929435,5634134767,6572020438,10731165789,4350433762,14510070047,2620274659,4573585837};
            generatingvectors[48592008053]={1,18560138411,9547098320,3785452840,12633112625,1079370434,2348350358,8251590624,5328275748,4759811440};
            generatingvectors[68719476767]={1,26244027730,14007824995,31591716930,26842308147,29732316964,19894805679,17981118397,29347359790,4732146810};
            return generatingvectors;

        }
    };
};
#endif
#ifndef QMC_CORE_CUDA_COMPUTE_KERNEL_H
#define QMC_CORE_CUDA_COMPUTE_KERNEL_H

#ifdef __CUDACC__

namespace integrators
{
    namespace core
    {
        namespace cuda
        {
            // TODO (V2) - investigate if using restricted pointers improves performance
            template <U M, typename T, typename D, typename I>
            __global__
            void compute_kernel(const U work_offset, const U work_this_iteration, const U total_work_packages, const U* z, const D* d, T* r, const U d_r_size_over_m, const U n, const U m, I* func)
            {
                U i = blockIdx.x*blockDim.x + threadIdx.x;
                if (i < work_this_iteration)
                {
                    for (U k = 0; k < m; k++)
                    {
                        for (U offset = work_offset + i; offset < n; offset += total_work_packages)
                        {
                            D wgt = 1.;
                            D mynull = 0;
                            D x[M];

                            for (U sDim = 0; sDim < func->number_of_integration_variables; sDim++)
                            {
                                x[sDim] = modf(integrators::math::mul_mod<D, D>(offset, z[sDim], n) / n + d[k*func->number_of_integration_variables + sDim], &mynull);
                            }

                            T point = (*func)(x);
                            
                            r[k*d_r_size_over_m + i] += wgt*point;
                        }
                    }
                }
            };

            template <U M, typename T, typename D, typename I>
            __global__
            void generate_samples_kernel(const U work_offset, const U work_this_iteration, const U* z, const D* d, T* r, const U n, I* func)
            {
                U i = blockIdx.x*blockDim.x + threadIdx.x;
                if (i < work_this_iteration)
                {
                    D mynull = 0;
                    D x[M];

                    for (U sDim = 0; sDim < func->number_of_integration_variables; sDim++)
                    {
                        x[sDim] = modf(integrators::math::mul_mod<D, D>(work_offset + i, z[sDim], n) / n + d[sDim], &mynull);
                    }

                    r[i] = (*func)(x);
                }
            };
        };
    };
};

#endif
#endif
#ifndef QMC_CORE_CUDA_COMPUTE_H
#define QMC_CORE_CUDA_COMPUTE_H

#ifdef __CUDACC__
#include <memory> // unique_ptr
#include <cuda_runtime_api.h> // cudadDeviceSynchronize

#ifndef QMC_CORE_CUDA_DETAIL_CUDA_MEMORY_H
#define QMC_CORE_CUDA_DETAIL_CUDA_MEMORY_H

#ifdef __CUDACC__
#include <type_traits> // remove_const
#include <cuda_runtime_api.h>

#ifndef QMC_CORE_CUDA_DETAIL_CUDA_SAFE_CALL_H
#define QMC_CORE_CUDA_DETAIL_CUDA_SAFE_CALL_H

#ifdef __CUDACC__
#include <stdexcept>
#include <exception>
#include <string>

#include <cuda_runtime_api.h>

namespace integrators
{
    namespace core
    {
        namespace cuda
        {
            namespace detail
            {
                struct cuda_error : public std::runtime_error { using std::runtime_error::runtime_error; };

                inline void cuda_safe_call(cudaError_t error, const char *file, int line)
                {
                    if (error != cudaSuccess)
                        throw cuda_error(std::string(cudaGetErrorString(error)) + ": " + std::string(file) + " line " + std::to_string(line));
                };
            };
        };
    };
};

#endif
#endif

#define QMC_CORE_CUDA_SAFE_CALL(err) { integrators::core::cuda::detail::cuda_safe_call((err), __FILE__, __LINE__); }

namespace integrators
{
    namespace core
    {
        namespace cuda
        {
            namespace detail
            {
                template<typename Tin>
                class cuda_memory
                {
                public:
                    using T = typename std::remove_const<Tin>::type;
                private:
                    T* memory;
                public:
                    operator T*() { return memory; }
                    cuda_memory(std::size_t s) { QMC_CORE_CUDA_SAFE_CALL(cudaMalloc(&memory, s*sizeof(T))); };
                    cuda_memory(int value, std::size_t s) { 
                        QMC_CORE_CUDA_SAFE_CALL(cudaMalloc(&memory, s*sizeof(T)));
                        QMC_CORE_CUDA_SAFE_CALL(cudaMemset(memory, value, s*sizeof(T)));
                    };
                    ~cuda_memory() { QMC_CORE_CUDA_SAFE_CALL(cudaFree(memory)); }
                };
            };
        };
    };
};

#undef QMC_CORE_CUDA_SAFE_CALL

#endif
#endif
// (Included Above): #include "detail/cuda_safe_call.hpp"

#define QMC_CORE_CUDA_SAFE_CALL(err) { integrators::core::cuda::detail::cuda_safe_call((err), __FILE__, __LINE__); }

namespace integrators
{
    namespace core
    {
        namespace cuda
        {
            template <U M, typename I, typename T, typename D, typename Q>
            void compute(
                         const Q& qmc,
                         const U i, const U work_this_iteration, const U total_work_packages,
                         std::unique_ptr<integrators::core::cuda::detail::cuda_memory<U>>& d_z,
                         std::unique_ptr<integrators::core::cuda::detail::cuda_memory<D>>& d_d,
                         std::unique_ptr<integrators::core::cuda::detail::cuda_memory<T>>& d_r,
                         const U d_r_size_over_m, const U n, const U m,
                         std::unique_ptr<integrators::core::cuda::detail::cuda_memory<I>>& d_func,
                         const int device
                         )
            {
                if (qmc.verbosity > 1) qmc.logger << "- (" << device << ") computing work_package " << i << ", work_this_iteration " << work_this_iteration << ", total_work_packages " << total_work_packages << std::endl;

                if(qmc.verbosity > 2) qmc.logger << "- (" << device << ") launching gpu kernel<<<" << qmc.cudablocks << "," << qmc.cudathreadsperblock << ">>>" << std::endl;

                integrators::core::cuda::compute_kernel<M><<< qmc.cudablocks, qmc.cudathreadsperblock >>>(i, work_this_iteration, total_work_packages,
                                                                                                       static_cast<U*>(*d_z),
                                                                                                       static_cast<D*>(*d_d),
                                                                                                       static_cast<T*>(*d_r),
                                                                                                       d_r_size_over_m, n, m,
                                                                                                       static_cast<I*>(*d_func)
                                                                                                       );
                QMC_CORE_CUDA_SAFE_CALL(cudaPeekAtLastError());
                QMC_CORE_CUDA_SAFE_CALL(cudaDeviceSynchronize());

            };

            template <U M, typename I, typename T, typename D, typename Q>
            void generate_samples(
                                     const Q& qmc,
                                     const U i_start, const U work_this_iteration,
                                     std::unique_ptr<integrators::core::cuda::detail::cuda_memory<U>>& d_z,
                                     std::unique_ptr<integrators::core::cuda::detail::cuda_memory<D>>& d_d,
                                     std::unique_ptr<integrators::core::cuda::detail::cuda_memory<T>>& d_r,
                                     const U n,
                                     std::unique_ptr<integrators::core::cuda::detail::cuda_memory<I>>& d_func,
                                     const int device
                                 )
            {
                if (qmc.verbosity > 1) qmc.logger << "- (" << device << ") computing samples " << i_start << ", work_this_iteration " << work_this_iteration << ", n " << n << std::endl;

                if(qmc.verbosity > 2) qmc.logger << "- (" << device << ") launching gpu kernel<<<" << qmc.cudablocks << "," << qmc.cudathreadsperblock << ">>>" << std::endl;

                integrators::core::cuda::generate_samples_kernel<M><<< qmc.cudablocks, qmc.cudathreadsperblock >>>(i_start, work_this_iteration,
                                                                                                                static_cast<U*>(*d_z),
                                                                                                                static_cast<D*>(*d_d),
                                                                                                                static_cast<T*>(*d_r),
                                                                                                                n,
                                                                                                                static_cast<I*>(*d_func)
                                                                                                                );
                QMC_CORE_CUDA_SAFE_CALL(cudaPeekAtLastError());
                QMC_CORE_CUDA_SAFE_CALL(cudaDeviceSynchronize());

            };
        };
    };
};

#undef QMC_CORE_CUDA_SAFE_CALL

#endif
#endif

#ifndef QMC_CORE_CUDA_SETUP_H
#define QMC_CORE_CUDA_SETUP_H

#ifdef __CUDACC__
#include <memory> // unique_ptr
#include <cassert> // assert
#include <cuda_runtime_api.h> // cudaMemcpy

// (Included Above): #include "detail/cuda_memory.hpp"
// (Included Above): #include "detail/cuda_safe_call.hpp"

#define QMC_CORE_CUDA_SAFE_CALL(err) { integrators::core::cuda::detail::cuda_safe_call((err), __FILE__, __LINE__); }

namespace integrators
{
    namespace core
    {
        namespace cuda
        {
            template <typename I, typename T, typename D>
            void setup_sample(
                                   std::unique_ptr<integrators::core::cuda::detail::cuda_memory<U>>& d_z, const std::vector<U>& z,
                                   std::unique_ptr<integrators::core::cuda::detail::cuda_memory<D>>& d_d, const std::vector<D>& d,
                                   std::unique_ptr<integrators::core::cuda::detail::cuda_memory<T>>& d_r, const U d_r_size_over_m,
                                   const T* r_element, const U r_size_over_m, const U m,
                                   std::unique_ptr<integrators::core::cuda::detail::cuda_memory<I>>& d_func, I& func,
                                   const int device, const U verbosity, const Logger& logger
                             )
            {
                // Set Device
                if (verbosity > 1) logger << "- (" << device << ") setting device" << std::endl;
                QMC_CORE_CUDA_SAFE_CALL(cudaSetDevice(device));
                if (verbosity > 1) logger << "- (" << device << ") device set" << std::endl;

                if(verbosity > 1) logger << "- (" << device << ") allocating d_func,d_z,d_d,d_r" << std::endl;
                d_func.reset( new integrators::core::cuda::detail::cuda_memory<I>(1) );
                d_z.reset( new integrators::core::cuda::detail::cuda_memory<U>(z.size()) );
                d_d.reset( new integrators::core::cuda::detail::cuda_memory<D>(d.size()) );
                d_r.reset( new integrators::core::cuda::detail::cuda_memory<T>(0,d_r_size_over_m*m) ); // allocate and set to 0
                if(verbosity > 1) logger << "- (" << device << ") allocated d_func,d_z,d_d,d_r" << std::endl;

                // copy func (initialize on new active device)
                if(verbosity > 1) logger << "- (" << device << ") initializing function on active device" << std::endl;
                I func_copy = func;
                if(verbosity > 1) logger << "- (" << device << ") initialized function on active device" << std::endl;

                if(verbosity > 1) logger << "- (" << device << ") copying d_func to device memory" << std::endl;
                QMC_CORE_CUDA_SAFE_CALL(cudaMemcpy(static_cast<typename std::remove_const<I>::type*>(*d_func), &func_copy, sizeof(I), cudaMemcpyHostToDevice));
                if(verbosity > 1) logger << "- (" << device << ") copied d_func to device memory" << std::endl;

                // Copy z,d,r,func to device
                if(verbosity > 1) logger << "- (" << device << ") copying z,d to device memory" << std::endl;
                QMC_CORE_CUDA_SAFE_CALL(cudaMemcpy(static_cast<U*>(*d_z), z.data(), z.size() * sizeof(U), cudaMemcpyHostToDevice));
                QMC_CORE_CUDA_SAFE_CALL(cudaMemcpy(static_cast<D*>(*d_d), d.data(), d.size() * sizeof(D), cudaMemcpyHostToDevice));
                // d_r not copied (initialised to 0 above)
                if(verbosity > 1) logger << "- (" << device << ") copied z,d to device memory" << std::endl;

                //        QMC_CORE_CUDA_SAFE_CALL(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1)); // TODO (V2) - investigate if this helps
                //        cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, MyKernel, 0, 0); // TODO (V2) - investigate if this helps - https://devblogs.nvidia.com/cuda-pro-tip-occupancy-api-simplifies-launch-configuration/
            };

            template <typename I, typename T, typename D>
            void setup_evaluate(
                                    std::unique_ptr<integrators::core::cuda::detail::cuda_memory<U>>& d_z, const std::vector<U>& z,
                                    std::unique_ptr<integrators::core::cuda::detail::cuda_memory<D>>& d_d, const std::vector<D>& d,
                                    std::unique_ptr<integrators::core::cuda::detail::cuda_memory<T>>& d_r, const U d_r_size,
                                    std::unique_ptr<integrators::core::cuda::detail::cuda_memory<I>>& d_func, I& func,
                                    const int device, const U verbosity, const Logger& logger
                                )
            {
                // Set Device
                if (verbosity > 1) logger << "- (" << device << ") setting device" << std::endl;
                QMC_CORE_CUDA_SAFE_CALL(cudaSetDevice(device));
                if (verbosity > 1) logger << "- (" << device << ") device set" << std::endl;

                if(verbosity > 1) logger << "- (" << device << ") allocating d_func,d_z,d_d,d_r" << std::endl;
                d_func.reset( new integrators::core::cuda::detail::cuda_memory<I>(1) );
                d_z.reset( new integrators::core::cuda::detail::cuda_memory<U>(z.size()) );
                d_d.reset( new integrators::core::cuda::detail::cuda_memory<D>(d.size()) );
                d_r.reset( new integrators::core::cuda::detail::cuda_memory<T>(d_r_size) );
                if(verbosity > 1) logger << "- (" << device << ") allocated d_func,d_z,d_d,d_r" << std::endl;

                // copy func (initialize on new active device)
                if(verbosity > 1) logger << "- (" << device << ") initializing function on active device" << std::endl;
                I func_copy = func;
                if(verbosity > 1) logger << "- (" << device << ") initialized function on active device" << std::endl;

                if(verbosity > 1) logger << "- (" << device << ") copying d_func to device memory" << std::endl;
                QMC_CORE_CUDA_SAFE_CALL(cudaMemcpy(static_cast<typename std::remove_const<I>::type*>(*d_func), &func_copy, sizeof(I), cudaMemcpyHostToDevice));
                if(verbosity > 1) logger << "- (" << device << ") copied d_func to device memory" << std::endl;

                // Copy z,d,func to device
                if(verbosity > 1) logger << "- (" << device << ") copying z,d,r to device memory" << std::endl;
                QMC_CORE_CUDA_SAFE_CALL(cudaMemcpy(static_cast<U*>(*d_z), z.data(), z.size() * sizeof(U), cudaMemcpyHostToDevice));
                QMC_CORE_CUDA_SAFE_CALL(cudaMemcpy(static_cast<D*>(*d_d), d.data(), d.size() * sizeof(D), cudaMemcpyHostToDevice));
                if(verbosity > 1) logger << "- (" << device << ") copied z,d to device memory" << std::endl;

                //        QMC_CORE_CUDA_SAFE_CALL(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1)); // TODO (V2) - investigate if this helps
                //        cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, MyKernel, 0, 0); // TODO (V2) - investigate if this helps - https://devblogs.nvidia.com/cuda-pro-tip-occupancy-api-simplifies-launch-configuration/
            };
        };
    };
};

#undef QMC_CORE_CUDA_SAFE_CALL

#endif
#endif
#ifndef QMC_CORE_CUDA_TEARDOWN_H
#define QMC_CORE_CUDA_TEARDOWN_H

#ifdef __CUDACC__
#include <memory> // unique_ptr
#include <cuda_runtime_api.h> // cudaMemcpy
#include <thrust/device_ptr.h> // thrust::device_ptr
#include <thrust/reduce.h> // thrust::reduce

// (Included Above): #include "detail/cuda_memory.hpp"
// (Included Above): #include "detail/cuda_safe_call.hpp"

#define QMC_CORE_CUDA_SAFE_CALL(err) { integrators::core::cuda::detail::cuda_safe_call((err), __FILE__, __LINE__); }

namespace integrators
{
    namespace core
    {
        namespace cuda
        {
            template <typename T>
            void teardown_sample(const std::unique_ptr<integrators::core::cuda::detail::cuda_memory<T>>& d_r, const U d_r_size_over_m, T* r_element, const U r_size_over_m, const U m, const int device, const U verbosity, const Logger& logger)
            {
                // Wrap raw device pointer into a thrust::device_ptr
                thrust::device_ptr<T> d_r_ptr(static_cast<T*>(*d_r.get()));
                // Reduce d_r on device and copy result to host
                for (U k = 0; k < m; k++)
                {
                    r_element[k*r_size_over_m] = thrust::reduce(d_r_ptr+k*d_r_size_over_m, d_r_ptr+k*d_r_size_over_m + d_r_size_over_m);
                }
                if (verbosity > 1) logger << "- (" << device << ") copied r to host memory" << std::endl;
            };

            template <typename T>
            void copy_back(const std::unique_ptr<integrators::core::cuda::detail::cuda_memory<T>>& d_r, const U d_r_size, T* r_element, const int device, const U verbosity, const Logger& logger)
            {
                // Copy r_element to host
                QMC_CORE_CUDA_SAFE_CALL(cudaMemcpy(r_element, static_cast<T*>(*d_r), d_r_size * sizeof(T), cudaMemcpyDeviceToHost));
                if (verbosity > 1) logger << "- (" << device << ") copied r_element to host memory" << std::endl;
            };
        };
    };
};

#undef QMC_CORE_CUDA_SAFE_CALL

#endif
#endif
#ifndef QMC_CORE_CUDA_GET_DEVICE_COUNT_H
#define QMC_CORE_CUDA_GET_DEVICE_COUNT_H

#ifdef __CUDACC__
#include <cassert> // assert
#include <cuda_runtime_api.h> // cudaMemcpy

// (Included Above): #include "detail/cuda_safe_call.hpp"

#define QMC_CORE_CUDA_SAFE_CALL(err) { integrators::core::cuda::detail::cuda_safe_call((err), __FILE__, __LINE__); }

namespace integrators
{
    namespace core
    {
        namespace cuda
        {
            inline int get_device_count()
            {
                int device_count;
                QMC_CORE_CUDA_SAFE_CALL(cudaGetDeviceCount(&device_count));
                assert(device_count >= 0);
                return device_count;
            };
        };
    };
};

#undef QMC_CORE_CUDA_SAFE_CALL

#endif
#endif
#ifndef QMC_CORE_GENERIC_COMPUTE_H
#define QMC_CORE_GENERIC_COMPUTE_H

#include <cmath>

// (Included Above): #include "../../math/mul_mod.hpp"

namespace integrators
{
    namespace core
    {
        namespace generic
        {
            template <typename T, typename D, typename I>
            void compute(const U i, const std::vector<U>& z, const std::vector<D>& d, T* r_element, const U r_size_over_m, const U total_work_packages, const U n, const U m, const bool batching, I& func)
            {
                using std::modf;
                
                U batchsize = 1;
                T* points = nullptr;
                if (batching)
                {
                    batchsize = (n / total_work_packages) + ((i < (n % total_work_packages)) ? 1 : 0);
                    points = new T[batchsize];
                }
                std::vector<D> x(func.number_of_integration_variables * batchsize, 0);

                for (U k = 0; k < m; k++)
                {
                    for( U offset = i; offset < n; offset += total_work_packages )
                    {
                        D mynull = 0;

                        for (U sDim = 0; sDim < func.number_of_integration_variables; sDim++)
                        {
                            #define QMC_MODF_CALL modf( integrators::math::mul_mod<D,D>(offset,z.at(sDim),n)/(static_cast<D>(n)) + d.at(k*func.number_of_integration_variables+sDim), &mynull)

                            static_assert(std::is_same<decltype(QMC_MODF_CALL),D>::value, "Downcast detected in integrators::core::generic::compute. Please implement \"D modf(D)\".");
                            x[sDim + (batching ? (func.number_of_integration_variables * (offset / total_work_packages)) : 0)] = QMC_MODF_CALL;

                            #undef QMC_MODF_CALL
                        }
                        if (!batching) {
                            D wgt = 1.;
                            T point = func(x.data());
                            r_element[k*r_size_over_m] += wgt*point;
                        }
                    }

                    if (batching)
                    {
                        func(x.data(), points, batchsize);
                        D wgt = 1.;
                        for ( U i = 0; i != batchsize; ++i)
                        {
                            r_element[k*r_size_over_m] += wgt*points[i];
                        }
                    }
                }
                
                if (batching)
                    delete[] points;
            }

            template <typename T, typename D, typename I>
            void generate_samples(const U i, const std::vector<U>& z, const std::vector<D>& d, T* r_element, const U n, I& func)
            {
                using std::modf;

                D mynull = 0;
                std::vector<D> x(func.number_of_integration_variables,0);

                for (U sDim = 0; sDim < func.number_of_integration_variables; sDim++)
                {
                    #define QMC_MODF_CALL modf( integrators::math::mul_mod<D,D>(i,z.at(sDim),n)/(static_cast<D>(n)) + d.at(sDim), &mynull)

                    static_assert(std::is_same<decltype(QMC_MODF_CALL),D>::value, "Downcast detected in integrators::core::generic::compute. Please implement \"D modf(D)\".");
                    x[sDim] = QMC_MODF_CALL;

                    #undef QMC_MODF_CALL
                }

                *r_element = func(x.data());
            }
        };
    };
};

#endif

// (Included Above): #include "core/has_batching.hpp"
#ifndef QMC_LEAST_SQUARES_H
#define QMC_LEAST_SQUARES_H

#include <algorithm> // std::max
#include <cassert> // assert
#include <cmath> // std::nan
#include <cstddef> // std::nullptr_t
#include <sstream> // std::ostringstream
#include <string> // std::to_string
#include <iostream> // std::endl
#include <iomanip> //  std::setw, std::setfill
#include <vector> // std::vector
#include <iterator> // std::distance
#include <stdexcept> // std::logic_error

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

namespace integrators
{
    namespace core
    {
        using fit_function_jacobian_wrapper_ptr = int (*)(const gsl_vector * , void *, gsl_matrix *);
        using fit_function_hessian_wrapper_ptr = int (*)(const gsl_vector *, const gsl_vector *, void *, gsl_vector *);

        template<typename D, typename F1, typename F2, typename F3>
        struct least_squares_wrapper_t {
            const F1 fit_function;
            const F2 fit_function_jacobian;
            const F3 fit_function_hessian;
            const std::vector<D> x;
            const std::vector<D> y;
        };

        template<typename D, typename F1, typename F2, typename F3>
        int fit_function_wrapper(const gsl_vector * parameters_vector, void *xyfunc_ptr, gsl_vector * f)
        {
            // Unpack data
            const std::vector<D>& x = reinterpret_cast<least_squares_wrapper_t<D,F1,F2,F3>*>(xyfunc_ptr)->x;
            const std::vector<D>& y = reinterpret_cast<least_squares_wrapper_t<D,F1,F2,F3>*>(xyfunc_ptr)->y;
            const F1& fit_function = reinterpret_cast<least_squares_wrapper_t<D,F1,F2,F3>*>(xyfunc_ptr)->fit_function;

            // Compute deviates of fit function for input points
            for (size_t i = 0; i < x.size(); i++)
            {
                gsl_vector_set(f, i, static_cast<double>(fit_function(x[i], gsl_vector_const_ptr(parameters_vector,0)) - y[i]) );
            }

            return GSL_SUCCESS;
        }

        template<typename D, typename F1, typename F2, typename F3>
        int fit_function_jacobian_wrapper(const gsl_vector * parameters_vector, void *xyfunc_ptr, gsl_matrix * J)
        {
            // Unpack data
            const std::vector<D>& x = reinterpret_cast<least_squares_wrapper_t<D,F1,F2,F3>*>(xyfunc_ptr)->x;
            const F2& fit_function_jacobian = reinterpret_cast<least_squares_wrapper_t<D,F1,F2,F3>*>(xyfunc_ptr)->fit_function_jacobian;

            for (size_t i = 0; i < x.size(); i++)
                for (size_t j = 0; j < fit_function_jacobian.num_parameters; j++)
                    gsl_matrix_set(J, i, j, static_cast<double>(fit_function_jacobian(x[i], gsl_vector_const_ptr(parameters_vector,0), j)) );

            return GSL_SUCCESS;
        }

        template<typename D, typename F1, typename F2, typename F3>
        int fit_function_hessian_wrapper(const gsl_vector * parameters_vector, const gsl_vector * v, void *xyfunc_ptr, gsl_vector * fvv)
        {
            // Unpack data
            const std::vector<D>& x = reinterpret_cast<least_squares_wrapper_t<D,F1,F2,F3>*>(xyfunc_ptr)->x;
            const F3& fit_function_hessian = reinterpret_cast<least_squares_wrapper_t<D,F1,F2,F3>*>(xyfunc_ptr)->fit_function_hessian;

            // Compute hessian of fit function
            for (size_t i = 0; i < x.size(); i++)
            {
                gsl_vector_set(fvv, i, static_cast<double>(fit_function_hessian(x[i], gsl_vector_const_ptr(v,0), gsl_vector_const_ptr(parameters_vector,0))) );
            }
            return GSL_SUCCESS;
        }

        struct callback_params_t {
            const U& verbosity;
            Logger& logger;
        };

        inline void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
        {
//            const U& verbosity = reinterpret_cast<callback_params_t*>(params)->verbosity;
            Logger& logger = reinterpret_cast<callback_params_t*>(params)->logger;

            gsl_vector *f = gsl_multifit_nlinear_residual(w);
            gsl_vector *x = gsl_multifit_nlinear_position(w);
            double avratio = gsl_multifit_nlinear_avratio(w);
            double rcond = std::nan("");

            // compute reciprocal condition number of J(x)
            if ( iter > 0 )
                gsl_multifit_nlinear_rcond(&rcond, w);

            const char separator       = ' ';
            const int name_width       = 11;
            const int small_name_width = 9;
            const int num_width        = 15;
            logger << std::left << std::setw(name_width) << std::setfill(separator) << "iter " + std::to_string(iter) + ": ";
            bool display_timing = logger.display_timing;
            logger.display_timing = false;
            for (size_t i = 0; i < x->size; i++)
            {
                logger << std::left << std::setw(small_name_width) << std::setfill(separator) << "p[" + std::to_string(i) + "] = ";
                logger << std::left << std::setw(num_width) << std::setfill(separator) << gsl_vector_get(x, i);
            }
            logger << std::left << std::setw(10) << std::setfill(separator) << "cond(J) = ";
            logger << std::left << std::setw(num_width) << std::setfill(separator) << 1.0 / rcond;
            logger << std::left << std::setw(10) << std::setfill(separator) << "|a|/|v| = ";
            logger << std::left << std::setw(num_width) << std::setfill(separator) << avratio;
            logger << std::left << std::setw(9) << std::setfill(separator) << "|f(x)| = ";
            logger << std::left << std::setw(num_width) << std::setfill(separator) << gsl_blas_dnrm2(f);
            logger << std::endl;
            logger.display_timing = display_timing;
        }

        template <typename D, typename F1, typename F2, typename F3>
        fit_function_jacobian_wrapper_ptr get_fit_function_jacobian_wrapper(const F2& fit_function_jacobian, const U& verbosity, Logger& logger)
        {
            if (verbosity > 1)
                logger << "using analytic jacobian" << std::endl;
            return fit_function_jacobian_wrapper<D,F1,F2,F3>;
        }
        template <typename D, typename F1, typename F2, typename F3>
        std::nullptr_t get_fit_function_jacobian_wrapper(const std::nullptr_t& fit_function_jacobian, const U& verbosity, Logger& logger)
        {
            if (verbosity > 1)
                logger << "using numeric jacobian" << std::endl;
            return nullptr;
        }

        template <typename D, typename F1, typename F2, typename F3>
        fit_function_hessian_wrapper_ptr get_fit_function_hessian_wrapper(const F3& fit_function_hessian, const U& verbosity, Logger& logger)
        {
            if (verbosity > 1)
                logger << "using analytic hessian" << std::endl;
            return fit_function_hessian_wrapper<D,F1,F2,F3>;
        }
        template <typename D, typename F1, typename F2, typename F3>
        std::nullptr_t get_fit_function_hessian_wrapper(const std::nullptr_t& fit_function_hessian, const U& verbosity, Logger& logger)
        {
            if (verbosity > 1)
                logger << "using numeric hessian" << std::endl;
            return nullptr;
        }

        template <typename D, typename F1, typename F2, typename F3, typename std::enable_if< F1::num_parameters != 0, int>::type = 0>
        std::vector<D> least_squares(F1& fit_function, F2& fit_function_jacobian, F3& fit_function_hessian, const std::vector<D>& x, const std::vector<D>& y, const U& verbosity, Logger& logger, const size_t maxiter, const double xtol, const double gtol, const double ftol, gsl_multifit_nlinear_parameters fitparametersgsl)
        {
            const size_t num_points = x.size();
            const size_t num_parameters = fit_function.num_parameters;

            assert(x.size() == y.size());
            assert(num_points > num_parameters + 1);

            least_squares_wrapper_t<D,F1,F2,F3> data = { fit_function, fit_function_jacobian, fit_function_hessian, x, y };

            const gsl_multifit_nlinear_type *method = gsl_multifit_nlinear_trust;
            gsl_multifit_nlinear_workspace *w;
            gsl_multifit_nlinear_fdf fdf;
            gsl_multifit_nlinear_parameters fdf_params = fitparametersgsl;
            gsl_vector *f;
            gsl_matrix *J;
            gsl_matrix *covar = gsl_matrix_alloc(num_parameters, num_parameters);

            // define the function to be minimized
            fdf.f = fit_function_wrapper<D,F1,F2,F3>;
            fdf.df = get_fit_function_jacobian_wrapper<D,F1,F2,F3>(fit_function_jacobian, verbosity, logger);
            fdf.fvv = get_fit_function_hessian_wrapper<D,F1,F2,F3>(fit_function_hessian, verbosity, logger);
            fdf.n = num_points;
            fdf.p = num_parameters;
            fdf.params = &data;
            
            // compute dx/dy of input points, which should be used as an additional weight in the evaluation of chisq
            std::vector<double> dxdy(x.size());
            double maxwgt = 0.;

            const size_t nsteps = 1; 
            for (size_t i = 0; i < x.size(); i++)
            {
                D dy = (i<nsteps) ? D(0) : -y[i-nsteps];  
                D dx = (i<nsteps) ? D(0) : -x[i-nsteps];
                if(i != x.size()-nsteps)
                {
                    dy += y[i+nsteps];
                    dx += x[i+nsteps];
                }
                else
                {
                    dy += D(1);
                    dx += D(1);
                }
                dxdy[i] = static_cast<double>(dx/dy);
                
                maxwgt=std::max(maxwgt,dxdy[i]);
            }
            
            // the gsl fit doesn't seem to work with weights>1 
            for(size_t i=0; i< x.size(); i++)
            {
                dxdy[i]/=maxwgt;
            }
            
            gsl_vector_view wgt = gsl_vector_view_array(dxdy.data(), dxdy.size());

            double chisq,chisq0;
            int status, info;

            // allocate workspace with parameters
            w = gsl_multifit_nlinear_alloc(method, &fdf_params, num_points, num_parameters);

            std::vector<std::vector<D>> fit_parameters;
            std::vector<double> fit_chisqs;
            fit_chisqs.reserve(fit_function.initial_parameters.size());
            for (size_t i = 0; i < fit_function.initial_parameters.size(); i++)
            {
                std::vector<double> initial_parameters(fit_function.initial_parameters.at(i).begin(),fit_function.initial_parameters.at(i).end()); // note: cast to double

                if( initial_parameters.size() != fit_function.num_parameters)
                    throw std::domain_error("least_squares called with incorrect number of initial_parameters (" + std::to_string(initial_parameters.size()) + "), expected " +  std::to_string(fit_function.num_parameters) + " parameters");

                if (verbosity > 0)
                {
                    logger << "-- running fit (run " << i << ")" << " --" << std::endl;
                    std::ostringstream initial_parameters_stream;
                    for(const auto& elem: initial_parameters)
                        initial_parameters_stream << elem << " ";
                    logger << "with initial_parameters " << initial_parameters_stream.str() << std::endl;
                }

                gsl_vector_view pv = gsl_vector_view_array(initial_parameters.data(), num_parameters);

                // initialize solver with starting point
                gsl_multifit_nlinear_winit(&pv.vector, &wgt.vector, &fdf, w);

                // compute initial cost function
                f = gsl_multifit_nlinear_residual(w);
                gsl_blas_ddot(f, f, &chisq0);

                // solve the system with a maximum of "maxiter" iterations
                callback_params_t callback_params{verbosity,logger};
                if (verbosity > 2)
                    status = gsl_multifit_nlinear_driver(maxiter, xtol, gtol, ftol, callback, &callback_params, &info, w);
                else
                    status = gsl_multifit_nlinear_driver(maxiter, xtol, gtol, ftol, nullptr, nullptr, &info, w);

                // compute covariance of best fit parameters
                J = gsl_multifit_nlinear_jac(w);
                gsl_multifit_nlinear_covar(J, 0., covar);

                // compute final cost
                gsl_blas_ddot(f, f, &chisq);

                // store fit parameters
                std::vector<D> this_fit_parameters;
                this_fit_parameters.reserve(num_parameters);
                for (size_t j = 0; j < num_parameters; j++)
                    this_fit_parameters.push_back( gsl_vector_get(w->x, j) );
                fit_parameters.push_back( this_fit_parameters );

                // Report output of fit function
                if (verbosity > 1)
                {
                    double dof = num_points - num_parameters - 1;
                    double c = std::max(1., sqrt(chisq/dof));

                    logger << "-- fit output (run " << i << ")" << " --" << std::endl;
                    logger << "summary from method "   << gsl_multifit_nlinear_name(w) << " " << gsl_multifit_nlinear_trs_name(w) << std::endl;
                    logger << "number of iterations: " << gsl_multifit_nlinear_niter(w) << std::endl;
                    logger << "function evaluations: " << fdf.nevalf << std::endl;
                    logger << "Jacobian evaluations: " << fdf.nevaldf << std::endl;
                    logger << "Hessian evaluations: " << fdf.nevalfvv << std::endl;
                    if (info == 0)
                        logger << "reason for stopping: " << "maximal number of iterations (maxiter)" << std::endl;
                    else if (info == 1)
                        logger << "reason for stopping: " << "small step size (xtol)" << std::endl;
                    else if (info == 2)
                        logger << "reason for stopping: " << "small gradient (gtol)" << std::endl;
                    else
                        logger << "reason for stopping: " << "unknown" << std::endl;
                    logger << "initial |f(x)| = "      << sqrt(chisq0) << std::endl;;
                    logger << "final   |f(x)| = "      << sqrt(chisq) << std::endl;
                    logger << "chisq/dof = "           << chisq/dof << std::endl;
                    for (size_t j = 0; j < num_parameters; j++)
                        logger << "fit_parameters[" << j << "] = " << this_fit_parameters.at(j) << " +/- " << c*sqrt(gsl_matrix_get(covar,j,j)*chisq/dof) << std::endl;
                    logger << "status = "              << gsl_strerror(status) << std::endl;
                    logger << "-----------" << std::endl;
                }

                fit_chisqs.push_back(chisq);

            }


            gsl_multifit_nlinear_free(w);
            gsl_matrix_free(covar);

            // get index of best fit (minimum chisq)
            const long best_fit_index = std::distance(fit_chisqs.begin(), std::min_element(fit_chisqs.begin(),fit_chisqs.end()));
            assert(best_fit_index >= 0);

            if (verbosity > 0)
            {
                if (verbosity>2) logger << "choosing fit run " << best_fit_index << std::endl;
                std::ostringstream final_parameters_stream;
                for(const auto& elem: fit_parameters.at(static_cast<size_t>(best_fit_index)))
                    final_parameters_stream << elem << " ";
                logger << "fit final_parameters " << final_parameters_stream.str() << std::endl;
                logger << "-----------" << std::endl;
            }

            return fit_parameters.at(static_cast<size_t>(best_fit_index));
        }
    
        template <typename D, typename F1, typename F2, typename F3, typename std::enable_if< F1::num_parameters == 0, int >::type = 0>
        std::vector<D> least_squares(F1& fit_function, F2& fit_function_jacobian, F3& fit_function_hessian, const std::vector<D>& x, const std::vector<D>& y, const U& verbosity, Logger& logger, const size_t maxiter, const double xtol, const double gtol, const double ftol, gsl_multifit_nlinear_parameters fitparametersgsl)
        {
            throw std::logic_error("least_squares called on function with no parameters to be fitted");
        }
    
    };
};

#endif
#ifndef QMC_CORE_REDUCE_H
#define QMC_CORE_REDUCE_H

// (Included Above): #include "../types/result.hpp"

namespace integrators
{
    namespace core
    {
        template <typename T>
        integrators::result<T> reduce(const std::vector<T>& r, const U n, const U m, std::vector<result<T>> & previous_iterations, const U& verbosity, const Logger& logger)
        {
            if (verbosity > 1)
            {
                logger << "-- qmc::reduce called --" << std::endl;
                for(const auto& previous_result : previous_iterations)
                {
                    logger << "previous_result: integral " << previous_result.integral << ", error " << previous_result.error << ", n " << previous_result.n << ", m " << previous_result.m << std::endl;
                }
            }

            T mean = {0.};
            T variance = {0.};
            U previous_m = 0;
            U previous_num_iterations = 0;
            U previous_num_evaluations = 0;
            if(!previous_iterations.empty())
            {
                result<T> & previous_res = previous_iterations.back();
                previous_num_iterations = previous_res.iterations;
                previous_num_evaluations = previous_res.evaluations;
                if(previous_res.n == n)
                {
                    if (verbosity>2) logger << "using additional shifts to improve previous iteration" << std::endl;
                    previous_m = previous_res.m;
                    mean = previous_res.integral*static_cast<T>(n);
                    variance = integrators::overloads::compute_variance_from_error(previous_res.error);
                    variance *= static_cast<T>(previous_res.m-1) * static_cast<T>(previous_res.m) * static_cast<T>(previous_res.n) * static_cast<T>(previous_res.n);
                    previous_iterations.pop_back();
                }
            }
            U r_size = r.size();
            for(U k = 0; k < m; k++)
            {
                T sum = {0.};
                T delta = {0.};
                for (U i = 0; i<r_size/m; i++)
                {
                    sum += r.at(k*r_size/m+i);
                }
                if (verbosity > 1) logger << "shift " << k+previous_m << " result: " << sum/static_cast<T>(n) << std::endl;
                // Compute Variance using online algorithm (Knuth, The Art of Computer Programming)
                delta = sum - mean;
                mean = mean + delta/(static_cast<T>(k+previous_m+1));
                variance = integrators::overloads::compute_variance(mean, variance, sum, delta);
            }
            T integral = mean/(static_cast<T>(n));
            variance = variance/( static_cast<T>(m+previous_m-1) * static_cast<T>(m+previous_m) * static_cast<T>(n) * static_cast<T>(n) ); // variance of the mean
            T error = integrators::overloads::compute_error(variance);
            previous_iterations.push_back({integral, error, n, m+previous_m, 1+previous_num_iterations, n*m+previous_num_evaluations});
            if (verbosity > 0)
                logger << "integral " << integral << ", error " << error << ", n " << n << ", m " << m+previous_m << std::endl;
            return {integral, error, n, m+previous_m, 1+previous_num_iterations, n*m+previous_num_evaluations};
        };
    };
};

#endif
#ifndef QMC_MEMBERS_H
#define QMC_MEMBERS_H

#include <cstddef> // size_t
#include <cmath> // modf, abs, sqrt, pow
#include <stdexcept> // domain_error, invalid_argument
#include <thread> // thread
#include <algorithm> // min, max
#include <type_traits> // make_signed, enable_if
#include <limits> // numeric_limits
#include <string> // to_string
#include <vector>
#include <iostream>
#include <iterator> // advance
#include <mutex>
#include <memory> // unique_ptr
#include <numeric> // partial_sum
#include <cassert> // assert
#include <chrono>

#include <gsl/gsl_multifit_nlinear.h>

namespace integrators
{
    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    void Qmc<T,D,M,P,F,G,H>::init_z(std::vector<U>& z, const U n, const U number_of_integration_variables) const
    {
        z = generatingvectors.at(n);
        if ( number_of_integration_variables > z.size() )
            throw std::domain_error("number_of_integration_variables > generating vector dimension. Please supply a generating vector table with a larger number of dimensions.");
        z.resize(number_of_integration_variables);
    };
    
    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    void Qmc<T,D,M,P,F,G,H>::init_d(std::vector<D>& d, const U m, const U number_of_integration_variables)
    {
        d.clear();
        d.reserve(m*number_of_integration_variables);
        for (U k = 0; k < m; k++)
            for (U sDim = 0; sDim < number_of_integration_variables; sDim++)
                d.push_back(uniform_distribution(randomgenerator));
    };
    
    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    void Qmc<T,D,M,P,F,G,H>::init_r(std::vector<T>& r, const U m, const U r_size_over_m) const
    {
        r.clear();
        r.resize(m * r_size_over_m, {0.});
    };
    
    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    template <typename I>
    void Qmc<T,D,M,P,F,G,H>::sample_worker(const U thread_id, U& work_queue, std::mutex& work_queue_mutex, const std::vector<U>& z, const std::vector<D>& d, std::vector<T>& r, const U total_work_packages, const U n, const U m, I& func, const int device, D& time_in_ns, U& points_computed) const
    {
        std::chrono::steady_clock::time_point time_before_compute = std::chrono::steady_clock::now();

        points_computed = 0;

        // Setup worker
#ifdef __CUDACC__
        // define device pointers (must be accessible in local scope of the entire function)
        U d_r_size = m*cudablocks*cudathreadsperblock;
        std::unique_ptr<integrators::core::cuda::detail::cuda_memory<I>> d_func;
        std::unique_ptr<integrators::core::cuda::detail::cuda_memory<U>> d_z;
        std::unique_ptr<integrators::core::cuda::detail::cuda_memory<D>> d_d;
        std::unique_ptr<integrators::core::cuda::detail::cuda_memory<T>> d_r;
#endif
        U i;
        U  work_this_iteration;
        if (device == -1) {
            work_this_iteration = 1;
        } else {
#ifdef __CUDACC__
            work_this_iteration = cudablocks*cudathreadsperblock;
            integrators::core::cuda::setup_sample(d_z, z, d_d, d, d_r, d_r_size/m, &r[thread_id], r.size()/m, m, d_func, func, device, verbosity, logger);
#else
            throw std::invalid_argument("qmc::sample called with device != -1 (CPU) but CUDA not supported by compiler, device: " + std::to_string(device));
#endif
        }

        bool work_remaining = true;
        while( work_remaining )
        {
            // Get work
            work_queue_mutex.lock();
            if (work_queue == 0)
            {
                work_remaining=false;
                i = 0;
            }
            else if (work_queue >= work_this_iteration)
            {
                work_queue-=work_this_iteration;
                i = work_queue;
            }
            else
            {
                work_this_iteration = work_queue;
                work_queue = 0;
                i = 0;
            }
            work_queue_mutex.unlock();
            
            if( !work_remaining )
                break;
            
            // Do work
            if (device == -1)
            {
                integrators::core::generic::compute(i, z, d, &r[thread_id], r.size()/m, total_work_packages, n, m, batching, func);
            }
            else
            {
#ifdef __CUDACC__
                integrators::core::cuda::compute<M>(*this, i, work_this_iteration, total_work_packages, d_z, d_d, d_r, d_r_size/m, n, m, d_func, device);
#endif
            }

            points_computed += work_this_iteration*m;

        }

        // Teardown worker
#ifdef __CUDACC__
        if (device != -1) {
            integrators::core::cuda::teardown_sample(d_r, d_r_size/m, &r[thread_id], r.size()/m, m, device, verbosity, logger);
        }
#endif

        std::chrono::steady_clock::time_point time_after_compute = std::chrono::steady_clock::now();
        time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_after_compute - time_before_compute).count();

    };
    
    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    template <typename I>
    result<T> Qmc<T,D,M,P,F,G,H>::sample(I& func, const U n, const U m, std::vector<result<T>> & previous_iterations)
    {
        std::vector<U> z;
        std::vector<D> d;
        std::vector<T> r;

        result<T> res;

        U points_per_package = std::min(maxnperpackage, n); // points to compute per thread per work_package
        U total_work_packages = n/points_per_package; // Set total number of work packages to be computed
        if( n%points_per_package != 0) total_work_packages++;

        U extra_threads = devices.size() - devices.count(-1);
        
        // Memory required for result vector
        U r_size_over_m = extra_threads; // non-cpu workers
        if (devices.count(-1) != 0 && cputhreads > 0)
        {
            r_size_over_m += cputhreads; // cpu-workers
        }

        U iterations = (m+maxmperpackage-1)/maxmperpackage;
        U shifts_per_iteration = std::min(m,maxmperpackage);
        for(U iteration = 0; iteration < iterations; iteration++)
        {
            U shifts = shifts_per_iteration;
            if ( iteration == iterations-1)
            {
                // last iteration => compute remaining shifts
                shifts = m%maxmperpackage == 0 ? std::min(m,maxmperpackage) : m%maxmperpackage;
            }

            // Generate z, d, r
            init_z(z, n, func.number_of_integration_variables);
            init_d(d, shifts, func.number_of_integration_variables);
            init_r(r, shifts, r_size_over_m);

            if (verbosity > 0)
            {
                logger << "-- qmc::sample called --" << std::endl;
                logger << "func.number_of_integration_variables " << func.number_of_integration_variables << std::endl;
                logger << "minn " << minn << std::endl;
                logger << "minm " << minm << std::endl;
                logger << "epsrel " << epsrel << std::endl;
                logger << "epsabs " << epsabs << std::endl;
                logger << "maxeval " << maxeval << std::endl;
                logger << "cputhreads " << cputhreads << std::endl;
                logger << "maxnperpackage " << maxnperpackage << std::endl;
                logger << "maxmperpackage " << maxmperpackage << std::endl;
                logger << "cudablocks " << cudablocks << std::endl;
                logger << "cudathreadsperblock " << cudathreadsperblock << std::endl;
                logger << "devices ";
                bool display_timing = logger.display_timing;
                logger.display_timing = false;
                for (const int& i : devices)
                    logger << i << " ";
                logger << std::endl;
                logger.display_timing = display_timing;
                logger << "n " << n << std::endl;
                logger << "m " << m << std::endl;
                logger << "shifts " << shifts << std::endl;
                logger << "iterations " << iterations << std::endl;
                logger << "total_work_packages " << total_work_packages << std::endl;
                logger << "points_per_package " << points_per_package << std::endl;
                logger << "r " << shifts << "*" << r_size_over_m << std::endl;
            }

            std::chrono::steady_clock::time_point time_before_compute = std::chrono::steady_clock::now();

            if ( cputhreads == 1 && devices.size() == 1 && devices.count(-1) == 1)
            {
                // Compute serially on cpu
                if (verbosity > 2) logger << "computing serially" << std::endl;
                for( U i=0; i < total_work_packages; i++)
                {
                    integrators::core::generic::compute(i, z, d, &r[0], r.size()/shifts, total_work_packages, n, shifts, batching, func);
                }
            }
            else
            {
                // Create threadpool
                if (verbosity > 2)
                {
                    logger << "distributing work" << std::endl;
                    if ( devices.count(-1) != 0)
                        logger << "creating " << std::to_string(cputhreads) << " cputhreads," << std::to_string(extra_threads) << " non-cpu threads" << std::endl;
                    else
                        logger << "creating " << std::to_string(extra_threads) << " non-cpu threads" << std::endl;
                }

                // Setup work queue
                std::mutex work_queue_mutex;
                U work_queue = total_work_packages;

                // Launch worker threads
                U thread_id = 0;
                U thread_number = 0;
                std::vector<std::thread> thread_pool;
                thread_pool.reserve(cputhreads+extra_threads);
                std::vector<D> time_in_ns_per_thread(cputhreads+extra_threads,D(0));
                std::vector<U> points_computed_per_thread(cputhreads+extra_threads,U(0));
                for (int device : devices)
                {
                    if( device != -1)
                    {
#ifdef __CUDACC__
                        thread_pool.push_back( std::thread( &Qmc<T,D,M,P,F,G,H>::sample_worker<I>, this, thread_id, std::ref(work_queue), std::ref(work_queue_mutex), std::cref(z), std::cref(d), std::ref(r), total_work_packages, n, shifts, std::ref(func), device, std::ref(time_in_ns_per_thread[thread_number]), std::ref(points_computed_per_thread[thread_number])  ) ); // Launch non-cpu workers
                        thread_id += 1;
                        thread_number += 1;
#else
                        throw std::invalid_argument("qmc::sample called with device != -1 (CPU) but CUDA not supported by compiler, device: " + std::to_string(device));
#endif
                    }
                }
                if( devices.count(-1) != 0 && cputhreads > 0)
                {
                    for ( U i=0; i < cputhreads; i++)
                    {
                        thread_pool.push_back( std::thread( &Qmc<T,D,M,P,F,G,H>::sample_worker<I>, this, thread_id, std::ref(work_queue), std::ref(work_queue_mutex), std::cref(z), std::cref(d), std::ref(r), total_work_packages, n, shifts, std::ref(func), -1, std::ref(time_in_ns_per_thread[thread_number]), std::ref(points_computed_per_thread[thread_number]) ) ); // Launch cpu workers
                        thread_id += 1;
                        thread_number += 1;
                    }
                }
                // Destroy threadpool
                for( std::thread& thread : thread_pool )
                    thread.join();
                thread_pool.clear();

                if(verbosity > 2)
                {
                    for( U i=0; i< extra_threads; i++)
                    {
                        logger << "(" << i << ") Million Function Evaluations/s: " << D(1000)*D(points_computed_per_thread[i])/D(time_in_ns_per_thread[i]) << " Mfeps (Approx)" << std::endl;
                    }
                    if( devices.count(-1) != 0 && cputhreads > 0)
                    {
                        D time_in_ns_on_cpu = 0;
                        U points_computed_on_cpu = 0;
                        for( U i=extra_threads; i< extra_threads+cputhreads; i++)
                        {
                            points_computed_on_cpu += points_computed_per_thread[i];
                            time_in_ns_on_cpu = std::max(time_in_ns_on_cpu,time_in_ns_per_thread[i]);
                        }
                        logger << "(-1) Million Function Evaluations/s: " << D(1000)*D(points_computed_on_cpu)/D(time_in_ns_on_cpu) << " Mfeps (Approx)" << std::endl;
                    }
                }
            }

            std::chrono::steady_clock::time_point time_after_compute = std::chrono::steady_clock::now();

            if(verbosity > 2)
            {
                D mfeps = D(1000)*D(n*shifts)/D(std::chrono::duration_cast<std::chrono::nanoseconds>(time_after_compute - time_before_compute).count()); // million function evaluations per second
                logger << "(Total) Million Function Evaluations/s: " << mfeps << " Mfeps" << std::endl;
            }
            res = integrators::core::reduce(r, n, shifts,  previous_iterations, verbosity, logger);
        }
        return res;
    };
    
    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    template <typename I>
    void Qmc<T,D,M,P,F,G,H>::evaluate_worker(const U thread_id, U& work_queue, std::mutex& work_queue_mutex, const std::vector<U>& z, const std::vector<D>& d, std::vector<T>& r, const U n, I& func, const int device, D& time_in_ns, U& points_computed) const
    {
        std::chrono::steady_clock::time_point time_before_compute = std::chrono::steady_clock::now();

        points_computed = 0;

        // Setup worker
#ifdef __CUDACC__
        // define device pointers (must be accessible in local scope of the entire function)
        U d_r_size = cudablocks*cudathreadsperblock;
        std::unique_ptr<integrators::core::cuda::detail::cuda_memory<I>> d_func;
        std::unique_ptr<integrators::core::cuda::detail::cuda_memory<U>> d_z;
        std::unique_ptr<integrators::core::cuda::detail::cuda_memory<D>> d_d;
        std::unique_ptr<integrators::core::cuda::detail::cuda_memory<T>> d_r;
#endif
        U i;
        U  work_this_iteration;
        if (device == -1) {
            work_this_iteration = 1;
        } else {
#ifdef __CUDACC__
            work_this_iteration = cudablocks*cudathreadsperblock;
            integrators::core::cuda::setup_evaluate(d_z, z, d_d, d, d_r, d_r_size, d_func, func, device, verbosity, logger);
#else
            throw std::invalid_argument("qmc::sample called with device != -1 (CPU) but CUDA not supported by compiler, device: " + std::to_string(device));
#endif
        }

        bool work_remaining = true;
        while( work_remaining )
        {
            // Get work
            work_queue_mutex.lock();
            if (work_queue == 0)
            {
                work_remaining=false;
                i = 0;
            }
            else if (work_queue >= work_this_iteration)
            {
                work_queue-=work_this_iteration;
                i = work_queue;
            }
            else
            {
                work_this_iteration = work_queue;
                work_queue = 0;
                i = 0;
            }
            work_queue_mutex.unlock();
            
            if( !work_remaining )
                break;
            
            // Do work
            if (device == -1)
            {
                integrators::core::generic::generate_samples(i, z, d, &r[i], n, func);
            }
            else
            {
#ifdef __CUDACC__
                integrators::core::cuda::generate_samples<M>(*this, i, work_this_iteration, d_z, d_d, d_r, n, d_func, device);
#endif
            }

            points_computed += work_this_iteration;


            // copy results to host
#ifdef __CUDACC__
            if (device != -1) {
                integrators::core::cuda::copy_back(d_r, work_this_iteration, &r[i], device, verbosity, logger);
            }
#endif
        }

        std::chrono::steady_clock::time_point time_after_compute = std::chrono::steady_clock::now();
        time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_after_compute - time_before_compute).count();

    };
    
    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    template <typename I>
    samples<T,D> Qmc<T,D,M,P,F,G,H>::evaluate(I& func)
    {
        if ( func.number_of_integration_variables < 1 )
            throw std::invalid_argument("qmc::evaluate called with func.number_of_integration_variables < 1. Check that your integrand depends on at least one variable of integration.");
        if ( func.number_of_integration_variables > M )
            throw std::invalid_argument("qmc::evaluate called with func.number_of_integration_variables > M. Please increase M (maximal number of integration variables).");
        if ( devices.size() == 1 && devices.count(-1) == 1 && cputhreads == 0 )
            throw std::domain_error("qmc::evaluate called with cputhreads = 0 and devices = {-1} (CPU). Please set cputhreads to a positive integer or specify at least one non-CPU device.");
#ifndef __CUDACC__
        if ( devices.size() != 1 || devices.count(-1) != 1)
            throw std::invalid_argument("qmc::evaluate called with devices != {-1} (CPU) but CUDA not supported by compiler.");
#endif

        // allocate memory
        samples<T,D> res;
        U& n = res.n;
        n = get_next_n(evaluateminn); // get next available n >= evaluateminn
        std::vector<U>& z = res.z;
        std::vector<D>& d = res.d;
        std::vector<T>& r = res.r;

        // initialize z, d, r
        init_z(z, n, func.number_of_integration_variables);
        init_d(d, 1, func.number_of_integration_variables);
        init_r(r, 1, n); // memory required for result vector

        U extra_threads = devices.size() - devices.count(-1);

        if (verbosity > 0)
        {
            logger << "-- qmc::evaluate called --" << std::endl;
            logger << "func.number_of_integration_variables " << func.number_of_integration_variables << std::endl;
            logger << "evaluateminn " << evaluateminn << std::endl;
            logger << "cputhreads " << cputhreads << std::endl;
            logger << "cudablocks " << cudablocks << std::endl;
            logger << "cudathreadsperblock " << cudathreadsperblock << std::endl;
            logger << "devices ";
            bool display_timing = logger.display_timing;
            logger.display_timing = false;
            for (const int& i : devices)
                logger << i << " ";
            logger << std::endl;
            logger.display_timing = display_timing;
            logger << "n " << n << std::endl;
        }

        std::chrono::steady_clock::time_point time_before_compute = std::chrono::steady_clock::now();

        if ( cputhreads == 1 && devices.size() == 1 && devices.count(-1) == 1)
        {
            // Compute serially on cpu
            if (verbosity > 2) logger << "computing serially" << std::endl;
            for( U i=0; i < n; i++)
            {
                integrators::core::generic::generate_samples(i, z, d, &r[i], n, func);
            }
        }
        else
        {
            // Create threadpool
            if (verbosity > 2)
            {
                logger << "distributing work" << std::endl;
                if ( devices.count(-1) != 0)
                    logger << "creating " << std::to_string(cputhreads) << " cputhreads," << std::to_string(extra_threads) << " non-cpu threads" << std::endl;
                else
                    logger << "creating " << std::to_string(extra_threads) << " non-cpu threads" << std::endl;
            }

            // Setup work queue
            std::mutex work_queue_mutex;
            U work_queue = n;

            // Launch worker threads
            U thread_id = 0;
            U thread_number = 0;
            std::vector<std::thread> thread_pool;
            thread_pool.reserve(cputhreads+extra_threads);
            std::vector<D> time_in_ns_per_thread(cputhreads+extra_threads,D(0));
            std::vector<U> points_computed_per_thread(cputhreads+extra_threads,U(0));
            for (int device : devices)
            {
                if( device != -1)
                {
#ifdef __CUDACC__
                    thread_pool.push_back( std::thread( &Qmc<T,D,M,P,F,G,H>::evaluate_worker<I>, this, thread_id, std::ref(work_queue), std::ref(work_queue_mutex), std::cref(z), std::cref(d), std::ref(r), n, std::ref(func), device, std::ref(time_in_ns_per_thread[thread_number]), std::ref(points_computed_per_thread[thread_number]) ) ); // Launch non-cpu workers
                    thread_id += cudablocks*cudathreadsperblock;
                    thread_number += 1;
#else
                    throw std::invalid_argument("qmc::sample called with device != -1 (CPU) but CUDA not supported by compiler, device: " + std::to_string(device));
#endif
                }
            }
            if( devices.count(-1) != 0 && cputhreads > 0)
            {
                for ( U i=0; i < cputhreads; i++)
                {
                    thread_pool.push_back( std::thread( &Qmc<T,D,M,P,F,G,H>::evaluate_worker<I>, this, thread_id, std::ref(work_queue), std::ref(work_queue_mutex), std::cref(z), std::cref(d), std::ref(r), n, std::ref(func), -1, std::ref(time_in_ns_per_thread[thread_number]), std::ref(points_computed_per_thread[thread_number]) ) ); // Launch cpu workers
                    thread_id += 1;
                    thread_number += 1;
                }
            }
            // Destroy threadpool
            for( std::thread& thread : thread_pool )
                thread.join();
            thread_pool.clear();

            if(verbosity > 2)
            {
                for( U i=0; i< extra_threads; i++)
                {
                    logger << "(" << i << ") Million Function Evaluations/s: " << D(1000)*D(points_computed_per_thread[i])/D(time_in_ns_per_thread[i]) << " Mfeps (Approx)" << std::endl;
                }
                if( devices.count(-1) != 0 && cputhreads > 0)
                {
                    D time_in_ns_on_cpu = 0;
                    U points_computed_on_cpu = 0;
                    for( U i=extra_threads; i< extra_threads+cputhreads; i++)
                    {
                        points_computed_on_cpu += points_computed_per_thread[i];
                        time_in_ns_on_cpu = std::max(time_in_ns_on_cpu,time_in_ns_per_thread[i]);
                    }
                    logger << "(-1) Million Function Evaluations/s: " << D(1000)*D(points_computed_on_cpu)/D(time_in_ns_on_cpu) << " Mfeps (Approx)" << std::endl;
                }
            }
        }

        std::chrono::steady_clock::time_point time_after_compute = std::chrono::steady_clock::now();

        if(verbosity > 2)
        {
            D mfeps = D(1000)*D(n)/D(std::chrono::duration_cast<std::chrono::nanoseconds>(time_after_compute - time_before_compute).count()); // million function evaluations per second
            logger << "(Total) Million Function Evaluations/s: " << mfeps << " Mfeps" << std::endl;
        }

        return res;
    };

    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    template <typename I>
    typename F<I,D,M>::transform_t Qmc<T,D,M,P,F,G,H>::fit(I& func)
    {
        using std::abs;

        typename F<I,D,M>::function_t fit_function;
        typename F<I,D,M>::jacobian_t fit_function_jacobian;
        typename F<I,D,M>::hessian_t fit_function_hessian;
        typename F<I,D,M>::transform_t fit_function_transform(func);

        if (fit_function.num_parameters <= 0) {
            return fit_function_transform;
        } else {
            std::vector<D> x,y;
            std::vector<std::vector<D>> fit_parameters;
            fit_parameters.reserve(func.number_of_integration_variables);

            // Generate data to be fitted
            integrators::samples<T,D> result = evaluate(func);

            // fit fit_function
            for (U sdim = 0; sdim < func.number_of_integration_variables; ++sdim)
            {
                // compute the x values
                std::vector<D> unordered_x;
                unordered_x.reserve(result.n);
                for (size_t i = 0; i < result.n; ++i)
                {
                    unordered_x.push_back( result.get_x(i, sdim) );
                }

                // sort by x value
                std::vector<size_t> sort_key = math::argsort(unordered_x);
                x.clear();
                y.clear();
                x.reserve( sort_key.size() );
                y.reserve( sort_key.size() );
                for (const auto& idx : sort_key)
                {
                    x.push_back( unordered_x.at(idx) );
                    y.push_back( abs(result.r.at(idx)) );
                }

                // compute cumulative sum
                std::partial_sum(y.begin(), y.end(), y.begin());
                for (auto& element : y)
                {
                    element /= y.back();
                }

                // reduce number of sampling points for fit
                std::vector<D> xx;
                std::vector<D> yy;
                for ( size_t i = fitstepsize/2; i<x.size(); i+=fitstepsize)
                {
                        xx.push_back(x.at(i));
                        yy.push_back(y.at(i));
                }

                // run a least squares fit
                fit_parameters.push_back( core::least_squares(fit_function,fit_function_jacobian, fit_function_hessian, yy,xx,verbosity,logger, fitmaxiter, fitxtol, fitgtol, fitftol, fitparametersgsl) );
            }

            for (size_t d = 0; d < fit_function_transform.number_of_integration_variables; ++d)
                for (size_t i = 0; i < fit_parameters.at(d).size(); ++i)
                    fit_function_transform.p[d][i] = fit_parameters.at(d).at(i);

            return fit_function_transform;
        }
    };

    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    void Qmc<T,D,M,P,F,G,H>::update(const result<T>& res, U& n, U& m) const
    {
        using std::pow;

        if (verbosity > 2) logger << "-- qmc::update called --" << std::endl;

        const D MAXIMUM_ERROR_RATIO = static_cast<D>(20);
        const D EXPECTED_SCALING = static_cast<D>(0.8); // assume error scales as n^(-expectedScaling)

        D error_ratio = std::min(integrators::overloads::compute_error_ratio(res, epsrel, epsabs, errormode),MAXIMUM_ERROR_RATIO);
        if (error_ratio < static_cast<D>(1))
        {
            if (verbosity > 2) logger << "error goal reached" << std::endl;
            return;
        }

        if(res.evaluations > maxeval)
        {
            if (verbosity > 2) logger << "maxeval reached" << std::endl;
            return;
        }

        U new_m = minm;
        #define QMC_POW_CALL pow(error_ratio,static_cast<D>(1)/EXPECTED_SCALING)
        static_assert(std::is_same<decltype(QMC_POW_CALL),D>::value, "Downcast detected in qmc::update(. Please implement \"D pow(D)\".");
        U new_n = get_next_n(static_cast<U>(static_cast<D>(n)*QMC_POW_CALL));
        #undef QMC_POW_CALL
        if ( new_n <= n or ( error_ratio*error_ratio - static_cast<D>(1) < static_cast<D>(new_n)/static_cast<D>(n)))
        {
            // n did not increase, or increasing m will be faster, increase m
            new_n = n;
            new_m = static_cast<U>(static_cast<D>(m)*error_ratio*error_ratio)+1-m;
            if (verbosity > 2)
                logger << "n did not increase, or increasing m will be faster, increasing m to " << new_m << "." << std::endl;
        }
        if ( maxeval < res.evaluations + new_n*new_m)
        {
            // Decrease new_n
            if ( verbosity > 2 )
                logger << "requested number of function evaluations greater than maxeval, reducing n." << std::endl;
            new_n = std::max(n, get_next_n((maxeval-res.evaluations)/new_m));
        }
        if ( n == new_n && maxeval < res.evaluations + new_n*new_m)
        {
            // Decrease new_m
            if ( verbosity > 2 )
                logger << "requested number of function evaluations greater than maxeval, reducing m." << std::endl;
            new_m = std::max(U(1),(maxeval-res.evaluations)/new_n);
        }
        n = new_n;
        m = new_m;
        if(verbosity > 1 ) logger << "updated n m " << n << " " << m << std::endl;

    };

    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    U Qmc<T,D,M,P,F,G,H>::get_next_n(U preferred_n) const
    {
        U n;
        if ( generatingvectors.lower_bound(preferred_n) == generatingvectors.end() )
        {
            n = generatingvectors.rbegin()->first;
            if (verbosity > 0)
                logger << "Qmc integrator does not have a generating vector with n larger than " << std::to_string(preferred_n) << ", using largest generating vector with size " << std::to_string(n) << "." << std::endl;
        } else {
            n = generatingvectors.lower_bound(preferred_n)->first;
        }

        // Check n satisfies requirements of mod_mul implementation
        if ( n >= std::numeric_limits<typename std::make_signed<U>::type>::max() )
            throw std::domain_error("Qmc integrator called with n larger than the largest finite value representable with the signed type corresponding to U. Please decrease minn or use a larger unsigned integer type for U.");
        if ( n >= std::pow(std::numeric_limits<D>::radix,std::numeric_limits<D>::digits-1) )
            throw std::domain_error("Qmc integrator called with n larger than the largest finite value representable by the mantissa.");

        return n;
    };

    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    template <typename I>
    result<T> Qmc<T,D,M,P,F,G,H>::integrate_no_fit_no_transform(I& func)
    {
        if ( func.number_of_integration_variables < 1 )
            throw std::invalid_argument("qmc::integrate called with func.number_of_integration_variables < 1. Check that your integrand depends on at least one variable of integration.");
        if ( func.number_of_integration_variables > M )
            throw std::invalid_argument("qmc::integrate called with func.number_of_integration_variables > M. Please increase M (maximal number of integration variables).");
        if ( minm < 2 )
            throw std::domain_error("qmc::integrate called with minm < 2. This algorithm can not be used with less than 2 random shifts. Please increase minm.");
        if ( maxmperpackage < 2 )
            throw std::domain_error("qmc::integrate called with maxmperpackage < 2. This algorithm can not be used with less than 2 concurrent random shifts. Please increase maxmperpackage.");
        if ( maxnperpackage == 0 )
            throw std::domain_error("qmc::integrate called with maxnperpackage = 0. Please set maxnperpackage to a positive integer.");
        if ( devices.size() == 1 && devices.count(-1) == 1 && cputhreads == 0 )
            throw std::domain_error("qmc::evaluate called with cputhreads = 0 and devices = {-1} (CPU). Please set cputhreads to a positive integer or specify at least one non-CPU device.");
#ifndef __CUDACC__
        if ( devices.size() != 1 || devices.count(-1) != 1)
            throw std::invalid_argument("qmc::integrate called with devices != {-1} (CPU) but CUDA not supported by compiler.");
#endif

        if (verbosity > 2)
            logger << "-- qmc::integrate called --" << std::endl;

        std::vector<result<T>> previous_iterations; // keep track of the different interations
        U n = get_next_n(minn); // get next available n >= minn
        U m = minm;
        result<T> res;
        do
        {
            if (verbosity > 1)
                logger << "iterating" << std::endl;
            res = sample(func,n,m, previous_iterations);
            if (verbosity > 1)
                logger << "result " << res.integral << " " << res.error << std::endl;
            update(res,n,m);
        } while  ( integrators::overloads::compute_error_ratio(res, epsrel, epsabs, errormode) > static_cast<D>(1) && res.evaluations < maxeval );
        if (verbosity > 0)
        {
            logger << "-- qmc::integrate returning result --" << std::endl;
            logger << "integral " << res.integral << std::endl;
            logger << "error " << res.error << std::endl;
            logger << "n " << res.n << std::endl;
            logger << "m " << res.m << std::endl;
            logger << "iterations " << res.iterations << std::endl;
            logger << "evaluations " << res.evaluations << std::endl;

        }
        return res;
    };

    /*
     * implementation of integrate
     */
    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    template <typename I>
    result<T> Qmc<T,D,M,P,F,G,H>::integrate(I& func)
    {
        typename F<I,D,M>::transform_t fitted_func = fit(func);
        P<typename F<I,D,M>::transform_t,D,M> transformed_fitted_func(fitted_func);
        return integrate_no_fit_no_transform(transformed_fitted_func);
    };

    template <typename T, typename D, U M, template<typename,typename,U> class P, template<typename,typename,U> class F, typename G, typename H>
    Qmc<T,D,M,P,F,G,H>::Qmc() :
    logger(std::cout), randomgenerator( G( std::random_device{}() ) ), minn(8191), minm(32), epsrel(0.01), epsabs(1e-7), maxeval(1000000), maxnperpackage(1), maxmperpackage(1024), errormode(integrators::ErrorMode::all), cputhreads(std::thread::hardware_concurrency()), cudablocks(1024), cudathreadsperblock(256), devices({-1}), generatingvectors(integrators::generatingvectors::cbcpt_dn1_100()), verbosity(0), batching(false), evaluateminn(100000), fitstepsize(10), fitmaxiter(40), fitxtol(3e-3), fitgtol(1e-8), fitftol(1e-8), fitparametersgsl({})
    {
        // Check U satisfies requirements of mod_mul implementation
        static_assert( std::numeric_limits<U>::is_modulo, "Qmc integrator constructed with a type U that is not modulo. Please use a different unsigned integer type for U.");
        static_assert( std::numeric_limits<D>::radix == 2, "Qmc integrator constructed with a type D that does not have radix == 2. Please use a different floating point type for D.");
        
        if ( cputhreads == 0 )
        {
            cputhreads = 1; // Correct cputhreads if hardware_concurrency is 0, i.e. not well defined or not computable
            if (verbosity > 1) logger << "Qmc increased cputhreads from 0 to 1." << std::endl;
        }

#ifdef __CUDACC__
        // Get available gpus and add to devices
        int device_count = integrators::core::cuda::get_device_count();
        for(int i = 0; i < device_count; i++)
            devices.insert(i);
#endif

        fitparametersgsl = gsl_multifit_nlinear_default_parameters();
        fitparametersgsl.trs = gsl_multifit_nlinear_trs_lmaccel;
    };
};

#endif

#endif
