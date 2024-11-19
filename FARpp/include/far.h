/*
 *    Copyright 2024 C.S.Brady & H.Ratcliffe

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
#ifndef FAR_H
#define FAR_H
#include <array>
#include <algorithm>
#include <any>
#include <bit>
#include <bitset>
#include <cassert>
#include <charconv>
#include <climits>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include "callableTraits.h"
#include "ascii_map.h"
#include "demangle.h"
#include "debug.h"
#include "LBrefWrapper.h"
#include "OMPDynamicShared.h"

#ifdef _OPENMP
#include <omp.h>
#define FAR_OMP_BARRIER _Pragma("omp barrier")
#define FAR_OMP_SINGLE _Pragma("omp single")
#define FAR_OMP_CRITICAL _Pragma("omp critical")
#else
#define FAR_OMP_BARRIER
#define FAR_OMP_SINGLE
#define FAR_OMP_CRITICAL
#endif

#ifndef FAR_SIMPLELENGTH
#define FAR_SIMPLELENGTH 1
#endif

#ifndef FAR_NON_ET_LENGTH
#ifndef __clang__
#define FAR_NON_ET_LENGTH 1
#else
#define FAR_NON_ET_LENGTH 1
#endif
#endif

#ifndef FAR_ET_LENGTH
#ifndef __clang__
#define FAR_ET_LENGTH 1
#else
#define FAR_ET_LENGTH 1
#endif
#endif

#ifndef FAR_SIGNED_INDEX_TYPE
#define FAR_SIGNED_INDEX_TYPE int64_t
#endif

#ifndef FAR_UNSIGNED_INDEX_TYPE
#define FAR_UNSIGNED_INDEX_TYPE size_t
#endif

#ifndef FAR_DEFAULT_INDEX_BUILDER
	#if defined(FAR_CSTYLE) || defined(FAR_USE_C_INDEX)
		#define FAR_DEFAULT_INDEX_BUILDER CIndexBuilder
	#else
		#define FAR_DEFAULT_INDEX_BUILDER FortranIndexBuilder
	#endif
#endif

#ifndef FAR_DEFAULT_LB
	#if defined(FAR_CSTYLE) || defined(FAR_USE_C_BOUNDARY)
		#define FAR_DEFAULT_LB 0
	#else
		#define FAR_DEFAULT_LB 1
	#endif
#endif


//If on any compiler that supports the __HAS_INT128 macro use it
//Otherwise leave it unset
#if defined(__SIZEOF_INT128__) && !defined FAR_INT128
#define FAR_INT128 __int128
#endif

#ifndef FAR_UNUSED
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_LLVM_COMPILER)
#define FAR_UNUSED [[maybe_unused]]
#else
#define FAR_UNUSED
#endif
#endif
#ifndef FAR_INLINE_PREFIX
#define FAR_INLINE_PREFIX __attribute__((flatten))
#endif
#define FAR_INLINE FAR_INLINE_PREFIX inline
#ifndef FAR_PURE_PREFIX
#define FAR_PURE_PREFIX __attribute__((pure,flatten))
#endif
#define FAR_PURE FAR_PURE_PREFIX inline
#ifndef FAR_FLATTEN_PREFIX
#define FAR_FLATTEN_PREFIX __attribute__((flatten))
#endif
#define FAR_FLATTEN FAR_FLATTEN_PREFIX

#ifndef FAR_FUNCTION_PREFIX
#define FAR_FUNCTION_PREFIX 
#endif

//Define the default global RNG engine
#ifndef FAR_DEFAULT_GLOBAL_RNG
#define FAR_DEFAULT_GLOBAL_RNG std::mt19937 //Use Mersenne twister RNG (good but slower)
//#define FAR_DEFAULT_GLOBAL_RNG std::minstd_rand //Use C minstd RNG (not good but fast)
#endif

#ifdef _OPENMP
#define FAR_NEST_LOCK
#endif

namespace far
{
	using namespace ompdynamicshared;
	static const FAR_SIGNED_INDEX_TYPE defaultLB = FAR_DEFAULT_LB;
	namespace{
		//inline static sharedVariable<bool> workshare=false;
		//inline bool workshare=false;
		inline std::map<int,bool> workshare;
		inline bool getTeamNum(){
			#ifdef _OPENMP
			return omp_get_team_num();
			#else
			return 0;
			#endif
		}
	};

	//Only the types of these variables are important
	constexpr auto defaultFunction = (int*)nullptr;
	constexpr auto nullFunction = (void*)nullptr;

	// Define a concept to check for begin and end methods
	template<typename T>
	concept isIterable = requires(T t) {
 	   { t.begin() } -> std::input_iterator;
  	  { t.end() } -> std::input_iterator;
	  { t.size() } -> std::integral;
	};

	using defaultFunctionType = int*;
	using nullFunctionType = void*;

	inline void beginWorkshare(){
#ifdef _OPENMP
#pragma omp single
		workshare[getTeamNum()] = true;
#else 
		workshare[0] = false;
#endif
	}
	inline void endWorkshare(){
#ifdef _OPENMP
#pragma omp single
		workshare[getTeamNum()]=false;
#endif
	}

	//Classes used to represent generic instances of types
	//Only used in compile time type resolution
	struct array_base
	{
	};
	struct array_operator : public array_base
	{
	};
	struct array_actual : public array_base
	{
	};

	//Enumerator for bounds check state
	enum bounds_check_state
	{
		bc_never = false,
#ifndef FAR_DEFAULT_BOUNDS_CHECK
		bc_default = false,
#else
		bc_default = true,
#endif
		bc_always = true
	};

	/**
	 * Enum for iteration actions
	 */
	enum iter_action
	{
		iaContinue = 0, /**< Continue with iterations */
		iaTerminate=1 /**< Terminate */
	};

	/**
	 * Enum for reduction directions
	 */
	enum reductionDir{
		rdNone=-2, /**< Perform no reduction */
		rdAll=-1 /**< Reduce in all directions */
	};

	/**
	 * Enum for index actions
	 */
	enum getIndexAction{
		giToDo = -1, /**< Need to auto generate an indexer */
		giExternal = -2 /**< An external indexer has been generated */
			// Any other value is the level at which an indexer has been found
	};

	/**
	 * Enum class for return types
	 */
	enum class ert{
		Auto=1, /**< Automatically determine the return type */
		Lazy=2, /**< Return a lazy array */
		Array=3, /**< Return a real array */
	};

	inline int getThreadCount() {
		int numThreads = 1;
#ifdef _OPENMP
		numThreads = omp_get_num_threads();
#endif
		return numThreads;
	}

	inline int getThreadID() {
		int threadId = 0;
#ifdef _OPENMP
		threadId = omp_get_thread_num();
#endif
		return threadId;
	}

	/**
	 * Function to access an element of an array using default lower bounds
	 */
	template<typename T_array, typename T_tuple>
		FAR_INLINE decltype(auto) callDefaultBase(const T_array &array, const T_tuple &tuple){
			return std::apply([&array](auto &&...args) -> decltype(auto)
					{  return array.defaultBase(args...); }, tuple);
		}

	/**
	 *Class representing a ranged,strided subsection of an array
	 */
	struct Range
	{
		FAR_SIGNED_INDEX_TYPE lb, ub, stride;
		Range() : lb(std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::min()), ub(std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::max()), stride(1) {}
		Range(FAR_SIGNED_INDEX_TYPE value) : lb(value), ub(value), stride(1){}
		Range(FAR_SIGNED_INDEX_TYPE lb, FAR_SIGNED_INDEX_TYPE ub) : lb(lb), ub(ub), stride(1) {}
		Range(FAR_SIGNED_INDEX_TYPE lb, FAR_SIGNED_INDEX_TYPE ub, FAR_SIGNED_INDEX_TYPE stride) : lb(lb), ub(ub), stride(stride) {}
	};

	/**
	 *Return a range from the lower bound to the specified upper bound
	 */
	inline Range fromLB(FAR_SIGNED_INDEX_TYPE ub, FAR_SIGNED_INDEX_TYPE stride = 1) { return Range(std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::min(), ub, stride); }
	/**
	 *Return a range from the specified upper bound to the specified lower bound, treating the stride as negative
	 */
	inline Range toLB(FAR_SIGNED_INDEX_TYPE ub, FAR_SIGNED_INDEX_TYPE stride = 1) { return Range(std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::min(), ub, -stride); }
	/**
	 *Return a range from a specified lower bound to the upper bound
	 */
	inline Range toUB(FAR_SIGNED_INDEX_TYPE lb, FAR_SIGNED_INDEX_TYPE stride = 1) { return Range(lb, std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::max(), stride); }
	/**
	 *Return a range from the upper bound to a specified lower bound, treating the stride as negative
	 */
	inline Range fromUB(FAR_SIGNED_INDEX_TYPE lb, FAR_SIGNED_INDEX_TYPE stride = 1) { return Range(lb, std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::max(), -stride); }
	/**
	 *Whole range from lower bound to upper bound
	 */
	inline Range LBToUB(FAR_SIGNED_INDEX_TYPE stride = 1) { return Range(std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::min(), std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::max(), stride); }
	/**
	 *Whole range reversed from upper bound to lower bound, treating the stride as negative
	 */
	inline Range UBToLB(FAR_SIGNED_INDEX_TYPE stride = 1) { return Range(std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::min(), std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::max(), -stride); }


	/**
	 * Check if a type is a fxed size char array
	 */
	template<typename T>
		struct isFixedSizeCharArray : std::false_type {constexpr static FAR_UNSIGNED_INDEX_TYPE size=0;};

	// Partial specialization for fixed size char array
	template<std::FAR_UNSIGNED_INDEX_TYPE N>
		struct isFixedSizeCharArray<char[N]> : std::true_type {constexpr static FAR_UNSIGNED_INDEX_TYPE size=N;};

	/**
	 *Check whether or not a type is a complex number
	 */
	template <typename T>
		struct is_complex_core : std::false_type
	{
	};
	template <typename T>
		struct is_complex_core<std::complex<T>> : std::true_type
		{
		};

	/**
	 * Check whether something is anything derived from std::basic_string
	 */
	template<typename T>
		struct is_stdstring : std::false_type {
			using type = T;
		};

	template<typename T>
		struct is_stdstring<std::basic_string<T>> : std::true_type {
			using type = std::add_pointer_t<T>;
		};

	template<typename T>
		constexpr bool is_stdstring_v = is_stdstring<T>::value;

	template<typename T>
		using is_stdstring_t = typename is_stdstring<T>::type;

	template<typename T>
		using refWrapifNeeded = std::conditional_t<std::is_reference_v<T>,LBrefWrapper<std::remove_reference_t<T>>,T>;

	/**
	 * Template to create a version of a function where the return value is passed as a parameter
	 */
	template<auto callable, typename T_params, auto reduce=defaultFunction>
		void wrapReturn(refWrapifNeeded<typename callableTraits<decltype(callable)>::type> &ret, T_params params){
			if constexpr(callableTraits<decltype(reduce)>::value){
				reduce(ret, std::apply(callable,params));
			} else {
				ret = std::apply(callable,params);
			}
		}

	/**
	 *Check whether or not a type is a logical type
	 */
	template <typename T>
		struct is_logical_core : std::false_type
	{
	};
	template <>
		struct is_logical_core<bool> : std::true_type
		{
		};

	/**
	 *Remove reference, pointer, and cv_qualification
	 */
	template <typename T>
		using clean = std::remove_pointer_t<std::decay_t<T>>;
	/**
	 *Get index type from an object
	 */
	template <typename T>
		using index_type = typename clean<T>::index_type;

	/**
	 * Check whether an array is std::array (used here for constant arrays
	 */
	template<typename T>
		struct isStdArray : std::false_type {};

	template<typename T, std::FAR_UNSIGNED_INDEX_TYPE N>
		struct isStdArray<std::array<T, N>> : std::true_type {};


	/**
	 *Check whether a given type is an array object
	 */
	template <typename T, typename = void>
		struct arrayInfo
		{
			private:
			static constexpr bool isLBWrapper = isLBrefWrapper<T>::value;
			using unwrappedType = isLBrefWrapper<std::decay_t<T>>::type;
			public:
			using qualifiedType = unwrappedType;
			using wrappedType = T;
			using type = std::decay_t<qualifiedType>;
			using referenceType = std::conditional_t<std::is_reference_v<qualifiedType>,type&,type>;
			using unrefType = std::remove_reference_t<qualifiedType>;
			//If T is a value type then the storage type is a value
			//If it is a pointer or reference then it is a pointer with appropriate constness
			using storageType = std::conditional_t<std::is_reference_v<T>,
						std::conditional_t<std::is_const_v<std::remove_reference_t<T>>,const clean<T> *,clean<T>*>,T>;
			using index = void;//std::remove_reference_t<T>;
			static constexpr int rank = 0;
			static constexpr int order = 0;
			static constexpr bool value = false;
			static constexpr bounds_check_state bounds_check = bounds_check_state::bc_default;
			static constexpr bool forceContiguous = false;
		};

	template <typename T>
		class arrayInfo<T, std::enable_if_t<std::is_base_of<array_base, std::decay_t<T>>::value>>
		{
			private:
				using uqcore = typename std::decay_t<T>::core_type;
				static constexpr bool isLBWrapper = isLBrefWrapper<uqcore>::value;
				using unwrappedType = isLBrefWrapper<std::decay_t<uqcore>>::type;
			public:
				using type = unwrappedType;
				using wrappedType = uqcore;
				using qualifiedType = typename clean<T>::core_type;
				using referenceType = std::conditional_t<std::is_reference_v<T>,type&,type>;
				using unrefType = typename std::remove_reference_t<typename clean<T>::core_type>;
				//Possible future idea - store lazy arrays by value and real arrays by reference?
				using storageType = std::conditional_t<std::is_const_v<std::remove_reference_t<T>>,const clean<T> *,clean<T>*>;
				using index = typename clean<T>::index_type;
				static constexpr int rank = clean<T>::core_rank;
				static constexpr int order = clean<T>::index_type::core_order;
				static constexpr bool value = true;
				static constexpr bounds_check_state bounds_check = clean<T>::core_bounds_check;
				static constexpr bool forceContiguous = clean<T>::coreForceContiguous;
		};

	template <typename... T>
	constexpr int maxRank = std::max({arrayInfo<T>::rank...});


	template<typename T>
	using arrayInfo_t = typename arrayInfo<T>::type;
	template<typename T>
	constexpr bool arrayInfo_v = arrayInfo<T>::value;
	template<typename T>
	constexpr int arrayInfo_r = arrayInfo<T>::rank;

	template<typename T>
		concept anyArray = arrayInfo<T>::value;

	template<typename T>
		concept notArray = !arrayInfo<T>::value;

	template<typename T>
		concept anyRealArray = arrayInfo<T>::value && std::is_floating_point_v<typename arrayInfo<T>::type>;

	template<typename T>
		concept anyIntegerArray = arrayInfo<T>::value && std::is_integral_v<typename arrayInfo<T>::type>;

	template<typename T, typename... T_inner>
		concept anyRankArray = std::disjunction_v<std::is_same<typename arrayInfo<T_inner>::type, typename arrayInfo<T>::type>..., std::bool_constant<arrayInfo<T_inner>::value>...>;

	template<typename T, typename... T_inner>
		concept anyRankConvertibleArray = std::disjunction_v<std::is_convertible<typename arrayInfo<T_inner>::type, typename arrayInfo<T>::type>..., std::bool_constant<arrayInfo<T_inner>::value>...>;

	template<typename T, int rank>
		concept anyTypeArray = arrayInfo<T>::value && (arrayInfo<T>::rank==rank);

	template<typename T, int rank, typename... T_inner>
		concept manyTypeArrayParameter = std::disjunction_v<std::is_same<typename arrayInfo<T_inner>::type, typename arrayInfo<T>::type>..., std::bool_constant<arrayInfo<T_inner>::value>...,std::bool_constant<arrayInfo<T_inner>::rank==rank>...>;;

	template<typename T,typename T_inner, int rank>
		concept arrayParameter = std::is_same_v<T_inner, typename arrayInfo<T>::type> && arrayInfo<T>::value && arrayInfo<T>::rank==rank;

	template<typename T, int rank>
		concept rankIntegerArray = arrayInfo<T>::value && (arrayInfo<T>::rank==rank && std::is_integral_v<typename arrayInfo<T>::type>);

	template<typename T, int rank>
		concept rankRealArray = arrayInfo<T>::value && (arrayInfo<T>::rank==rank && std::is_floating_point_v<typename arrayInfo<T>::type>);


	template<typename T, typename... T_inner>
		concept elementalParameter = std::disjunction_v<std::is_same<typename arrayInfo<T_inner>::type, typename arrayInfo<T>::type>...>;

	template<typename... Args>
		constexpr FAR_UNSIGNED_INDEX_TYPE countRanges() {
			return (... + std::is_same_v<Args, Range>);
		}

	template <typename... Ts>
		struct first_arrayInfo;

	template <typename T, typename... Ts>
		struct first_arrayInfo<T, Ts...> {
			using info = std::conditional_t<arrayInfo<T>::value, arrayInfo<T>, typename first_arrayInfo<Ts...>::info>;
		};

	template <typename T>
		struct first_arrayInfo<T> {
			using info = arrayInfo<T>;
		};
	//namespace {

	/**
	 * Return the first index object from a parameter pack
	 * Return as a pointer and return nullptr if no index object is found
	 */
	template <typename T, typename... Ts>
		constexpr auto getFirstIndex(T &&t, Ts &&... ts) {
			if constexpr (arrayInfo<T>::value) {
				return &t;
			} else if constexpr (sizeof...(Ts) > 0) {
				return getFirstIndex(std::forward<Ts>(ts)...);
			}
		}

	/** 
	 * Create an N element tuple of a given type
	 * use as N_ary_tuple_type<int,4> to create std::tuple<int,int,int,int>
	 */
	template <typename T, int N, typename... REST>
		struct N_ary_tuple_type
		{
			typedef typename N_ary_tuple_type<T, N-1, T, REST...>::type type;
		};
	template <typename T, typename... REST>
		struct N_ary_tuple_type<T, 0, REST...>
		{
			typedef std::tuple<REST...> type;
		};

		template <typename T, int N>
		using N_ary_tuple_type_t = typename N_ary_tuple_type<T, N>::type;


	/** 
	 * Create an N level pointer
	 * use as N_ary_pointer<int,4> to create int****
	 */
	template <typename T, int N>
		struct N_ary_pointer
		{
			using type = std::add_pointer_t<typename N_ary_pointer<T, N-1>::type>;
		};
	template <typename T>
		struct N_ary_pointer<T, 1>
		{
			using type = std::add_pointer_t<T>;
		};

		template <typename T, int N>
		using N_ary_tuple_type_t = typename N_ary_tuple_type<T, N>::type;

		template<int Ti_rank>
			struct indexInfo {
				constexpr static int rank = Ti_rank;
				using type = typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type;
				using reduced_type = std::conditional_t<(rank>1),typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,std::max(rank-1,0)>::type,std::tuple<double>>;
				type &indices; //Index with default base
				type &SBindices; //Index with specific base
				int direction;
				reduced_type reduced_indices{};
				indexInfo( type &idc) : indices(idc), SBindices(idc), direction(-1) {;}
				indexInfo( type &idc, int direction) : indices(idc), SBindices(idc), direction(direction) {;}
				indexInfo( type &idc,  type &sbidc, int direction) : indices(idc), SBindices(sbidc), direction(direction) {;}
				indexInfo( type &idc,  reduced_type &rti, int direction) : indices(idc),  SBindices(idc),
				direction(direction), reduced_indices(rti) {;}
				indexInfo( type &idc,  type &sbidc, reduced_type &rti, int direction) : indices(idc),  SBindices(sbidc),
				direction(direction), reduced_indices(rti) {;}
			};

			template<typename T>
			struct isReductionIndexInfo : std::false_type {
				static constexpr int rank = 0;
			};

			template<int Ti_rank>
			struct isReductionIndexInfo<indexInfo<Ti_rank>> : std::true_type {
				static constexpr int rank = Ti_rank;
			};

			//This is a helper struct to find which parameters are arrays in the set of parameters
			//that the elemental function is called with, but are not arrays in the inner function
			template<int level, typename T, auto callable>
				struct is_unmatched_array{
					using callable_params = typename callableTraits<decltype(callable)>::params;
					static constexpr bool value = !std::is_convertible_v<clean<T>,clean<std::tuple_element_t<level,callable_params>>> 
						&& arrayInfo<T>::value;
				};

			/**
			 * Function to build either a reference to an element of an array,
			 * a reference to an entire variable (array or otherwise), or a copy of a value
			 * if the input parameter is an rvalue
			 */
			template<auto callable, int level, typename T_info, typename T>
				static decltype(auto) buildRef(T_info& info, T& current){
					auto& itpl = info.indices;
					auto& direction = info.direction;
					constexpr bool deref = [](){if constexpr(level<0 && arrayInfo<T>::value){
						return true; //Manually specifying an array
					} else if constexpr(level<0){
						return false; //Manually specifying a non-array
					} else if constexpr(is_unmatched_array<level,T,callable>::value){
						return true;//Automatically detecting an array
					} else { 
						return false;} //Automatically detecting an non-array
					}();
					if constexpr(deref){
						auto &&v = ([&current, &itpl, direction]() -> decltype(auto){
								if constexpr(arrayInfo<T>::rank==std::tuple_size_v<std::decay_t<decltype(itpl)>>){
									FAR_UNUSED auto capture = direction;//Silence warnings from ifx
								return callDefaultBase(current,itpl);
								} else {
								typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,arrayInfo<T>::rank>::type cutTuple;
								reduceTuple(itpl, cutTuple, direction);
								return callDefaultBase(current,cutTuple);
								}
								})();
						if constexpr(std::is_rvalue_reference_v<decltype(v)>){
							return std::decay_t<decltype(v)>(v);
						} else {
							return std::ref(v);
						}
					} else {
						if constexpr(std::is_rvalue_reference_v<decltype(current)>){
							return std::decay_t<decltype(current)>(current);
						} else {
							return std::ref(current);
						}
					}
				}


			/**
			 * Function to build the tuple of parameters that is needed to call the inner function
			 */
			template<auto callable, int level, typename T_info, typename T, typename... T_others>
				static auto buildCallTuple(T_info& rri, T&& current, T_others&&... others){
					using callable_params = typename callableTraits<decltype(callable)>::params;
					auto &itpl = rri.indices;
					auto &direction = rri.direction;
					auto val = buildRef<callable, level>(rri, current);
					constexpr int rank = std::tuple_size_v<std::decay_t<decltype(itpl)>>;
					if constexpr(level < std::tuple_size_v<callable_params>-1){
						if constexpr(isReductionIndexInfo<std::decay_t<std::tuple_element_t<level, callable_params>>>::value){
							return std::tuple_cat(std::make_tuple(std::ref(rri)),buildCallTuple<callable,level+1>(
										rri, current, others...));
						} else {
							if constexpr(sizeof...(others)>0) {
								return std::tuple_cat(std::make_tuple(val),buildCallTuple<callable,level+1>(rri, others...));
							} else {
								//Deal with case where indexInfo is last parameter - use a dummy for current
								return std::tuple_cat(std::make_tuple(val),buildCallTuple<callable,level+1>(rri, direction, 1));
							}
						}
					} else {
						if constexpr(std::is_same_v<std::decay_t<std::tuple_element_t<level, callable_params>>,
								indexInfo<rank>>){
							return std::make_tuple(std::ref(rri));
						} else {
							return std::make_tuple(val);
						}
					}
				}

	//}

	/**
	 * Check if a callable has a return type and all parameters are either values or const references
	 */
	template<auto callable,int level=std::tuple_size_v<typename callableTraits<decltype(callable)>::params>-1>
	struct paramsCanChange{
		using param = std::tuple_element_t<level,typename callableTraits<decltype(callable)>::params>;
		static constexpr int q=1;
		static constexpr bool isref = std::is_reference_v<param>;
		static constexpr bool isconst = std::is_const_v<std::remove_reference_t<param>>;
		static constexpr bool value = (!isReductionIndexInfo<std::decay_t<param>>::value && //Is not a reduction index info
			(
				(std::is_reference_v<param> && !std::is_const_v<std::remove_reference_t<param>>) // Is non const reference
		 	   || (std::is_pointer_v<param> && !std::is_const_v<std::remove_pointer_t<param>>) // Is non const pointer
			))
			|| paramsCanChange<callable,level-1>::value; //Check the next parameter
	};

	//Specialise for the end of the tuple
	template<auto callable>
	struct paramsCanChange<callable,-1>{
		using param = void;
		static constexpr int q=0;
		static constexpr bool value = false;
		static constexpr bool isref = false;
		static constexpr bool isconst = false;
	};

	template<auto callable, int level=0, typename ...T_other>
		struct arraySetInfo;


	template<auto callable, int level, typename T, typename... T_others>
		struct arraySetInfo<callable, level, T, T_others...> {
			template<auto callable2, int level2, typename... T_other2>
				friend struct arraySetInfo;
			private:
			using callableType = std::decay_t<std::tuple_element_t<level,typename callableTraits<decltype(callable)>::params>>;
			static constexpr bool isRRI = isReductionIndexInfo<callableType>::value;
			using nextType = std::conditional_t<isRRI,
				arraySetInfo<callable, level+1, T, T_others...>,
				arraySetInfo<callable, level+1, T_others...>
				>;
			static constexpr bool hasRRI = isRRI || nextType::hasRRI;
			static constexpr bool isUnmatchedThisLevel = (arrayInfo<T>::value && 
				!std::is_same_v<std::decay_t<T>, callableType>) || hasRRI;
			static constexpr int maxUnmatchedRank = isUnmatchedThisLevel ? std::max(arrayInfo<T>::rank,nextType::maxUnmatchedRank) : nextType::maxUnmatchedRank;
			static constexpr int maxRank = std::max(arrayInfo<T>::rank,nextType::maxRank);
			static constexpr int minUnmatchedRank = isUnmatchedThisLevel ? (nextType::minUnmatchedRank>0?std::min(arrayInfo<T>::rank,nextType::minUnmatchedRank):arrayInfo<T>::rank) : nextType::minUnmatchedRank;
			static constexpr bool isGoodThisLevel = ((arrayInfo<T>::value && 
				std::is_convertible_v<std::decay_t<T>, callableType>) //Is matched array
				|| (arrayInfo<T>::value && std::is_convertible_v<typename arrayInfo<T>::type, callableType>)) //Is convertible array
				|| (std::is_convertible_v<T, callableType>) //Is convertible scalar 
				|| hasRRI; //Is a reduction index info
			public:
			static constexpr bool anyUnmatched = isUnmatchedThisLevel || nextType::anyUnmatched;
			static constexpr bool isGood = (maxUnmatchedRank == minUnmatchedRank) && isGoodThisLevel && nextType::isGood;
			static constexpr int foundLevel = isUnmatchedThisLevel ? level : nextType::foundLevel;
			static constexpr int badLevel = !isGoodThisLevel?level:nextType::badLevel;
			using info = std::conditional_t<isUnmatchedThisLevel && arrayInfo<T>::rank == (hasRRI?maxRank:maxUnmatchedRank) , arrayInfo<T>, typename nextType::info>;	
		};

	template<auto callable, int level, typename T>
		struct arraySetInfo<callable, level, T> {
			template<auto callable2, int level2, typename... T_other>
				friend struct arraySetInfo;
			private:
			using callableType = std::decay_t<std::tuple_element_t<level,typename callableTraits<decltype(callable)>::params>>;
			static constexpr bool isRRI = isReductionIndexInfo<callableType>::value;
			static constexpr bool hasRRI = isRRI;
			static constexpr bool isUnmatchedThisLevel = (arrayInfo<T>::value && 
				!std::is_same_v<std::decay_t<T>, callableType>) || isRRI;
			static constexpr int maxUnmatchedRank = isUnmatchedThisLevel ? arrayInfo<T>::rank : 0;
			static constexpr int maxRank = arrayInfo<T>::rank;
			static constexpr int minUnmatchedRank = isUnmatchedThisLevel ? arrayInfo<T>::rank : 0;
			static constexpr bool isGoodThisLevel = ((arrayInfo<T>::value && 
				std::is_same_v<std::decay_t<T>, callableType>) //Is matched array
				|| (arrayInfo<T>::value && std::is_convertible_v<typename arrayInfo<T>::type, callableType>)) //Is convertible array
				|| (std::is_convertible_v<T, callableType>) //Is convertible scalar 
				|| isRRI; //Is a reduction index info
			public:
			static constexpr bool anyUnmatched = isUnmatchedThisLevel;
			static constexpr bool isGood = isGoodThisLevel;
			static constexpr int foundLevel = isUnmatchedThisLevel ? level : -1;
			static constexpr int badLevel = !isGoodThisLevel?level:-1;
			using info = std::conditional_t<isUnmatchedThisLevel, arrayInfo<T>, arrayInfo<int>>;
		};
	

	//Get the largest and smallest ranks
	template <typename... Ts>
		struct arrayRankRange;

	template <typename T, typename... Ts>
		struct arrayRankRange<T, Ts...> {
			constexpr static int minRank = std::min(arrayInfo<T>::rank, arrayRankRange<Ts...>::minRank);
			constexpr static int maxRank = std::max(arrayInfo<T>::rank, arrayRankRange<Ts...>::maxRank);
		};

	template <typename T>
		struct arrayRankRange<T> {
			constexpr static int minRank = arrayInfo<T>::rank;
			constexpr static int maxRank = arrayInfo<T>::rank;
		};


		/**
		 * @brief Convert a parameter pack to a tuple of the same types and values
		 * @tparam order The order in which the parameters are packed into the tuple
		 * @tparam T The type of the first parameter
		 * @tparam Args The types of the remaining parameters
		 * @param t The first parameter
		 * @param args The remaining parameters
		 * @return A tuple containing the parameters
		 */
		template<int order, typename T, typename... Args>
		FAR_PURE auto packToTuple(T && t, Args&&... args)
		{
			if constexpr(sizeof...(args)==0){
				return std::make_tuple(t);
			} else {
				if constexpr(order==0){
					return std::tuple_cat(std::make_tuple(t),packToTuple<order>(std::forward<Args>(args)...));
				} else {
					return std::tuple_cat(packToTuple<order>(std::forward<Args>(args)...),std::make_tuple(t));
				}
			}
		}


	/**
	 * @brief Function implementing the Fortran ishftc function for a single value of any integer type
	 * @tparam T The type of the input value
	 * @tparam T_shift The type of the shift value
	 * @tparam T_size The type of the size value
	 * @param i The input value
	 * @param shift The shift value
	 * @param size The size value
	 * @return The shifted value
	 */
	template<typename T, typename T_shift, typename T_size>
		inline T ishftc_core(T i, T_shift shift, T_size size = sizeof(arrayInfo<T>::type)*CHAR_BIT) {
			//I would be unsurprised to find that there is a cleverer solution than this
			//but this works
			using uT = std::make_unsigned_t<T>;
			uT val = static_cast<uT>(i);
			uT shift_mask = size == sizeof(T)*CHAR_BIT ? ~0 : (1 << size) - 1;
			uT unshift_mask = ~shift_mask;
			return static_cast<T>((((shift>=0)?((val << shift) | ((val) >> (size - shift))):((val) >> std::abs(shift)) | (val << (size - std::abs(shift))))& shift_mask) | (val & unshift_mask));
		}

	/**
	 * @brief Function implementing the Fortran shiftl function for a single value of any integer type
	 * @tparam T The type of the input value
	 * @tparam T_shift The type of the shift value
	 * @param i The input value
	 * @param shift The shift value
	 * @return The shifted value
	 */
	template<typename T, typename T_shift>
		inline T shiftl_core(T i, T_shift shift){
			using unsigned_T = std::make_unsigned_t<T>;
			//If the source is unsigned get a reference to it
			//Otherwise copy the data into an unsigned value and get the reference to that
			unsigned_T rValue;

			if constexpr(std::is_unsigned_v<T>){
				rValue = i;
			} else {
				std::memcpy(&rValue,&i,sizeof(T));
			}
			
			unsigned_T result = (rValue<<shift);//Doing this shift with an unsigned integer automatically gives the right behaviour

			//If the input type is unsigned then return it as unsigned
			//Otherwise copy the unsigned value into a signed value and return it
			if constexpr(std::is_unsigned_v<T>){
				return result;
			} else {
				T finalResult;
				std::memcpy(&finalResult,&result,sizeof(T));
				return finalResult;
			}
		}

	/**
	 * @brief Function implementing the Fortran shiftr function for a single value of any integer type
	 * @tparam T The type of the input value
	 * @tparam T_shift The type of the shift value
	 * @param i The input value
	 * @param shift The shift value
	 * @return The shifted value
	 */
	template<typename T, typename T_shift>
		inline T shiftr_core(T i, T_shift shift){
			using unsigned_T = std::make_unsigned_t<T>;
			//If the source is unsigned get a reference to it
			//Otherwise copy the data into an unsigned value and get the reference to that
			unsigned_T rValue;
			if constexpr(std::is_unsigned_v<T>){
				rValue = i;
			} else {
				std::memcpy(&rValue,&i,sizeof(T));
			}
			
			unsigned_T result = (rValue>>shift);//Doing this shift with an unsigned integer automatically gives the right behaviour

			//If the input type is unsigned then return it as unsigned
			//Otherwise copy the unsigned value into a signed value and return it
			if constexpr(std::is_unsigned_v<T>){
				return result;
			} else {
				T finalResult;
				std::memcpy(&finalResult,&result,sizeof(T));
				return finalResult;
			}
		}

	/**
	 * @brief Check if a type is a pointer even if it is wrapped in a reference
	 * @tparam T The type to check
	 */
	template<typename T>
		static const bool is_true_pointer_v = std::is_pointer_v<std::decay_t<T>>;

	/**
	 * @brief Check if at least one of two types is an array
	 * @tparam T_lhs The first type
	 * @tparam T_rhs The second type
	 */
	template<typename T_lhs, typename T_rhs>
		static const bool oneSideArray = std::is_base_of_v<array_base, clean<T_lhs>> || std::is_base_of_v<array_base, clean<T_rhs>>;

	/**
	 * @brief Find the type of the an index for one of two types
	 * @tparam T_lhs The first type
	 * @tparam T_rhs The second type
	 */
	template<typename T_lhs, typename T_rhs>
		using firstIndex = std::conditional_t<arrayInfo<T_lhs>::value, typename arrayInfo<T_lhs>::index, typename arrayInfo<T_rhs>::index>;

	/**
	 * @brief Check if a type is an range object or an integral type
	 * @tparam T The type to check
	 */
	template <typename T>
		struct rangeOrIntegral{
			static const bool value = std::is_same_v<T,Range> || std::is_integral_v<T>;
		};

	template <typename... Ts>
		constexpr bool is_integral_all = std::conjunction<std::is_integral<typename arrayInfo<std::decay_t<Ts>>::type>...>::value;
	template <typename... Ts>
		constexpr bool is_integral_any = std::disjunction<std::is_integral<typename arrayInfo<std::decay_t<Ts>>::type>...>::value;
	template<typename T>
		constexpr bool is_integral = is_integral_any<T>;

	template <typename... Ts>
		constexpr bool is_floating_point_all = std::conjunction<std::is_floating_point<typename arrayInfo<std::decay_t<Ts>>::type>...>::value;

	template <typename... Ts>
		constexpr bool is_floating_point_any = std::disjunction<std::is_floating_point<typename arrayInfo<std::decay_t<Ts>>::type>...>::value;

	template<typename T>
		constexpr bool is_floating_point = is_floating_point_any<T>;


	template <typename... Ts>
		constexpr bool is_arithmetic_all = std::conjunction<std::is_arithmetic<typename arrayInfo<std::decay_t<Ts>>::type>...>::value;

	template <typename... Ts>
		constexpr bool is_arithmetic_any = std::disjunction<std::is_arithmetic<typename arrayInfo<std::decay_t<Ts>>::type>...>::value;

	template<typename T>
		constexpr bool is_arithmetic = is_arithmetic_any<T>;


	template <typename... Ts>
		constexpr bool is_logical_all = std::conjunction<is_logical_core<typename arrayInfo<std::decay_t<Ts>>::type>...>::value;

	template <typename... Ts>
		constexpr bool is_logical_any = std::disjunction<is_logical_core<typename arrayInfo<std::decay_t<Ts>>::type>...>::value;

	template<typename T>
		constexpr bool is_logical = is_logical_any<T>;


	template <typename... Ts>
		constexpr bool is_complex_all = std::conjunction<is_complex_core<typename arrayInfo<std::decay_t<Ts>>::type>...>::value;

	template <typename... Ts>
		constexpr bool is_complex_any = std::disjunction<is_complex_core<typename arrayInfo<std::decay_t<Ts>>::type>...>::value;

	template<typename T>
		constexpr bool is_complex = is_complex_any<T>;


	template <typename... Ts>
		constexpr bool validSubscript_all = std::conjunction<rangeOrIntegral<typename arrayInfo<std::decay_t<Ts>>::type>...>::value;
	template <typename... Ts>
		constexpr bool validSubscript_any = std::disjunction<rangeOrIntegral<typename arrayInfo<std::decay_t<Ts>>::type>...>::value;
	template <typename T>
		constexpr bool validSubscript = validSubscript_all<T>;

	template<int rank, typename... Ts>
		constexpr bool is_correct_subscript = (sizeof...(Ts)==rank);// && validSubscript_all<Ts...>;



	template <typename... Ts>
		constexpr bool is_same = std::conjunction<
		std::is_same<typename arrayInfo<Ts>::type,typename arrayInfo<typename std::tuple_element<0,std::tuple<Ts...>>::type>::type>...>::value;

	template <typename... Ts>
		constexpr bool is_convertible = std::conjunction<
		std::is_convertible<typename arrayInfo<Ts>::type,typename arrayInfo<typename std::tuple_element<0,std::tuple<Ts...>>::type>::type>...>::value;


	template <typename T1, typename T2>
		constexpr bool isConformable = std::is_same_v<typename arrayInfo<T1>::type, typename arrayInfo<T2>::type> && (arrayInfo<T1>::rank == arrayInfo<T2>::rank);

	template <typename T, typename U>
		using larger_type = typename std::conditional<(sizeof(T) >= sizeof(U)), T, U>::type;

	template <typename T>
		using arrayInfo_st = typename arrayInfo<T>::storageType;

	template <typename T>
		using arrayInfo_t = typename arrayInfo<T>::type;

	/**
	 *Check whether a parameter pack contains any Range objects. Used to disambiguate
	 *slice and dereference operators
	 */
	template <typename... Ts>
		constexpr bool containsRange = std::disjunction<std::is_same<Ts, Range>...>::value;

	template <typename T, typename = void>
		struct is_integral_convertible : std::false_type {};

	template <typename T>
		struct is_integral_convertible<T, std::enable_if_t<std::is_convertible_v<T, int>>> : std::true_type {};

	template<typename T>
		constexpr inline bool is_integral_convertible_v = is_integral_convertible<T>::value;

	template <typename T, typename = void>
		struct is_real_convertible : std::false_type {};

	template <typename T>
		struct is_real_convertible<T, std::enable_if_t<std::is_convertible_v<T, float>>> : std::true_type {};

	template<typename T>
		constexpr inline bool is_real_convertible_v = is_real_convertible<T>::value;

	template <typename T, typename = void>
		struct is_bool_convertible : std::false_type {};

	template <typename T>
		struct is_bool_convertible<T, std::enable_if_t<std::is_convertible_v<T, bool>>> : std::true_type {};

	template <typename T>
		constexpr inline bool is_bool_convertible_v = is_bool_convertible<T>::value;

	template <typename T, typename = void>
		struct is_complex_convertible : std::false_type {};

	template <typename T>
		struct is_complex_convertible<T, std::enable_if_t<std::is_convertible_v<T, std::complex<float>> || 
		std::is_convertible_v<T, std::complex<double>>>> : std::true_type {};

	template <typename T>
		constexpr inline bool is_complex_convertible_v = is_complex_convertible<T>::value;

	template <typename T>
		struct is_string : std::false_type {};

	template <>
		struct is_string<std::string> : std::true_type {};

	template <>
		struct is_string<std::wstring> : std::true_type {};

	template <>
		struct is_string<std::u16string> : std::true_type {};

	template <>
		struct is_string<std::u32string> : std::true_type {};

	template <>
		struct is_string<char*> : std::true_type {};

	template <>
		struct is_string<const char*> : std::true_type {};

	template <>
		struct is_string<wchar_t*> : std::true_type {};

	template <>
		struct is_string<const wchar_t*> : std::true_type {};

	template<FAR_SIGNED_INDEX_TYPE N>
		struct is_string<char[N]> : std::true_type {};

	template<FAR_SIGNED_INDEX_TYPE N>
		struct is_string<const char[N]> : std::true_type {};

	template<FAR_SIGNED_INDEX_TYPE N>
		struct is_string<wchar_t[N]> : std::true_type {};

	template<FAR_SIGNED_INDEX_TYPE N>
		struct is_string<const wchar_t[N]> : std::true_type {};



	template <typename T>
		constexpr inline bool is_string_v = is_string<T>::value;

	template <typename T>
		struct is_string_convertible {
			static std::ostringstream os;
			static const bool value = std::is_convertible<T, decltype(os << std::declval<T>())>::value;
		};

	template <typename T>
		constexpr inline bool is_string_convertible_v = is_string_convertible<T>::value;

#define STRICT_TYPE

#ifdef STRICT_TYPE
	//Strict type checking is for exactly matching types
#define AERRCHECKREAL(T, FUNC, PARAM) static_assert(std::is_floating_point_v<typename arrayInfo<T>::type>,#PARAM" must be floating point in "#FUNC);
#define AERRCHECKINTEGER(T, FUNC, PARAM) static_assert(std::is_integral_v<typename arrayInfo<T>::type>,#PARAM" must be integer in "#FUNC);
#define AERRCHECKCOMPLEX(T, FUNC, PARAM) static_assert(is_complex<typename arrayInfo<T>::type>,#PARAM " must be complex in " #FUNC);
#define AERRCHECKREALINTEGER(T, FUNC, PARAM) static_assert(std::is_integral_v<typename arrayInfo<T>::type> || std::is_floating_point_v<typename arrayInfo<T>::type>,#PARAM" must be floating point or integer in "#FUNC);
#define AERRCHECKREALINTEGERCOMPLEX(T,FUNC,PARAM) static_assert(std::is_integral_v<typename arrayInfo<T>::type> || std::is_floating_point_v<typename arrayInfo<T>::type> || is_complex<typename arrayInfo<T>::type>,#PARAM " must be floating point, integer or complex in " #FUNC);
#define AERRCHECKREALCOMPLEX(T,FUNC,PARAM) static_assert(std::is_floating_point_v<typename arrayInfo<T>::type> || is_complex<typename arrayInfo<T>::type>,#PARAM" must be floating point or complex in "#FUNC);
#define AERRCHECKREALINTEGERSTRING(T, FUNC, PARAM) static_assert(std::is_integral_v<typename arrayInfo<T>::type> || std::is_floating_point_v<typename arrayInfo<T>::type> || is_string_v<T>,#PARAM" must be floating point, string or integer in "#FUNC);
#define AERRCHECKTRUTHY(T,FUNC,PARAM) static_assert(std::is_same_v<bool, typename arrayInfo<T>::type>,#PARAM" must be boolean in "#FUNC);
#define AERRCHECKSTRINGY(T,FUNC,PARAM) static_assert(is_string_v<std::decay_t<typename arrayInfo<T>::type>>,#PARAM" must be string in "#FUNC);

#elif !defined(NO_TYPES)
	//Convertible type checking checks if the types are capable of being converted to
	//ones that are known to work
#define AERRCHECKREAL(T, FUNC, PARAM) static_assert(is_real_convertible_v<typename arrayInfo<T>::type>,#PARAM" must be floating point convertible in "#FUNC);
#define AERRCHECKINTEGER(T, FUNC, PARAM) static_assert(is_integral_convertible_v<typename arrayInfo<T>::type>,#PARAM "must be integer convertible in "#FUNC);
#define AERRCHECKCOMPLEX(T, FUNC, PARAM) static_assert(is_complex_convertible_v<typename arrayInfo<T>::type>,#PARAM "must be complex convertible in "#FUNC);
#define AERRCHECKREALINTEGER(T, FUNC, PARAM) static_assert(is_integral_convertible_v<typename arrayInfo<T>::type> || is_real_convertible_v<typename arrayInfo<T>::type>,#PARAM" must be floating point or integer convertible in "#FUNC);
#define AERRCHECKREALINTEGERCOMPLEX(T,FUNC,PARAM) static_assert(is_integral_convertible_v<typename arrayInfo<T>::type> || is_real_convertible_v<typename arrayInfo<T>::type> || is_complex_convertible_v<typename arrayInfo<T>::type>,#PARAM" must be floating point, integer or complex convertible in "#FUNC);
#define AERRCHECKREALCOMPLEX(T,FUNC,PARAM) static_assert(std::is_floating_point_v<typename arrayInfo<T>::type> || is_complex<typename arrayInfo<T>::type>,#PARAM" must be floating point or complex in "#FUNC);
#define AERRCHECKREALINTEGERSTRING(T, FUNC, PARAM) static_assert(is_integral_convertible_v<typename arrayInfo<T>::type> || is_real_convertible_v<typename arrayInfo<T>::type> || is_string_convertible_v<T>,#PARAM" must be floating point, string or integer convertible in "#FUNC);
#define AERRCHECKTRUTHY(T,FUNC,PARAM) static_assert(std::is_same_v<bool,arrayInfo<T>::type>,#PARAM" must be boolean in" #FUNC);
#define AERRCHECKSTRINGY(T,FUNC,PARAM) static_assert(is_string_convertible_v<std::remove_reference_t<T>>,#PARAM" must be string in "#FUNC);
#else
	//Don't check at all
#define AERRCHECKREAL(T, FUNC, PARAM)
#define AERRCHECKINTEGER(T, FUNC, PARAM)
#define AERRCHECKCOMPLEX(T, FUNC, PARAM)
#define AERRCHECKREALINTEGER(T, FUNC, PARAM)
#define AERRCHECKREALINTEGERCOMPLEX(T,FUNC, PARAM)
#define AERRCHECKREALCOMPLEX(T,FUNC, PARAM)
#define AERRCHECKREALINTEGERSTRING(T,FUNC,PARAM)
#define AERRCHECKTRUTHY(T,FUNC, PARAM)
#define AERRCHECKSTRINGY(T,FUNC,PARAM)
#endif

#if (defined(STRICT_TYPE) || defined(CONVERTIBLE_TYPES)) && !defined(NO_TYPES)
	//Only check on strict or convertible type checking type checking
#define ERRCHECKREAL(T, FUNC, PARAM) AERRCHECKREAL(T, FUNC, PARAM)
#define ERRCHECKINTEGER(T, FUNC, PARAM) AERRCHECKINTEGER(T, FUNC, PARAM)
#define ERRCHECKCOMPLEX(T, FUNC, PARAM) AERRCHECKCOMPLEX(T, FUNC, PARAM)
#define ERRCHECKREALINTEGER(T, FUNC, PARAM) AERRCHECKREALINTEGER(T, FUNC, PARAM)
#define ERRCHECKREALINTEGERCOMPLEX(T,FUNC, PARAM) AERRCHECKREALINTEGERCOMPLEX(T,FUNC,PARAM)
#define ERRCHECKREALCOMPLEX(T,FUNC, PARAM) AERRCHECKREALCOMPLEX(T,FUNC,PARAM)
#define ERRCHECKREALINTEGERSTRING(T,FUNC,PARAM) AERRCHECKREALINTEGERCOMPLEX(T,FUNC,PARAM)
#define ERRCHECKTRUTHY(T,FUNC, PARAM) AERRCHECKTRUTHY(T,FUNC,PARAM)
#define ERRCHECKSTRINGY(T,FUNC,PARAM) AERRCHECKSTRINGY(T,FUNC,PARAM)

#else

	//Don't check at all
#define ERRCHECKREAL(T, FUNC, PARAM)
#define ERRCHECKINTEGER(T, FUNC, PARAM)
#define ERRCHECKCOMPLEX(T, FUNC, PARAM)
#define ERRCHECKREALINTEGER(T, FUNC, PARAM)
#define ERRCHECKREALINTEGERCOMPLEX(T,FUNC, PARAM)
#define ERRCHECKREALCOMPLEX(T,FUNC, PARAM)
#define ERRCHECKREALINTEGERSTRING(T,FUNC,PARAM)
#define ERRCHECKTRUTHY(T,FUNC, PARAM)
#define ERRCHECKTRUTHY(T,FUNC,PARAM)

#endif

	/**
	 * @brief assign a single value to all elements of a tuple
	 * @tparam level The level of the tuple to assign to
	 * @tparam T1 The type of the tuple
	 * @tparam T2 The type of the value
	 * @param tuple The tuple to assign to
	 * @param value The value to assign
	 */
	template<int level=0, typename T1, typename T2>
		inline	void assignTuple(T1& tuple, const T2 value){
			std::get<level>(tuple)=value;
			if constexpr(level < std::tuple_size_v<T1>-1) assignTuple<level+1>(tuple,value);
		}

	/**
	 * @brief assign a single value to a specific element of a tuple
	 * @tparam level The level of the tuple to assign to
	 * @tparam T1 The type of the tuple
	 * @tparam T2 The type of the value
	 * @param tuple The tuple to assign to
	 * @param target The element to assign to
	 * @param value The value to assign
	 */
	template<int level=0, typename T1, typename T2>
		inline	void assignTupleLevel(T1& tuple, const int target, const T2 value){
			if (level == target) {
				std::get<level>(tuple)=value;
			}
			if constexpr (level<std::tuple_size_v<T1>-1) assignTupleLevel<level+1>(tuple,target,value);
		}

	/**
	 * @brief assign a tuple from values in a rank 1 FAR++ array
	 * @tparam level The level of the tuple to assign to
	 * @tparam T_Tuple The type of the tuple
	 * @tparam T_Array The type of the array
	 * @param tuple The tuple to assign to
	 * @param array The array to assign from
	 */
	template<int level=0,typename T_Tuple=void, typename T_Array=void>
		inline void assignTupleArray(T_Tuple & tuple, const T_Array & array){
			std::get<level>(tuple) = array.defaultBase(level+defaultLB);
			if constexpr (level<std::tuple_size_v<T_Tuple>-1) assignTupleArray<level+1>(tuple,array);
		}

	/**
	 * @brief assign the elements of a tuple to a rank 1 FAR++ array
	 * @tparam level The level of the tuple to assign to
	 * @tparam T_Tuple The type of the tuple
	 * @tparam T_Array The type of the array
	 * @param tuple The tuple to assign to
	 * @param array The array to assign from
	 */
	template<int level=0,typename T_Tuple=void, typename T_Array=void>
		inline void assignArrayTuple(T_Array & array, const T_Tuple &tuple){
			array.defaultBase(level+defaultLB) =  std::get<level>(tuple);
			if constexpr (level<std::tuple_size_v<T_Tuple>-1) assignArrayTuple<level+1>(array,tuple);
		}

	/**
	 * Type to remove only rvalue reference types
	 */
	template<typename T>
		struct remove_rvalue_ref {
			typedef T type; //Keep everything other than rvalue references unchanged
		};

	/**
	 * Type to remove only rvalue reference types
	 */
	template<typename T>
		struct remove_rvalue_ref<T&&> {
			typedef const T& type; //Convert rvalue references to const lvalue references
		};

	/**
	 * Equivalent to std::forward_as_tuple, but rvalue are converted to values
	 */
	template<typename... Args>
		inline constexpr auto forward_as_tuple_no_rvalue(Args&&... args) {
			return std::tuple<typename remove_rvalue_ref<Args>::type...>(std::forward<Args>(args)...);
		}


	/**
	 * @brief Function to create a destination tuple copied from a source tuple
	 * with one element left out
	 * @tparam level_s The level of the source tuple
	 * @tparam level_d The level of the destination tuple
	 * @tparam T1 The type of the source tuple
	 * @tparam T2 The type of the destination tuple
	 * @param src The source tuple
	 * @param dest The destination tuple
	 * @param reduction The element to leave out
	 */
	template<int level_s=0, int level_d=0, typename T1, typename T2>
		inline void reduceTuple(const T1& src, T2 & dest, const int reduction){
			if constexpr(level_s == std::tuple_size<T1>() || 
					level_d == std::tuple_size<T2>()){
				return;
			} else {
				if (level_s!=reduction){
					std::get<level_d>(dest) = std::get<level_s>(src);
					reduceTuple<level_s+1,level_d+1>(src,dest,reduction);
				} else{
					reduceTuple<level_s+1,level_d>(src,dest,reduction);
				}
			}
		}

	/**
	 * @brief function to get a specific element of a tuple
	 * @tparam level The level of the tuple to get
	 * @tparam T1 The type of the tuple
	 * @param tuple The tuple to get from
	 * @param target The element to get
	 */
	template<int level=0, typename T1=double>
		inline decltype(auto) getTupleLevel(T1& tuple, const int target){
			if (level == target) {
				return std::get<level>(tuple);
			}
			if constexpr (level<std::tuple_size_v<T1>-1) 
			{
				return getTupleLevel<level+1>(tuple,target);
			} else {
				static_assert(level>=std::tuple_size_v<T1>-1,"Tuple level requested out of range");
				return std::get<0>(tuple); //This line is just to suppress the warning about not returning
			}
		}

	/**
	* @brief Compare two tuples for greater than element wise from the last element
	* @tparam T_tuple The type of the tuple
	* @tparam level The level of the tuple to compare
	* @param lhs The left hand side tuple
	* @param rhs The right hand side tuple
	* @return True if lhs > rhs element wise starting from the last element
	*/
	template<typename T_tuple, int level=0>
		inline bool indexTupleGT(const T_tuple &lhs, const T_tuple &rhs){
			if constexpr(std::tuple_size_v<T_tuple> == 0){
				return false;
			} else {
				bool depth;
				if constexpr(level < std::tuple_size_v<T_tuple>-1){
					depth = indexTupleGT<T_tuple,level+1>(lhs,rhs);
					return depth && std::get<level>(lhs) > std::get<level>(rhs);
				} else {
					return std::get<level>(lhs) > std::get<level>(rhs);
				}
			}
		}

	/**
	* @brief Compare two tuples for less than element wise from the last element
	* @tparam T_tuple The type of the tuple
	* @tparam level The level of the tuple to compare
	* @param lhs The left hand side tuple
	* @param rhs The right hand side tuple
	* @return True if lhs < rhs element wise starting from the last element
	*/
	template<typename T_tuple, int level=0>
		inline bool indexTupleLT(const T_tuple &lhs, const T_tuple &rhs){
			if constexpr(std::tuple_size_v<T_tuple> == 0){
				return false;
			} else {
				bool depth;
				if constexpr(level < std::tuple_size_v<T_tuple>-1){
					depth = indexTupleGT<T_tuple,level+1>(lhs,rhs);
					return depth && std::get<level>(lhs) < std::get<level>(rhs);
				} else {
					return std::get<level>(lhs) < std::get<level>(rhs);
				}
			}
		}

	/**
	 * A version of max where you can get a pointer to a specific instance
	 */
	template<typename T_first, typename... T_others>
		inline auto pmax(const T_first &first, const T_others&... others){
			if constexpr(sizeof...(others)>0){
				decltype(auto) v = pmax(others...);
				return first>v?first:v;
			} else {
				return first;
			}
		}

	/**
	 * A version of min where you can get a pointer to a specific instance
	 */
	template<typename T_first, typename... T_others>
		inline auto pmin(const T_first &first, const T_others&... others){
			if constexpr(sizeof...(others)>0){
				auto v = pmin(others...);
				return first<v?first:v;
			} else {
				return first;
			}
		}

	/**
	 * @brief Get either a specific element of an array or return
	 * the value if data is not an array. Use the default based lower bounds
	 * of the array
	 * @tparam T_info The type of the index
	 * @tparam T_data The type of the data
	 * @param data The data to get from
	 * @param info The index to get
	 * @return The value at the index
	 * @note This function strips off any LBrefWrapper
	 */
	template<typename T_info, typename T_data>
		FAR_INLINE decltype(auto) getItem(T_data &&data, const T_info &info){
			constexpr int arrayRank = arrayInfo<T_data>::rank;
			if constexpr(arrayRank>0){
				if constexpr(std::is_same_v<std::decay_t<T_info>,typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,arrayRank>::type>){
					using TRet = std::decay_t<decltype(std::apply([&data](auto... params) -> decltype(auto){return data.defaultBase(params...);},info))>;
					//Called with a simple tuple
					if constexpr(isLBrefWrapper<TRet>::value){
						return std::apply([&data](auto... params) -> decltype(auto){return data.defaultBaseNoUnpack(params...);},info).get();
					} else {
						return std::apply([&data](auto... params) -> decltype(auto){return data.defaultBaseNoUnpack(params...);},info);
					}
				} else {
					//Called with an info object
					constexpr int infoRank = T_info::rank;
					if constexpr(infoRank==arrayRank) {
						using TRet = std::decay_t<decltype(std::apply([&data](auto... params) -> decltype(auto){return data.defaultBase(params...);},info.indices))>;
						if constexpr(isLBrefWrapper<TRet>::value){
							return std::apply([&data](auto... params) -> decltype(auto){return data.defaultBaseNoUnpack(params...);},info.indices).get();
						} else {
							return std::apply([&data](auto... params) -> decltype(auto){return data.defaultBaseNoUnpack(params...);},info.indices);
						}
					} else {
						using TRet = std::decay_t<decltype(std::apply([&data](auto... params) -> decltype(auto){return data.defaultBase(params...);},info.reduced_indices))>;
						if constexpr(isLBrefWrapper<TRet>::value){
							return std::apply([&data](auto... params) -> decltype(auto){return data.defaultBaseNoUnpack(params...);},info.reduced_indices).get();
						} else {
							return std::apply([&data](auto... params) -> decltype(auto){return data.defaultBaseNoUnpack(params...);},info.reduced_indices);
						}
					}
					//Deliberately no other else here - code will fail to compile if these conditions are not met
				}
			} else {
					return std::forward<T_data>(data);
			}
		}


		/**
		 * @brief Get either a specific element of an array or return
		 * the value if data is not an array. Use the default based lower bounds
		 * of the array
		 * @tparam T_info The type of the index
		 * @tparam T_data The type of the data
		 * @param data The data to get from
		 * @param info The index to get
		 * @return The value at the index
		 * @note This function does not strip off any LBrefWrapper
		 */
		template<typename T_info, typename T_data>
			FAR_INLINE decltype(auto) getItemNoUnpack(T_data &&data, const T_info &info){
				constexpr int arrayRank = arrayInfo<T_data>::rank;
				if constexpr(arrayRank>0){
					if constexpr(std::is_same_v<std::decay_t<T_info>,typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,arrayRank>::type>){
						//Called with a simple tuple
						return std::apply([&data](auto... params) -> decltype(auto){return data.defaultBaseNoUnpack(params...);},info);
					} else {
						//Called with an info object
						constexpr int infoRank= T_info::rank;
						if constexpr(infoRank==arrayRank) {
							return std::apply([&data](auto... params)->decltype(auto){return data.defaultBaseNoUnpack(params...);},info.indices);
						} else {
							return std::apply([&data](auto... params)->decltype(auto){return data.defaultBaseNoUnpack(params...);},info.reduced_indices);
						}
						//Deliberately no other else here - code will fail to compile if these conditions are not met
					}
				} else {
						return std::forward<T_data>(data);
				}
			}


			/**
			 * @brief Mapping function that maps Arrays to elements. Uses getItem to get the actual elements. If current is a value it is returned as is
			 * @tparam T_info The type of the index
			 * @tparam T_current The type of the current element
			 * @tparam T_others The type of the other elements
			 * @param info The index to get
			 * @param current The current element
			 * @param others The other elements
			 * @return The mapped elements as a tuple
			 */
			template<typename T_info, typename T_current, typename... T_others>
				FAR_INLINE decltype(auto) minimalNMap(T_info info, T_current &&current, T_others&&... others){
					if constexpr(sizeof...(others)>0){
						return std::tuple_cat(forward_as_tuple_no_rvalue(getItem(current,info)),minimalNMap(info,others...));
					} else {
						return forward_as_tuple_no_rvalue(getItem(current,info));
					}
				}

			/**
			 * @brief Mapping function that maps Arrays to elements. Uses getItemNoUnpack to get the actual elements. If current is a value it is returned as is
			 * @tparam T_info The type of the index
			 * @tparam T_current The type of the current element
			 * @tparam T_others The type of the other elements
			 * @param info The index to get
			 * @param current The current element
			 * @param others The other elements
			 * @return The mapped elements as a tuple
			 */
			template<typename T_info, typename T_current, typename... T_others>
				FAR_INLINE decltype(auto) minimalNMapNoRefUnwrap(T_info info, T_current &&current, T_others&&... others){
					if constexpr(sizeof...(others)>0){
						return std::tuple_cat(forward_as_tuple_no_rvalue(getItemNoUnpack(current,info)),minimalNMapNoRefUnwrap(info,others...));
					} else {
						return forward_as_tuple_no_rvalue(getItemNoUnpack(current,info));
					}
				}


			/** Class representing information about a shift to be applied to a tuple
			*/
			template<int I_rank,int S_rank>
				struct offset_info{
					static const int rank = I_rank;
					static const int source_rank = S_rank;
					FAR_SIGNED_INDEX_TYPE offset[rank];
					FAR_SIGNED_INDEX_TYPE stride[source_rank];
					bool index_bound[source_rank]={};
					FAR_SIGNED_INDEX_TYPE bound_index[source_rank];
				};

			/**
			 * Shift a tuple of indices by a runtime offset in each direction
			 * Works for both positive and negative strides due to subshift
			 * (Is there a better solution?)
			 */
			template<typename T_offset_info, typename T_dest, typename T_src, int dlevel=0, int slevel=0>
				inline void shiftPartialOffsetTuple(const T_offset_info &oi, T_dest &offset_tuple, T_src &src_tuple){
					if (!oi.index_bound[slevel]){
						char subshift = oi.stride[dlevel]>0?1:0;
						if constexpr(dlevel < T_offset_info::rank) {
							std::get<slevel>(src_tuple) = oi.offset[dlevel] + (std::get<dlevel>(offset_tuple)) * oi.stride[dlevel] -
								(oi.stride[dlevel]-subshift);
						}
						if constexpr(slevel!=T_offset_info::source_rank-1 
							/*&& dlevel !=T_offset_info::rank-1*/) shiftPartialOffsetTuple<T_offset_info, T_dest, T_src, dlevel+1, slevel+1>(oi, offset_tuple, src_tuple);
					} else {
						std::get<slevel>(src_tuple) = oi.bound_index[slevel];
						if constexpr(slevel!=T_offset_info::source_rank-1) shiftPartialOffsetTuple<T_offset_info, T_dest, T_src, dlevel, slevel+1>(oi, offset_tuple, src_tuple);
					}
				}

			/**
			 * Shift a tuple of indices by a runtime offset in each direction
			 * Works for both positive and negative strides due to subshift
			 * (Is there a better solution?)
			 */
			template<typename T_offset_info, typename T_offset, int level=0>
				void shiftOffsetTuple(const T_offset_info &oi, T_offset &offset_tuple){
					char subshift = oi.stride[level]>0?1:0;
					std::get<level>(offset_tuple) = oi.offset[level] + (std::get<level>(offset_tuple)) * oi.stride[level] -
						(oi.stride[level]-subshift);
					if constexpr(level!=T_offset_info::rank-1) shiftOffsetTuple<T_offset_info, T_offset, level+1>(oi, offset_tuple);
				}

			/**
			 * Get the first element of an array. Used in MAX and MIN operations
			 */
			template <typename T_src, int level = 0>
				typename T_src::core_type getFirstValue(const T_src &src, typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, T_src::core_rank>::type tuple_src = {})
				{
					if constexpr (level != T_src::core_rank)
					{
						std::get<level>(tuple_src) = src.getLB(level+1);
						return getFirstValue<T_src, level + 1>(src, tuple_src);
					}
					else
					{
						return std::apply(src, tuple_src);
					}
				}

			/**
			 * @brief Routine to check if two arrays have memory that overlaps potentially destructively
			 * Should be improved to check more than just overlapi
			 * @tparam T1 The type of the first array
			 * @tparam T2 The type of the second array
			 * @tparam level The level of the array to check
			 * @param t1 The first array
			 * @param t2 The second array
			 * @param index_shifted True if the indices are shifted
			 */
			template <typename T1, typename T2, int level = 0>
				bool checkMemoryOverlap(const T1 &t1, const T2 &t2, bool index_shifted, 
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, T1::core_rank>::type tuple_e1_t1 = {},
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, T1::core_rank>::type tuple_e2_t1 = {},
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, T1::core_rank>::type order_t1 = {},
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, T2::core_rank>::type tuple_e1_t2 = {},
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, T2::core_rank>::type tuple_e2_t2 = {},
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, T2::core_rank>::type order_t2 = {})
				{

					if constexpr (level < T1::core_rank || level < T2::core_rank)
					{
						if constexpr(level < T1::core_rank){
							std::get<level>(tuple_e1_t1) = t1.getLB(level+1);
							std::get<level>(tuple_e2_t1) = t1.getUB(level+1);
							std::get<level>(order_t1) = t1.getStride(level+1);
						}
						if constexpr(level < T2::core_rank){
							std::get<level>(tuple_e1_t2) = t2.getLB(level+1);
							std::get<level>(tuple_e2_t2) = t2.getUB(level+1);
							std::get<level>(order_t2) = t2.getStride(level+1);
						}
						return checkMemoryOverlap<T1, T2, level + 1>(t1, t2, index_shifted, tuple_e1_t1, tuple_e2_t1, order_t1, tuple_e1_t2, tuple_e2_t2, order_t2);
					}
					else
					{
						FAR_UNSIGNED_INDEX_TYPE start1 = reinterpret_cast<FAR_UNSIGNED_INDEX_TYPE>(&std::apply(t1, tuple_e1_t1));
						FAR_UNSIGNED_INDEX_TYPE end1 = reinterpret_cast<FAR_UNSIGNED_INDEX_TYPE>(&std::apply(t1, tuple_e2_t1));
						FAR_UNSIGNED_INDEX_TYPE start2 = reinterpret_cast<FAR_UNSIGNED_INDEX_TYPE>(&std::apply(t2, tuple_e1_t2));
						FAR_UNSIGNED_INDEX_TYPE end2 = reinterpret_cast<FAR_UNSIGNED_INDEX_TYPE>(&std::apply(t2, tuple_e2_t2));
						bool swapped = false;
						if (start1 > end1)
						{
							std::swap(start1, end1);
							swapped = !swapped;
						}
						if (start2 > end2)
						{
							std::swap(start2, end2);
							swapped = !swapped;
						}

						return !((start1 > end2 || start2 > end1)				   // If this is not true memory does not overlap at all
								|| (start1 == start2 && end1 == end2 && !swapped && !index_shifted) // Full overlap of source and desination with no reverse ordering
								);
					}
				}

			template <typename T, int rank, bounds_check_state bounds_check, template <int, bounds_check_state> class indexer, bool forceContiguous>
				class Array;

			/**
			 * @brief Class to build an index for a contiguous array
			 * @tparam rank The rank of the array
			 * @tparam bounds_check The bounds check state
			 * @tparam order Fortran order or C order
			 */
			template <int rank, bounds_check_state bounds_check = bounds_check_state::bc_default, int order=0>
				class contiguousIndexBuilder
				{
					template<int rank2, bounds_check_state bounds_check2, int order2>
						friend class contiguousIndexBuilder;
					private:
					public:
					static constexpr int core_rank = rank;
					static constexpr int core_order = order;
					static constexpr bounds_check_state core_bounds = bounds_check;
					// n = size of strides for index calculation
					// n_subsize = size of array reported to user
					FAR_SIGNED_INDEX_TYPE n[rank] = {}, n_subsize[rank] = {};
					// lb = lower bounds
					// ub = upper bounds
					// stride = user specified stride
					FAR_SIGNED_INDEX_TYPE lb[rank] = {}, ub[rank] = {}, stride[rank] = {};
					// size = Total number of allocated elements
					// slice_size = Number of elements in specified slice
					FAR_SIGNED_INDEX_TYPE size = 0, slice_size = 0;
					// offset = Index offset for lower bounds (returns (lb,lb,lb,...) to 0)
					// slice_offset = Index offset for slice start (actually offsets into array)
					FAR_SIGNED_INDEX_TYPE offset = 0, slice_offset = 0;
					//contiguous = is the index to a contiguous chunk of memory
					bool contiguous = true;
					// Function to calculate all of the sizes etc. when an array is sized
					// Each parameter must be an integer or a Range object. If an integer then sizes are given from
					// 1 to N (or 0 to N if DEFAULT_ZERO_LB is set). If a Range object then the lower and upper bounds
					// From the range are used. Stride in the range object is ignored
					/**
					 *Build the size of an array
					 */
					template <int level=0, typename current, typename... rest>
						void build_size(current c, rest... r)
						{
							constexpr int rlevel = order==0?level:rank-level-1;
							if constexpr (level == 0)
							{
								size = 1;
								slice_size = 1;
								offset = 0;
								slice_offset = 0;
							}
							stride[level] = 1;
							if constexpr (std::is_integral_v<current>)
							{
								lb[rlevel] = defaultLB;
								n[rlevel] = static_cast<FAR_UNSIGNED_INDEX_TYPE>(c);
								ub[rlevel] = lb[rlevel] + n[rlevel] - 1;
							}
							else
							{
								lb[rlevel] = c.lb;
								ub[rlevel] = c.ub;
								n[rlevel] = c.ub - c.lb + 1;
							}
							n_subsize[rlevel] = n[rlevel];
							size *= n[rlevel];
							slice_size *= n_subsize[rlevel];
							if constexpr(sizeof...(rest)==0)
							{
								offset -= lb[rlevel];
							}
							else
							{
								build_size<level + 1>(r...);
								offset = offset * n[rlevel] - lb[rlevel];
							}
						}

					/**
					 * Build just the offsets of an array. Used when an array's bounds are shifted
					 * Sets n, n_subsize, size and slice_size from ub and lb
					 */
					template <int level=0>
						void build_offsets()
						{
							constexpr int rlevel = order==0?level:rank-level-1;
							if constexpr (level == 0)
							{
								size = 1;
								slice_size = 1;
								offset = 0;
								slice_offset = 0;
								contiguous = true;
							}
							n[rlevel] = ub[rlevel] - lb[rlevel] + 1;
							n_subsize[rlevel] = n[rlevel];
							size *= n[rlevel];
							slice_size *= n_subsize[rlevel];
							if constexpr (rlevel == (order==0?rank-1:0))
							{
								offset -= lb[rlevel];
							}
							else
							{
								build_offsets<level + 1>();
								offset = offset * n[rlevel] - lb[rlevel];
							}
						}

					template<int level, typename T_other, typename current, typename... others>
					FAR_INLINE void build_slice_offset(T_other & other, current in, others... r)
					{
						Range c;
						constexpr FAR_UNSIGNED_INDEX_TYPE rlevel = (order==0?level:rank-level-1);
						constexpr bool is_range = !std::is_integral_v<current>;
						if constexpr (!is_range)
						{
							c = Range(in, in, 1);
						}
						else
						{
							c = in;
							c.lb = c.lb == std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::min() ? lb[rlevel] : c.lb;
							c.ub = c.ub == std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::max() ? ub[rlevel] : c.ub;
#ifndef NO_REVERSE_NEGATIVE_STRIDE
							if (c.stride<0) std::swap(c.lb,c.ub);
#endif
						}
						if constexpr (sizeof...(others)>0) build_slice_offset<level+1>(other, r...);
						slice_offset = slice_offset * other.n[rlevel] + ((c.stride>0?c.lb:c.ub) - other.lb[rlevel]) * other.stride[rlevel];	
					}

					/**
					 * Build the information needed for a slice of an array
					 */
					template <int level=0, int slevel=0, bool firstcall=true, typename T_other, typename current, typename... rest>
						FAR_INLINE void build_slice(T_other &other, current in, rest... r)
						{
							if constexpr( level == rank) {
								build_slice_offset<slevel>(other,in,r...);
							} else {

								Range c;
								constexpr FAR_UNSIGNED_INDEX_TYPE rlevel = (order==0?level:rank-level-1);
								constexpr FAR_UNSIGNED_INDEX_TYPE rslevel = (T_other::core_order==0?slevel:T_other::core_rank-slevel-1);
								constexpr bool is_range = !std::is_integral_v<current>;
								constexpr bool should_advance = is_range;
								if constexpr (!is_range)
								{
									c = Range(in, in, 1);
								}
								else
								{
									c = in;
									c.lb = c.lb == std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::min() ? other.lb[rslevel] : c.lb;
									c.ub = c.ub == std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::max() ? other.ub[rslevel] : c.ub;
	#ifndef NO_REVERSE_NEGATIVE_STRIDE
									if (c.stride<0) std::swap(c.lb,c.ub);
	#endif
								}

								if constexpr (firstcall)
								{
									size=1;
									slice_size = 1;
									slice_offset = 0;
									offset = 0;
									for (int i=0;i<rank;++i) {
										n[i]=1;
										n_subsize[i]=1;
										stride[i]=1;
									}
									contiguous = true;
								}
								if constexpr(is_range){
									n_subsize[rlevel] = std::ceil(double(c.ub - c.lb + 1) / std::abs((double)c.stride));
									lb[rlevel] = defaultLB;
									ub[rlevel] = lb[rlevel] + n_subsize[rlevel] - 1;
									//Number of elements in slice is always just based on *this* slice other information 	
									//has already been folded in
								}
								n[rlevel] *= other.n[rslevel]; // Number of elements in the whole dimension
								size *= n_subsize[rlevel];
								slice_size *= n_subsize[rlevel];
								stride[rlevel] *= c.stride * other.stride[rslevel] * (is_range?1:(other.n[rslevel]*other.stride[rslevel]*c.stride));

								if constexpr(sizeof...(rest)>0){
									build_slice<level + (should_advance?1:0),slevel+1,false>(other, r...);
								}
								if constexpr(is_range){
									offset = offset * n[rlevel] - stride[rlevel] * lb[rlevel];
								}
								FAR_SIGNED_INDEX_TYPE upper = (c.stride>0?c.lb:c.ub);
								//upper = c.lb;
								//The purpose of slice offset is to map 1 in this array to the start of the
								//slice in the other array
								slice_offset = slice_offset * other.n[rslevel] + (upper - other.lb[rslevel]) * other.stride[rslevel];

								//May still be contiguous if
								//1) Other end being sliced is contiguous
								//2) Stride is uniform (AND)
								//(
								//3) Slice is in the last dimension(that is not sliced out) (OR)
								//4) Slice is the whole length of a dimension (OR)
								//5) Slice is a single element
								// )
								contiguous = contiguous && other.contiguous && (((stride[rlevel] == 1 && 
										((c.lb == other.lb[rslevel] && c.ub == other.ub[rslevel]) || c.lb==c.ub))
										) || (rlevel == rank-1 && stride[rlevel] == 1));
							}
						}


					/**
					 * Build the offset for a slice in a lazy array
					 */
					template<int dlevel=0, int slevel=0, typename T_off=void, typename T_other=void, typename T_offset_info=void, typename... T_ranges>
						void buildLazySliceOffsets(T_offset_info &oi,  T_other && otherIndex, T_off && current, T_ranges&&... others){

							if constexpr(slevel==0) {
								slice_size=1;
								size=1;
								slice_offset=0;
							}
							Range c;
							if constexpr(std::is_same_v<std::decay_t<T_off>,Range>){
								lb[dlevel]=otherIndex.lb[slevel];
								ub[dlevel]=otherIndex.ub[slevel];
								c = current;
								if (c.lb == std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::min()) c.lb = lb[dlevel];
								if (c.ub == std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::max()) c.ub = ub[dlevel];
								stride[dlevel]=c.stride;
								n_subsize[dlevel] = std::ceil(double(c.ub - c.lb + 1) / std::abs((double)stride[dlevel]));
								n[dlevel] = n_subsize[dlevel];
								size *= n[dlevel];
								slice_size *= n_subsize[dlevel];
								oi.offset[dlevel]=c.stride>0?c.lb-lb[dlevel]:c.ub;
								oi.stride[dlevel]=c.stride;

								lb[dlevel]=defaultLB;
								ub[dlevel]=lb[dlevel] + n_subsize[dlevel]-1;
								if constexpr(sizeof...(others)>0) buildLazySliceOffsets<dlevel+1, slevel+1>(oi, otherIndex, others...);
							} else {
								oi.index_bound[slevel] = true;
								oi.bound_index[slevel] = current-otherIndex.lb[slevel]+defaultLB;
								if constexpr(sizeof...(others)>0) buildLazySliceOffsets<dlevel, slevel+1>(oi, otherIndex, others...);
							}
						}

					/**
					 * Build the actual offset into the rank 1 array
					 */
					template <int level=0, typename current, typename... rest>
						FAR_PURE FAR_SIGNED_INDEX_TYPE buildIndex(current c, rest... r) const
						{
							constexpr FAR_SIGNED_INDEX_TYPE rlevel = (order==0?level:rank-level-1);
							FAR_SIGNED_INDEX_TYPE item = static_cast<FAR_SIGNED_INDEX_TYPE>(c) * stride[rlevel];
							if constexpr (bounds_check)
							{
								if (c < lb[rlevel] || c > ub[rlevel])
								{
									std::stringstream s;
									s << "In index " << rlevel + 1 << ", " << c << " is outside the range of " << lb[rlevel] << " to " << ub[rlevel];
									throw std::out_of_range(s.str());
								}
							}
							if constexpr (level < rank - 1)
							{
								item += buildIndex<level + 1>(r...) * n[rlevel];
							}
							if constexpr(level==0) item += offset;
							return item;
						}

					/**
					 * Build offset into the rank 1 array if the array starts at default index
					 */
					template <int level=0, typename current, typename... rest>
						void build_default_index(FAR_SIGNED_INDEX_TYPE &offset, current c, rest... r) const
						{
							constexpr FAR_SIGNED_INDEX_TYPE rlevel = (order==0?level:rank-level-1);
							if constexpr (bounds_check)
							{
								if (c < defaultLB || c >= n[rlevel]+defaultLB)
								{
									std::stringstream s;
									s << "In index " << rlevel + 1 << ", " << c << " is outside the range of " << defaultLB << " to " << defaultLB+n_subsize[rlevel];
									throw std::out_of_range(s.str());
								}
							}
							if constexpr (sizeof...(r)>0)
							{
								build_default_index<level + 1>(offset, r...);
								offset *= n[rlevel];
							}
							offset += static_cast<FAR_SIGNED_INDEX_TYPE>((c - defaultLB) * stride[rlevel]);
						}

					//Function to build a tuple of indices from a single index
					template <int current=0>
						auto reverseIndex(FAR_UNSIGNED_INDEX_TYPE index) const
						{
							constexpr FAR_SIGNED_INDEX_TYPE rlevel = (order==0?current:rank-current-1);
							FAR_UNSIGNED_INDEX_TYPE rsize = n_subsize[rlevel];
							if constexpr (current==rank-1){
								return std::make_tuple(static_cast<FAR_SIGNED_INDEX_TYPE>(index%rsize + lb[rlevel]));
							}
							else{
								if constexpr(order==0){
									return std::tuple_cat(std::make_tuple(static_cast<FAR_SIGNED_INDEX_TYPE>(index%rsize + lb[rlevel])),reverseIndex<current+1>(index/rsize));
								} else {
									return std::tuple_cat(reverseIndex<current+1>(index/rsize),std::make_tuple(static_cast<FAR_SIGNED_INDEX_TYPE>(index%rsize + lb[rlevel])));
								}
							}
						}

					/**
					 * Calculate the shift in an index due to a given circular shift to the index
					 * This allows shifts in all directions simultaneously*/
					template <int level = 0>
						void cshiftIndex(typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, rank>::type &index, 
								typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, rank>::type shift) const
						{
							std::get<level>(index) = ((std::get<level>(index)-lb[level]) + 
									n_subsize[level] + std::get<level>(shift)) % n_subsize[level] + lb[level];
							if constexpr (level < rank-1)
								cshiftIndex<level + 1>(index,shift);
						}

					/**
					 * Calculate the shift in an index due to a given end off shift to the index
					 * This allows shifts in all directions simultaneously*/
					template <int level = 0>
						bool eoshiftIndex(typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, rank>::type &index,
								typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, rank>::type shift) const
						{
							std::get<level>(index) += std::get<level>(shift);
							bool oob = std::get<level>(index) < lb[level] || std::get<level>(index) > ub[level];
							if constexpr (level < rank-1){
								return oob || eoshiftIndex<level + 1>(index,shift);
							} else {
								return oob;
							}
						}

					public:
					/**
					 * Construct index from list of either sizes or ranges.
					 * If ranges are used then custom upper and lower bounds are used
					 * other wise all array start from 1 in each rank. Strides are ignored
					 * if specified
					 */
					template <typename... Args,std::enable_if_t<is_correct_subscript<rank,Args...>, int> = 0>
						contiguousIndexBuilder(Args... args) { 
							std::apply([this](auto... args) { this->build_size(args...); }, packToTuple<order>(args...));
						}

					/**
					 * Constructor from array of sizes. Used internally
					 */
					contiguousIndexBuilder(FAR_UNSIGNED_INDEX_TYPE sizes[rank])
					{
						size = 1;
						for (int i = 0; i < rank; ++i)
						{
							stride[i] = 1;
							lb[i] = defaultLB;
							ub[i] = lb[i] + sizes[i] - 1;
						}
						build_offsets<0>();
					}

					/**
					 * Construct from std::array of sizes. Used internally
					 */

					template<typename T>
						contiguousIndexBuilder(std::array<T, rank> sizes)
						{
							size = 1;
							for (int i = 0; i < rank; ++i)
							{
								stride[i] = 1;
								lb[i] = defaultLB;
								ub[i] = lb[i] + sizes[i] - 1;
							}
							slice_size = size;
							build_offsets<0>();
						}

					/**
					 * Constructor from FAR++ array of sizes. Used internally
					 */
					template<typename T, std::enable_if_t<arrayInfo<T>::value,int>* =nullptr>
						contiguousIndexBuilder(T sizes)
						{
							size = 1;
							for (int i = defaultLB; i < rank+defaultLB; ++i)
							{
								stride[i] = 1;
								lb[i] = defaultLB;
								ub[i] = lb[i] + sizes[i] - 1;
							}
							slice_size=size;
							build_offsets<0>();
						}

					/**
					 * Default constructor
					 */
					contiguousIndexBuilder() = default;

					/**
					 * Copy constructor
					 */
					template <bounds_check_state bounds_other, int order_other>
						contiguousIndexBuilder(const contiguousIndexBuilder<rank, bounds_other, order_other> &other)
						{
							for (int i = 0; i < rank; ++i)
							{
								n[i] = other.n[i];
								n_subsize[i] = other.n_subsize[i];
								stride[i] = other.stride[i];
								lb[i] = other.lb[i];
								ub[i] = other.ub[i];
							}
							size = other.size;
							slice_size = other.slice_size;
							offset = other.offset;
							slice_offset = other.slice_offset;
							contiguous = other.contiguous;
						}

					/**
					 * Copy constructor
					 */
					contiguousIndexBuilder(const contiguousIndexBuilder &other)
					{
						for (int i = 0; i < rank; ++i)
						{
							n[i] = other.n[i];
							n_subsize[i] = other.n_subsize[i];
							stride[i] = other.stride[i];
							lb[i] = other.lb[i];
							ub[i] = other.ub[i];
						}
						size = other.size;
						slice_size = other.slice_size;
						offset = other.offset;
						slice_offset = other.slice_offset;
						contiguous = other.contiguous;
					}

					/**
					 * Copy explicitly
					 */
					template <bounds_check_state bounds_other, int order_other>
						void copyFrom(const contiguousIndexBuilder<rank, bounds_other, order_other> &other)
						{
							size = 1;
							for (int i = 0; i < rank; ++i)
							{
								n[i] = other.n[i];
								n_subsize[i] = other.n_subsize[i];
								stride[i] = other.stride[i];
								lb[i] = other.lb[i];
								ub[i] = other.ub[i];
							}
							size = other.size;
							slice_size = other.slice_size;
							offset = other.offset;
							slice_offset = other.slice_offset;
							contiguous = other.contiguous;
						}

					/**
					 * Copy explicitly and normalise
					 */
					template <bounds_check_state bounds_other, int order_other>
						void normaliseCopy(const contiguousIndexBuilder<rank, bounds_other, order_other> &other)
						{
							for (int i = 0; i < rank; ++i)
							{
								n[i] = other.n_subsize[i];
								n_subsize[i] = n[i];
								stride[i] = 1;
								lb[i] = other.lb[i];
								ub[i] = other.ub[i];
							}
							build_offsets();
						}

					/**
					 * Move constructor
					 */
					template <bounds_check_state bounds_other, int order_other>
						contiguousIndexBuilder(contiguousIndexBuilder<rank, bounds_other, order_other> &&other) noexcept
						{
							for (int i = 0; i < rank; ++i)
							{
								n[i] = other.n[i];
								n_subsize[i] = other.n_subsize[i];
								stride[i] = other.stride[i];
								lb[i] = other.lb[i];
								ub[i] = other.ub[i];
							}
							size = other.size;
							slice_size = other.slice_size;
							offset = other.offset;
							slice_offset = other.slice_offset;
							contiguous = other.contiguous;
						}

					/**
					 * Move constructor
					 */
					contiguousIndexBuilder(contiguousIndexBuilder&&other) noexcept
						{
							for (int i = 0; i < rank; ++i)
							{
								n[i] = other.n[i];
								n_subsize[i] = other.n_subsize[i];
								stride[i] = other.stride[i];
								lb[i] = other.lb[i];
								ub[i] = other.ub[i];
							}
							size = other.size;
							slice_size = other.slice_size;
							offset = other.offset;
							slice_offset = other.slice_offset;
							contiguous = other.contiguous;
						}	

					/**
					 * Copy assignment operator
					 */
					contiguousIndexBuilder &operator=(const contiguousIndexBuilder &other)
					{
						// Same logic applies for assignment as for construction
						for (int i = 0; i < rank; ++i)
						{
							n[i] = other.n[i];
							n_subsize[i] = other.n_subsize[i];
							stride[i] = other.stride[i];
							lb[i] = other.lb[i];
							ub[i] = other.ub[i];
						}
						slice_size = other.size;
						size = other.size;
						offset = other.offset;
						slice_offset = other.slice_offset;
						return *this;
					}

					/**
					 * Copy assignment operator
					 * from an array with another bounds check
					 */
					template <bounds_check_state bounds_other, int order_other>
						contiguousIndexBuilder &operator=(const contiguousIndexBuilder<rank, bounds_other, order_other> &other)
						{
							// Same logic applies for assignment as for construction
							for (int i = 0; i < rank; ++i)
							{
								n[i] = other.n[i];
								n_subsize[i] = other.n_subsize[i];
								stride[i] = other.stride[i];
								lb[i] = other.lb[i];
								ub[i] = other.ub[i];
							}
							slice_size = other.size;
							size = other.size;
							offset = other.offset;
							slice_offset = other.slice_offset;
							contiguous = other.contiguous;
							return *this;
						}

					/**
					 * Function to copy an indexer assuming that the resulting array
					 * will point to the same memory as the original
					 */
					template <typename T>
						void pointCopy(const T &other)
						{
							// Same logic applies for assignment as for construction
							for (int i = 0; i < rank; ++i)
							{
								n[i] = other.n[i];
								n_subsize[i] = other.n_subsize[i];
								stride[i] = other.stride[i];
								lb[i] = other.lb[i];
								ub[i] = other.ub[i];
							}
							size = other.size;
							slice_size = other.slice_size;
							offset = other.offset;
							slice_offset = other.slice_offset;
							contiguous = other.contiguous;
						}

					/**
					 * Move assignment operator
					 */
					template <bounds_check_state bounds_other, int order_other>
						contiguousIndexBuilder &operator=(contiguousIndexBuilder<rank, bounds_other, order_other> &&other) noexcept
						{
							for (int i = 0; i < rank; ++i)
							{
								n[i] = other.n[i];
								n_subsize[i] = other.n_subsize[i];
								stride[i] = other.stride[i];
								lb[i] = other.lb[i];
								ub[i] = other.ub[i];
							}
							size = other.size;
							slice_size = other.size;
							offset = other.offset;
							slice_offset = other.slice_offset;
							return *this;
						}

					/**
					 * Function to set an indexer to have the same number of elements but be lb_in based
					 */
					void rebase(FAR_SIGNED_INDEX_TYPE lb_in)
					{
						for (int i = 0; i < rank; ++i)
						{
							n[i] = n_subsize[i];
							lb[i] = lb_in;
							ub[i] = lb[i] + n[i] - 1;
						}
						build_offsets<0>();
					}

					/**
					 * Function to take an indexer that has been created possibly from a slice
					 * And convert the indexer to be a non-slice indexer
					 */
					void normalise()
					{
						for (int i = 0; i < rank; ++i)
						{
							stride[i] = 1;
							n[i] = n_subsize[i];
						}
						build_offsets<0>();
					}

					/**
					 * Set the size of an indexer to be a given size. Can take integers
					 * or Range objects. See constructor
					 */
					template <typename... Args, std::enable_if_t<is_correct_subscript<rank,Args...>, int> = 0>
						void setSize(Args... args) { build_size<0>(args...); }

					/**
					 * Set the size of an indexer to be a given size
					 * using an array of sizes*/
					void setSize(FAR_UNSIGNED_INDEX_TYPE sizes[rank])
					{
						size = 1;
						for (int i = 0; i < rank; ++i)
						{
							stride[i] = 1;
							n[i] = sizes[i];
							n_subsize[i] = n[i];
							lb[i] = defaultLB;
							ub[i] = lb[i] + sizes[i] - 1;
							size *= n[i];
						}
						slice_size = size;
						build_offsets<0>();
					}

					/**
					 * Set the size of an indexer to be a given size
					 * using an std::array of sizes*/
					template<typename T, FAR_UNSIGNED_INDEX_TYPE s>
						void setSize(const std::array<T,s> &sizes)
						{
							size = 1;
							offset = 0;
							for (int i = 0; i < rank; ++i)
							{
								stride[i]=1;
								n[i] = sizes[i];
								n_subsize[i] = n[i];
								lb[i] = defaultLB;
								ub[i] = lb[i] + sizes[i] - 1;
								size *= n[i];
							}
							slice_size = size;
							build_offsets<0>();
						}

					/**
					 * Function to set the upper and lower bounds of an array
					 */
					void setBounds(FAR_SIGNED_INDEX_TYPE LB[rank], FAR_SIGNED_INDEX_TYPE UB[rank])
					{
						size = 1;
						for (int i = 0; i < rank; ++i)
						{
							lb[i] = LB[i];
							ub[i] = UB[i];
						}
						build_offsets<0>();
					}

					/**
					 * Get the total number of elements
					 * This is the number of elements available to the user
					 * so is smaller than the number of elements in memory
					 * if the array is sliced
					 */
					FAR_UNSIGNED_INDEX_TYPE getSize() const { return slice_size; }
					/** 
					 * Return the slice size because a user wants to know
					 * what the available indices are, not the real number
					 * of elements used for a stride
					 */
					FAR_UNSIGNED_INDEX_TYPE getRankSize(int irank) const { return n_subsize[irank-1]; }
					/**
					 * Get rank
					 */
					FAR_UNSIGNED_INDEX_TYPE getRank() const {return rank;}
					/**
					 * Get the lower bound for a specified direction
					 */
					FAR_SIGNED_INDEX_TYPE getLB(int irank) const { return lb[irank-1]; }
					/**
					 * Get the upper bound for a specified direction
					 */
					FAR_SIGNED_INDEX_TYPE getUB(int irank) const { return ub[irank-1]; }
					/**
					 * Get the stride for a specified direction
					 */
					FAR_SIGNED_INDEX_TYPE getStride(int irank) const { return stride[irank-1]; }
					/**
					 * Get the contiguous status for the indexer
					 */
					bool getContiguous() const {return contiguous;}

					/**
					 * Round bracket operator. Gets the actual offset
					 */
					template <typename... Args,std::enable_if_t<is_correct_subscript<rank,Args...>, int> = 0>
						FAR_PURE FAR_SIGNED_INDEX_TYPE operator()(Args... args) const
						{
							FAR_SIGNED_INDEX_TYPE value = std::apply([this](auto... args){return this->buildIndex(args...);},packToTuple<order>(args...));
							return value;
						}
					/**
					 * Get the actual offset assuming that the directions are all default based
					 */
					template <typename... Args>
						FAR_PURE FAR_SIGNED_INDEX_TYPE defaultBase(Args... args) const
						{
							static_assert(sizeof...(args) == rank, "Number of arguments must match rank");
							FAR_SIGNED_INDEX_TYPE offset = 0;
							std::apply([&offset,this](auto... args){(this->build_default_index(offset,args...));},packToTuple<order>(args...));
							return offset;
						}

					/**
					 * Square bracket operator. Works as expected, but is very slow because it involves a forwards and then
					 * reverse index calculation
					 */
					__attribute__((pure)) auto operator[](FAR_SIGNED_INDEX_TYPE index) const {
						//auto a = reverseIndex(index);
						return std::apply(*this,reverseIndex(index)); 
					}


					/**
					 * Iterate over an parameter pack of items applying appropriate functions to them in parallel
					 * Function only exists if code compiled in OpenMP mode as a sanity check
					 */
#ifdef _OPENMP
					template <auto callable, auto interposer = defaultFunction, bool reduction=false, bool backwards=false, int level=0, FAR_SIGNED_INDEX_TYPE far_et_length = FAR_ET_LENGTH, FAR_SIGNED_INDEX_TYPE far_non_et_length = FAR_NON_ET_LENGTH, FAR_SIGNED_INDEX_TYPE simple_length = FAR_SIMPLELENGTH, bool parallel = true, typename T_tuple = int, typename... T_params>
						FAR_INLINE iter_action forEachParallel(T_tuple && tupleBase, T_tuple && tupleLB, const int reduction_direction=rdNone, T_params&&... others) const
						{
							constexpr FAR_SIGNED_INDEX_TYPE rlevel = std::max((order==0?rank-level-1:level),0);
							FAR_SIGNED_INDEX_TYPE indexesPerThread = n_subsize[rlevel]/omp_get_num_threads();
							FAR_SIGNED_INDEX_TYPE tlb = defaultLB + omp_get_thread_num() * indexesPerThread;
							FAR_SIGNED_INDEX_TYPE tub = omp_get_thread_num()!=omp_get_num_threads()-1?defaultLB + (omp_get_thread_num()+1) * indexesPerThread:defaultLB + n_subsize[rlevel];
							FAR_OMP_BARRIER
								if constexpr(backwards) {
									for (FAR_SIGNED_INDEX_TYPE i=tub-1;i>=tlb;--i){
										std::get<rlevel>(tupleBase)=i;
										std::get<rlevel>(tupleLB)=i+lb[rlevel]-defaultLB;
										forEachCore<callable, interposer, reduction, backwards, level+1, far_et_length, far_non_et_length, simple_length, parallel>
											(tupleBase, tupleLB, reduction_direction, std::forward<T_params>(others)...);
									}
								} else {
									for (FAR_SIGNED_INDEX_TYPE i=tlb;i<tub;++i){
										std::get<rlevel>(tupleBase)=i;
										std::get<rlevel>(tupleLB)=i+lb[rlevel]-defaultLB;
										forEachCore<callable, interposer, reduction, backwards, level+1, far_et_length, far_non_et_length, simple_length, parallel>
											(tupleBase, tupleLB, reduction_direction, std::forward<T_params>(others)...);
									}
								}
							FAR_OMP_BARRIER
								return iaContinue;
						}
#endif

					/** 
					 * Iterate over an parameter pack of items applying appropriate functions to them
					 */
					template <auto callable, auto interposer = defaultFunction, bool reduction=false, bool backwards=false, int level=0, FAR_SIGNED_INDEX_TYPE et_length = FAR_ET_LENGTH, FAR_SIGNED_INDEX_TYPE non_et_length = FAR_NON_ET_LENGTH, FAR_SIGNED_INDEX_TYPE simple_length = FAR_SIMPLELENGTH, bool parallel = true, typename T_tuple = int, typename... T_params>
						FAR_INLINE iter_action forEachCore(T_tuple && tupleBase, T_tuple && tupleLB, const int reduction_direction=rdNone, T_params&&... others) const
						{
							constexpr FAR_SIGNED_INDEX_TYPE rlevel = std::max((order==0?rank-level-1:level),0);
							constexpr bool dcat = true;
#ifdef _OPENMP
							if constexpr((level == (order==0?0:rank-1)) && parallel){
								if (workshare[getTeamNum()]){
									forEachParallel<callable, interposer, reduction, backwards, level, et_length, non_et_length, simple_length, parallel>
										(tupleBase, tupleLB, reduction_direction, std::forward<T_params>(others)...);
									return iaContinue;
								}
							}
#endif
							constexpr bool can_early_terminate=reduction && std::is_same_v<typename callableTraits<decltype(callable)>::type,iter_action> && dcat;
							if constexpr(level<rank){
								//This is ugly. These values are obtained by experiment
								//If you don't do this then early termination really kills performance
								//Unrolling the loops isn't bad, the fact that gcc and clang require different
								//unrolls is the problem
								constexpr FAR_SIGNED_INDEX_TYPE veclength=((rlevel==(rank-1))?(can_early_terminate?et_length:non_et_length):simple_length);
								FAR_SIGNED_INDEX_TYPE llb = defaultLB;
								FAR_SIGNED_INDEX_TYPE lub = defaultLB+n_subsize[rlevel]-1;
								constexpr FAR_SIGNED_INDEX_TYPE delta = backwards?-1ll:1ll;								
								bool terminate[veclength]={};
								FAR_SIGNED_INDEX_TYPE p1lb = !backwards?llb:lub;
								FAR_SIGNED_INDEX_TYPE p1ub = !backwards?lub-(veclength-1):llb+(veclength-1);
								FAR_SIGNED_INDEX_TYPE i=p1lb;
								//This is the code for arbitrary vector size

								for (i=p1lb;[&p1ub](FAR_SIGNED_INDEX_TYPE &i){return backwards?i>=p1ub:i<=p1ub;}(i);i+=delta*veclength){
									for (FAR_SIGNED_INDEX_TYPE j=0;j<veclength;++j){
										std::get<rlevel>(tupleBase)=i+j*delta;
										std::get<rlevel>(tupleLB)=i+lb[rlevel]-defaultLB+j*delta;
										terminate[j]=
											(forEachCore<callable, interposer, reduction, backwards, level+1, et_length, non_et_length, simple_length, parallel>
												 (tupleBase, tupleLB, reduction_direction, std::forward<T_params>(others)...) == iaTerminate);
									}
									if constexpr(can_early_terminate){
										bool term=false;
										for (int i=0;i<veclength;++i) {
											term |= terminate[i];
											terminate[i]=false;
										}
										if (term) {
											if (reduction_direction == rdAll || rlevel == reduction_direction) {
												return iaTerminate;
											}
										}
									}
								}
								if constexpr (veclength >1){
									FAR_SIGNED_INDEX_TYPE p2lb = i;
									FAR_SIGNED_INDEX_TYPE p2ub = !backwards?lub:llb;
									FAR_SIGNED_INDEX_TYPE t_offset = !backwards?i:p2ub;
									auto l = [&p2ub](FAR_SIGNED_INDEX_TYPE j){return backwards?j>=p2ub:j<=p2ub;};
									for (FAR_SIGNED_INDEX_TYPE j=p2lb;l(j);j+=delta){
										std::get<rlevel>(tupleBase)=j;
										std::get<rlevel>(tupleLB)=j+lb[rlevel]-defaultLB;
										terminate[j-t_offset]=
											(forEachCore<callable, interposer, reduction, backwards, level+1, et_length, non_et_length, simple_length, parallel>
											 (tupleBase, tupleLB, reduction_direction, std::forward<T_params>(others)...)
											 == iaTerminate);
									}
									if constexpr(can_early_terminate){
										bool term=false;
										for (int i=0;i<veclength;++i) term |= terminate[i];
										if (term) {
											if (reduction_direction == rdAll || rlevel == reduction_direction) {
												return iaTerminate;
											}
										}
									}
								}
								return iaContinue;
							} else {
								indexInfo<rank> rri(tupleBase, tupleLB, reduction_direction);
								if constexpr(reduction){if (reduction_direction!=rdNone) reduceTuple(rri.indices, rri.reduced_indices,reduction_direction);}
								if constexpr(callableTraits<decltype(interposer)>::value){
									if constexpr(can_early_terminate) {
										return std::apply(callable, interposer(rri, std::forward<T_params>(others)...));
									} else {
										decltype(auto) a =interposer(rri, std::forward<T_params>(others)...);
										std::apply(callable,a);
										return iaContinue;
									}
								} else {
									if constexpr(can_early_terminate) {
										return callable(rri, std::forward<T_params>(others)...);
									} else {
										callable(rri, std::forward<T_params>(others)...);
										return iaContinue;
									}
								}
							}
						}

					template <auto callable, typename T_matrix_a, typename T_matrix_b>
						static auto matmul(const T_matrix_a & A, const T_matrix_b & B){

							using T1 = typename arrayInfo<T_matrix_a>::type;
							constexpr int rank1 =  arrayInfo<T_matrix_a>::rank; 
							FAR_UNUSED constexpr bounds_check_state b1 = arrayInfo<T_matrix_a>::bounds_check;

							using T2 = typename arrayInfo<T_matrix_b>::type;
							constexpr int rank2 =  arrayInfo<T_matrix_b>::rank;
							FAR_UNUSED constexpr bounds_check_state b2 = arrayInfo<T_matrix_b>::bounds_check;

							if constexpr(rank1==2 && rank2==1) //Matrix vector multiply
							{
								using dest_type = std::invoke_result_t<decltype(std::decay_t<decltype(B)>:: template buildVariant<T2, rank2, b2 ,false>)>;
								dest_type C;
								C.mold(B);
								for (FAR_UNSIGNED_INDEX_TYPE i=defaultLB;i<defaultLB+C.getRankSize(1);++i) C.defaultBase(i)={};
								for (FAR_UNSIGNED_INDEX_TYPE j=defaultLB;j<defaultLB+A.getRankSize(2);++j){
									for (FAR_UNSIGNED_INDEX_TYPE i=defaultLB;i<defaultLB+A.getRankSize(1);++i){
										callable(A.defaultBase(i,j), B.defaultBase(j), C.defaultBase(i));
									}
								}
								return C;
							}
							else if constexpr(rank1==1 && rank2==2) //Vector matrix multiply
							{
								using dest_type = std::invoke_result_t<decltype(std::decay_t<decltype(A)>:: template buildVariant<T1, rank1, b1 ,false>)>;
								dest_type C;
								C.mold(A);
								for (FAR_UNSIGNED_INDEX_TYPE j=defaultLB;j<defaultLB+B.getRankSize(2);++j){
									C.defaultBase(j)={};
									for (FAR_UNSIGNED_INDEX_TYPE i=defaultLB;i<defaultLB+B.getRankSize(1);++i){
										callable(A.defaultBase(i),B.defaultBase(i,j),C.defaultBase(j));
									}
								}
								return C;
							}
							else if constexpr(rank1==2 && rank2==2) //Matrix matrix multiply
							{
								using dest_type = std::invoke_result_t<decltype(std::decay_t<decltype(B)>:: template buildVariant<T1, rank1, b1 ,false>)>;
								dest_type C(A.getRankSize(1),B.getRankSize(2));
								for (FAR_UNSIGNED_INDEX_TYPE j=defaultLB;j<defaultLB+B.getRankSize(2);++j){
									for (FAR_UNSIGNED_INDEX_TYPE i=defaultLB;i<defaultLB+A.getRankSize(1);++i){
										C.defaultBase(i,j)={};
										for (FAR_UNSIGNED_INDEX_TYPE inner =defaultLB;inner<defaultLB+A.getRankSize(2);++inner){
											callable(A.defaultBase(i,inner),B.defaultBase(inner,j),C.defaultBase(i,j));

										}
									}
								}
								return C;
							}
						}

					template<int direction=1, int level=0, bool default_base=true>
						FAR_INLINE bool advanceTuple(typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type &index, FAR_SIGNED_INDEX_TYPE shift=1) const{
							const int rlevel = order==0?level:rank-level-1;
							FAR_SIGNED_INDEX_TYPE UB,LB,sz;
							if constexpr(default_base){
								LB=defaultLB;
								UB=defaultLB+getRankSize(rlevel+1)-1;
							} else {
								LB=getLB(rlevel+1);
								UB=getUB(rlevel+1);
							}
							sz = getRankSize(rlevel+1);
							FAR_SIGNED_INDEX_TYPE &val = std::get<level>(index);
							val+=shift*direction;
							//Deal with val > UB
							FAR_SIGNED_INDEX_TYPE diff;
							if constexpr(direction>0) {
								diff = val - UB;
								if (diff > 0){
									if constexpr (level==rank-1){ //Not rlevel because we are checking for finishing the last level. Equivalent would be rlevel=(order==0)?rank-1:0;
										//Have gone beyond the end of the array
										//Set finished and doine
										val = UB+1;//One after last value
										return true;
									} else {
										FAR_SIGNED_INDEX_TYPE next_delta = diff/sz+1;
										bool finished = advanceTuple<direction, level+1>(index, direction*next_delta);
										val=finished?LB:val-next_delta * sz;
										return finished;
									}
								}
							} else {
								//Deal with val < LB
								diff = val - LB;
								if (diff < 0) {
									if constexpr (level==rank-1){
										//Have gone beyond the end of the array
										//Set finished and done
										val = LB-1;//One before first value
										return true;
									} else {
										FAR_SIGNED_INDEX_TYPE next_delta = diff/(FAR_SIGNED_INDEX_TYPE)sz+1;
										bool finished = advanceTuple<direction, level+1>(index, -direction*next_delta);
										val=finished?UB:val+next_delta * sz;
										return finished;
									}
								}
							}
							return false;
						}

					template<typename T_index, int level=0>
						FAR_INLINE void advanceTuple2(T_index &index, FAR_SIGNED_INDEX_TYPE delta) const{
							FAR_SIGNED_INDEX_TYPE &val = std::get<level>(index);
							FAR_SIGNED_INDEX_TYPE N = [this](){
								if constexpr(level!=rank-1){
									return n_subsize[level];
								} else {
									return n_subsize[level]+1;
								}
							}();
							FAR_SIGNED_INDEX_TYPE newVal = val + delta - defaultLB;
							if constexpr(level<rank-1){
								//Bit twiddling hack gets the sign bit from the result of the modulo division.
								//Subtracts one if the result is negative
								FAR_SIGNED_INDEX_TYPE carry = newVal / N - ((newVal % N >> (sizeof(FAR_SIGNED_INDEX_TYPE) * CHAR_BIT - 1)) & 1);
								advanceTuple2<T_index,level+1>(index, carry);
							}
							//Map val back onto 0->N-1 
							val = (newVal % N + N) % N + defaultLB;
						}

				};

			//END contiguousIndexBuilder

			template <int rank, bounds_check_state bounds_check = bounds_check_state::bc_default>
				using FortranIndexBuilder = contiguousIndexBuilder<rank, bounds_check, 0>;

			template <int rank, bounds_check_state bounds = bounds_check_state::bc_default>
				using CIndexBuilder = contiguousIndexBuilder<rank, bounds, 1>;


			template<typename T, typename... T_others>
				auto getIndexer(T&& current, T_others&&... others){
					constexpr int maxRank =  arrayRankRange<T,T_others...>::maxRank;
					if constexpr(arrayInfo<T>::value && arrayInfo<T>::rank == maxRank){
						return &current.getIndexer();
					} else if constexpr(sizeof...(others)>0) {
						return getIndexer(others...);
					} else {
						return nullFunction;
					}
				}

			auto& fortFile(int n);
			//Produce an array the same size as a given array with one dimension removed
			template <typename T_result = void, typename T=double>
				auto reduce_rank(T &&source, int dir);

			template<typename T_dest, typename T_src>
				FAR_INLINE void forEachCopyTight(T_dest &&dest, T_src &&src){
					using T_dest_item = typename arrayInfo<T_dest>::wrappedType;
					using T_src_item = typename arrayInfo<T_src>::wrappedType;
					static_assert(std::is_assignable_v<T_dest_item,T_src_item> || std::is_constructible_v<T_dest_item,T_src_item>,"Cannot copy from source to destination");
					auto index = getIndexer(dest,src);
					constexpr int rank = std::remove_pointer_t<decltype(index)>::core_rank;
					typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
					//Simple assignment is much more effective with unit strides, so override it here
					index->template forEachCore<
						[](T_dest_item &dest, std::conditional_t<std::is_pointer_v<T_src_item>,T_src_item, const T_src_item&>src){
							if constexpr(std::is_assignable_v<T_dest_item,T_src_item>){
								dest = src;
							} else {
								dest.~T_dest_item();
								new (&dest) T_dest_item(src);
							}
						}
						,minimalNMapNoRefUnwrap<indexInfo<rank>,T_dest,T_src>, false, false, 0
							,2,FAR_SIMPLELENGTH,FAR_SIMPLELENGTH //This line sets use of unit strides
							>(t1,t2,rdNone,std::forward<T_dest>(dest), std::forward<T_src>(src));
				}

			template<typename T_dest, typename T_src>
				FAR_INLINE void forEachCopyTightSerial(T_dest &&dest, T_src &&src){
					using T_dest_item = typename arrayInfo<T_dest>::wrappedType;
					using T_src_item = typename arrayInfo<T_src>::wrappedType;
					auto index = getIndexer(dest,src);
					constexpr int rank = std::remove_pointer_t<decltype(index)>::core_rank;
					typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
					//Simple assignment is much more effective with unit strides, so override it here
					index->template forEachCore<
						[](T_dest_item &dest, std::conditional_t<std::is_pointer_v<T_src_item>,T_src_item,T_src_item&>src){
								if constexpr(std::is_assignable_v<T_dest_item,T_src_item>){
									dest = src;
								} else {
									dest.~T_dest_item();
									new (&dest) T_dest_item(src);
								}
						},
						minimalNMapNoRefUnwrap<indexInfo<rank>,T_dest,T_src>, false, false, 0
							,1,1,1,false //This line sets use of unit strides and disables parallelism
							>(t1,t2,rdNone,std::forward<T_dest>(dest), std::forward<T_src>(src));
				}

			template<auto callable, typename... T_params>
				FAR_INLINE void forEachTight(T_params&&... params){
					auto indexer = getIndexer(params...);
					constexpr int rank = std::decay_t<decltype(*indexer)>::core_rank;
					typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
					indexer->template forEachCore<
						callable,
						minimalNMap<indexInfo<rank>,T_params...>, false, false, 0
							,1,1,1 //This line sets use of unit strides
							>(t1,t2,rdNone,std::forward<T_params>(params)...);
				}

			template<auto callable, bool backwards = false, typename... T_params>
				FAR_INLINE void forEach(T_params&&... params){
					auto indexer = getIndexer(params...);
					if constexpr(!std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>){
						constexpr int rank = std::decay_t<decltype(*indexer)>::core_rank;
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
						//printPack(params...);
						indexer->template forEachCore<callable, minimalNMap<indexInfo<rank>,T_params...>,false,backwards>(t1,t2,rdNone,std::forward<T_params>(params)...);
					} else {
						//Applying forEach to scalar
						callable(std::forward<T_params>(params)...);
					}
				}
			template<auto callable, auto mapper, bool backwards = false, typename... T_params>
				FAR_INLINE void forEachMapper(T_params&&... params){
					auto indexer = getIndexer(params...);
					if constexpr(!std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>){
						constexpr int rank = std::decay_t<decltype(*indexer)>::core_rank;
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
						indexer->template forEachCore<callable, mapper,false,backwards>(t1,t2,rdNone,params...);
					} else {
						static_assert(!std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>,"Unable to find indexer in forEach");
					}
				}

			template<auto callable, bool backwards = false, typename... T_params>
				FAR_INLINE void forEachSerial(T_params&&... params){
					auto indexer = getIndexer(params...);
					if constexpr(!std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>){
						constexpr int rank = std::decay_t<decltype(*indexer)>::core_rank;
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
						indexer->template forEachCore<
							callable,
							minimalNMap<indexInfo<rank>,T_params...>, false, backwards, 0
								,1,1,1,false //This line sets use of unit strides
								>(t1,t2,rdNone,std::forward<T_params>(params)...);
					} else {
						static_assert(!std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>,"Unable to find indexer in forEachSerial");
					}
				}

			template<auto callable, auto mapper, bool backwards = false, typename... T_params>
				FAR_INLINE void forEachSerialMapper(T_params&&... params){
					auto indexer = getIndexer(params...);
					if constexpr(!std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>){
						constexpr int rank = std::decay_t<decltype(*indexer)>::core_rank;
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
						indexer->template forEachCore<
							callable,
							mapper, false, backwards, 0
								,1,1,1,false //This line sets use of unit strides
								>(t1,t2,rdNone,std::forward<T_params>(params)...);
					} else {
						static_assert(!std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>,"Unable to find indexer in forEachSerial");
					}
				}				

			template<bool backwards = false, typename T_runtime, typename... T_params>
				FAR_INLINE void forEachRuntime(T_runtime &runtime, T_params&&... params){
					auto indexer = getIndexer(params...);
					if constexpr(!std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>){
						constexpr int rank = std::decay_t<decltype(*indexer)>::core_rank;
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
						auto callable = [](T_runtime &runtime, T_params&&... params){runtime(params...);};
						indexer->template forEachCore<callable, minimalNMap<indexInfo<rank>,T_params...>,false,backwards>(t1,t2,rdNone,runtime,params...);
					} else {
						static_assert(std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>,"Unable to find indexer in forEachRuntime");
					} 
				}

			template<auto callable, bool backwards = false, typename... T_params>
				FAR_INLINE void forEachReduce(int direction, T_params&&... params){
					auto indexer = getIndexer(params...);
					if constexpr(!std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>){
						constexpr int rank = std::decay_t<decltype(*indexer)>::core_rank;
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
						indexer->template forEachCore<callable, minimalNMap<indexInfo<rank>,T_params...>,true,backwards>(t1,t2,direction-1,params...);
					} else {
						static_assert(std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>,"Unable to find indexer in forEachReduce");
					} 
				}

			template<auto callable, auto mapper, bool backwards = false, typename... T_params>
				FAR_INLINE void forEachReduceMapper(int direction, T_params&&... params){
					auto indexer = getIndexer(params...);
					if constexpr(!std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>){
						constexpr int rank = std::decay_t<decltype(*indexer)>::core_rank;
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
						indexer->template forEachCore<callable, mapper,true,backwards>(t1,t2,direction-1,params...);
					} else {
						static_assert(std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>,"Unable to find indexer in forEachReduce");
					} 
				}

			template<auto callable, bool backwards = false, typename... T_params>
				FAR_INLINE void forEachReduceSerial(int direction, T_params&&... params){
					auto indexer = getIndexer(params...);
					if constexpr(!std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>){
						constexpr int rank = std::decay_t<decltype(*indexer)>::core_rank;
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
						indexer->template forEachCore<callable, minimalNMap<indexInfo<rank>,T_params...>,true,backwards,
							0,FAR_ET_LENGTH,FAR_NON_ET_LENGTH,FAR_SIMPLELENGTH,false>(t1,t2,direction-1,params...);
					} else {
						static_assert(std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>,"Unable to find indexer in forEachReduceSerial");
					}
				}

			template<auto callable, auto mapper, bool backwards = false, typename... T_params>
				FAR_INLINE void forEachReduceSerialMapper(int direction, T_params&&... params){
					auto indexer = getIndexer(params...);
					if constexpr(!std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>){
						constexpr int rank = std::decay_t<decltype(*indexer)>::core_rank;
						typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
						indexer->template forEachCore<callable, mapper,true,backwards,
							0,FAR_ET_LENGTH,FAR_NON_ET_LENGTH,FAR_SIMPLELENGTH,false>(t1,t2,direction-1,params...);
					} else {
						static_assert(std::is_same_v<std::decay_t<decltype(indexer)>,nullFunctionType>,"Unable to find indexer in forEachReduceSerial");
					}
				}

			/*
			 * Backend for reduction functions
			 * Deals with parallel and serial reductions
			 */
			template<auto parallelCombine, 
				auto serialCombine=parallelCombine, 
				auto initialise=nullFunction,
				auto mapper=defaultFunction,
				bool back=false,
				typename T_dest = typename std::decay_t<std::tuple_element_t<std::is_same_v<decltype(mapper),nullFunctionType>?2:1,typename callableTraits<decltype(parallelCombine)>::params>>,
				typename T_array, typename... T_others>
					FAR_INLINE auto reduction(T_array&& array, T_others&&... others){

						if constexpr (arrayInfo<std::decay_t<T_array>>::rank == 0){
							return array;
						} else {

#ifdef FAR_PARALLEL_REDUCE
							using dtype = T_dest;
							if (workshare[getTeamNum()] && omp_get_num_threads()>1){
								sharedVariable<dtype> acc;
								dtype dest;
								if constexpr(std::is_same_v<decltype(initialise),nullFunctionType>){
									//Null initialisation means, just don't
									//i.e. default construct
								} else if constexpr(std::is_same_v<decltype(initialise),defaultFunctionType>){
									//Default initialisation means value initialise
									dest={};
								} else {
									//Otherwise call the initialisation function
									initialise(dest);
								}
								if constexpr(!std::is_same_v<decltype(mapper),defaultFunctionType>){
									forEachMapper<parallelCombine,mapper,back>(array,dest,others...);
								} else {
									forEach<parallelCombine,back>(array,dest,others...);
								}
								dtype &result = acc.getData();
#pragma omp single
								{
									if constexpr(std::is_same_v<decltype(initialise),nullFunctionType>){
									}
									else if constexpr(std::is_same_v<decltype(initialise),defaultFunctionType>){
										result ={};
									} else {
										initialise(result);
									}
								}
#pragma omp critical
								{
									serialCombine(dest,result);
								}
								FAR_OMP_BARRIER
									return result;
							} else
#endif
							{
								return reductionSerial<parallelCombine,initialise,mapper,back,T_dest>(array,others...);
							}
						}
					}

			template<auto parallelCombine, 
				auto initialise=nullFunction,
				auto mapper=defaultFunction,
				bool back=false,
				typename T_dest = typename std::decay_t<std::tuple_element_t<std::is_same_v<decltype(mapper),nullFunctionType>?
					2:1,typename callableTraits<decltype(parallelCombine)>::params>>,
				typename T_array, typename... T_others>
					FAR_INLINE auto reductionSerial(T_array&& array, T_others&&... others){
						using dtype = T_dest;
						dtype dest;
						if constexpr(std::is_same_v<decltype(initialise),nullFunctionType>){
							//Do nothing
						} else if constexpr(std::is_same_v<decltype(initialise),defaultFunctionType>){
							//Value initialise
							dest = {};
						} else {
							//Call the initialisation function
							initialise(dest);
						}
						if constexpr(!std::is_same_v<decltype(mapper),defaultFunctionType>){
							forEachMapper<parallelCombine,mapper,back>(array,dest,others...);
						} else {
							forEach<parallelCombine,back>(array,dest,others...);
						}
						return dest;
					}

			/**
			 * Backend for reduction functions with a direction
			 * Deals with parallel and serial reductions
			 */
			template<auto parallelCombine, auto serialCombine=parallelCombine, 
				auto initialise=defaultFunction,
				auto mapper = defaultFunction,
				bool back=false,
				typename T_dest = typename std::decay_t<std::tuple_element_t<std::is_same_v<decltype(mapper),nullFunctionType>?2:1,typename callableTraits<decltype(parallelCombine)>::params>>,
				typename T_array, typename... T_others>
					FAR_INLINE auto reductionWithDirection(int direction, T_array&& array, T_others&&... others){

#ifdef FAR_PARALLEL_REDUCE
						using type = typename arrayInfo<std::decay_t<T_array>>::type;
						using dtype = std::conditional_t<std::is_void_v<T_dest>,type, typename arrayInfo<T_dest>::type>;
						if (workshare[getTeamNum()] && omp_get_num_threads()>1){
							//Fast parallel reduction creates copies for every thread of the reduction array
							using rtype = decltype(reduce_rank<dtype>(array,direction));
							perThreadVariable<rtype> mt;
							rtype &inter= mt.getData();
							sharedVariable<rtype> sres;
							rtype &result = sres.getData();
							result.setShared(true);
							inter.mold(reduce_rank(array,direction));
							result.mold(inter);
							if constexpr(std::is_same_v<decltype(initialise),defaultFunctionType>){
								forEach<[](type &src){src={};}>(result);
							} else {
								forEach<initialise>(result);
							}
							if constexpr(std::is_same_v<decltype(initialise),defaultFunctionType>){
								forEachSerial<[](type &src){src={};}>(inter);
							} else {
								forEachSerial<initialise>(inter);
							}
							if constexpr(std::is_same_v<decltype(mapper),defaultFunctionType>){
								forEachReduce<parallelCombine,back>(direction,array,inter,others...);
							} else {
								forEachReduceMapper<parallelCombine,mapper,back>(direction,array,inter,others...);
							}
							FAR_OMP_BARRIER
								rtype *all = mt.getAllData();
							for (int i=0;i<omp_get_num_threads();++i){
								if constexpr(std::is_base_of_v<array_base,std::decay_t<std::tuple_element_t<0,typename callableTraits<decltype(serialCombine)>::params>>>){
									serialCombine(all[i],result);
								} else {
									forEach<serialCombine>(all[i],result);
								}
							}
							FAR_OMP_BARRIER
								return result;
						} else
#endif
							//If we are either in serial or not in a workshare directive then the
							//current thread must do all of the work
						{
							return reductionSerialWithDirection<parallelCombine,initialise,mapper,back,T_dest>(direction,array,others...);

						}
					}

			/**
			 * Backend for reduction functions with a direction
			 * Deals with parallel and serial reductions
			 */
			template<auto parallelCombine, auto serialCombine=parallelCombine, 
				auto initialise=defaultFunction,
				auto initialise_meta = defaultFunction,
				auto mapper = defaultFunction,
				bool back=false,
				typename T_dest = typename std::decay_t<std::tuple_element_t<std::is_same_v<decltype(mapper),nullFunctionType>?2:1,typename callableTraits<decltype(parallelCombine)>::params>>,
				typename T_meta = typename std::decay_t<std::tuple_element_t<std::is_same_v<decltype(mapper),nullFunctionType>?3:2,typename callableTraits<decltype(parallelCombine)>::params>>,
				typename T_array, typename... T_others>
					FAR_INLINE auto reductionWithDirectionAndMeta(int direction, T_array&& array, T_others&&... others){
#ifdef FAR_PARALLEL_REDUCE
						using type = typename arrayInfo<std::decay_t<T_array>>::type;
						using meta_type = typename arrayInfo<std::decay_t<T_meta>>::type;
						using dtype = std::conditional_t<std::is_void_v<T_dest>,type, typename arrayInfo<T_dest>::type>;
						if (workshare[getTeamNum()] && omp_get_num_threads()>1){
							//Fast parallel reduction creates copies for every thread of the reduction array
							using rtype = decltype(reduce_rank<dtype>(array,direction));
							using mtype = decltype(reduce_rank<meta_type>(array,direction));
							perThreadVariable<rtype> mt;
							//Inter is intermediate result of the reduction
							rtype &inter= mt.getData();

							perThreadVariable<mtype> mmt;
							//Inter is intermediate metadata of the reduction
							mtype &minter= mmt.getData();

							//This is the shared value on all threads of the final result of the reduction
							sharedVariable<rtype> sres;
							rtype &result = sres.getData();
							result.setShared(true);

							//This is the shared value on all threads of the final metadata of the reduction
							sharedVariable<mtype> smeta;
							mtype &resultMeta = smeta.getData();
							resultMeta.setShared(true);

							//Allocate the arrays
							inter.mold(reduce_rank(array,direction));
							minter.mold(inter);
							result.mold(inter);
							resultMeta.mold(inter);

							//Initialise the final result array
							if constexpr(std::is_same_v<decltype(initialise),defaultFunctionType>){
								forEach<[](type &src){src={};}>(result);
							} else {
								forEach<initialise>(result);
							}

							//Initialise the final metadata array
							if constexpr(std::is_same_v<decltype(initialise_meta),defaultFunctionType>){
								forEach<[](meta_type &src){src={};}>(resultMeta);
							} else {
								forEach<initialise_meta>(resultMeta);
							}

							//Initialise the intermediate result array
							if constexpr(std::is_same_v<decltype(initialise),defaultFunctionType>){
								forEachSerial<[](type &src){src={};}>(inter);
							} else {
								forEachSerial<initialise>(inter);
							}
							//Initialise the intermediate metadata
							if constexpr(std::is_same_v<decltype(initialise_meta),defaultFunctionType>){
								forEachSerial<[](mtype &src){src={};}>(minter);
							} else {
								forEachSerial<initialise_meta>(minter);
							}
							//Perform the local reduction on this thread
							if constexpr(std::is_same_v<decltype(mapper),defaultFunctionType>){
								forEachReduce<parallelCombine,back>(direction,array,inter, minter, others...);
							} else {
								forEachReduceMapper<parallelCombine,mapper,back>(direction,array,inter, minter, others...);
							}
							FAR_OMP_BARRIER
								rtype *all = mt.getAllData();
							mtype *mall = mmt.getAllData();
							for (int i=0;i<omp_get_num_threads();++i){
								if constexpr(std::is_base_of_v<array_base,std::decay_t<std::tuple_element_t<0,typename callableTraits<decltype(serialCombine)>::params>>>){
									serialCombine(all[i], mall[i], result, resultMeta);
								} else {
									forEach<serialCombine>(all[i], mall[i], result, resultMeta);
								}
							}
							FAR_OMP_BARRIER
								return result;
						} else
#endif
							//If we are either in serial or not in a workshare directive then the
							//current thread must do all of the work
						{
							return reductionSerialWithDirectionAndMeta<parallelCombine,initialise,initialise_meta,mapper,back,T_dest>(direction,array,others...);

						}
					}

			template<auto parallelCombine, 
				auto initialise=defaultFunction,
				auto mapper = defaultFunction,
				bool back=false,
				typename T_dest = typename std::decay_t<std::tuple_element_t<std::is_same_v<decltype(mapper),nullFunctionType>?2:1,typename callableTraits<decltype(parallelCombine)>::params>>,
				typename T_array, typename... T_others>
					FAR_INLINE auto reductionSerialWithDirection(int direction, T_array&& array, T_others&&... others){
						using type = typename arrayInfo<std::decay_t<T_array>>::type;
						using dtype = std::conditional_t<std::is_void_v<T_dest>,type, typename arrayInfo<T_dest>::type>;
						auto result = reduce_rank<dtype>(array, direction);
						if constexpr(std::is_same_v<decltype(initialise),defaultFunctionType>){
							forEach<[](type &src){src={};}>(result);
						} else {
							forEach<initialise>(result);
						}

						if constexpr(std::is_same_v<decltype(mapper),defaultFunctionType>){
							forEachReduceSerial<parallelCombine,back>(direction,array,result,others...);
						} else {
							forEachReduceSerialMapper<parallelCombine,mapper,back>(direction,array,result,others...);
						}
						return result;
					}


			template<auto parallelCombine, 
				auto initialise=defaultFunction,
				auto initialise_meta = defaultFunction,
				auto mapper = defaultFunction,
				bool back=false,
				typename T_dest = typename std::decay_t<std::tuple_element_t<std::is_same_v<decltype(mapper),nullFunctionType>?2:1,typename callableTraits<decltype(parallelCombine)>::params>>,
				typename T_meta = typename std::decay_t<std::tuple_element_t<std::is_same_v<decltype(mapper),nullFunctionType>?3:2,typename callableTraits<decltype(parallelCombine)>::params>>,
				typename T_array, typename... T_others>
					FAR_INLINE auto reductionSerialWithDirectionAndMeta(int direction, T_array&& array, T_others&&... others){
						using type = typename arrayInfo<std::decay_t<T_array>>::type;
						using mtype = typename arrayInfo<std::decay_t<T_meta>>::type;
						using dtype = std::conditional_t<std::is_void_v<T_dest>,type, typename arrayInfo<T_dest>::type>;
						auto result = reduce_rank<dtype>(array, direction);
						auto meta = reduce_rank<mtype>(array, direction);
						if constexpr(std::is_same_v<decltype(initialise),defaultFunctionType>){
							forEach<[](type &src){src={};}>(result);
						} else {
							forEach<initialise>(result);
						}
						if constexpr(std::is_same_v<decltype(initialise_meta),defaultFunctionType>){
							forEach<[](mtype &src){src={};}>(meta);
						} else {
							forEach<initialise_meta>(meta);
						}
						if constexpr(std::is_same_v<decltype(mapper),defaultFunctionType>){
							forEachReduceSerial<parallelCombine,back>(direction,array,result ,meta, others...);
						} else {
							forEachReduceSerialMapper<parallelCombine,mapper,back>(direction,array,result,meta,others...);
						}
						return result;
					}

			template <typename T, int rank, bounds_check_state bounds_check = bounds_check_state::bc_default, template <int, bounds_check_state> class indexer = FAR_DEFAULT_INDEX_BUILDER, bool forceContiguous=false>
				class Array;

#define ITERATORFAR_INLINE __attribute__((flatten)) inline

			template <typename T_core, FAR_SIGNED_INDEX_TYPE mult=1>
				class array_iterator{
					friend T_core;
					static const int rank = T_core::core_rank;
					using return_type = T_core::core_type;
					private:

					FAR_SIGNED_INDEX_TYPE offset=-1;

					template<int level>
						void init_lb(){
							offset=0;
						}

					template<int level>
						void init_ub(){
							offset = core.getSize()-1;
						}

					const T_core &core;
					public:

					using value_type = std::remove_reference_t<return_type>;
					using difference_type = FAR_SIGNED_INDEX_TYPE;
					using pointer = std::remove_reference_t<return_type>*;
					using reference = std::remove_reference_t<return_type>&;
					using iterator_category = std::random_access_iterator_tag;

					ITERATORFAR_INLINE decltype(auto) operator*(){
						return core.getItem(offset);
					}

					array_iterator(const array_iterator &other) = default;
					array_iterator & operator=(const array_iterator &other){
						offset = other.offset;
						return *this;
					}

					ITERATORFAR_INLINE decltype(auto) operator*() const{
						return core.getItem(offset);
					}
					ITERATORFAR_INLINE decltype(auto) operator->(){
						return core->getItem(offset);
					}
					ITERATORFAR_INLINE auto operator+=(FAR_SIGNED_INDEX_TYPE shift) noexcept {
						offset+=shift*mult;
						return *this;
					}
					ITERATORFAR_INLINE auto operator++(int){
						auto a=*this;
						(*this)+=1;
						return a;
					}
					ITERATORFAR_INLINE auto operator++(){
						(*this)+=1;
						return *this;
					}

					ITERATORFAR_INLINE auto operator-=(FAR_SIGNED_INDEX_TYPE shift){
						offset-=shift*mult;
						return *this;
					}

					ITERATORFAR_INLINE auto operator--(int){
						auto a=*this;
						(*this)-=1;
						return a;
					}
					ITERATORFAR_INLINE auto operator--(){
						(*this)-=1;
						return *this;
					}

					ITERATORFAR_INLINE auto operator-(array_iterator& other){
						return offset - other.offset;
					}

					ITERATORFAR_INLINE array_iterator(const T_core &src) : core(src){init_lb<0>();}
					ITERATORFAR_INLINE array_iterator(const T_core &src, bool backwards) : core(src){
						if (!backwards) {
							init_lb<0>();
						} else {
							init_ub<0>();
						}
					}

					ITERATORFAR_INLINE bool operator==(array_iterator &other){
						return offset == other.offset;
					}

					ITERATORFAR_INLINE bool operator!=(array_iterator &other){
						return !(*this==other);
					}

					ITERATORFAR_INLINE bool operator<(array_iterator &other){
						return offset < other.offset;
					}

					ITERATORFAR_INLINE bool operator<=(array_iterator &other){
						return offset <= other.offset;
					}

					ITERATORFAR_INLINE bool operator>(array_iterator &other){
						return offset > other.offset;
					}

					ITERATORFAR_INLINE bool operator>=(array_iterator &other){
						return offset > other.offset;
					}

					ITERATORFAR_INLINE auto operator+(FAR_SIGNED_INDEX_TYPE offset){
						auto a = *this;
						a+=offset * mult;
						return a;
					}

					ITERATORFAR_INLINE auto operator-(FAR_SIGNED_INDEX_TYPE offset){
						auto a = *this;
						a-=offset * mult;
						return a;
					}

					ITERATORFAR_INLINE decltype(auto) operator[](FAR_SIGNED_INDEX_TYPE offset){
						auto a = *this;
						a+=offset *  mult;
						return *a;
					}

				};

			/** 
			 * Class representing an elemental function
			 */
			template <auto callable, auto combination=defaultFunction, bool partialReduce=false, auto initializer=defaultFunction>
				struct elementalFunction{
					private:
						using callable_params = typename callableTraits<decltype(callable)>::params;
						using rtype = typename callableTraits<decltype(callable)>::type;
						static constexpr bool isReducing = callableTraits<decltype(combination)>::value;
						static constexpr bool customInit = callableTraits<decltype(initializer)>::value;

						int reduction=0;
						bool backwards=false;

						//This is a helper function to find information about the first unmatched array in the set of parameters
						//that the elemental function are called with
						template <int level, auto callable2, typename... Ts>
							struct first_base;

						template <int level, auto callable2, typename T, typename... Ts>
							struct first_base<level, callable2, T, Ts...> {
								constexpr static bool isRRI = isReductionIndexInfo<T>::value;
								constexpr static bool hasRRI = isRRI || first_base<level+1, callable2, Ts...>::hasRRI;
								//Point of condition is to detect an "elemental" parameter, so an array specified in the outer function
								//that is not an array in the inner function
								constexpr static bool condition = is_unmatched_array<level,T,callable2>::value;
								//Maxrank and minrank contain the largest and smallest ranks for non "true" parameters
								static constexpr int maxrank = condition?std::max(arrayInfo<T>::rank,first_base<level+1, callable2, Ts...>::maxrank)
									:first_base<level+1, callable2, Ts...>::maxrank;
								static constexpr int minrank = condition?std::min(arrayInfo<T>::rank,first_base<level+1, callable2, Ts...>::minrank)
									:first_base<level+1, callable2, Ts...>::minrank;
								static constexpr int truerank = isRRI?(isReductionIndexInfo<T>::rank):hasRRI?first_base<level+1, callable2, Ts...>::truerank:maxrank;
								//This is the level of the first elemental parameter
								static constexpr int alevel = condition?level:first_base<level+1, callable2, Ts...>::alevel;
								using array_type = std::conditional_t<condition, std::decay_t<T>, typename first_base<level+1, callable2, Ts...>::array_type>;
								//Datatype of first elemental parameter
								using type = std::conditional_t<condition, typename arrayInfo<T>::type, typename first_base<level+1, callable2, Ts...>::type>;
								using index = std::conditional_t<condition, typename arrayInfo<T>::index, typename first_base<level+1, callable2, Ts...>::index>;
							};

						template <int level, auto callable2, typename T>
							struct first_base<level, callable2, T> {
								constexpr static bool isRRI = isReductionIndexInfo<T>::value;
								constexpr static bool hasRRI = isRRI;
								constexpr static bool condition = is_unmatched_array<level,T,callable2>::value;
								static constexpr int maxrank = condition?arrayInfo<T>::rank:0;
								static constexpr int minrank = condition?arrayInfo<T>::rank:std::numeric_limits<int>::max();
								static constexpr int truerank = isRRI?(isReductionIndexInfo<T>::rank):maxrank;
								static constexpr int alevel = condition?level:-1;
								using array_type = std::decay_t<T>;
								using type = typename arrayInfo<T>::type;
								using index = typename arrayInfo<T>::index;
							};

						template<int level, typename T_dest, typename T_cur, typename... T_others>
							void sizeArray(T_dest&dest, T_cur&& current, T_others&&... params){
								if constexpr(is_unmatched_array<level,T_cur,callable>::value){
									if constexpr(isReducing){
										dest = reduce_rank<typename arrayInfo<T_dest>::type>(current,reduction);
									} else {
										dest.mold(current);
									}
								} else if constexpr(sizeof...(T_others)>0){
									sizeArray<level+1>(dest, params...);
								}
							}

						/**
						 * Function to find the first index object for an array that is unmatched in the inner
						 * array. Should not be called if there are no arrays to find values for
						 */
						template<int level=0, int truerank, typename T_cur=int, typename... T_others>
							auto& getIndex(T_cur&& current, T_others&&... params){
								if constexpr(is_unmatched_array<level,T_cur,callable>::value &&
										arrayInfo<T_cur>::rank == truerank){
									return current.getIndexer();
								} else if constexpr(sizeof...(T_others)>0){
									return getIndex<level+1,truerank>(params...);
								}
							}

					public:
						/**
						 * Function to set which axis reduction should occur over if partial reduction is occuring
						 */
						void setReductionAxis(int reduction=1){
							this->reduction = reduction-1;
						}

						/**
						 * Funxtion to set whether the elements of the arrays should be iterated through backwards
						 */
						void setBackwards(bool backwards=false){this->backwards=backwards;}

						/**
						 * Round brackets operator. Meat of the elemental function. Takes the parameters that this is called with
						 * and determines what to do with them. Variables that are non-arrays are passed through unchanged, as are
						 * arrays that are also arrays in the inner function. Arrays that are passed in here that are not matched
						 * by arrays in the inner function are iterated over and the inner function called for each element.
						 * All arrays must be the same size and rank, unless reduction is active in which case arrays of rank one less
						 * than the largest array in the parameters are accessed with the reduction access removed from their indices
						 */
						template<typename... T>
							decltype(auto) operator()(T&&... params){
								constexpr bool returns = !std::is_void_v<typename callableTraits<decltype(callable)>::type>;
								//This branch triggers if the function is being called non-elementally
								if constexpr (first_base<0,callable,T...>::alevel==-1){
									auto inter = std::make_tuple();
									if constexpr(isReducing){
										using rt = typename arrayInfo<std::decay_t<std::tuple_element_t<0,typename callableTraits<decltype(combination)>::params>>>::type;
										rt result;
										if constexpr(customInit){
											initializer(result);
										} else {
											result = {};
										}
										if constexpr(std::tuple_size_v<typename callableTraits<decltype(combination)>::params> == 2){
											combination(result, std::apply(callable,buildCallTuple<callable, 0>(indexInfo<0>(inter,reduction), params...)));
										} else {
											auto m = std::make_tuple();
											combination(result, std::apply(callable,buildCallTuple<callable, 0>(indexInfo<0>(inter,reduction), params...)),
													indexInfo<0>(m,reduction));
										}
										return result;
									} else {
										auto ii=indexInfo<0>(inter,reduction);
										return std::apply(callable,buildCallTuple<callable, 0>(ii, params...));
									}
								} else {
									constexpr bool localFullReduce = isReducing && (first_base<0,callable,T...>::maxrank == 1 || !partialReduce);
									constexpr int rank = first_base<0,callable,T...>::truerank;
									if constexpr(returns){
										using rt = std::conditional_t<isReducing,
													typename arrayInfo<std::decay_t<std::tuple_element_t<0,typename callableTraits<decltype(combination)>::params>>>::type,
													rtype>;
										//Elemental function

										std::conditional_t<localFullReduce, rt, 
											std::invoke_result_t<decltype(first_base<0, callable, T...>::array_type:: template buildVariant<
													rt, isReducing?rank-1:rank>)>> result;

										auto &index = getIndex<0,first_base<0,callable,T...>::truerank>(params...);
										// Deal with reduction
										if constexpr(!localFullReduce) sizeArray<0>(result,params...);
										if constexpr(isReducing) {
											if constexpr(customInit){
												if constexpr(localFullReduce) {
													//Full reduction with initializerr
													initializer(result);
												} else {
													//Partial reduction with initializer
													for (auto& el:result) {initializer(el);}
												}
											} else {
												if constexpr(localFullReduce){
													//Full reduction with value initialization
													result = {};
												} else {
													//Partial reduction with value initialization
													for (auto& el:result) {el={};}
												}
											}
										}
										typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
										using ptype = typename callableTraits<decltype(callable)>::params;
										if constexpr(!isReducing){
											auto l = [](indexInfo<rank> &rii, decltype(result)& result, T&&... others){
												return std::make_tuple(buildRef<callable,-1>(rii, result),
														buildCallTuple<callable, 0,indexInfo<rank>,T...>(rii,std::forward<T>(others)...));
											};
											if (backwards){
												index.template forEachCore<wrapReturn<callable,ptype>,l,false,true>(t1, t2, reduction, result, std::forward<T>(params)...);
											} else {
												index.template forEachCore<wrapReturn<callable,ptype>,l,false,false>(t1, t2, reduction, result, std::forward<T>(params)...);
											}
										} else {
											auto l = [](indexInfo<rank> &rii, decltype(result)& result, T&&... others){
												return std::tuple_cat(std::make_tuple(buildRef<callable,-1>(rii, result)),
														buildCallTuple<callable, 0,indexInfo<rank>,T...>(rii,std::forward<T>(others)...));
											};
											if (backwards){
												index.template forEachCore<wrapReturn<callable,ptype,combination>,l,true,true>(t1, t2, reduction, result, 
														std::forward<T>(params)...);
											} else {
												index.template forEachCore<wrapReturn<callable,ptype,combination>,l,true,false>(t1, t2, reduction, result, 
														std::forward<T>(params)...);
											}
										}
										return result;
									} else {
										//Elemental subroutine
										auto &index = getIndex<0,first_base<0,callable,T...>::truerank>(params...);
										typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
										if (backwards){
											index.template forEachCore<callable,buildCallTuple<callable, 0,indexInfo<rank>,T...>,false,true>(t1,t2,reduction,
													std::forward<T>(params)...);
										} else {
											index.template forEachCore<callable,buildCallTuple<callable, 0,indexInfo<rank>,T...>,false,false>
												(t1,t2,reduction,std::forward<T>(params)...);
										}
									}
								}
							}
				};

			/**
			 * Class representing a "lazy" array. Implements the same operations as a real array
			 * but only calculates values when needed rather than storing them in memory
			 * Note that it is templated on an INSTANCE of a callable, NOT a class
			 */
			template <auto callable, typename lindexer, typename... types>
				struct lazyArray : public array_operator
			{
				template <typename T2, int rank2, bounds_check_state bounds_check2, template <int, bounds_check_state> class indexer2, bool forceContiguous2>
					friend class Array;
				public:
				/**
				 * Tuple type representing the types of all of the types passed to the template for the type
				 * if any of them are of type array or array_operator they are stored as pointers
				 */
				using op_parameters = std::tuple<arrayInfo_st<types>...>;
				using op_bare = std::tuple<types...>;
				/** 
				 * Retun type of the callable object that this was created with
				 */
				using return_type = typename callableTraits<decltype(callable)>::type;
				/**
				 * Types of the callable object parameters
				 */
				using callable_params = typename callableTraits<decltype(callable)>::params;
				/**
				 * Type of the indexer
				 */
				using index_type = lindexer;

				/**
				 * Rank of the lazy array
				 */
				static const int core_rank = lindexer::core_rank;

				public:
				/** Tuple containing the bound operands to this lazy array
				 * NOTE : May be rvalue references
				 */
				op_parameters operands;

				bool shouldDelete[std::tuple_size_v<op_parameters>]={};
				bool rValue[std::tuple_size_v<op_parameters>]={};

				/**
				 * Indexer object used to access elements
				 */
				lindexer index;

				bool index_shifted=false;

				/**
				 * Build an index for this object from another index
				 * Rebase it to defaultLB because all lazy arrays are defaultLB based
				 */

				template <typename T_index>
					void buildIndex(T_index& other_index)
					{
						index = other_index;
						index.normalise();
						index.rebase(defaultLB);
					}

				template<int level=0, typename T_src>
					void copyOperands(const T_src&& src){
						std::get<level>(operands) = std::get<level>(src.operands);
						src.shouldDelete[level]=false;
						if constexpr(level<std::tuple_size_v<decltype(operands)>-1) copyOperands<level+1>(src);
					}

				template<int level=0, typename T_src>
					void copyOperands(const T_src& src){
						if (!src.rValue[level]) {
							std::get<level>(operands) = std::get<level>(src.operands);
						} else {
							if constexpr(std::is_pointer_v<std::tuple_element_t<level, op_parameters>>){
								using tt = clean<std::tuple_element_t<level, op_parameters>>;
								std::get<level>(operands) = new tt(*std::get<level>(src.operands));
							}
						}
						if constexpr(level<std::tuple_size_v<decltype(operands)>-1) copyOperands<level+1>(src);
					}

				template<int level=0, typename T_src>
					void moveOperands(T_src&& src){
						if constexpr(level<std::tuple_size_v<decltype(operands)>){
							std::get<level>(operands) = std::get<level>(src.operands);
							src.shouldDelete[level]=false;
							if constexpr(level<std::tuple_size_v<decltype(operands)>-1) moveOperands<level+1>(src);
						}
					}

				template<int level=0>
					void deleteOperands(){
						//If operand is not a pointer, then can't even try to delete it
						if constexpr(std::tuple_size_v<decltype(operands)> > level){
							if constexpr(std::is_pointer_v<std::tuple_element_t<level, op_parameters>>){
								if (shouldDelete[level]) delete std::get<level>(operands);
							}
							if constexpr(level<std::tuple_size_v<decltype(operands)>-1) deleteOperands<level+1>();
						}
					}

				template<int level=0>
					auto& stabilise(){
						if constexpr(std::is_base_of_v<array_operator,clean<std::tuple_element_t<level,op_parameters>>>){
							std::get<level>(operands)->stabilise();
						}
						//Have to have this here
						if constexpr(std::is_pointer_v<std::tuple_element_t<level, op_parameters>>){
							if (rValue[level]){
								auto *tmp = std::get<level>(operands);
								using n = clean<std::tuple_element_t<level,op_parameters>>;
								shouldDelete[level]=true;
								std::get<level>(operands) = new n(*tmp);
							}
						}
						if constexpr(level<std::tuple_size_v<decltype(operands)>-1) stabilise<level+1>();
						return *this;
					}

				/**
				 * Store the operands to the lazy array
				 * Recursive function that finds the first operand that is an array and copies
				 * its indexer to get sizes and then rebases it to one
				 */
				template <int level, typename c, typename... r>
					void storeOperands_inner(int index_index, c &&p1, r &&...others)
					{
						constexpr bool unmatchedArray = 
							arrayInfo<clean<std::tuple_element_t<level,op_parameters>>>::value;
							/*&& !std::is_same_v<std::decay_t<std::remove_pointer_t<decltype(p1)>>,
							std::decay_t<std::tuple_element_t<level,callable_params>>>;*/
						// If type is an array or array_operator then store the index from it if one is not specified
						if constexpr(unmatchedArray)
						{
							// The first array type that you find build your index ranges etc. from them
							// This checks for the first array having the right index type
							// There is a problem if none of them have the right index type, but this shouldn't happen!
							if constexpr(std::is_same_v<typename clean<decltype(p1)>::index_type,lindexer>){
								if constexpr(!is_true_pointer_v<typename std::decay_t<decltype(p1)>>){
									if (index_index==level || index_index==giToDo)
									{
										buildIndex(p1.getIndexer());
										index_index=level;
									}
								} else {
									if (index_index==level || index_index==giToDo)
									{
										buildIndex(p1->getIndexer());
										index_index=level;
									}
								}
							}
						}
						//This checks whether a type that has been stored is actually a pointer or not
						//And how to store the bound data to it if it is a pointer or not
						if constexpr (is_true_pointer_v<std::tuple_element_t<level,op_parameters>>) {
							if constexpr(is_true_pointer_v<decltype(p1)>){
								std::get<level>(operands) = p1; //Storage is pointer and bound data is pointer
							} else {
								if constexpr(std::is_rvalue_reference_v<decltype(p1)>) rValue[level]=true;
								std::get<level>(operands) = &p1; //Storage is pointer but bound data is not, so take address
							}
						} else {
							if constexpr(is_true_pointer_v<decltype(p1)>){
								std::get<level>(operands) = *p1; //Storage is value but bound data is pointer so dereference
							} else {
								std::get<level>(operands) = p1;  //Storage is value and bound data is value so assign
							}
						}
						if constexpr (sizeof...(r) > 0)
						{
							storeOperands_inner<level + 1>(index_index, std::forward<r>(others)...);
						}
			}

				template <typename ...T_params>
					void storeOperands(int index_index, T_params &&...params){
					if (sizeof...(params) == 0) return;
					storeOperands_inner<0>(index_index, std::forward<T_params>(params)...);
					}


			/**
			 * Create the tuple used to call the bound callable
			 */
			template <int level=0,int param=0>
				FAR_INLINE auto buildCallTuple(typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, core_rank>::type &indices) const
				{
					if constexpr(level == std::tuple_size<callable_params>()){// || param == std::tuple_size<op_parameters>()) {
						//Terminate the regression with an empty tuple
						return std::tuple<>();
					} else {
						//Get index info - callable is going to do its own indexing
						using riiType = indexInfo<core_rank>;
						if constexpr(std::is_same_v<std::decay_t<std::tuple_element_t<level,callable_params>>,
								riiType>){
							return std::tuple_cat(std::make_tuple(riiType(indices)),buildCallTuple<level + 1,param>(indices));
						}else{
							// Parameter type is an array or array operator, so subscript it and return a reference to the value
							// Unless the matching parameter to the function is also an array in which case pass it whole
							if constexpr(arrayInfo<clean<decltype(std::get<param>(operands))>>::value // is it an array
									&& !std::is_same_v<std::decay_t<std::remove_pointer_t<std::decay_t<decltype(std::get<param>(operands))>>>,
									std::decay_t<std::tuple_element_t<level,callable_params>>> //Is it not the same type 
									)
							{
								//This checks whether a parameter is a pointer or not
								if constexpr(is_true_pointer_v<decltype(std::get<param>(operands))>){
									auto &value = *std::get<param>(operands);
									// Interposer lambda
									auto calltuple = [&value](auto &&...args) -> decltype(auto)
									{ 
										return value.defaultBase(args...);
									};
									decltype(auto) a = std::apply(calltuple, indices);
									if constexpr(std::is_reference_v<decltype(a)>){
										auto b = std::tuple_cat(std::make_tuple(std::ref(a)),
												buildCallTuple<level + 1,param+1>(indices));
										return b;
									} else {
										auto b = std::tuple_cat(std::make_tuple(a),
												buildCallTuple<level + 1,param+1>(indices));
										return b;
									}
								} else {
									auto &value = std::get<param>(operands);
									// Interposer lambda
									auto calltuple = [&value](auto &&...args) -> decltype(auto)
									{
										return value.defaultBase(args...);
									};
									decltype(auto) a = std::apply(calltuple, indices);
									if constexpr(std::is_reference_v<decltype(a)>){
										auto b = std::tuple_cat(std::make_tuple(std::ref(a)),
												buildCallTuple<level + 1,param+1>(indices));
										return b;
									} else {
										auto b = std::tuple_cat(std::make_tuple(a),
												buildCallTuple<level + 1,param+1>(indices));
										return b;
									}
								}

							}
							else
							{//Parameter is a real parameter (not an array in any way)
								if constexpr(!is_true_pointer_v<std::decay_t<std::tuple_element_t<param,callable_params>>>){
									if constexpr(!is_true_pointer_v<decltype(std::get<param>(operands))>){
										return std::tuple_cat(std::make_tuple(std::ref(std::get<param>(operands))),buildCallTuple<level + 1,param+1>(indices));
									} else {
										return std::tuple_cat(std::make_tuple(std::ref(*std::get<param>(operands))),buildCallTuple<level + 1,param+1>(indices));
									}

								} else {
									return std::tuple_cat(std::make_tuple(std::get<param>(operands)),buildCallTuple<level + 1,param+1>(indices));
								}
							}
						}
					}

					}

					/**
					 * Build overlap between this and a specified other array by checking operands
					 */
					template <int level, typename T_other>
						bool buildOverlap(T_other &&other, bool index_shifted) const
						{
							if constexpr (is_true_pointer_v<std::decay_t<decltype(std::get<level>(operands))>> &&
									arrayInfo<std::remove_pointer_t<std::tuple_element_t<level,op_parameters>>>::value)
							{
								if ((std::get<level>(operands))->getOverlap(other,index_shifted || this->index_shifted))
									return true;
							}
							if constexpr (level < std::tuple_size<op_parameters>() - 1)
							{
								return buildOverlap<level + 1>(other, index_shifted);
							}
							return false;
						}

					public:
					/**
					 * Fulfil the contract to have these types
					 */
					using core_type = std::decay_t<return_type>;
					static const bounds_check_state core_bounds_check = bc_never;//lindexer::core_bounds_check;
					static const bool coreForceContiguous = false;

					void setIndexShifted(){this->index_shifted=true;}

					lazyArray()=default;

					lazyArray(const lazyArray &other){
						index = other.index;
						copyOperands(other);
						for (FAR_UNSIGNED_INDEX_TYPE i=0;i<std::tuple_size_v<op_parameters>;++i){
							rValue[i] = other.rValue[i];
							shouldDelete[i] = other.shouldDelete[i];
						}
					};

					lazyArray(lazyArray &&other){
						index = other.index;
						moveOperands(other);
					};

					~lazyArray(){
						deleteOperands();
					}

					/**
					 * Construct a lazy array from a parameter pack of operands
					 */
					//template <typename... params>
					template <typename... params, std::enable_if_t<std::conjunction_v<std::negation<std::is_same<lazyArray, std::decay_t<params>>>...>, int> = 0>
						FAR_INLINE lazyArray(params &&...to_store)
						{
							if constexpr(sizeof...(to_store)>0)
								storeOperands(giToDo, std::forward<params>(to_store)...);
						}

					/**
					 * Construct a lazy array from a parameter pack of operands
					 */
					template <typename... params, std::enable_if_t<std::conjunction_v<std::negation<std::is_same<lazyArray, std::decay_t<params>>>...>, int> = 0>
						FAR_INLINE lazyArray(lindexer &index, params &&...to_store)
						{
							buildIndex(index);
							if constexpr(sizeof...(to_store)>0)
								storeOperands(giExternal, std::forward<params>(to_store)...);
						}

					using value_type = return_type;
					using iterator = array_iterator<lazyArray>;
					using const_iterator = array_iterator<lazyArray>;
					constexpr static int isVoid = std::is_void_v<return_type>;
					//using reference = std::conditional<!std::is_void_v<value_type>,std::decay_t<value_type>&,void>;
					//using const_reference = const std::decay_t<value_type>&;
					using difference_type = std::ptrdiff_t;
					using pointer = std::decay_t<value_type>*;
					using const_pointer = const std::decay_t<value_type>*;

					auto begin() const {return array_iterator<typename std::decay_t<decltype(*this)>>(*this);}
					auto end() const {array_iterator<typename std::decay_t<decltype(*this)>> i(*this);i. template init_ub<0>(); i++; return i;}

					auto rbegin() const {array_iterator<typename std::decay_t<decltype(*this)>,-1> i(*this); i. template init_ub<0>(); return i;}
					auto rend() const {array_iterator<typename std::decay_t<decltype(*this)>,-1> i(*this);i. template init_lb<0>(); i++; return i;}

					/**
					 * Function to get an operand by index
					 */
					template<int operand>
						auto& getOperand(){
							if constexpr(!is_true_pointer_v<typename std::tuple_element<operand,op_parameters>::type>) {
								return std::get<operand>(operands);
							} else {
								return *std::get<operand>(operands);
							}
						}

					/** 
					 * Operator to get a calculated value
					 */
					template <typename... T_params, std::enable_if_t<is_correct_subscript<core_rank,T_params...>, int>* = nullptr,
						std::enable_if_t<!containsRange<T_params...>, int> = 0>
						FAR_PURE return_type operator()(T_params &&...params) const
						{
							auto params_i64 = typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,core_rank>::type(static_cast<FAR_SIGNED_INDEX_TYPE>(params)...);
							auto a = buildCallTuple<0>(params_i64);
							return_type val = std::apply(callable,a);
							return val;
						}

					template<int level=0, typename T_target=return_type, int rank_target=core_rank, bounds_check_state bc_target=bc_default, bool forceContigTarget = false>
						static auto coreBuildVariant(){
							if constexpr(arrayInfo<std::remove_pointer_t<std::tuple_element_t<level,op_parameters>>>::value){
								return std::remove_pointer_t<std::tuple_element_t<level,op_parameters>>:: template buildVariant<T_target,rank_target,bc_target,forceContigTarget>();
							} else if constexpr(level<std::tuple_size<op_parameters>()-1){
								return coreBuildVariant<level+1,T_target,rank_target,bc_target, forceContigTarget>();
							}
						}


					template<typename T_target=return_type, int rank_target=core_rank, bounds_check_state bc_target=bc_default, bool forceContigTarget = false>
						static auto buildVariant(){
							return coreBuildVariant<0,T_target,rank_target,bc_target, forceContigTarget>();
						}


					template<int level=0, int rank_target=core_rank, bounds_check_state bc_target=bc_default, int ordering_target=0>
						static auto coreBuildIndexVariant(){
							if constexpr(arrayInfo<std::remove_pointer_t<std::tuple_element_t<level,op_parameters>>>::value){
								return std::remove_pointer_t<std::tuple_element_t<level,op_parameters>>:: template buildIndexVariant<rank_target,bc_target,ordering_target>();
							} else if constexpr(level<std::tuple_size<op_parameters>()-1){
								return coreBuildVariant<level+1, rank_target, bc_target, ordering_target>();
							}
						}


					template<int rank_target=core_rank, bounds_check_state bc_target=bc_default, int ordering_target=0>
						static auto buildIndexVariant(){
							return coreBuildIndexVariant<0, rank_target, bc_target, ordering_target>();
						}

					/**
					 * Assign to the references behind the lazy array if the lazy array returns a reference
					 */
					template<typename To>
					void operator=(To&& other){
						static_assert(std::is_reference_v<return_type>,"Assignment to lazy array must be to a reference");
						using other_core = arrayInfo<To>::type;
						forEach<[](return_type a, const other_core &other) {a=other;}>(*this, std::forward<To>(other));
					}

					/**
					 * Operator to get a calculated value using an array as the indices
					 */
					template <typename T_params, std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_params>>, int> = 0>
						FAR_PURE return_type operator()(T_params && param_array) const
						{
							static_assert(std::is_integral_v<typename arrayInfo<T_params>::type>,"Index array must be of integers");
							static_assert(arrayInfo<T_params>::rank==1,"Index array must be rank 1");
							typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, core_rank>::type indices;
							assignTupleArray(indices,param_array);
							return std::apply(*this,indices);
						}


					/** 
					 * Operator to get a calculated value with a guaranteed default lower bound
					 */
					template <typename... T_params, std::enable_if_t<is_correct_subscript<core_rank,T_params...>, int> = 0>
						FAR_PURE return_type defaultBase(T_params &&...params) const
						{
							auto params_i64 = typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,core_rank>::type(static_cast<FAR_SIGNED_INDEX_TYPE>(params)...);
							auto a = buildCallTuple<0>(params_i64);
							return_type val = std::apply(callable,a);
							return val;
						}

					/** 
					 * This function is a duplicate of defaultBase because it has to be here to be equivalent to 
					 * the same function in the array class
					 */
					template <typename... T_params, std::enable_if_t<is_correct_subscript<core_rank,T_params...>, int> = 0>
						FAR_PURE return_type defaultBaseNoUnpack(T_params &&...params) const
						{
							auto params_i64 = typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,core_rank>::type(static_cast<FAR_SIGNED_INDEX_TYPE>(params)...);
							auto a = buildCallTuple<0>(params_i64);
							return_type val = std::apply(callable,a);
							return val;
						}


					/**
					 * Round brackets operator with ranges. Returns a lazy array slice
					 */
					template<typename... T_indices, std::enable_if_t<is_correct_subscript<core_rank, T_indices...>, int> = 0, std::enable_if_t<containsRange<T_indices...>, int> = 0>
						auto operator()(const T_indices... indices){
							constexpr int rrank = countRanges<T_indices...>();
							if constexpr(rrank == core_rank){
								auto l = [](indexInfo<core_rank> info,
										const std::decay_t<decltype(*this)> &op, offset_info<rrank,core_rank> oi) -> auto
								{
									shiftOffsetTuple(oi, info.indices);
									return std::apply(op,info.indices);
								};

								offset_info<rrank,core_rank> o;
								auto ret = lazyArray<l, lindexer, decltype(*this), offset_info<rrank,core_rank>>(*this,o);
								auto& val = ret. template getOperand<1>();
								ret.getIndexer().buildLazySliceOffsets(val, this->getIndexer(), indices...);
								return ret;
							} else {
								using dest_type = std::invoke_result_t<decltype(std::decay_t<decltype(*this)>:: template buildVariant<int, rrank, lindexer::core_bounds,false>)>;
								using dest_itype = std::decay_t<typename dest_type::index_type>;
								auto l = [](indexInfo<rrank> info,
										const std::decay_t<decltype(*this)> op, offset_info<rrank,core_rank> oi) -> auto
								{
									typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, core_rank>::type call_indices;
									shiftPartialOffsetTuple(oi, info.indices, call_indices);
									return std::apply(op,call_indices);
								};
								offset_info<rrank, core_rank> o;
								auto ret = lazyArray<l, dest_itype, decltype(*this), offset_info<rrank,core_rank>>(*this,o);
								auto& val = ret.template getOperand<1>();
								ret.getIndexer().buildLazySliceOffsets(val, this->getIndexer(), indices...);
								return ret;
							}
						}

					/**
					 * Operator to get a slice using an arrays as the upper and lower bounds
					 */
					template <typename T_LB, typename T_UB, std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_LB>>, int> = 0,
									 std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_UB>>, int> = 0>
										 FAR_PURE auto operator()(T_LB && LB, T_UB && UB)
										 {
											 static_assert(std::is_integral_v<typename arrayInfo<T_LB>::type>,"Lower bound array must be of integers");
											 static_assert(std::is_integral_v<typename arrayInfo<T_UB>::type>,"Upper bound array must be of integers");
											 static_assert(arrayInfo<T_LB>::rank==1,"Lower bound array must be rank 1");
											 static_assert(arrayInfo<T_UB>::rank==1,"Upper bound array must be rank 1");
											 typename N_ary_tuple_type<Range, core_rank>::type indices;
											 for (int i = defaultLB; i<defaultLB+LB.getSize(); ++i){
												 assignTupleLevel<0>(indices,i-defaultLB,Range(LB.defaultBase(i), UB.defaultBase(i)));
											 }
											 return std::apply(*this,indices);
										 }

					/**
					 * Operator to get a slice using an arrays as the upper and lower bounds and stride
					 */
					template <typename T_LB, typename T_UB, typename T_Stride, std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_LB>>, int> = 0,
									 std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_UB>>, int> = 0,
									 std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_Stride>>, int> = 0
										 >
										 FAR_PURE auto operator()(T_LB && LB, T_UB && UB, T_Stride &&Stride)
										 { 
											 static_assert(std::is_integral_v<typename arrayInfo<T_LB>::type>,"Lower bound array must be of integers");
											 static_assert(std::is_integral_v<typename arrayInfo<T_UB>::type>,"Upper bound array must be of integers");
											 static_assert(std::is_integral_v<typename arrayInfo<T_Stride>::type>,"Stride array must be of integers");
											 static_assert(arrayInfo<T_LB>::rank==1,"Lower bound array must be rank 1");
											 static_assert(arrayInfo<T_UB>::rank==1,"Upper bound array must be rank 1");
											 static_assert(arrayInfo<T_Stride>::rank==1,"Stride array must be rank 1");
											 typename N_ary_tuple_type<Range, core_rank>::type indices;
											 for (int i = defaultLB; i<defaultLB+LB.getSize(); ++i){
												 assignTupleLevel<0>(indices,i-defaultLB,Range(LB.defaultBase(i), UB.defaultBase(i), Stride.defaultBase(i)));
											 }
											 return std::apply(*this,indices);
										 }

							template<auto ptr>
							auto selectComponent(){
								auto l = [](std::decay_t<core_type> d) {
									return d.*ptr;
								};
								return makeLazy<l>(*this);
							}



					/**
					 * Square bracket operator. Very slow compared to other
					 * accessors
					 */
					return_type getItem(FAR_UNSIGNED_INDEX_TYPE index)
					{
						return std::apply(*this, this->index.reverseIndex(index));
					}

					/**
					 * Const square bracket operator.
					 */
					const return_type getItem(FAR_UNSIGNED_INDEX_TYPE index) const
					{
						return std::apply(*this, this->index.reverseIndex(index));
					}

					/**
					 * Public function to check if there is an overlap with a specified target
					 */
					template <typename T_target>
						bool getOverlap(T_target &&target, bool index_shifted) const
						{
							return buildOverlap<0>(target, index_shifted);
						}
					/**
					 * Get lower bounds
					 */
					FAR_SIGNED_INDEX_TYPE getLB(int irank) const { return index.getLB(irank); }
					FAR_SIGNED_INDEX_TYPE LB(int irank) const { return index.getLB(irank); }
					/**
					 * Get upper bounds
					 */
					FAR_SIGNED_INDEX_TYPE getUB(int irank) const { return index.getUB(irank); }
					FAR_SIGNED_INDEX_TYPE UB(int irank) const { return index.getUB(irank); }
					/**
					 * Get the total size of the array
					 */
					FAR_UNSIGNED_INDEX_TYPE getSize() const { return index.getSize(); }
					FAR_UNSIGNED_INDEX_TYPE size() const { return index.getSize(); }
					/**
					 * Get the size of the array in a given direction
					 */
					FAR_UNSIGNED_INDEX_TYPE getRankSize(int irank) const { return index.getRankSize(irank); }
					FAR_UNSIGNED_INDEX_TYPE rankSize(int irank) const { return index.getRankSize(irank); }
					/**
					 * Get contiguous state
					 * lazyArrays are never contiguous because they do not exist in memory at all
					 */
					bool getContiguous() const {return false;}
					/**
					 * Get allocated state
					 * Always true for lazyArrays since they always produce data
					 */
					bool getAllocated() const {return true;}
					/**
					 * Get associated state
					 * Always false for lazyArrays since they behave like a pointer in many sense
					 * i.e. you cannot move from a lazyArray
					 */
					bool getAssociated() const {return true;}
					/**
					 * Get access to the indexer object
					 */
					auto &getIndexer() { return this->index; }
					const auto &getIndexer() const { return this->index; }

					core_type* data() const {
						#ifndef FAR_BAD_DATA_IS_NULL
						throw std::runtime_error("Cannot get data from a lazy array");
						#else
						return nullptr;
						#endif
					}
				};

			template <auto callable, typename lindexer=void, typename...types>
				static auto makeLazy(types&&... params){
					using paramInfo =  arraySetInfo<callable,0,types...>;
					static_assert(paramInfo::isGood,"Parameters to function are inconsistent with function signature");
					using indexer = std::conditional_t<std::is_void_v<lindexer>,typename paramInfo::info::index,lindexer>;
					if constexpr(!std::is_void_v<indexer> && paramInfo::anyUnmatched){
						return lazyArray<callable,indexer,types...>(std::forward<types>(params)...);
					} else {
						return callable(std::forward<types>(params)...);
					}
				}

			template <typename T_callable, typename lindexer=void, typename...types>
				static auto makeLazyRuntime(T_callable &callable, types&&... params){
					using indexer = std::conditional_t<std::is_void_v<lindexer>,typename first_arrayInfo<types...>::info::index,lindexer>;
					if constexpr(!std::is_void_v<indexer>){
						auto wrapper = [](T_callable &callable, types&&... params){
							return callable(std::forward<types>(params)...);
						};
						return lazyArray<wrapper,indexer,decltype(callable),types...>(callable, std::forward<types>(params)...);
					} else {
						return callable(std::forward<types>(params)...);
					}
				}

			template<auto callable, ert returnType=ert::Auto, typename... types>
				auto makeElemental(types&&... params){
					using indexer = typename first_arrayInfo<types...>::info::index;
					constexpr bool useLazy = (returnType==ert::Lazy) || (returnType==ert::Auto && !paramsCanChange<callable>::value);
					if constexpr(!std::is_void_v<indexer>){
						if constexpr(!useLazy){
							//Parameters can change so have to use elemental function
							static auto func = elementalFunction<callable>();
							return func(std::forward<types>(params)...);
						} else {
							return makeLazy<callable>(std::forward<types>(params)...);
						}
					} else {
						std::tuple<> dummy;
						indexInfo<0> rri(dummy);
						return std::apply(callable,buildCallTuple<callable,0>(rri, std::forward<types>(params)...));
					}
				}

			/**
			 * Class representing a Fortran like array
			 * Templated on the type of the array, the rank of the array, whether bounds checking should
			 * be on or off and the type of the indexer - Fortran type by default
			 */
			template <typename T, int rank, bounds_check_state bounds_check, template <int, bounds_check_state> class indexer, bool forceContiguous>
				class Array : public array_actual
			{

				public:
					template <typename T2, int rank2, bounds_check_state bounds_check2, template <int, bounds_check_state> class indexer2, bool forceContiguous2>
						friend class Array;

					template <typename T2, FAR_SIGNED_INDEX_TYPE rank2>
						friend class array_iterator;

					/**
					 * Store the type of the stored data, the type of the indexer,
					 * the rank and the bounds check state for later retrieval.
					 * Required by contract by other bits of the code
					 */
					using valueType = std::conditional_t<std::is_reference_v<T>,LBrefWrapper<std::remove_reference_t<T>>,T>;
					using core_type = valueType;
					using index_type = indexer<rank, bounds_check>;
					static constexpr FAR_UNSIGNED_INDEX_TYPE core_rank = rank;
					static constexpr bounds_check_state core_bounds_check = bounds_check;
					static constexpr bool coreForceContiguous = forceContiguous;

#ifdef _OPENMP
#ifdef LOCK
					omp_lock_t lock;
#elif defined(FAR_NEST_LOCK)
					omp_nest_lock_t lock;
#endif
#endif

					bool isPointer = false;
					bool copyPointer = false;
					bool shared = false;

				public:

					void setShared(bool value){shared=value;}

					/**
					 * Raw data item. May be shifted by slicing
					 */
					valueType *rdata;
					/**
					 * Base of the raw data. Will NOT be shifted by slicing
					 */
					valueType *bdata;

					/**
					 * Is reallocation an error?
					 */

					bool locked = false;

					indexer<rank, bounds_check> index; // Indexer item

					int procsAtConstruct;

					/**
					 * Logical for whether the array owns the memory
					 * in rdata and bdata. If true will free the memory
					 * on destruction
					 */
					bool shouldDelete = false;

					bool pointerCopy = false;

					/**
					 * Check if the memory of this array overlaps with the memory of another array
					 */
					template <typename T_target>
						bool getOverlap(T_target &&target, bool index_shifted) const { return checkMemoryOverlap(*this, target, index_shifted); }

					template<typename T_target=valueType, int rank_target=rank, bounds_check_state bc_target=bounds_check, bool forceContigTarget = false>
						static auto buildVariant(){
							return Array<T_target, rank_target, bc_target, indexer, forceContigTarget>();
						}

					template<int rank_target=rank, bounds_check_state bc_target=bounds_check>
						static auto buildIndexVariant(){
							return indexer<rank_target,bc_target>();
						}
					/**
					 * If this array owns its memory free it
					 */
					void free_core_inner()
					{
#ifdef LOCK
						omp_set_lock(&lock);
#elif defined(FAR_NEST_LOCK)
						omp_set_nest_lock(&lock);
#elif defined(_OPENMP)
#pragma omp critical(free_core)
#endif
						{
							if (shouldDelete && rdata)
							{
#ifndef SIMPLE_MEMORY
								delete[] rdata;
#else
								std::free(rdata);
#endif
							}
							rdata = nullptr;
							bdata = nullptr;
							shouldDelete = false;
						}
#ifdef LOCK
						omp_unset_lock(&lock);
#elif defined(FAR_NEST_LOCK)
						omp_unset_nest_lock(&lock);
#endif
					}

					/**
					 * Dellocate memory for this array
					 * deal correctly with an OpenMP shared array
					 */
					void free_core()
					{
						free_core_inner();
						return;
					}

					/**
					 * Allocate memory for this array using the sizes from its indexer
					 * for a thread local array
					 */
					void allocate_core_inner()
					{
						free_core();
#ifdef LOCK
						omp_set_lock(&lock);
#elif defined(FAR_NEST_LOCK)
						omp_set_nest_lock(&lock);
#elif defined(_OPENMP)
#pragma omp critical(allocate_core)
						{
#endif
							if (index.getSize()==0) {
								rdata=nullptr;
								bdata=nullptr;
								shouldDelete=false;
							} else {
#ifndef SIMPLE_MEMORY
								rdata = new valueType[index.getSize()];
#else
								rdata = (valueType*)std::malloc(sizeof(valueType)*index.getSize());
								for (FAR_UNSIGNED_INDEX_TYPE i=0;i<index.getSize();++i){
									new(rdata+i) valueType();
								}
#endif
								bdata = rdata;
								shouldDelete = true;
							}
#ifdef LOCK
							omp_unset_lock(&lock);
#elif defined(FAR_NEST_LOCK)
							omp_unset_nest_lock(&lock);
#endif
						}

						/**
						 * Allocate memory for this array using the sizes from its indexer
						 * for a shared OpenMP array
						 */
						void allocate_core()
						{
							allocate_core_inner();
							return;
						}


						/**
						 * Get the indexer for this array
						 */
						auto &getIndexer() { return this->index; }
						/**
						 * Const version of get indexer
						 */
						auto &getIndexer() const { return this->index; }

						/**
						 * Function to copy another array or lazy array to this array
						 * If there is no overlap of source and destination then just copy into existing memory
						 * If there is an overlap then create a temporary to hold the result
						 * and either swap the temporary for this array's memory if this array is not
						 * a pointer or copy again if it is
						 */
						template <typename T_other>
							void safeCopy(T_other &&other)
							{
								//If rank>1 only needed because different indices doesn't matter for rank 1
								if constexpr((!std::is_same_v<typename arrayInfo<T_other>::index,decltype(index)>) && (rank>1)){
#ifdef _OPENMP
										if (isThreadPrivate()){
											//If source is thread private then have to copy the
											//whole of the source to this array
											forEachCopyTightSerial(*this,std::forward<T_other>(other));
										} else
#endif
										{
											//Otherwise can do a parallel copy
											forEachCopyTight(*this, std::forward<T_other>(other));
										}
									return;
								}
								using other_core = std::remove_reference_t<typename arrayInfo<T_other>::type>;
								//Check for overlap between operands
								//This only happens for overlapping slices of the same array
								//except in some very, very weird edge cases that shouldn't happen
								//There is code to deal with them, but it is poorly tested
								if (other.getOverlap(*this,false))
								{
									Array A;
									A.mold(other);
									forEachCopyTight(A,other);

									if constexpr(std::is_same_v<other_core,valueType>){
										//If assignee is either a pointer or has a different indexer
										//then have to copy back
										if (!std::is_same_v<decltype(index),typename arrayInfo<T_other>::index> || getAssociated() || other.getAssociated())
										{
											forEachCopyTight(*this,A);
										}
										else
										{
											// We are not a pointer and have the same indexer so just switch the array
											free_core();
											A.shouldDelete=false;
											this->rdata=A.rdata;
											this->bdata=A.bdata;
										}
									} else {
										//Destination is of different type to source so have to copy
										//This will only trigger in odd cases because of same type requirement
										forEachCopyTight(*this,A);
									}
								} else {
									// It is safe to do a simple copy (same type and source and destination contiguous)
									if (!locked && std::is_same_v<other_core,core_type> && std::is_trivially_copyable_v<core_type>
									 		&& getContiguous() && other.getContiguous() 
									 		&& !(workshare[getTeamNum()] || isThreadPrivate())){
										if constexpr(std::is_base_of_v<array_actual,std::decay_t<T_other>> 
										&& std::is_trivially_copyable_v<core_type> && std::is_same_v<other_core,core_type>){
											std::memcpy(rdata,other.data(),getSize()*sizeof(valueType)); 
										}
									} else {
#ifdef _OPENMP
										if (isThreadPrivate()){
											//If source is thread private then have to copy the
											//whole of the source to this array
											forEachCopyTightSerial(*this,std::forward<T_other>(other));
										} else
#endif
										{
											//Otherwise can do a parallel copy
											forEachCopyTight(*this, std::forward<T_other>(other));
										}
									}
								}
								FAR_OMP_BARRIER
							}

						/**
						 * Routine to determine if the shape of this array is different to another
						 */
						template <typename T_src>
						bool needIndexUpdate(T_src &&src) const
						{
							if (getAssociated()) return false;
							for (int i=1;i<=rank;++i){
								if (src.getRankSize(i)!=getRankSize(i)){
									return true;
								}
							}
							return false;
						}


						/**
						 * Routine to reallocate memory if needed during assignment
						 */
						template <typename T_src>
							void allocateIfNeeded(T_src &&src)
							{
								// If we are an associated pointer nothing we *can* do
								if (getAssociated()) return;
								// If we are allocated and the right size then nothing *to* do
								bool needToAllocate = !rdata;
								if (rdata && !needToAllocate){
									needToAllocate = needIndexUpdate(src);
								}
								if (!needToAllocate) return;
								if (locked) throw std::runtime_error("Allocating memory for a locked array is an error");
								//Now need to update the indexer whatever the source
								//Not referring to the same memory so copy and normalise the index
								FAR_UNSIGNED_INDEX_TYPE sz = rdata?getSize():0;
								index.normaliseCopy(src.getIndexer());
								//Only need to actually reallocate if the sizes are actually different
								if (src.getSize()!=sz) 
									allocate_core();
							}
						/**
						 * Function to store the number of processors when it is created
						 */
						void OMPInit(){
#ifndef _OPENMP
							procsAtConstruct=1;
#else
#ifdef LOCK
							omp_init_lock(&lock);
#elif defined(FAR_NEST_LOCK)
							omp_init_nest_lock(&lock);
#endif
							procsAtConstruct = omp_get_num_threads();
#endif
						}
						//Function to test if a variable is shared or thread private
						bool isThreadPrivate() const{
#ifndef _OPENMP
							return false;
#else
							return (omp_get_num_threads() == procsAtConstruct && !shared);
#endif
						}

						public:
						/**
						 * Constructor from another array or lazy array
						 */
						template <typename T_src, typename = std::enable_if_t<std::is_base_of<array_base, clean<T_src>>::value>>
							Array(T_src const &src) : rdata(nullptr), bdata(nullptr), shouldDelete(false)
						{
							OMPInit();
							*this = src;
						}

						using value_type = core_type;
						using iterator = array_iterator<Array>;
						using const_iterator = array_iterator<Array>;
						using reference = std::decay_t<value_type>&;
						using const_reference = const std::decay_t<value_type>&;
						using difference_type = std::ptrdiff_t;
						using pointer = std::decay_t<value_type>*;
						using const_pointer = const std::decay_t<value_type>*;

						auto begin() const {return array_iterator<typename std::decay_t<decltype(*this)>>(*this);}
						auto end() const {array_iterator<typename std::decay_t<decltype(*this)>> i(*this,true); i++; return i;}

						auto rbegin() const {array_iterator<typename std::decay_t<decltype(*this)>,-1> i(*this,true); return i;}
						auto rend() const {array_iterator<typename std::decay_t<decltype(*this)>,-1>i(*this); i++; return i;}

						/**
						 * Constructor from another array (copy constructor)
						 */
						Array(const Array &src) : rdata(nullptr), bdata(nullptr), shouldDelete(false)
						{
							OMPInit();
							if (src.copyPointer) {this->pointTo(src); return;}
							this->index = src.getIndexer();
							this->index.normalise();
							if (src.rdata) {
								allocate_core();
								safeCopy(std::forward<decltype(src)>(src));
							}
							this->isPointer = src.isPointer;
							this->shared = src.shared;
						}

						/**
						 * Construct from std::array. Only available if rank==1
						 */
						template<typename T_other, FAR_UNSIGNED_INDEX_TYPE sz_other>
							Array(const std::array<T_other,sz_other> &other){
								OMPInit();
#ifdef LOCK
								omp_set_lock(&lock);
#elif defined(FAR_NEST_LOCK)
								omp_set_nest_lock(&lock);
#elif defined(_OPENMP)
#pragma omp critical(fromarray)
#endif
								{
									static_assert(rank==1,"Can only initialise from std::array for rank 1");
									index = index_type(sz_other);
									allocate_core();
									for (FAR_UNSIGNED_INDEX_TYPE i=0;i<sz_other;++i){
										defaultBase(i+defaultLB) = other[i];
									}
								}
#ifdef LOCK
								omp_unset_lock(&lock);
#elif defined(FAR_NEST_LOCK)
								omp_unset_nest_lock(&lock);
#endif
							}

						/**
						 * Move constructor
						 */
						Array(Array &&other) noexcept : rdata(other.rdata), bdata(other.bdata), index(other.index), shouldDelete(other.shouldDelete)
						{
							OMPInit();
#ifdef LOCK
							omp_set_lock(&lock);
#elif defined(FAR_NEST_LOCK)
							omp_set_nest_lock(&lock);
#elif defined(_OPENMP)
#pragma omp critical(move)
#endif
							{
								other.rdata = nullptr;
								other.bdata = nullptr;
								other.shouldDelete = false;
								this->isPointer = other.isPointer;
								this->shared = other.shared;
							}
#ifdef LOCK
							omp_unset_lock(&lock);
#elif defined(FAR_NEST_LOCK)
							omp_unset_nest_lock(&lock);
#endif

						}

						/**
						 * Construct from a raw pointer and sizes/Ranges in each direction. Don't use this without a good reason
						 */
						template <typename... Args>
							explicit Array(valueType *rdata, Args... sizes) : rdata(rdata), bdata(rdata), index(sizes...), shouldDelete(false)
						{
							OMPInit();
							static_assert(sizeof...(sizes) == rank, "Number of arguments must match rank");
						}
						/**
						 * Construct from a list of sizes or Ranges. If a simple size is provided then
						 * that rank runs from 1->size, if a range is specified then it runs from lb to ub
						 */
						template <typename... Args,std::enable_if_t<validSubscript_all<Args...>, int> = 0>
							explicit Array(Args... sizes) : rdata(nullptr), bdata(nullptr), index(sizes...), shouldDelete(false)
						{
							static_assert(sizeof...(sizes) == rank, "Number of size parameters must match rank");
							OMPInit();
							allocate_core();
						}

						/**
						 * Construct from an array of sizes rather than a parameter pack
						 */
						explicit Array(FAR_UNSIGNED_INDEX_TYPE sizes[rank]) : rdata(nullptr), bdata(nullptr), index(sizes), shouldDelete(false)
						{
							OMPInit();
							allocate_core();
						}

						/**
						 * Construct from an std::array of sizes rather than a parameter pack
						 */
						explicit Array(std::array<FAR_SIGNED_INDEX_TYPE,rank> &sizes) : rdata(nullptr), bdata(nullptr), index(sizes), shouldDelete(false)
						{
							OMPInit();
							allocate_core();
						}


						/**
						 * If this is a rank 1 array then construct from std::initializer_list<valueType>
						 */

						template <typename T_other, std::enable_if_t<rank==1 && std::is_convertible<T_other,valueType>::value, int> = 0>
							Array(std::initializer_list<T_other> other) : rdata(nullptr), bdata(nullptr), shouldDelete(false){
								OMPInit();
								index = index_type(other.size());
								allocate_core();
								FAR_UNSIGNED_INDEX_TYPE i=0;
								for (auto &val : other){
									rdata[i] = val;
									++i;
								}
							}

						/**
						 * If this is a rank 2 array then construct from std::initializer_list<std::initializer_list<valueType>>
						 */
						template <typename T_other, std::enable_if_t<rank==2 && std::is_convertible<T_other,valueType>::value, int> = 0>
							Array(std::initializer_list<std::initializer_list<T_other>> other) : rdata(nullptr), bdata(nullptr), shouldDelete(false){
								OMPInit();
								index = index_type(other.begin()->size(),other.size());
								allocate_core();
								FAR_UNSIGNED_INDEX_TYPE i=0;
								for (auto &row : other){
									for (auto &val : row){
										rdata[i] = val;
										++i;
									}
								}
							}

						/**
						 * If this is a rank 3 array then construct from std::initializer_list<std::initializer_list<valueType>>
						 */
						template <typename T_other, std::enable_if_t<rank==3 && std::is_convertible<T_other,valueType>::value, int> = 0>
							Array(std::initializer_list<std::initializer_list<std::initializer_list<T_other>>> other) : rdata(nullptr), bdata(nullptr), shouldDelete(false){
								OMPInit();
								index = index_type(other.begin()->begin()->size(), other.begin()->size(), other.size());
								allocate_core();
								FAR_UNSIGNED_INDEX_TYPE i=0;
								for (auto &plane : other){
									for (auto &row : plane){
										for (auto &val : row){
											rdata[i] = val;
											++i;
										}
									}
								}
							}

						/** 
						 * Default constructor
						 */
						Array() : rdata(nullptr), bdata(nullptr), shouldDelete(false) {OMPInit();}

						/**
						 * Destructor
						 */
						~Array()
						{
							free_core_inner();
							rdata = nullptr;
							bdata = nullptr;
#ifdef LOCK
							omp_destroy_lock(&lock);
#elif defined(FAR_NEST_LOCK)
							omp_destroy_nest_lock(&lock);
#endif
						}

					/**
					 * Bind a new array to the wrapper
					 */
					template <typename... Args>
						void bind(valueType *rdata, Args... sizes)
						{
							free_core();
							this->rdata = rdata;
							this->bdata = rdata;
							index.setSize(sizes...);
							shouldDelete = false;
							isPointer = true;
						}

					/**
					 * Bind a new array to the wrapper, taking ownership of the memory
					 */
					template <typename... Args>
						void takeOwnership(valueType *rdata, Args... sizes)
						{
							free_core();
							this->rdata = rdata;
							this->bdata = rdata;
							index.setSize(sizes...);
							shouldDelete = true;
							isPointer = false;
						}

						void lockAllocation(){
							locked = true;
						}

						void unlockAllocation(){
							locked = false;
						}

						/**
						 * Make this array a pointer pointing to another array
						 * All FortArray Arrays can be pointers
						 */
						template<typename Ta>
						void pointTo(const Ta &other)
						{
							using stype = decltype(*this);
							static_assert(std::is_base_of<array_base, Ta>::value, "Can only point to another array");
							static_assert(arrayInfo<Ta>::rank==arrayInfo<stype>::rank, "Can only point to an array of the same rank");
							static_assert(std::is_same_v<typename arrayInfo<Ta>::type, typename arrayInfo<stype>::type>, "Can only point to an array of the same type");
							static_assert(arrayInfo<Ta>::order==arrayInfo<stype>::order, "Can only point to an array of the same order");
							if constexpr(forceContiguous){
								if (!other.getContiguous())
									throw std::runtime_error("Taking contiguous pointer to non contiguous array");
							}
							free_core();
#ifdef LOCK
							omp_set_lock(&lock);
#elif defined(FAR_NEST_LOCK)
							omp_set_nest_lock(&lock);
#elif defined(_OPENMP)
#pragma omp critical(pointto)
#endif
							{
								this->shared = !other.isThreadPrivate();
								this->rdata = other.rdata;
								this->bdata = other.bdata;
								index.pointCopy(other.index);
								shouldDelete = false;
								isPointer = true;
							}
#ifdef LOCK
							omp_unset_lock(&lock);
#elif defined(FAR_NEST_LOCK)
							omp_unset_nest_lock(&lock);
#endif
						}

						/**
						 * Make this array a pointer pointing to another array
						 * All FortArray Arrays can be pointers
						 */
						template<typename T_other>
							void pointSimple(const T_other &other)
							{
								if constexpr(forceContiguous){
									if (!other.getContiguous())
										throw std::runtime_error("Taking contiguous pointer to non contiguous array");
								}
								free_core();
#ifdef LOCK
								omp_set_lock(&lock);
#elif defined(FAR_NEST_LOCK)
								omp_set_nest_lock(&lock);
#elif defined(_OPENMP)
#pragma omp critical(pointto)
#endif
								{
									this->shared = !other.isThreadPrivate();
									this->rdata = other.rdata;
									this->bdata = other.bdata;
									//Don't copy the index, it will be sorted out later
									shouldDelete = false;
									isPointer = true;
								}
#ifdef LOCK
								omp_unset_lock(&lock);
#elif defined(FAR_NEST_LOCK)
								omp_unset_nest_lock(&lock);
#endif
							}

						void pointTo(){
							free_core();
							isPointer = false;
							this->rdata = nullptr;
							this->bdata = nullptr;
							shouldDelete = false;
							FAR_UNSIGNED_INDEX_TYPE sizes[rank];
							for (int i=0;i<rank;++i){
								sizes[i]=0;
							}
							index.setSize(sizes);
						}

						/**
						 * Make this array point to another array remapping bounds as it does so
						 */
						template<typename T_other, typename... T_ranges, std::enable_if_t<is_correct_subscript<rank,T_ranges...>, int> = 0>
							void pointTo(const T_other &other, T_ranges... ranges){
								using stype = decltype(*this);
								static_assert(std::is_base_of<array_base, T_other>::value, "Can only point to another array");
								static_assert(std::is_same_v<typename arrayInfo<T_other>::type, typename arrayInfo<stype>::type>, "Can only point to an array of the same type");
								static_assert(arrayInfo<T_other>::order==arrayInfo<stype>::order, "Can only point to an array of the same order");

								if (!other.getContiguous())
									throw std::runtime_error("Can only remap array when source is contiguous");

								free_core();

#ifdef LOCK
								omp_set_lock(&lock);
#elif defined(FAR_NEST_LOCK)
								omp_set_nest_lock(&lock);
#elif defined(_OPENMP)
#pragma omp critical(pointtoremap)
#endif
								{
									this->shared = !other.isThreadPrivate();
									this->rdata = other.rdata;
									this->bdata = other.bdata;
									index = index_type(ranges...);
									shouldDelete = false;
									isPointer = true;
								}
#ifdef LOCK
								omp_unset_lock(&lock);
#elif defined(FAR_NEST_LOCK)
								omp_unset_nest_lock(&lock);
#endif

							}

						/**
						 * Return another array that is a slice of the current array
						 * Return is a real array, not a lazy array
						 */
						template <typename... Args, std::enable_if_t<is_correct_subscript<rank,Args...>, int> = 0>
							auto slice(Args... slice_info)
							{
								constexpr static int rrank = countRanges<Args...>();
								using type = arrayInfo<decltype(*this)>::type;
								using dest_type = std::invoke_result_t<decltype(std::decay_t<decltype(*this)>:: template buildVariant<type, rrank, bounds_check, false>)>;
								dest_type ret;
#ifdef LOCK
								omp_set_lock(&lock);
#elif defined(FAR_NEST_LOCK)
								omp_set_nest_lock(&lock);
#elif defined(_OPENMP)
#pragma omp critical(pointtoremap)
#endif
								{
									ret.shared = !isThreadPrivate();
									decltype(*this) &demo = *this;
									ret.pointSimple(demo);
									constexpr int order = arrayInfo<decltype(*this)>::order;
									std::apply([&ret,this](auto... args){ret.index.build_slice(index, args...);},packToTuple<order>(slice_info...));
									ret.rdata += ret.index.slice_offset;
									// DO NOT SHIFT BDATA

#ifdef LOCK
									omp_unset_lock(&lock);
#elif defined(FAR_NEST_LOCK)
									omp_unset_nest_lock(&lock);
#endif
								}
								return ret;
							}

							/**
							 * Function to access a member of a type inside an array
							 */
							template<auto ptr>
							auto selectComponent(){
								auto l = [](core_type &d) -> decltype(auto){
									return d.*ptr;
								};
								return makeLazy<l>(*this);
							}

							template<auto ptr, ert returnType=ert::Auto, typename... T_others>
							auto callMethod(T_others&&... params){

								auto l = [](core_type &d, arrayInfo<T_others>::referenceType... params) -> decltype(auto){
									return (d.*ptr)(std::forward<typename arrayInfo<T_others>::referenceType>(params)...);
								};
								constexpr bool noReturn = std::is_void_v<typename callableTraits<decltype(l)>::type>;
								constexpr bool canChange = paramsCanChange<l>::value;
								constexpr bool returnsReference = std::is_reference_v<typename callableTraits<decltype(l)>::type>;
								constexpr bool autoLazy = ((!canChange || returnsReference) && !noReturn);
								constexpr bool useLazy = (returnType==ert::Lazy) || (returnType==ert::Auto && autoLazy);
								if constexpr(useLazy){
									if constexpr(noReturn){
										makeElemental<l,ert::Lazy>(std::forward<decltype(*this)&>(*this), std::forward<T_others>(params)...);
									} else {
										return makeElemental<l,ert::Lazy>(std::forward<decltype(*this)&>(*this), std::forward<T_others>(params)...);
									}
								} else {
									if constexpr(noReturn){
										makeElemental<l,ert::Array>(std::forward<decltype(*this)&>(*this), std::forward<T_others>(params)...);
									} else {
										return makeElemental<l,ert::Array>(std::forward<decltype(*this)&>(*this), std::forward<T_others>(params)...);
									}
								}
							}


						/**
						 * Allocate the array to a new size
						 * Each element can be either a number or a Range object
						 */
						template <typename... Args>
							void allocate(Args... sizes)
							{
								if (getAllocated()) throw std::runtime_error("Attempting to allocate an already allocated array");
								index.setSize(sizes...);
								allocate_core();
							}

						/**
						 * Allocate from std::array. Only available if rank==1
						 */
						template<typename T_other, FAR_UNSIGNED_INDEX_TYPE sz_other>
							void allocate(const std::array<T_other,sz_other> &other){
								if (getAllocated()) throw std::runtime_error("Attempting to allocate an already allocated array");
								index.setSize(other);
								allocate_core();
							}

						/**
						 * Allocate to a size with an array of sizes
						 */
						void allocate(FAR_UNSIGNED_INDEX_TYPE sizes[rank])
						{
							if (getAllocated()) throw std::runtime_error("Attempting to allocate an already allocated array");
							index.setSize(sizes);
							allocate_core();
						}

						/**
						 * Allocate from arrays of upper and lower bounds
						 */
						void allocate(FAR_SIGNED_INDEX_TYPE LB[rank], FAR_SIGNED_INDEX_TYPE UB[rank])
						{
							if (getAllocated()) throw std::runtime_error("Attempting to allocate an already allocated array");
							index.setBounds(LB, UB);
							allocate_core();
						}

						/**
						 * Deallocate this array
						 */
						void deallocate()
						{
							if (!getAllocated()) throw std::runtime_error("Attempting to deallocate an unallocated array");
							free_core_inner();
							isPointer = false;
						}

						/**
						 * Allocate if not allocated (debug function
						 */
						void renew()
						{
							if (rdata)
								return;
							allocate_core();
						}

						/**
						 * Allocate an array to be the same size as another array, but not have the same contents
						 */
						template<typename T_other>
							void mold(const T_other &other)
							{
#ifdef LOCK
								omp_set_lock(&lock);
#elif defined(FAR_NEST_LOCK)
								omp_set_nest_lock(&lock);
#endif
								{
									free_core();
									index.copyFrom(other.index);
									index.normalise();
									rdata = nullptr;
									bdata = nullptr;
								}
#ifdef LOCK
								omp_unset_lock(&lock);
#elif defined(FAR_NEST_LOCK)
								omp_unset_nest_lock(&lock);
#endif
								allocate_core();
								
							}

						/**
						 * Allocate an array to be the same size as another array, but not have the same contents
						 */
						auto mold()
						{
							Array result;
							result.index = this->index;
							result.index.normalise();
							result.allocate_core();
							return result;
						}

						/**
						 * Get the size of a given rank
						 */
						FAR_UNSIGNED_INDEX_TYPE getRankSize(int r) const { return index.getRankSize(r); }
						/**
						 * Get the size of a given rank
						 */
						FAR_UNSIGNED_INDEX_TYPE rankSize(int r) const { return index.getRankSize(r); }
						/**
						 * Get the lower bound of a given rank
						 */
						FAR_SIGNED_INDEX_TYPE getLB(int r) const { return index.getLB(r); }
						/**
						 * Get the lower bound of a given rank
						 */
						FAR_SIGNED_INDEX_TYPE LB(int r) const { return index.getLB(r); }
						/**
						 * Get the upper bound of a given rank
						 */
						FAR_SIGNED_INDEX_TYPE getUB(int r) const { return index.getUB(r); }
						/**
						 * Get the upper bound of a given rank
						 */
						FAR_SIGNED_INDEX_TYPE UB(int r) const { return index.getUB(r); }
						/**
						 * Get the stride of a given rank
						 */
						FAR_SIGNED_INDEX_TYPE getStride(int r) const { return index.getStride(r); }
						/**
						 * Get the stride of a given rank
						 */
						FAR_SIGNED_INDEX_TYPE stride(int r) const { return index.getStride(r); }
						/**
						 * Get the size of the array
						 */
						FAR_UNSIGNED_INDEX_TYPE getSize() const { return index.getSize(); }
						/**
						 * Get the size of the array
						 */
						FAR_UNSIGNED_INDEX_TYPE size() const { return index.getSize(); }
						/**
						 * Get the contiguous status of the array
						 */
						bool getContiguous() const { return index.getContiguous(); }
						/**
						 * Get the contiguous status of the array
						 */
						bool contiguous() const { return index.getContiguous(); }
						/**
						 * Check if this array has memory associated with it
						 */
						bool getAllocated() const { return this->rdata; }
						/**
						 * Check if this array has memory associated with it
						 */
						bool allocated() const { return this->rdata; }
						/**
						 * Check if this array is a pointer
						 */
						bool getAssociated() const { return isPointer; }
						/**
						 * Check if this array is a pointer
						 */
						bool associated() const { return isPointer; }
						/**
						 * Check if this array is a pointer to a specific other array
						 */
						template <typename T_target>
							bool getAssociated(T_target &&target) { return (this->bdata == target.bdata && this->getSize() == target.getSize()); }
						/**
						 * Check if this array is a pointer to a specific other array
						 */
						template <typename T_target>
							bool associated(T_target &&target) { return (this->bdata == target.bdata && this->getSize() == target.getSize()); }

						void detach()
						{
							if (!this->shouldDelete)
							{
								Array newData;
								newData.mold(*this);
								newData.safeCopy(*this);

								index.normalise();
								newData.shouldDelete=false;
#ifdef LOCK
								omp_set_lock(&lock);
#elif defined(FAR_NEST_LOCK)
								omp_set_nest_lock(&lock);
#endif
								{
									rdata = newData.rdata;
									bdata = newData.bdata;
									shouldDelete = true;
								}
#ifdef LOCK
								omp_unset_lock(&lock);
#elif defined(FAR_NEST_LOCK)
								omp_unset_nest_lock(&lock);
#endif

							}
						}

						/**
						 * Copy assignment
						 */
						Array &operator=(const Array &other)
						{
							// Copy assignment
							if (this != &other)
							{
								if (other.copyPointer) {this->pointTo(other); return *this;}
								if (!this->getAssociated()) allocateIfNeeded(other);
								safeCopy(other);
							}
							return *this;
						}

						/**
						 * Assign from an array operator
						 */
						template <typename T_src, std::enable_if_t<std::is_base_of<array_base, T_src>::value, int> = 0>
							FAR_INLINE Array &operator=(T_src const &other)
							{
								allocateIfNeeded(other);
								safeCopy(other);
								return *this;
							}

						/**
						 * Assign a single value
						 */
						template <typename T_src, std::enable_if_t<!std::is_base_of<array_base, T_src>::value &&
							std::is_assignable_v<core_type&,std::decay_t<T_src>>, int>* = nullptr>
							FAR_INLINE Array &operator=(T_src &&value)
							{
								// Assignment!
								if (getContiguous() && !workshare[getTeamNum()] && std::is_same_v<core_type, std::decay_t<T_src>> && false){
									std::fill(rdata,rdata+getSize(),value);
								} else{
									forEachCopyTight(*this, std::forward<T_src>(value));
								}
								return *this;
							}

						/**
						 * If rank is 1 then assign from std::initializer_list<T_other>, so long as T_other is
						 * convertible to valueType
						 */

						template <typename T_other, std::enable_if_t<rank==1 && std::is_convertible<T_other,valueType>::value, int> = 0>
							Array &operator=(std::initializer_list<T_other> other){
								if (getSize()!=other.size()){
									index = index_type(other.size());
									free_core();
									allocate_core();
								}
								FAR_UNSIGNED_INDEX_TYPE index=0;
								for (auto &val : other){
									this->getItem(index) = val;
									++index;
								}
								return *this;
							}
						/**
						 * If rank is 2 then assign from std::initializer_list<std::initializer_list<T_other>>, so long as T_other is
						 * convertible to valueType
						 */
						template <typename T_other, std::enable_if_t<rank==2 && std::is_convertible<T_other,valueType>::value, int> = 0>
							Array &operator=(std::initializer_list<std::initializer_list<T_other>> other){
								if (getSize()!=other.size()*other.begin()->size()){
									index = index_type(other.begin()->size(),other.size());
									free_core();
									allocate_core();
								}
								FAR_UNSIGNED_INDEX_TYPE index=0;
								for (auto &row : other){
									for (auto &val : row){
										this->getItem(index) = val;
										++index;
									}
								}
								return *this;
							}

						/**
						 * If rank is 3 then assign from std::initializer_list<std::initializer_list<T_other>>, so long as T_other is
						 * convertible to valueType
						 */
						template <typename T_other, std::enable_if_t<rank==3 && std::is_convertible<T_other,valueType>::value, int> = 0>
							Array &operator=(std::initializer_list<std::initializer_list<std::initializer_list<T_other>>> other){
								if (getSize()!=other.size()*other.begin()->size()*other.begin()->begin()->size()){
									index = index_type(other.begin()->begin()->size(), other.begin()->size(), other.size());
									free_core();
									allocate_core();
								}
								FAR_UNSIGNED_INDEX_TYPE index=0;
								for (auto &plane : other){
									for (auto &row : plane){
										for (auto &val : row){
											this->getItem(index) = val;
											++index;
										}
									}
								}
								return *this;
							}

						//Move assignment operator
						Array &operator=(Array &&other) noexcept
						{
							if (other.copyPointer) {this->pointTo(other); return *this;}
							bool allPrivate = isThreadPrivate() && other.isThreadPrivate();
							bool nonePrivate = !isThreadPrivate() && !other.isThreadPrivate();
							bool xorPrivate = !allPrivate && !nonePrivate;
							if (this != &other)
							{
								//Logic here in parallel is
								//1) If source is thread private and destination is shared
								//then this is only correct code in a workshare environment so MUST do copy
								//to get the parts from all threads
								//2) If source is shared and destination is thread private then MUST also copy
								//because the result may well also be being copied by other threads
								//Don't need to explicitly do a single threaded copy though because this case
								//is only valid OUTSIDE a workshare environment
								if (getAssociated() || other.getAssociated() || xorPrivate || locked)
								{
									// If this array is associated with another array
									// or the source is a pointer then you have to copy the data
									// even though this is move assignment operator
									if (!getAssociated()) allocateIfNeeded(other);
									safeCopy(other);
								}
								else
								{
									// If this array is either empty or not associated then
									// Can actually perform a move
									if (allPrivate){
										//Source and destination are both thread local, move is fine
										//Shouldn't need any kind of security
										{
											free_core();
											if (needIndexUpdate(other)){
												index=other.index;
											}
											rdata = other.rdata;
											bdata = other.bdata;
											shouldDelete = other.shouldDelete;
											other.shouldDelete = false;
											other.rdata = nullptr;
											other.bdata = nullptr;
										}
									} else {
										//Remaining case is everything both source and destination are shared arrays
										//Can move so long as only one thread does so
										FAR_OMP_SINGLE
										{
											free_core();
											if (needIndexUpdate(other)){
												index=other.index;
											}
											rdata = other.rdata;
											bdata = other.bdata;
											shouldDelete = other.shouldDelete;
											other.shouldDelete = false;
											other.rdata = nullptr;
											other.bdata = nullptr;
										}
									}
								}
							}
							return *this;
						}


						//Move from operator
						template<typename T_other, std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_other>> &&
							!std::is_same_v<std::decay_t<T_other>,std::decay_t<Array>> //&&
							//			std::is_rvalue_reference_v<T_other>
							, int> = 0>
							Array &operator=(T_other &&other)
							{
								//Must copy if
								//1) RHS is not an rvalue (i.e. must be left intact)
								//2) Indexer is not the same on RHS or LHS
								//3) Type is not the same on RHS and LHS
								//4) Either source or destination are pointers or slices
								if constexpr(!std::is_rvalue_reference_v<decltype(other)> ||
										!std::is_same_v<index_type,typename arrayInfo<T_other>::index> ||
										!std::is_same_v<core_type,typename arrayInfo<T_other>::type>
										) {
									// If this array is associated with another array
									// or the source is a pointer then you have to copy the data
									// even though this is move assignment operator
									allocateIfNeeded(std::forward<T_other>(other));
									safeCopy(std::forward<T_other>(other));
								} else {
									// If this array is either empty or not associated then
									// Can actually perform a move
									// The second branch here shouldn't ever trigger because lazy arrays are associated
									// but here just to prevent a nasty surprise when expanding in future
									if constexpr(std::is_base_of_v<array_actual,std::decay_t<T_other>> &&
										std::is_same_v<core_type,typename arrayInfo<T_other>::type>){
										free_core();
										if (needIndexUpdate(other)){
											index=other.index;
										}
										rdata = other.rdata;
										bdata = other.bdata;
										shouldDelete = other.shouldDelete;
										other.shouldDelete = false;
										other.rdata = nullptr;
										other.bdata = nullptr;
									} else {
										allocateIfNeeded(std::forward<T_other>(other));
										safeCopy(std::forward<T_other>(other));
									}
								}
								return *this;
							}


						/**
						 * Round bracket operator with no ranges. Returns reference to data
						 */
						template <typename... Args, std::enable_if_t<!containsRange<Args...>, int> = 0>
							FAR_PURE decltype(auto) operator()(Args... args)
							{
								if constexpr (isLBrefWrapper<std::remove_reference_t<decltype(rdata)>>::value){
									return rdata[index(args...)].get();
								} else {
									return rdata[index(args...)];
								}
							}


						/**
						 * Operator to get a calculated value using an array as the indices
						 */
						template <typename T_params, std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_params>>, int> = 0>
							FAR_PURE auto& operator()(T_params && param_array)
							{
								static_assert(std::is_integral_v<typename arrayInfo<std::decay_t<T_params>>::type>,"Index array must be of integers");
								static_assert(arrayInfo<std::decay_t<T_params>>::rank==1,"Index array must be rank 1");
								typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, core_rank>::type indices;
								assignTupleArray(indices,param_array);
								return std::apply(*this,indices);
							}

						/**
						 * Operator to get a slice using an arrays as the upper and lower bounds
						 */
						template <typename T_LB, typename T_UB, std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_LB>>, int> = 0,
										 std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_UB>>, int> = 0>
											 FAR_PURE auto operator()(T_LB && LB, T_UB && UB) const
											 {
												 static_assert(std::is_integral_v<typename arrayInfo<T_LB>::type>,"Lower bound array must be of integers");
												 static_assert(std::is_integral_v<typename arrayInfo<T_UB>::type>,"Upper bound array must be of integers");
												 static_assert(arrayInfo<T_LB>::rank==1,"Lower bound array must be rank 1");
												 static_assert(arrayInfo<T_UB>::rank==1,"Upper bound array must be rank 1");
												 typename N_ary_tuple_type<Range, core_rank>::type indices;
												 for (int i = defaultLB; i<defaultLB+LB.getSize(); ++i){
													 assignTupleLevel<0>(indices,i-defaultLB,Range(LB.defaultBase(i), UB.defaultBase(i)));
												 }
												 return std::apply(*this,indices);
											 }

						/**
						 * Operator to get a slice using an arrays as the upper and lower bounds and stride
						 */
						template <typename T_LB, typename T_UB, typename T_Stride, 
										 std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_LB>>, int> = 0,
										 std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_UB>>, int> = 0,
										 std::enable_if_t<std::is_base_of_v<array_base,std::decay_t<T_Stride>>, int> = 0
											 >
											 FAR_PURE auto operator()(T_LB && LB, T_UB && UB, T_Stride && Stride) const
											 {
												 static_assert(std::is_integral_v<typename arrayInfo<T_LB>::type>,"Lower bound array must be of integers");
												 static_assert(std::is_integral_v<typename arrayInfo<T_UB>::type>,"Upper bound array must be of integers");
												 static_assert(std::is_integral_v<typename arrayInfo<T_UB>::type>,"Stride array must be of integers");
												 static_assert(arrayInfo<T_LB>::rank==1,"Lower bound array must be rank 1");
												 static_assert(arrayInfo<T_UB>::rank==1,"Upper bound array must be rank 1");
												 static_assert(arrayInfo<T_UB>::rank==1,"Stride array must be rank 1");
												 typename N_ary_tuple_type<Range, core_rank>::type indices;
												 for (int i = defaultLB; i<defaultLB+LB.getSize(); ++i){
													 assignTupleLevel<0>(indices,i-defaultLB,Range(LB.defaultBase(i), UB.defaultBase(i), Stride.defaultBase(i)));
												 }
												 return std::apply(*this,indices);
											 }

						/**
						 * Round bracket operator with no ranges. Returns reference to data. Const version
						 */
						template <typename... Args, std::enable_if_t<!containsRange<Args...>, int> = 0>
							FAR_PURE decltype(auto) operator()(Args... args) const
							{
								if constexpr (isLBrefWrapper<std::remove_reference_t<decltype(rdata)>>::value){
									return rdata[index(args...)].get();
								} else {
									return rdata[index(args...)];
								}
							}

						/**
						 * Round bracket operator with ranges. Return slice of array
						 */
						template <typename... Args, std::enable_if_t<containsRange<Args...>, int> = 0>
							FAR_PURE auto operator()(Args... args)
							{
								return slice(args...);
							}

						/**
						 * Round bracket operator with a one base. Only takes values
						 */
						template <typename... Args>
							FAR_PURE decltype(auto) defaultBase(Args... args)
							{
								if constexpr (isLBrefWrapper<std::remove_reference_t<decltype(rdata)>>::value){
									return rdata[index.defaultBase(args...)].get();
								} else {
									return rdata[index.defaultBase(args...)];
								}

							}

						/**
						 * Const one base round bracket operator, const version
						 */
						template <typename... Args>
							FAR_PURE decltype(auto) defaultBase(Args... args) const
							{
								if constexpr (isLBrefWrapper<std::remove_reference_t<decltype(rdata)>>::value){
									return rdata[index.defaultBase(args...)].get();
								} else {
									return rdata[index.defaultBase(args...)];
								}
							}

						/**
						 * Round bracket operator with a one base. Only takes values
						 */
						template <typename... Args>
							FAR_PURE decltype(auto) defaultBaseNoUnpack(Args... args)
							{
								return rdata[index.defaultBase(args...)];
							}

						/**
						 * Const one base round bracket operator, const version
						 */
						template <typename... Args>
							FAR_PURE decltype(auto) defaultBaseNoUnpack(Args... args) const
							{
								return rdata[index.defaultBase(args...)];
							}

						/**
						 * At method. Always range checks C++ style
						 */
						template <typename... Args>
							valueType &at(Args... args)
							{
								FAR_UNSIGNED_INDEX_TYPE offset = 0;
								offset = index(args...);
								if (offset < 0 || offset >= index.getSize())
									throw std::out_of_range("Array::at");
								return rdata[offset];
							}


						/**
						 * Square bracket operator. Very slow compared to other
						 * accessors
						 */
						valueType& getItem(FAR_UNSIGNED_INDEX_TYPE index)
						{
							if constexpr(forceContiguous){
								return rdata[index];
							} else {
								if (getContiguous() && !isPointer){
									return rdata[index];
								} else {
									return rdata[this->index[index]];
								}
							}
						}

						/**
						 * Const square bracket operator.
						 */
						valueType& getItem(FAR_UNSIGNED_INDEX_TYPE index) const
						{
							if constexpr(forceContiguous){
								return rdata[index];
							} else {
								if (getContiguous() && !isPointer){
									return rdata[index];
								} else {
									return rdata[this->index[index]];
								}
							}
						}

						/**
						 * Return the raw data, including any slice offset
						 */
						valueType *data() const { return this->rdata; }

						/**
						 * Return the raw data
						 * For a sliced array DOES NOT return the slice
						 * but returns the whole item
						 * Mainly for debugging
						 */
						valueType *baseData() const {return this->bdata; }

					};

					template <int rank, bounds_check_state bounds = bounds_check_state::bc_default>
						using defaultarray = Array<double, rank, bounds>;

					/**
					 * Alias wrapping a double array
					 */
					template <int rank, bounds_check_state bounds = bounds_check_state::bc_default>
						using doubleArray = Array<double, rank, bounds>;

					/**
					 * Alias wrapping a float array
					 */
					template <int rank, bounds_check_state bounds = bounds_check_state::bc_default>
						using floatArray = Array<float, rank, bounds>;

					/**
					 * Alias wrapping an int32 array
					 */
					template <int rank, bounds_check_state bounds = bounds_check_state::bc_default>
						using int32Array = Array<int32_t, rank, bounds>;

					/**
					 * Alias wrapping an int64 array
					 */
					template <int rank, bounds_check_state bounds = bounds_check_state::bc_default>
						using int64Array = Array<FAR_SIGNED_INDEX_TYPE, rank, bounds>;

					template <typename T, int rank, bounds_check_state bounds = bounds_check_state::bc_default>
						using contiguousArray = Array<T,rank,bounds,FortranIndexBuilder,true>;

					template <typename T, int rank, bounds_check_state bounds = bounds_check_state::bc_default>
						using FortranArray = Array<T, rank, bounds, FortranIndexBuilder>;

					template <typename T, int rank, bounds_check_state bounds = bounds_check_state::bc_default>
						using CArray = Array<T, rank, bounds, CIndexBuilder>;


					//Fortran intrinsic functions replicated here
					/**
					 * Fortran ABS(https://fortranwiki.org/fortran/show/abs)
					 */
					template <typename T>
						auto abs(T &&a){
							ERRCHECKREALINTEGERCOMPLEX(T,"abs","a")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::abs(val); };
							return makeLazy<l>(std::forward<T>(a));
						}

					/**
					 * Fortran ACHAR(https://fortranwiki.org/fortran/show/achar)
					 */
					template <typename kind=char, typename T=char>
						auto achar(T &&i)
						{
							ERRCHECKINTEGER(T,"achar","i")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> kind
								{  return (kind)ascii_to_char[val]; };
							return makeLazy<l>(std::forward<T>(i));
						}

					/**
					 * Fortran ACOS(https://fortranwiki.org/fortran/show/acos)
					 */
					template <typename T>
						auto acos(T &&x)
						{
							ERRCHECKREALCOMPLEX(T,"acos","x");
							auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
							{  return std::acos(val); };
							return makeLazy<l>(std::forward<T>(x));
						}
					/** 
					 * Fortran ACOSH(https://fortranwiki.org/fortran/show/acosh)
					 */
					template <typename T>
						auto acosh(T &&x)
						{
							ERRCHECKREALCOMPLEX(T,"acosh","x");
							auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
							{  return std::acosh(val); };
							return makeLazy<l>(std::forward<T>(x));
						}


					/**
					 * Fortran ACOSPI (Arccosine on [-1,1], F2023)
					 */
					template <typename T>
						auto acospi(T &&x)
						{
							ERRCHECKREALCOMPLEX(T,"acospi","x");
							auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
							{  return std::acos(val)/M_PI; };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran ACOSD (Arccosine in degrees, F2023)
					 */
					template <typename T>
						auto acosd(T &&x)
						{
							ERRCHECKREALCOMPLEX(T,"acosd","x");
							auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
							{  return std::acos(val)*180.0/M_PI; };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran AIMAG(https://fortranwiki.org/fortran/show/aimag)
					 */
					template <typename T>
						auto aimag(T &&z)
						{
							AERRCHECKCOMPLEX(T,"aimag","z");
							auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
							{  return std::imag(val); };
							return makeLazy<l>(std::forward<T>(z));
						}

					/** 
					 * Fortran AINT(https://fortranwiki.org/fortran/show/aint)
					 */
					template <typename kind = int, typename T = int>
						auto aint(T &&a)
						{
							ERRCHECKREAL(T,"aint","a")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> kind
								{  return (kind)std::trunc(val); };
							return makeLazy<l>(std::forward<T>(a));
						}

					/**
					 * Fortran ALL(https://fortranwiki.org/fortran/show/all)
					 */
					template <typename T>
						auto all(T &&source)
						{
							AERRCHECKTRUTHY(T,"all","mask")
								return reduction<[](const bool &src, bool &dest){dest=dest&&src;return dest?iaContinue:iaTerminate;},
											 [](const bool &src, bool &dest){dest=dest&&src;},
											 [](bool &dest){dest=true;}
							>(std::forward<T>(source));
						}

					/**
					 * Fortran ALL with dim(https://fortranwiki.org/fortran/show/all)
					 */
					template <typename T>
						auto all(T &&source, int direction)
						{
							AERRCHECKTRUTHY(T,"all","mask")
								static_assert(arrayInfo<std::decay_t<T>>::value,"Source type must be an array in all with direction");
							return reductionWithDirection<[](const bool &src, bool &dest){dest=dest&&src;return dest?iaContinue:iaTerminate;},
										 [](const bool &src, bool &dest){dest=dest&&src;},
										 [](bool &dest){dest=true;}
							>(direction, std::forward<T>(source));	
						}

					/**
					 * Fortran ALL(https://fortranwiki.org/fortran/show/all)
					 * with mask
					 * This actually isn't something in Fortran, but is included here for testing purposes
					 */
					template <typename T, typename T_mask>
						auto all_with_mask(T &&source, T_mask &&mask)
						{
							AERRCHECKTRUTHY(T,"all","mask")
								return reduction<[](const bool &src, bool &dest, const bool& mask){dest=dest&&(mask?src:dest);return dest?iaContinue:iaTerminate;},
											 [](const bool &src, bool &dest){dest=dest&&src;},
											 [](bool &dest){dest=true;}
							>(std::forward<T>(source),std::forward<T_mask>(mask));
						}

					/**
					 * Fortran ALL with dim and mask(https://fortranwiki.org/fortran/show/all)
					 */
					template <typename T, typename T_mask>
						auto all_with_mask(T &&source, int direction, T_mask &&mask)
						{
							AERRCHECKTRUTHY(T,"all","mask")
								static_assert(arrayInfo<std::decay_t<T>>::value,"Source type must be an array in all with direction");
							return reductionWithDirection<[](const bool &src, bool &dest, const bool& mask){dest=dest&&(mask?src:dest);return dest?iaContinue:iaTerminate;},
										 [](const bool &src, bool &dest){dest=dest&&src;},
										 [](bool &dest){dest=true;}
							>(direction, std::forward<T>(source), std::forward<T_mask>(mask));
						}

					//Anonymous namespace for functions to do allocation
					namespace {
						/**
						 * Forward declaration of allocate middle function
						 */
						template<bool deallocate = false, typename T_item = int, typename... T_others>//, typename = std::enable_if_t<std::is_base_of_v<array, std::decay_t<T_item>>>>
							void alloc_middle(T_item & item, T_others&... others);

						/**
						 * Build a tuple of parameters than can be passed to the allocate method
						 * Helper function for allocate
						 */
						template<bool deallocate, typename T_item, typename T_params, typename T_next, typename... T_others,
							typename = std::enable_if_t<std::is_integral_v<std::decay_t<T_next>> || std::is_same_v<Range,std::decay_t<T_next>>>>
								void build_alloc(T_item && item, T_params&& params, T_next && next, T_others&&... others){
									auto next_tuple = std::tuple_cat(params,std::make_tuple(next));
									if constexpr(std::tuple_size<decltype(next_tuple)>() == std::remove_reference_t<decltype(item)>::core_rank){
										std::apply([&item](auto&&... params){item.allocate(params...);}, next_tuple);
										if constexpr(sizeof...(others)>0)
											alloc_middle<deallocate>(others...);
									} else {
										build_alloc<deallocate>(item, next_tuple, others...);
									}
								}
						/**
						 * Middle function to separate the object to be allocated and the allocation parameters
						 * Helper function for ALLOCATE
						 */
						template<bool deallocate, typename T_item, typename... T_others>//, typename = std::enable_if_t<std::is_base_of_v<array, std::decay_t<T_item>>>>
							void alloc_middle(T_item & item, T_others&... others){
								if constexpr(deallocate) {
									if (!item.getAllocated()) throw std::runtime_error("Trying to reallocate an unallocated array");
									item.deallocate();
								}
								build_alloc<deallocate>(item,std::tuple<>(),others...);
							}
					}

					/**
					 * Fortran ALLOCATE
					 * Allocate function. Parameters should be of the form
					 * Allocatable, RankSpec, RankSpec, Rankspec, Allocatable, ....
					 */
					template<typename... T>
						void allocate(T&& ...items){
							alloc_middle(items...);
						}


					/*
					 * Fortran ALLOCATED (https://fortranwiki.org/fortran/show/allocated)
					 */
					template <typename T>
						bool allocated(T &&source)
						{
							return source.getAllocated();
						}

					namespace{
						template<typename kind, typename param>
							kind templatedRound(param val){
								using signed_kind = std::make_signed_t<kind>;
								if constexpr(std::is_same_v<signed_kind,int>){
									return std::round(val);
								} else if constexpr(std::is_same_v<signed_kind,long>){
									return std::lround(val);
								} else if constexpr(std::is_same_v<signed_kind,long long>){
									return std::llround(val);
								} else if constexpr(std::is_same_v<signed_kind,float>){
									return std::roundf(val);
								} else if constexpr(std::is_same_v<signed_kind,double>){
									return std::round(val);
								} else {
									throw std::runtime_error("Invalid type in round");
								}
							}
					}

					/*
					 * Fortran ANINT(https://fortranwiki.org/fortran/show/aint)
					 */
					template <typename kind = int, typename T = int>
						auto anint(T &&a)
						{
							ERRCHECKREAL(T,"anint","a")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> kind
								{  return templatedRound<kind>(val); };
							return makeLazy<l>(std::forward<T>(a));
						}

					/*
					 * Fortran ANY(https://fortranwiki.org/fortran/show/any)
					 */
					template <typename T>
						auto any(T &&source)
						{
							AERRCHECKTRUTHY(T,"any","mask")
								return reduction<[](const bool &src, bool &dest){dest=dest||src;return dest?iaContinue:iaTerminate;},
											 [](const bool &src, bool &dest){dest=dest||src;},
											 [](bool &dest){dest=false;}
							>(std::forward<T>(source));
						}

					/*
					 * Fortran ANY(https://fortranwiki.org/fortran/show/any)
					 * with mask
					 */
					template <typename T, typename T_mask>
						auto any_with_mask(T &&source, T_mask &&mask)
						{
							AERRCHECKTRUTHY(T,"any_with_mask","source")
								AERRCHECKTRUTHY(T_mask,"any_with_mask","mask")
								return reduction<[](const bool &src, bool &dest, const bool &mask)
								{dest=dest||mask?src:false;return dest?iaContinue:iaTerminate;},
								[](const bool &src, bool &dest){dest=dest||src;},
								[](bool &dest){dest=false;}
							>(std::forward<T>(source), std::forward<T_mask>(mask));
						}

					/**
					 * Fortran ANY(https://fortranwiki.org/fortran/show/any)
					 * With dim
					 */
					template <typename T>
						auto any(T &&source, int dim)
						{
							AERRCHECKTRUTHY(T,"any","mask")
								return reductionWithDirection<[](const bool &src, bool &dest){dest=dest||src;return dest?iaTerminate:iaContinue;},
											 [](const bool &src, bool &dest){dest=dest||src;},
											 [](bool &dest){dest=false;}
							>(dim, std::forward<T>(source));
						}

					/**
					 * Fortran ANY(https://fortranwiki.org/fortran/show/any)
					 * With dim and mask
					 */
					template <typename T, typename T_mask>
						auto any_with_mask(T &&source, int dim, T_mask &&mask)
						{
							AERRCHECKTRUTHY(T,"any_with_mask","source")
								AERRCHECKTRUTHY(T_mask,"any_with_mask","mask")
								return reductionWithDirection<[](const bool &src, bool &dest, const bool &mask)
								{dest=dest||mask?src:false;return dest?iaContinue:iaTerminate;},
								[](const bool &src, bool &dest){dest=dest||src;},
								[](bool &dest){dest=false;}
							>(dim, std::forward<T>(source), std::forward<T_mask>(mask));
						}

					/**
					 * Fortran ASIN(https://fortranwiki.org/fortran/show/asin)
					 */
					template <typename T>
						auto asin(T &&x)
						{
							AERRCHECKREALCOMPLEX(T,"asin","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::asin(val); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran ASIND(asin with result in degrees F2023)
					 */
					template <typename T>
						auto asind(T &&x)
						{
							AERRCHECKREALCOMPLEX(T,"asind","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::asin(val) * 180.0/M_PI; };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran ASINH(https://fortranwiki.org/fortran/show/asinh)
					 */
					template <typename T>
						auto asinh(T &&x)
						{ 
							AERRCHECKREALCOMPLEX(T,"asinh","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::asinh(val); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran ASINPI(asin with result in degrees F2023)
					 */
					template <typename T>
						auto asinpi(T &&x)
						{
							AERRCHECKREALCOMPLEX(T,"asinpi","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::asin(val)/M_PI; };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran ASSOCIATED(https://fortranwiki.org/fortran/show/associated)
					 */
					template <typename T>
						bool associated(T &&pointer)
						{
							return pointer.getAssociated();
						}

					/**
					 * Fortran ASSOCIATED with target(https://fortranwiki.org/fortran/show/associated)
					 */
					template <typename T,typename T_target>
						bool associated(T &&pointer, T_target&&target)
						{ 
							return pointer.getAssociated(target);
						}

					/**
					 * Fortran ATAN(https://fortranwiki.org/fortran/show/atan)
					 */
					template <typename T>
						auto atan(T &&x)
						{
							AERRCHECKREALCOMPLEX(T,"atan","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::atan(val); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran ATAN with two parameters(https://fortranwiki.org/fortran/show/atan)
					 */
					template <typename T1, typename T2>
						auto atan(T1 &&y, T2 &&x)
						{
							AERRCHECKREALCOMPLEX(T1,"atan","y")
								AERRCHECKREALCOMPLEX(T2,"atan","x")
								auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &y, typename arrayInfo<std::decay_t<T2>>::type const &x) -> decltype(auto)
								{return std::atan2(y,x); };
							return makeLazy<l>(std::forward<T1>(y),std::forward<T2>(x));
						}

					/**
					 * Fortran ATAN2(https://fortranwiki.org/fortran/show/atan2)
					 */
					template <typename T1, typename T2>
						auto atan2(T1 &&y, T2 &&x)
						{
							AERRCHECKREALCOMPLEX(T1,"atan2","y")
								AERRCHECKREALCOMPLEX(T2,"atan2","x")
								auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &y, typename arrayInfo<std::decay_t<T2>>::type const &x) -> decltype(auto)
								{  return std::atan2(y,x); };
							return makeLazy<l>(std::forward<T1>(y),std::forward<T2>(x));
						}

					/**
					 * Fortran ATAN2D(ATAN2 with parameters in degrees)
					 */
					template <typename T1, typename T2>
						auto atan2d(T1 &&y, T2 &&x)
						{
							AERRCHECKREALCOMPLEX(T1,"atan2d","y")
								AERRCHECKREALCOMPLEX(T2,"atan2d","x")
								auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &y, typename arrayInfo<std::decay_t<T2>>::type const &x) -> decltype(auto)
								{  return std::atan2(y,x) * 180.0/M_PI; };
							return makeLazy<l>(std::forward<T1>(y),std::forward<T2>(x));
						}

					/**
					 * Fortran ATAN2PI(ATAN2 with parameters on [-1,1])
					 */
					template <typename T1, typename T2>
						auto atan2pi(T1 &&y, T2 &&x)
						{
							AERRCHECKREALCOMPLEX(T1,"atan2d","y")
								AERRCHECKREALCOMPLEX(T2,"atan2d","x")
								auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &y, typename arrayInfo<std::decay_t<T2>>::type const &x) -> decltype(auto)
								{  return std::atan2(y,x)/M_PI; };
							return makeLazy<l>(std::forward<T1>(y),std::forward<T2>(x));
						}

					/**
					 * Fortran ATAND(ATAN with result in degrees)
					 */
					template <typename T>
						auto atand(T &&x)
						{
							AERRCHECKREALCOMPLEX(T,"atand","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::atan(val) * 180.0/M_PI; };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran ATAND with two parameters with result in degres
					 */
					template <typename T1, typename T2>
						auto atand(T1 &&y, T2 &&x)
						{
							AERRCHECKREALCOMPLEX(T1,"atand","y")
								AERRCHECKREALCOMPLEX(T2,"atand","x")
								auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &y, typename arrayInfo<std::decay_t<T2>>::type const &x) -> decltype(auto)
								{  return std::atan2(y,x) * 180.0/M_PI; };

							return makeLazy<l>(std::forward<T1>(y),std::forward<T2>(x));

						}

					/**
					 * Fortran ATANH(https://fortranwiki.org/fortran/show/atanh)
					 */
					template <typename T>
						auto atanh(T &&x)
						{
							AERRCHECKREALCOMPLEX(T,"atanh","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::atanh(val); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran ATANPI(ATAN with result on [-1,1])
					 */
					template <typename T>
						auto atanpi(T &&x)
						{
							AERRCHECKREALCOMPLEX(T,"atanpi","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::atan(val)/M_PI; };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran ATANPI with two parameters with result on[-1,1]
					 */
					template <typename T1, typename T2>
						auto atanpi(T1 &&y, T2 &&x)
						{
							AERRCHECKREALCOMPLEX(T1,"atanpi","y")
								AERRCHECKREALCOMPLEX(T2,"atanpi","x")
								auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &y, typename arrayInfo<std::decay_t<T2>>::type const &x) -> decltype(auto)
								{  return std::atan2(y,x)/M_PI; };

							return makeLazy<l>(std::forward<T1>(y),std::forward<T2>(x));

						}

					/**
					 * Fortran BESSEL_J0(https://fortranwiki.org/fortran/show/bessel_j0)
					 */
					template <typename T>
						auto bessel_j0(T &&x)
						{
							AERRCHECKREAL(T,"bessel_j0","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::cyl_bessel_j(0,val); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran BESSEL_J1(https://fortranwiki.org/fortran/show/bessel_j1)
					 */
					template <typename T>
						auto bessel_j1(T &&x)
						{
							AERRCHECKREAL(T,"bessel_j1","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::cyl_bessel_j(1,val); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran BESSEL_JN(https://fortranwiki.org/fortran/show/bessel_jn)
					 */
					template <typename T1, typename T2>
						auto bessel_jn(T1 &&n, T2 &&x)
						{
							AERRCHECKINTEGER(T1,"bessel_jn","n")
								AERRCHECKREAL(T2,"bessel_jn","x")
								auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &n, 
										typename arrayInfo<std::decay_t<T2>>::type const &x) -> decltype(auto)
								{  return std::cyl_bessel_j(n,x); };
							return makeLazy<l>(std::forward<T1>(n),std::forward<T2>(x));
						}

					/**
					 * Fortran BESSEL_JN(https://fortranwiki.org/fortran/show/bessel_jn)
					 */
					template <typename T1, typename T2, typename T3>
						auto bessel_jn(T1 &&n1, T2 &&n2, T3 &&x)
						{
							static_assert(!arrayInfo<T1>::value,"n1 must not be an array in bessel_jn");
							static_assert(!arrayInfo<T2>::value,"n2 must not be an array in bessel_jn");
							static_assert(!arrayInfo<T3>::value,"x must not be an array in bessel_jn");

							AERRCHECKINTEGER(T1,"bessel_jn","n1")
								AERRCHECKINTEGER(T2,"bessel_jn","n2")
								AERRCHECKREAL(T3,"bessel_jn","x")
								Array<std::decay_t<T3>,1> result (n2-n1+1);
							for (int i=1;i<=n2-n1+1;++i){
								result(i) = std::cyl_bessel_j(i-1+n1,x);
							}
							return result;
						}

					/**
					 * Fortran BESSEL_Y0(https://fortranwiki.org/fortran/show/bessel_y0)
					 */
					template <typename T>
						auto bessel_y0(T &&x)
						{
							AERRCHECKREAL(T,"bessel_y0","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::cyl_neumann(0,val); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran BESSEL_Y1(https://fortranwiki.org/fortran/show/bessel_y1)
					 */
					template <typename T>
						auto bessel_y1(T &&x)
						{
							AERRCHECKREAL(T,"bessel_y1","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
								{  return std::cyl_neumann(1,val); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran BESSEL_YN(https://fortranwiki.org/fortran/show/bessel_yn)
					 */
					template <typename T1, typename T2>
						auto bessel_yn(T1 &&n, T2 &&x)
						{
							AERRCHECKINTEGER(T1,"bessel_yn","n")
								AERRCHECKREAL(T2,"bessel_yn","x")
								auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &n, typename arrayInfo<std::decay_t<T2>>::type const &x) -> decltype(auto)
								{  return std::cyl_neumann(n,x); };
							return makeLazy<l>(std::forward<T1>(n),std::forward<T2>(x));
						}

					/**
					 * Fortran BESSEL_YN(https://fortranwiki.org/fortran/show/bessel_jn)
					 */
					template <typename T1, typename T2, typename T3>
						auto bessel_yn(T1 &&n1, T2 &&n2, T3 &&x)
						{
							static_assert(!arrayInfo<T1>::value,"n1 must not be an array in bessel_jn");
							static_assert(!arrayInfo<T2>::value,"n2 must not be an array in bessel_jn");
							static_assert(!arrayInfo<T3>::value,"x must not be an array in bessel_jn");

							AERRCHECKINTEGER(T1,"bessel_yn","n1")
								AERRCHECKINTEGER(T2,"bessel_yn","n2")
								AERRCHECKREAL(T3,"bessel_yn","x")
								Array<std::decay_t<T3>,1> result (n2-n1+1);
							for (int i=1;i<=n2-n1+1;++i){
								result(i) = std::cyl_neumann(i-1+n1,x);
							}
							return result;
						}

					/**
					 * Fortran BGE(https://fortranwiki.org/fortran/show/bge)
					 */
					template <typename T1, typename T2>
						auto bge(T1 &&i, T2 &&j)
						{
							AERRCHECKINTEGER(T1,"bge","i")
								AERRCHECKINTEGER(T2,"bge","j")
								using At1 = typename arrayInfo<T1>::type;
							using At2 = typename arrayInfo<T2>::type;
							using Unsigned1 = typename std::make_unsigned<typename arrayInfo<T1>::type>::type;
							using Unsigned2 = typename std::make_unsigned<typename arrayInfo<T2>::type>::type;
							//Lambda for comparison
							auto l =[](At1 &u1, At2 &u2){return std::bit_cast<Unsigned1>(u1)>=std::bit_cast<Unsigned2>(u2);};
							return makeLazy<l>(std::forward<T1>(i),std::forward<T2>(j));
						}

					/**
					 * Fortran BGT(https://fortranwiki.org/fortran/show/bgt)
					 */
					template <typename T1, typename T2>
						auto bgt(T1 &&i, T2 &&j)
						{
							AERRCHECKINTEGER(T1,"bgt","i")
								AERRCHECKINTEGER(T2,"bgt","j")
								using At1 = typename arrayInfo<T1>::type;
							using At2 = typename arrayInfo<T2>::type;
							using Unsigned1 = typename std::make_unsigned<typename arrayInfo<T1>::type>::type;
							using Unsigned2 = typename std::make_unsigned<typename arrayInfo<T2>::type>::type;
							//Lambda for comparison
							auto l =[](At1 &u1, At2 &u2){return std::bit_cast<Unsigned1>(u1)>std::bit_cast<Unsigned2>(u2);};
							return makeLazy<l>(std::forward<T1>(i),std::forward<T2>(j));

						}

					/**
					 * Fortran BIT_SIZE(https://fortranwiki.org/fortran/show/bit_size)
					 */
					template<typename T>
						auto bit_size(FAR_UNUSED T&& i){
							return CHAR_BIT * sizeof(typename arrayInfo<T>::type);
						}

					/**
					 * Fortran BLE(https://fortranwiki.org/fortran/show/ble)
					 */
					template <typename T1, typename T2>
						auto ble(T1 &&i, T2 &&j)
						{
							AERRCHECKINTEGER(T1,"ble","i")
								AERRCHECKINTEGER(T2,"ble","j")
								using At1 = typename arrayInfo<T1>::type;
							using At2 = typename arrayInfo<T2>::type;
							using Unsigned1 = typename std::make_unsigned<typename arrayInfo<T1>::type>::type;
							using Unsigned2 = typename std::make_unsigned<typename arrayInfo<T2>::type>::type;
							//Lambda for comparison
							auto l =[](At1 &u1, At2 &u2){return std::bit_cast<Unsigned1>(u1)<=std::bit_cast<Unsigned2>(u2);};
							return makeLazy<l>(std::forward<T1>(i),std::forward<T2>(j));
						}

					/**
					 * Fortran BLT(https://fortranwiki.org/fortran/show/blt)
					 */
					template <typename T1, typename T2>
						auto blt(T1 &&i, T2 &&j)
						{
							AERRCHECKINTEGER(T1,"blt","i")
								AERRCHECKINTEGER(T2,"blt","j")
								using At1 = typename arrayInfo<T1>::type;
							using At2 = typename arrayInfo<T2>::type;
							using Unsigned1 = typename std::make_unsigned<typename arrayInfo<T1>::type>::type;
							using Unsigned2 = typename std::make_unsigned<typename arrayInfo<T2>::type>::type;
							//Lambda for comparison
							auto l =[](At1 &u1, At2 &u2){return std::bit_cast<Unsigned1>(u1)<std::bit_cast<Unsigned2>(u2);};
							return makeLazy<l>(std::forward<T1>(i),std::forward<T2>(j));
						}

					/**
					 * Fortran BTEST(https://fortranwiki.org/fortran/show/btest)
					 */
					template <typename T1, typename T2>
						auto btest(T1 &&i, T2 &&pos)
						{
							AERRCHECKINTEGER(T1,"btest","i")
								AERRCHECKINTEGER(T2,"btest","pos")
								auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &i, typename arrayInfo<std::decay_t<T2>>::type const &pos) -> decltype(auto)
								{  return i & 1<<pos; };
							return makeLazy<l>(std::forward<T1>(i),std::forward<T2>(pos));
						}

					/**
					 * Fortran CEILING(https://fortranwiki.org/fortran/show/ceiling)
					 */
					template <typename itype=int,typename T=double>
						auto ceiling(T &&a)
						{
							AERRCHECKREAL(T,"ceiling","a")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &a) -> itype
								{  return std::ceil(a); };
							return makeLazy<l>(std::forward<T>(a));
						}

					/**
					 * Fortran CHAR(https://fortranwiki.org/fortran/show/char)
					 */
					template <typename ctype=char,typename T=char>
						auto f_char(T &&a)
						{ ERRCHECKINTEGER(T,"f_char","a");
							auto l = [](typename arrayInfo<std::decay_t<T>>::type const &a) -> ctype
							{  return static_cast<ctype>(a); };
							return makeLazy<l>(std::forward<T>(a));
						}

					/**
					 * Fortran CMPLX(https://fortranwiki.org/fortran/show/cmplx)
					 */
					template <typename ctype=double,typename T=double>
						auto cmplx(T &&x)
						{
							ERRCHECKREALINTEGERCOMPLEX(T,"cmplx","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> std::complex<ctype>
								{  return std::complex<ctype>(x); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran CMPLX(https://fortranwiki.org/fortran/show/cmplx)
					 * real and complex part
					 */
					template <typename citype=double, typename T1=double, typename T2=double>
						auto cmplx(T1 &&x, T2 &&y)
						{
							ERRCHECKREALINTEGER(T1,"cmplx","x")
								ERRCHECKREALINTEGER(T2,"cmplx","y")
								using ctype = std::complex<citype>;
							auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &x, typename arrayInfo<std::decay_t<T2>>::type const &y) -> ctype
							{  return ctype(x,y); };
							return makeLazy<l>(std::forward<T1>(x),std::forward<T2>(y));
						}

					/**
					 * Fortran CONJG(https://fortranwiki.org/fortran/show/conjg)
					 */
					template <typename T>
						auto conjg(T &&z)
						{
							AERRCHECKCOMPLEX(T,"conjg","z")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &z) -> decltype(auto)
								{  return std::conj(z); };
							return makeLazy<l>(std::forward<T>(z));
						}

					/**
					 * Fortran COS(https://fortranwiki.org/fortran/show/cos)
					 */
					template <typename T>
						auto cos(T &&x)
						{
							ERRCHECKREALINTEGERCOMPLEX(T,"cos","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
								{  return std::cos(x); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran COSD(cos with answer in degrees)
					 */
					template <typename T>
						auto cosd(T &&x)
						{ 
							ERRCHECKREALINTEGERCOMPLEX(T,"cosd","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
								{  return std::cos(x*M_PI/180.0); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran COSH(https://fortranwiki.org/fortran/show/cosh)
					 */
					template <typename T>
						auto cosh(T &&x)
						{
							ERRCHECKREALINTEGERCOMPLEX(T,"cosh","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
								{  return std::cosh(x); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran COSPI(cos with argument on [-1,1])
					 */
					template <typename T>
						auto cospi(T &&x)
						{
							ERRCHECKREALINTEGERCOMPLEX(T,"cospi","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
								{  return std::cos(x*M_PI); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran COUNT(https://fortranwiki.org/fortran/show/count)
					 */
					template <typename kind=FAR_UNSIGNED_INDEX_TYPE, typename T=double>
						auto count(T &&source)
						{
							AERRCHECKTRUTHY(T,"count","mask")
								using srcType = typename arrayInfo<T>::type;
							return reduction<[](const srcType &src, kind &dest){dest+=src?kind(1):kind(0);},
										 [](const kind &src, kind &dest){dest=dest+src;},
										 [](kind &dest){dest=0;}
							>(std::forward<T>(source));
						}


					/**
					 * Fortran COUNT(https://fortranwiki.org/fortran/show/count)
					 * with dir
					 */
					template <typename kind=FAR_UNSIGNED_INDEX_TYPE, typename T=double>
						auto count(T &&source, int dir)
						{
							AERRCHECKTRUTHY(T,"count","mask")
								using srcType = typename arrayInfo<T>::type;
							return reductionWithDirection<[](const srcType &src, kind &dest){dest+=src?1ll:0ll;},
										 [](const kind &src, kind &dest){dest=dest+src;},
										 [](kind &dest){dest=0;}
							>(dir,std::forward<T>(source));
						}


					/**
					 * Fortran CSHIFT(https://fortranwiki.org/fortran/show/cshift)
					 */
					template <typename T, typename T_shift>
						auto cshift(T &&array, T_shift &&shift, int dim=1)
						{
							static_assert(arrayInfo<std::decay_t<T>>::value,"cshift can only be applied to arrays");
							static_assert(arrayInfo<T_shift>::rank == arrayInfo<T>::rank -1 
									|| arrayInfo<T_shift>::rank == 0, "In cshift, shift must either be scalar or rank 1 less than array");
							using index = typename arrayInfo<std::decay_t<T>>::index;
							const int rank = arrayInfo<std::decay_t<T>>::rank;
							//Note! Internally dim is already shifted down to zero base!
							auto l = [](indexInfo<rank> info,
									const T &op, const T_shift& shift, const int dim) -> decltype(auto)
							{
								typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, rank>::type tpl;
								typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, rank-1>::type shift_tpl;
								if constexpr(arrayInfo<std::decay_t<T_shift>>::value){
									reduceTuple<0,0>(info.indices, shift_tpl, dim);
									int lshift=callDefaultBase(shift,shift_tpl);
									assignTupleLevel<0> (tpl, dim, lshift);
								} else {
									assignTupleLevel<0> (tpl, dim, shift);
								}
								op.getIndexer().cshiftIndex(info.indices,tpl);
								return std::apply(op,info.indices);
							};

							//This line shifts dim to zero base
							auto a = lazyArray<l, index, std::decay_t<T>, T_shift, int>(std::forward<T>(array), 
									std::forward<T_shift>(shift), dim-1);
							a.setIndexShifted();
							return a;
						}

					/**
					 * Fortran DBLE(https://fortranwiki.org/fortran/show/dble)
					 */
					template <typename conv=double, typename T=double>
						auto dble(T &&a)
						{
							auto l = [](typename arrayInfo<std::decay_t<T>>::type const &a) -> conv
							{  
								if constexpr(is_complex<typename arrayInfo<T>::type>){
									return static_cast<conv>(a.real());
								} else {
									return static_cast<conv>(a);
								}
							};
							return makeLazy<l>(std::forward<T>(a));
						}

					/**
					 * Fortran DEALLOCATE statement
					 */
					template<typename T_first, typename... T_others>
						void deallocate(T_first && first, T_others&&... others){
							first.deallocate();
							if constexpr(sizeof...(others)>0) deallocate(others...);
						}

					/**
					 * Dereference operator. Converts an array of pointers into a lazyArray that returns references
					 * to the pointed to values
					 */
					template<typename T, std::enable_if_t<std::is_pointer_v<typename arrayInfo<T>::type>,int> = 0>
						auto dereference(T &ptr){
							auto l = [](typename arrayInfo<T>::type ptr) -> decltype(auto){
								return *ptr;
							};
							return makeLazy<l>(ptr);
						}

					/**
					 * Fortran DIGITS(https://fortranwiki.org/fortran/show/digits)
					 */
					template <typename T>
						auto digits(FAR_UNUSED T &&a)
						{
							if constexpr (!is_complex<T>){
								return std::numeric_limits<typename arrayInfo<T>::type>::digits;
							} else {
								return std::numeric_limits<typename arrayInfo<T>::type::value_type>::digits;
							}
						}


					/**
					 * Fortran DIM(https://fortranwiki.org/fortran/show/dim)
					 */
					template <typename T1, typename T2>
						auto dim(T1 &&x, T2 &&y)
						{
							auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &x, 
									typename arrayInfo<std::decay_t<T2>>::type const &y) -> decltype(auto)
							{  return x-y>0?x-y:0; };

							return makeLazy<l>(std::forward<T1>(x),std::forward<T2>(y));
						}


					/**
					 * Fortran DOT_PRODUCT(https://fortranwiki.org/fortran/show/dot_product)
					 * numeric case
					 */
					template <typename T1, typename T2, std::enable_if_t<is_arithmetic_all<T1,T2> && !is_logical_all<T1,T2>, int> = 0>
						auto dot_product(T1 &&vector_a, T2 &&vector_b)
						{
							static_assert(arrayInfo<T1>::rank == 1 && arrayInfo<T2>::rank==1,"Both arrays must be rank 1 in dot_product");
							typename arrayInfo<T1>::type result={};
							forEachSerial<[](const typename arrayInfo<T1>::type &vector_a, const typename arrayInfo<T2>::type &vector_b, 
									typename arrayInfo<T1>::type &result){
								result+=(vector_a*vector_b); 
								}>(std::forward<T1>(vector_a),std::forward<T2>(vector_b),result);
							return result;
						}

					/**
					 * Fortran DOT_PRODUCT(https://fortranwiki.org/fortran/show/dot_product)
					 * complex case
					 */
					template <typename T1, typename T2, std::enable_if_t<is_complex_all<T1,T2>, int> = 0>
						auto dot_product(T1 &&vector_a, T2 &&vector_b)
						{
							static_assert(arrayInfo<T1>::rank == 1 && arrayInfo<T2>::rank==1,"Both arrays must be rank 1 in dot_product");
							typename arrayInfo<T1>::type result{};
							forEach<[](const typename arrayInfo<T1>::type &vector_a, const typename arrayInfo<T2>::type &vector_b,
									typename arrayInfo<T1>::type &result){
								result+=(std::conj(vector_a)*vector_b); }>(std::forward<T1>(vector_a),std::forward<T2>(vector_b),result);
							return result;
						}

					/**
					 * Fortran DOT_PRODUCT(https://fortranwiki.org/fortran/show/dot_product)
					 * logical case
					 */
					template <typename T1, typename T2, std::enable_if_t<is_logical_all<T1,T2>, int> = 0>
						auto dot_product(T1 &&vector_a, T2 &&vector_b)
						{
							static_assert(arrayInfo<T1>::rank == 1 && arrayInfo<T2>::rank==1,"Both arrays must be rank 1 in dot_product");
							typename arrayInfo<T1>::type result{};
							forEach<[](const typename arrayInfo<T1>::type &vector_a, const typename arrayInfo<T2>::type &vector_b,
									typename arrayInfo<T1>::type &result){
								result=result || (vector_a&&vector_b); }>(std::forward<T1>(vector_a),std::forward<T2>(vector_b),result);
							return result;
						}

#ifndef STRICT_TYPE

					/**
					 * Fortran DOT_PRODUCT(https://fortranwiki.org/fortran/show/dot_product)
					 * all other cases. Only exists with strict typing turned off
					 */
					template <typename T1, typename T2, std::enable_if_t<!(is_logical_all<T1,T2> || is_complex_all<T1,T2> || is_arithmetic_all<T1,T2>), int> = 0>
						auto dot_product(T1 &&vector_a, T2 &&vector_b)
						{
							static_assert(arrayInfo<T1>::rank == 1 && arrayInfo<T2>::rank==1,"Both arrays must be rank 1 in dot_product");
							typename arrayInfo<T1>::type result{};
							forEach<[](const typename arrayInfo<T1>::type &vector_a, const typename arrayInfo<T2>::type &vector_b,
									typename arrayInfo<T1>::type &result){
								result=result + (vector_a*vector_b); }>(std::forward<T1>(vector_a),std::forward<T2>(vector_b),result);
							return result;
						}

#endif

					/**
					 * Fortran DPROD(https://fortranwiki.org/fortran/show/dprod)
					 */
					template <typename T1, typename T2>
						auto dprod(T1 &&x, T2 &&y)
						{
							AERRCHECKREAL(T1,"dprod","x")
								AERRCHECKREAL(T2,"dprod","y")
								auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &x, typename arrayInfo<std::decay_t<T2>>::type const &y) -> double
								{  return (double)x*(double)y; };

							return makeLazy<l>(std::forward<T1>(x),std::forward<T2>(y));
						}

					/**
					 * Fortran DSHIFTL(https://fortranwiki.org/fortran/show/dshiftl)
					 */
					template <typename T1, typename T2, typename T3>
						auto dshiftl(T1 &&i, T2 &&j, T3 &&shift)
						{
							AERRCHECKINTEGER(T1,"dshiftl","x")
								AERRCHECKINTEGER(T2,"dshiftl","y")
								AERRCHECKINTEGER(T3,"dshiftl","shift")
								assert(shift>0);
							auto l = [](typename arrayInfo<T1>::type const &i, typename arrayInfo<T2>::type const &j, int shift) -> decltype(auto)
							{
								using type = larger_type<std::decay_t<decltype(i)>,std::decay_t<decltype(j)>>;
								using utype = std::make_unsigned_t<type>;
								utype ip=i;
								utype jp=j;
								utype left = ip << shift;
								utype right = jp >> (CHAR_BIT * sizeof(type) - shift);
								return type(left|right);
							};
							return makeLazy<l>(std::forward<T1>(i), std::forward<T2>(j), std::forward<T3>(shift));
						}

					/**
					 * Fortran DSHIFTR(https://fortranwiki.org/fortran/show/dshiftr)
					 */
					template <typename T1, typename T2, typename T3>
						auto dshiftr(T1 &&i, T2 &&j, T3 &&shift)
						{
							assert(shift>0);
							auto l = [](typename arrayInfo<T1>::type const &i, typename arrayInfo<T2>::type const &j, int shift) -> decltype(auto)
							{
								using type = larger_type<std::decay_t<decltype(i)>,std::decay_t<decltype(j)>>;
								using utype = std::make_unsigned_t<type>;
								utype ip=i;
								utype jp=j;
								utype right = jp >> shift;
								utype left = ip << (CHAR_BIT * sizeof(type) - shift);
								return type(left | right);
							};

							return makeLazy<l>(std::forward<T1>(i),std::forward<T2>(j),std::forward<T3>(shift));
						}

					/**
					 * Fortran EOSHIFT(https://fortranwiki.org/fortran/show/eoshift)
					 */
					template <typename T, typename T_shift>
						auto eoshift(T &&array, T_shift &&shift, int dim=1)
						{
							static_assert(arrayInfo<std::decay_t<T>>::value,"eoshift can only be applied to arrays");
							static_assert(!arrayInfo<std::decay_t<T_shift>>::value
									|| arrayInfo<std::decay_t<T_shift>>::rank == arrayInfo<std::decay_t<T>>::rank -1,
									"eoshift shift must be scalar or have rank 1 less than array");
							using index = typename arrayInfo<std::decay_t<T>>::index;
							const typename arrayInfo<std::decay_t<T>>::type boundary = {};
							const int rank = arrayInfo<std::decay_t<T>>::rank;
							//Note! Internally dim is already shifted down to zero base!
							auto l = [](indexInfo<rank> info,
									const T &op, const T_shift& shift, const typename std::remove_reference_t<arrayInfo<std::decay_t<T>>>::type boundary, 
									const int dim) -> typename arrayInfo<std::decay_t<T>>::type
							{
								typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, rank>::type tpl;
								typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, rank-1>::type shift_tpl;
								if constexpr(arrayInfo<std::decay_t<T_shift>>::value){
									reduceTuple<0,0>(info.indices, shift_tpl, dim);
									int lshift=callDefaultBase(shift,shift_tpl);
									assignTupleLevel<0> (tpl, dim, lshift);
								} else {
									assignTupleLevel<0> (tpl, dim, shift);
								}
								bool oob = op.getIndexer().eoshiftIndex(info.indices,tpl);
								return oob?boundary:std::apply(op,info.indices);
							};

							//This line shifts dim to zero base
							auto a = lazyArray<l, index, std::decay_t<T>, T_shift, typename arrayInfo<std::decay_t<T>>::type, int>
								(std::forward<T>(array), std::forward<T_shift>(shift), boundary, dim-1);
							a.setIndexShifted();
							return a;
						}

					/**
					 * Fortran EOSHIFT(https://fortranwiki.org/fortran/show/eoshift)
					 */
					template <typename T, typename T_shift, typename T_boundary>
						auto eoshift_with_boundary(T &&array, T_shift &&shift, T_boundary&& boundary, int dim=1)
						{
							static_assert(arrayInfo<std::decay_t<T>>::value,"eoshift can only be applied to arrays");
							static_assert(!arrayInfo<std::decay_t<T_shift>>::value
									|| arrayInfo<std::decay_t<T_shift>>::rank == arrayInfo<std::decay_t<T>>::rank -1,
									"eoshift shift must be scalar or have rank 1 less than array");
							using index = typename arrayInfo<std::decay_t<T>>::index;
							const int rank = arrayInfo<std::decay_t<T>>::rank;
							//Note! Internally dim is already shifted down to zero base!
							auto l = [](indexInfo<rank> info,
									const T &op, const T_shift& shift, const T_boundary& boundary,
									const int dim) -> typename arrayInfo<std::decay_t<T>>::type
							{
								typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, rank>::type tpl;
								typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE, rank-1>::type shift_tpl;
								if constexpr(arrayInfo<std::decay_t<T_shift>>::value || arrayInfo<std::decay_t<T_boundary>>::value){
									reduceTuple<0,0>(info.indices, shift_tpl, dim);
								}
								if constexpr(arrayInfo<std::decay_t<T_shift>>::value){
									int lshift=callDefaultBase(shift,shift_tpl);
									assignTupleLevel<0> (tpl, dim, lshift);
								} else {
									assignTupleLevel<0> (tpl, dim, shift);
								}
								bool oob = op.getIndexer().eoshiftIndex(info.indices,tpl);
								if constexpr(arrayInfo<T_boundary>::value){
									return oob?callDefaultBase(boundary,shift_tpl):std::apply(op,info.indices);
								} else {
									return oob?boundary:std::apply(op,info.indices);
								}
							};

							//This line shifts dim to zero base
							auto a = lazyArray<l, index, std::decay_t<T>, std::decay_t<T_shift>, std::decay_t<T_boundary>, int>
								(std::forward<T>(array), std::forward<T_shift>(shift), std::forward<T_boundary>(boundary), dim-1);
							a.setIndexShifted();
							return a;
						}

					/**
					 * Fortran EPSILON(https://fortranwiki.org/fortran/show/epsilon)
					 */
					template <typename T, typename = std::enable_if<is_floating_point<T>>>
						auto epsilon(FAR_UNUSED T &&x)
						{
							return std::numeric_limits<typename arrayInfo<T>::type>::epsilon;
						}

					/**
					 * Fortran ERF(https://fortranwiki.org/fortran/show/erf)
					 */
					template <typename T>
						auto erf(T &&x)
						{
							AERRCHECKREAL(T,"erf","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
								{  return std::erf(x); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran ERFC(https://fortranwiki.org/fortran/show/erfc)
					 */
					template <typename T>
						auto erfc(T &&x)
						{
							AERRCHECKREAL(T,"erfc","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
								{  return std::erfc(x); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran ERFC_SCALED(https://fortranwiki.org/fortran/show/erfc_scaled)
					 */
					template <typename T>
						auto erfc_scaled(T &&x)
						{
							AERRCHECKREAL(T,"erfc_scaled","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
								{  return std::exp(x*x) * std::erfc(x); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran EXP(https://fortranwiki.org/fortran/show/exp)
					 */
					template <typename T>
						auto exp(T &&x)
						{
							AERRCHECKREAL(T,"exp","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
								{  return std::exp(x); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran EXPONENT(https://fortranwiki.org/fortran/show/exponent)
					 */
					template <typename T>
						auto exponent(T &&x)
						{
							AERRCHECKREAL(T,"exponent","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
								{
									//This is unnecessary since C++ pretty much assumes that radix is 2, but Fortran
									//doesn't and C++ might not in the future
									constexpr int radix = std::numeric_limits<typename arrayInfo<T>::type>::radix;
									int exponent;
									std::frexp(x, &exponent);
									//If support C++26 use std::log2(radix) as constexpr
									//Use preprocessor to determine if C++26 is supported
									#if __cplusplus > 202303L
										constexpr typename arrayInfo<T>::type log_radix = std::log2(radix);
									#else
										typename arrayInfo<T>::type log_radix = std::log2(radix);
									#endif
									return static_cast<int>(exponent / log_radix);
								};
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran EXTENDS_TYPE_OF(https://fortranwiki.org/fortran/show/extends_type_of)
					 * There is no equivalent of a Fortran unlimited polymorphic 
					 */
					template <typename T, typename T_mold>
						std::enable_if<!std::is_same_v<T,std::any>,bool> extends_type_of(FAR_UNUSED T &&a, FAR_UNUSED T_mold &&mold)
						{
							return std::is_base_of_v<T_mold, T>;
						}

					/**
					 * Fortran EXTENDS_TYPE_OF(https://fortranwiki.org/fortran/show/extends_type_of)
					 */
					template <typename T, typename T_mold, typename = std::enable_if<std::is_same_v<T,std::any>>>
						std::enable_if<std::is_same_v<T,std::any>,bool> extends_type_of(T &&a, FAR_UNUSED T_mold &&mold)
						{
							try{
								auto &ref = any_cast<T_mold>(a);
								return true;
							} catch(const std::bad_any_cast &e){
								return false;
							}
						}


					/**
					 * Fortran FINDLOC
					 * Array, value, kind and back
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T=void, typename T_value = void, bool natural=false>
						auto findloc(T &&array, T_value &&value, bool back = false){
							ERRCHECKREALINTEGERSTRING(T,"findloc","array");
							static_assert(!arrayInfo<T_value>::value,"Value must be scalar in findloc");
							constexpr int rank = arrayInfo<T>::rank;
							using rtype = indexInfo<rank>;
							bool active;
							Array<kind,1> loc(rank);
							loc=defaultLB-1;
							active=true;

							//Find the first valid element
							auto findFirst=[](const rtype& rri, const T& value, decltype(loc)& loc, const T_value& target, decltype(active)& active){

								if (active && getItem(value,rri)==target){
									active = false;
									if constexpr(natural){
										assignArrayTuple(loc,rri.SBindices);
									} else {
										assignArrayTuple(loc,rri.indices);
									}
									return iaTerminate;
								}
								return iaContinue;
							};
							if (back){
								auto combine = [](const Array<kind,1> &src, Array<kind,1> &dest){
									for (FAR_UNSIGNED_INDEX_TYPE i=src.getSize();i>=1;--i){
										if (src(i)==0) return;
										if (src(i)<dest(i) && dest(i)!=0) return;
									}
									//Remember can't use assignment here because
									//Won't work with OpenMP
									for (FAR_UNSIGNED_INDEX_TYPE i=defaultLB;i<defaultLB+src.getSize();++i){
										dest(i)=src(i);
									}
								};
								for (FAR_UNSIGNED_INDEX_TYPE i=defaultLB;i<defaultLB+rank;++i) loc(i)=0;
								loc = reduction<findFirst,
										combine,
										[](Array<kind,1> &dest){dest.allocate(rank);for (int i=defaultLB;i<defaultLB+rank;++i) dest(i)=0;},
										nullFunction, true
											>(std::forward<T>(array),std::forward<T_value>(value),active);
							} else {
								auto combine = [](const Array<kind,1> &src, Array<kind,1> &dest){
									for (FAR_UNSIGNED_INDEX_TYPE i=src.getSize();i>=1;--i){
										if (src(i)==0) return;
										if (src(i)>dest(i) && dest(i)!=0) return;
									}
									//Remember can't use assignment here because
									//Won't work with OpenMP
									for (FAR_UNSIGNED_INDEX_TYPE i=defaultLB;i<defaultLB+src.getSize();++i){
										dest(i)=src(i);
									}
								};
								for (FAR_UNSIGNED_INDEX_TYPE i=defaultLB;i<defaultLB+rank;++i) loc(i)=0;
								loc = reduction<findFirst,
										combine,
										[](Array<kind,1> &dest){dest.allocate(rank);for (int i=defaultLB;i<defaultLB+rank;++i) dest(i)=0;},
										nullFunction
											>(std::forward<T>(array),std::forward<T_value>(value),active);
							}
							return loc;
						}

					/**
					 * Fortran FINDLOC
					 * Array, value, mask, kind and back
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T=void, typename T_value = void, typename T_mask=void, bool natural=false>
						auto findloc_with_mask(T &&array, T_value &&value, T_mask &&mask, bool back = false){
							ERRCHECKREALINTEGERSTRING(T,"findloc","array");
							AERRCHECKTRUTHY(T_mask,"findloc","mask");
							static_assert(!arrayInfo<T_value>::value,"Value must be scalar in findloc");
							constexpr int rank = arrayInfo<T>::rank;
							using rtype = indexInfo<rank>;
							bool active;
							Array<kind,1> loc(rank);
							loc=defaultLB-1;
							active=true;

							//Find the first valid element
							auto findFirst=[](const rtype& rri, const T& value, decltype(loc)& loc, const T_value& target, decltype(active)& active, const T_mask& mask){

								if (active && getItem(mask,rri) && getItem(value,rri)==target){
									active = false;
									if constexpr(natural){
										assignArrayTuple(loc,rri.SBindices);
									} else {
										assignArrayTuple(loc,rri.indices);
									}
									return iaTerminate;
								}
								return iaContinue;
							};
							if (back){
								auto combine = [](const Array<kind,1> &src, Array<kind,1> &dest){
									for (FAR_UNSIGNED_INDEX_TYPE i=src.getSize();i>=1;--i){
										if (src(i)==0) return;
										if (src(i)<dest(i) && dest(i)!=0) return;
									}
									//Remember can't use assignment here because
									//Won't work with OpenMP
									for (FAR_UNSIGNED_INDEX_TYPE i=defaultLB;i<defaultLB+src.getSize();++i){
										dest(i)=src(i);
									}
								};
								for (FAR_UNSIGNED_INDEX_TYPE i=defaultLB;i<defaultLB+rank;++i) loc(i)=0;
								loc = reduction<findFirst,
										combine,
										[](Array<kind,1> &dest){dest.allocate(rank);for (int i=defaultLB;i<defaultLB+rank;++i) dest(i)=0;},
										nullFunction, true
											>(std::forward<T>(array),std::forward<T_value>(value),active,std::forward<T_mask>(mask));
							} else {
								auto combine = [](const Array<kind,1> &src, Array<kind,1> &dest){
									for (FAR_UNSIGNED_INDEX_TYPE i=src.getSize();i>=1;--i){
										if (src(i)==0) return;
										if (src(i)>dest(i) && dest(i)!=0) return;
									}
									//Remember can't use assignment here because
									//Won't work with OpenMP
									for (FAR_UNSIGNED_INDEX_TYPE i=defaultLB;i<defaultLB+src.getSize();++i){
										dest(i)=src(i);
									}
								};
								for (FAR_UNSIGNED_INDEX_TYPE i=defaultLB;i<defaultLB+rank;++i) loc(i)=0;
								loc = reduction<findFirst,
										combine,
										[](Array<kind,1> &dest){dest.allocate(rank);for (int i=defaultLB;i<defaultLB+rank;++i) dest(i)=0;},
										nullFunction
											>(std::forward<T>(array),std::forward<T_value>(value),active,std::forward<T_mask>(mask));
							}
							return loc;
						}



					/**
					 * Fortran FINDLOC
					 * Array, mask, kind and back
					 * with direction
					 */

					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T=void, typename T_value, typename T_mask=void, bool natural=false>
						auto findloc_with_mask(T &&array, T_value &&value, T_mask&& mask, int direction, bool back=false){
							ERRCHECKREALINTEGERSTRING(T,"findloc","array");
							AERRCHECKTRUTHY(T_mask,"findloc","mask");
							static_assert(!arrayInfo<T_value>::value,"value must be scalar in findloc");
							constexpr int rank = arrayInfo<T>::rank;
							using rtype = indexInfo<rank>;
							int rdirection = direction-1;
							auto active = reduce_rank<bool>(array,rdirection);
							auto loc = reduce_rank<kind>(array,rdirection);
							active = true;
							//Find the first valid element
							auto findFirst=[](const rtype& rri, const T& value, decltype(loc)& loc, const T_value& target, const T_mask &mask, decltype(active)& active){
								if (getItem(mask,rri) && getItem(active,rri) && getItem(value,rri)==target){
									getItem(active,rri) = false;
									if constexpr(natural){
										getItem(loc,rri) = getTupleLevel(rri.SBindices,rri.direction);
									} else {
										getItem(loc,rri) = getTupleLevel(rri.indices,rri.direction);
									}
									return iaTerminate;
								}
								return iaContinue;
							};
							if (back){
								auto combine = [](const Array<kind,rank-1> &src, Array<kind,rank-1> &dest){
									auto elementCombine = [](const kind &src, kind &dest){
										if (src==0) return;
										if (src>dest && dest!=0) return;
										dest=src;
									};
									forEachSerial<elementCombine,true>(src,dest);
								};
								for (int i=defaultLB;i<defaultLB+rank;++i) loc(i)=0;
								loc = reductionWithDirection<findFirst,
										combine,
										[](kind &val){val=defaultLB-1;},
										nullFunction, true
											>(direction, std::forward<T>(array), std::forward<T_value>(value), std::forward<T_mask>(mask), active);

							} else {
								auto combine = [](const Array<kind,rank-1> &src, Array<kind,rank-1> &dest){
									auto elementCombine = [](const kind &src, kind &dest){
										if (src==0) return;
										if (src<dest && dest!=0) return;
										dest=src;
									};
									forEachSerial<elementCombine,true>(src,dest);
								};
								for (int i=defaultLB;i<defaultLB+rank;++i) loc(i)=0;
								loc = reductionWithDirection<findFirst,
										combine,
										[](kind &val){val=defaultLB-1;},
										nullFunction
											>(direction, std::forward<T>(array), std::forward<T_value>(value), std::forward<T_mask>(mask), active);
							}
							return loc;				
						}


					/**
					 * Fortran FINDLOC
					 * Array, kind and back
					 * with direction
					 */

					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T=void, typename T_value, bool natural=false>
						auto findloc(T &&array, T_value &&value, int direction, bool back=false){
							ERRCHECKREALINTEGERSTRING(T,"findloc","array");
							static_assert(!arrayInfo<T_value>::value,"value must be scalar in findloc");
							constexpr int rank = arrayInfo<T>::rank;
							using rtype = indexInfo<rank>;
							int rdirection = direction-1;
							auto active = reduce_rank<bool>(array,rdirection);
							auto loc = reduce_rank<kind>(array,rdirection);
							active = true;
							//Find the first valid element
							auto findFirst=[](const rtype& rri, const T& value, decltype(loc)& loc, const T_value& target, decltype(active)& active){
								if (getItem(active,rri) && getItem(value,rri)==target){
									getItem(active,rri) = false;
									if constexpr(natural){
										getItem(loc,rri) = getTupleLevel(rri.SBindices,rri.direction);
									} else {
										getItem(loc,rri) = getTupleLevel(rri.indices,rri.direction);
									}
									return iaTerminate;
								}
								return iaContinue;
							};
							if (back){
								auto combine = [](const Array<kind,rank-1> &src, Array<kind,rank-1> &dest){
									auto elementCombine = [](const kind &src, kind &dest){
										if (src==0) return;
										if (src>dest && dest!=0) return;
										//#pragma omp critical
										dest=src;
									};
									forEachSerial<elementCombine,true>(src,dest);
								};
								for (int i=defaultLB;i<defaultLB+rank;++i) loc(i)=0;
								loc = reductionWithDirection<findFirst,
										combine,
										[](kind&val){val=0;},
										nullFunction, true
											>(direction, std::forward<T>(array),std::forward<T_value>(value),active);

							} else {
								auto combine = [](const Array<kind,rank-1> &src, Array<kind,rank-1> &dest){
									auto elementCombine = [](const kind &src, kind &dest){
										if (src==0) return;
										if (src<dest && dest!=0) return;
										dest=src;
									};
									forEachSerial<elementCombine,true>(src,dest);
								};
								for (int i=defaultLB;i<defaultLB+rank;++i) loc(i)=0;
								loc = reductionWithDirection<findFirst,
										combine,
										[](kind &val){val=defaultLB-1;},
										nullFunction
											>(direction, std::forward<T>(array),std::forward<T_value>(value),active);
							}
							return loc;
						}


					/**
					 * Fortran FLOAT(https://fortranwiki.org/fortran/show/float)
					 */
					template <typename T, typename = std::enable_if<is_integral<T>>>
						auto f_float(T &&a)
						{
							AERRCHECKREAL(T,"f_float","a")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &a) -> float
								{ return static_cast<float>(a); };
							return makeLazy<l>(std::forward<T>(a));
						}

					/**
					 * Fortran FLOOR(https://fortranwiki.org/fortran/show/floor)
					 */
					template <typename kind=double, typename T=void, typename = std::enable_if<is_integral<T>>>
						auto floor(T &&a)
						{
							AERRCHECKREAL(T,"floor","a")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &a) -> float
								{ return static_cast<kind>(std::floor(a)); };
							return makeLazy<l>(std::forward<T>(a));
						}

					/**
					 * far::Array equivalent of std::for_each
					 * NB: This only works outside workshare mode if the generator has state
					 */
					template <typename T, typename T_generator, std::enable_if_t<(std::is_invocable_v<T_generator, typename arrayInfo<T>::type&>), int>* = nullptr>
						void for_each(T &&array, T_generator generator)
						{
							using core_type = typename arrayInfo<T>::type;
							forEachTight<([](core_type &dest, T_generator generator)
									{  generator(dest); })>(std::forward<T>(array),std::forward<T_generator>(generator));
						}

					/**
					 * Fortran FRACTION(https://fortranwiki.org/fortran/show/fraction)
					 */
					template <typename T, typename = std::enable_if<is_floating_point<T>>>
						auto fraction(T &&x)
						{
							AERRCHECKREAL(T,"fraction","x")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
								{ 
									int exponent;
									typename arrayInfo<T>::type fraction = std::frexp(x, &exponent);
									return fraction;
								};
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran GAMMA(https://fortranwiki.org/fortran/show/gamma)
					 */
					template <typename T>
						auto f_gamma(T &&x)
						{
							ERRCHECKREAL(T,"f_gamma","x");
							auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> decltype(auto)
							{  return std::tgamma(val); };
							return makeLazy<l>(std::forward<T>(x));
						}
					/**
					 * far::Array equivalent of std::generate
					 */
					template <typename T, typename T_generator, std::enable_if_t<(std::is_invocable_v<T_generator>), int>* = nullptr>
						void generate(T &&array, T_generator generator)
						{
							using core_type = typename arrayInfo<T>::type;
							forEachSerial<([](core_type &dest, T_generator generator)
									{  dest = generator(); })>(std::forward<T>(array),std::forward<T_generator>(generator));
						}

					/**
					 * Fortran HUGE(https://fortranwiki.org/fortran/show/huge)
					 */
					template <typename T>
						constexpr auto huge(FAR_UNUSED T &&x)
						{
							return std::numeric_limits<typename arrayInfo<T>::type>::max();
						}

					/**
					 * Fortran HYPOT(https://fortranwiki.org/fortran/show/hypot)
					 */
					template <typename T1, typename T2>
						auto hypot(T1 &&x, T2 &&y)
						{
							AERRCHECKREAL(T1,"hypot","x")
								AERRCHECKREAL(T2,"hypot","y")
								auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &x, typename arrayInfo<std::decay_t<T2>>::type const &y) 
								-> decltype(auto)
								{
									return std::hypot(x,y);
								};
							return makeLazy<l>(std::forward<T1>(x),std::forward<T2>(y));
						}

					/**
					 * Fortran IACHAR(https://fortranwiki.org/fortran/show/iachar)
					 */
					template <typename kind=char, typename T=char>
						auto iachar(T &&source)
						{
							AERRCHECKINTEGER(T,"iachar","c")
								auto l = [](typename arrayInfo<std::decay_t<T>>::type const &val) -> kind
								{  return char_to_ascii[val]; };
							return makeLazy<l>(std::forward<T>(source));
						}

					/**
					 * Fortran IALL(https://fortranwiki.org/fortran/show/iall)
					 */
					template <typename T>
						auto iall(T &&source)
						{
							AERRCHECKINTEGER(T,"iall","array")					
								using type=arrayInfo<T>::type;
							return reduction<[](const type &src, type &dest){dest&=src;},
										 defaultFunction,
										 [](type &dest){dest=~type(0);}
							>(std::forward<T>(source));
						}

					/**
					 * Fortran IALL(https://fortranwiki.org/fortran/show/iall)
					 * with mask
					 */
					template <typename T, typename T_mask>
						auto iall_with_mask(T &&source, T_mask &&mask)
						{
							AERRCHECKINTEGER(T,"iall","array")
								AERRCHECKINTEGER(T_mask,"iall","mask")
								using type=arrayInfo<T>::type;
							using maskType = arrayInfo<T_mask>::type;
							return reduction<[](const type &src, type &dest, const maskType &mask){dest&=mask?src:~type(0);},
										 defaultFunction,
										 [](type &dest){dest=~type(0);}
							>(std::forward<T>(source),std::forward<T_mask>(mask));
						}

					/**
					 * Fortran IALL(https://fortranwiki.org/fortran/show/iall)
					 * With dim
					 */
					template <typename T>
						auto iall(T &&source, int dim)
						{
							AERRCHECKINTEGER(T,"iall","array")
								// Produce an array with the summed direction removed
								using type=arrayInfo<T>::type;
							return reductionWithDirection<[](const type &src, type &dest){dest&=src;},
										 defaultFunction,
										 [](type &dest){dest=~type(0);}
							>(dim, std::forward<T>(source));
						}

					/**
					 * Fortran IALL(https://fortranwiki.org/fortran/show/iall)
					 * With dim and mask
					 */
					template <typename T, typename T_mask>
						auto iall_with_mask(T &&source, int dim, T_mask &&mask)
						{
							AERRCHECKINTEGER(T,"iall","array")
								AERRCHECKTRUTHY(T_mask,"iall","mask")
								// Produce an array with the summed direction removed
								using type=arrayInfo<T>::type;
							using maskType = arrayInfo<T_mask>::type;
							return reductionWithDirection<[](const type &src, type &dest, const maskType &mask){
								dest&=mask?src:~type(0);},
										 defaultFunction,
										 [](type &dest){dest=~type(0);}
							>(dim,std::forward<T>(source),std::forward<T_mask>(mask));
						}

					/**
					 * Fortran IAND(https://fortranwiki.org/fortran/show/iand)
					 */
					template <typename T_i, typename T_j>
						auto iand(T_i &&i, T_j &&j)
						{
							AERRCHECKINTEGER(T_i,"iand","i")
								AERRCHECKINTEGER(T_j,"iand","j")
								auto l = [](typename arrayInfo<std::decay_t<T_i>>::type const &i,
										typename arrayInfo<std::decay_t<T_j>>::type const &j) -> decltype(auto)
								{  return i&j; };

							return makeLazy<l>(std::forward<T_i>(i),std::forward<T_j>(j));
						}

					/**
					 * Fortran IANY(https://fortranwiki.org/fortran/show/iany)
					 */
					template <typename T>
						auto iany(T &&source)
						{
							AERRCHECKINTEGER(T,"iany","array")
								using type=arrayInfo<T>::type;
							return reduction<[](const type &src, type &dest){dest|=src;},
										 defaultFunction,
										 [](type &dest){dest=type(0);}
							>(std::forward<T>(source));
						}

					/**
					 * Fortran IANY(https://fortranwiki.org/fortran/show/iany)
					 * with mask
					 */
					template <typename T, typename T_mask>
						auto iany_with_mask(T &&source, T_mask &&mask)
						{
							AERRCHECKINTEGER(T,"iany","array")
							AERRCHECKTRUTHY(T_mask,"iany","mask")
							using type=arrayInfo<T>::type;
							using maskType = arrayInfo<T_mask>::type;
							return reduction<[](const type &src, type &dest, const maskType &mask){dest|=mask?src:type(0);},
										 defaultFunction,
										 [](type &dest){dest=type(0);}
							>(std::forward<T>(source),std::forward<T_mask>(mask));

						}

					/**
					 * Fortran IANY(https://fortranwiki.org/fortran/show/iany)
					 * With dim
					 */
					template <typename T>
						auto iany(T &&source, int dim)
						{
							AERRCHECKINTEGER(T,"iany","array")
							// Produce an array with the summed direction removed
							using type=arrayInfo<T>::type;
							return reductionWithDirection<[](const type &src, type &dest){dest|=src;},
										 defaultFunction,
										 [](type &dest){dest=type(0);}
							>(dim, std::forward<T>(source));
						}

					/**
					 * Fortran IANY(https://fortranwiki.org/fortran/show/iany)
					 * With dim and mask
					 */
					template <typename T, typename T_mask>
						auto iany_with_mask(T &&source, int dim, T_mask &&mask)
						{
							AERRCHECKINTEGER(T,"iany","array")
							AERRCHECKTRUTHY(T_mask,"iany","mask")
							using type=arrayInfo<T>::type;
							using maskType = arrayInfo<T_mask>::type;
							return reductionWithDirection<[](const type &src, type &dest, const maskType &mask){dest|=mask?src:type(0);},
										 defaultFunction,
										 [](type &dest){dest=type(0);}
							>(dim,std::forward<T>(source),std::forward<T_mask>(mask));
						}

					/**
					 * Fortran IBCLR(https://fortranwiki.org/fortran/show/ibclr)
					 */
					template <typename T_i, typename T_pos>
						auto ibclr(T_i &&i, T_pos &&pos)
						{
							AERRCHECKINTEGER(T_i,"ibclr","i");
							AERRCHECKINTEGER(T_pos,"ibclr","pos");
							using rtype = typename arrayInfo<std::decay_t<T_i>>::type;
							auto l = [](typename arrayInfo<std::decay_t<T_i>>::type const &i,
									typename arrayInfo<std::decay_t<T_pos>>::type const &pos) -> rtype
							{  return i & ~((typename std::make_unsigned_t<typename arrayInfo<std::decay_t<T_i>>::type>(1)) << pos);};
							return makeLazy<l>(std::forward<T_i>(i),std::forward<T_pos>(pos));
						}

					/**
					 * Fortran IBITS(https://fortranwiki.org/fortran/show/ibits)
					 */
					template <typename T_i, typename T_pos, typename T_len>
						auto ibits(T_i &&i, T_pos &&pos, T_len &&len)
						{
							AERRCHECKINTEGER(T_i,"ibits","i");
							AERRCHECKINTEGER(T_pos,"ibits","pos");
							AERRCHECKINTEGER(T_len,"ibits","len");
							using rtype = typename arrayInfo<std::decay_t<T_i>>::type;
							auto l = [](typename arrayInfo<std::decay_t<T_i>>::type const &i,
									typename arrayInfo<std::decay_t<T_pos>>::type const &pos, 
									typename arrayInfo<std::decay_t<T_len>>::type const &len) -> rtype
							{  return (i&((typename std::make_unsigned_t<typename arrayInfo<std::decay_t<T_i>>::type>(1) << len) - 1) << pos)>>pos; };
							return makeLazy<l>(std::forward<T_i>(i),std::forward<T_pos>(pos),std::forward<T_len>(len));
						}

					/**
					 * Fortran IBSET(https://fortranwiki.org/fortran/show/ibset)
					 */
					template <typename T_i, typename T_pos>
						auto ibset(T_i &&i, T_pos &&pos)
						{
							AERRCHECKINTEGER(T_i,"ibset","i");
							AERRCHECKINTEGER(T_pos,"ibset","pos");
							using rtype = typename arrayInfo<std::decay_t<T_i>>::type;
							auto l = [](typename arrayInfo<std::decay_t<T_i>>::type const &i,
									typename arrayInfo<std::decay_t<T_pos>>::type const &pos) -> rtype
							{  return i | ((typename std::make_unsigned_t<typename arrayInfo<std::decay_t<T_i>>::type>(1)) << pos); };
							return makeLazy<l>(std::forward<T_i>(i),std::forward<T_pos>(pos));
						}

					/**
					 * Fortran ICHAR(https://fortranwiki.org/fortran/show/ichar)
					 */
					template <typename kind=char, typename T_c=char>
						auto ichar(T_c &&c)
						{
							AERRCHECKINTEGER(T_c,"ichar","c");
							auto l = [](typename arrayInfo<std::decay_t<T_c>>::type const &val) -> kind
							{  return static_cast<kind>(val); };
							return makeLazy<l>(std::forward<T_c>(c));
						}


					/**
					 * Fortran IEOR(https://fortranwiki.org/fortran/show/ieor)
					 */
					template <typename T1, typename T2, typename = std::enable_if<std::is_integral_v<T1> && std::is_same_v<T1,T2>>>
						auto ieor(T1 &&i, T2 &&j)
						{
							AERRCHECKINTEGER(T1,"ieor","i");
							AERRCHECKINTEGER(T2,"ieor","j");
							using rtype = typename arrayInfo<std::decay_t<T1>>::type;
							auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &i, typename arrayInfo<std::decay_t<T2>>::type const &j)
								-> rtype
								{
									return i^j;
								};
							return makeLazy<l>(std::forward<T1>(i),std::forward<T2>(j));
						}


					/**
					 * Fortran INT(https://fortranwiki.org/fortran/show/int)
					 */
					template <typename kind=int, typename T_a=int>
						auto f_int(T_a &&a)
						{
							ERRCHECKREALINTEGERCOMPLEX(T_a,"f_int","a");
							auto l = [](typename arrayInfo<std::decay_t<T_a>>::type const &a) -> kind
							{  
								if constexpr(is_complex<typename arrayInfo<T_a>::type>){
									return static_cast<kind>(a.real());
								} else {
									return static_cast<kind>(a);
								}
							};
							return makeLazy<l>(std::forward<T_a>(a));
						}


					/**
					 * Fortran IOR(https://fortranwiki.org/fortran/show/ieor)
					 */
					template <typename T1, typename T2, typename = std::enable_if<std::is_integral_v<T1> && std::is_same_v<T1,T2>>>
						auto ior(T1 &&i, T2 &&j)
						{
							AERRCHECKINTEGER(T1,"ior","i");
							AERRCHECKINTEGER(T2,"ior","j");
							using rtype = typename arrayInfo<std::decay_t<T1>>::type;
							auto l = [](typename arrayInfo<std::decay_t<T1>>::type const &i, typename arrayInfo<std::decay_t<T2>>::type const &j)
								-> rtype
								{
									return i|j;
								};
							return makeLazy<l>(std::forward<T1>(i),std::forward<T2>(j));
						}

					/**
					 * C++ iota equivalent
					 */
					template <typename T_array, typename T_value>
						void iota(T_array &&array, T_value value)
						{
							forEachSerial<([](arrayInfo<T_array>::type& item, T_value &value)
									{item = value++;})>(std::forward<T_array>(array),std::forward<T_value>(value));
						}

					/**
					 * Fortran IPARITY(https://fortranwiki.org/fortran/show/iparity)
					 */
					template <typename T_array>
						auto iparity(T_array &&array)
						{
							AERRCHECKINTEGER(T_array,"iparity","array");
							using type=arrayInfo<T_array>::type;
							return reduction<[](const type &src, type &dest){dest^=src;},
										 defaultFunction,
										 [](type &dest){dest=type(0);}
							>(std::forward<T_array>(array));

						}

					/**
					 * Fortran IPARITY(https://fortranwiki.org/fortran/show/iparity)
					 * with mask
					 */
					template <typename T_array, typename T_mask>
						auto iparity_with_mask(T_array &&array, T_mask &&mask)
						{
							AERRCHECKINTEGER(T_array,"iparity","array");
							AERRCHECKTRUTHY(T_mask,"iparity","mask");
							using type=arrayInfo<T_array>::type;
							using maskType = arrayInfo<T_mask>::type;
							return reduction<[](const type &src, type &dest, const maskType &mask){dest^=mask?src:type(0);},
										 defaultFunction,
										 [](type &dest){dest=type(0);}
							>(std::forward<T_array>(array),std::forward<T_mask>(mask));

						}

					/**
					 * Fortran IPARITY(https://fortranwiki.org/fortran/show/iparity)
					 * With dim
					 */
					template <typename T_array>
						auto iparity(T_array &&array, int dim)
						{
							AERRCHECKINTEGER(T_array,"iparity","array");
							// Produce an array with the summed direction removed
							using type=arrayInfo<T_array>::type;
							return reductionWithDirection<[](const type &src, type &dest){dest^=src;},
										 defaultFunction,
										 [](type &dest){dest=type(0);}
							>(dim, std::forward<T_array>(array));
						}

					/**
					 * Fortran IPARITY(https://fortranwiki.org/fortran/show/iparity)
					 * With dim and mask
					 */
					template <typename T_array, typename T_mask>
						auto iparity_with_mask(T_array &&array, int dim, T_mask &&mask)
						{
							AERRCHECKINTEGER(T_array,"iparity","array")
								AERRCHECKTRUTHY(T_mask,"iparity","mask")
								using type=arrayInfo<T_array>::type;
							using maskType = arrayInfo<T_mask>::type;
							return reductionWithDirection<[](const type &src, type &dest, const maskType &mask){dest^=mask?src:type(0);},
										 defaultFunction,
										 [](type &dest){dest=type(0);}
							>(dim,std::forward<T_array>(array),std::forward<T_mask>(mask));

						}

					/**
					 * Fortran ISHFT(https://fortranwiki.org/fortran/show/ishft)
					 */
					template <typename T_i, typename T_shift>
						auto ishft(T_i &&i, T_shift&& shift)
						{
							AERRCHECKINTEGER(T_i,"ishft","i");
							AERRCHECKINTEGER(T_shift,"ishft","shift");
							using rtype = typename arrayInfo<std::decay_t<T_i>>::type;
							auto l = [](typename arrayInfo<std::decay_t<T_i>>::type const &i,
									typename arrayInfo<std::decay_t<T_shift>>::type const &shift) -> rtype
							{  return shift>0?i<<shift:i>>std::abs(shift);};
							return makeLazy<l>(std::forward<T_i>(i),std::forward<T_shift>(shift));
						}


					/**
					 * Fortran ISHFTC(https://fortranwiki.org/fortran/show/ishftc)
					 */
					template <typename T_i, typename T_shift, typename T_size=FAR_SIGNED_INDEX_TYPE>
						auto ishftc(T_i &&i, T_shift&& shift, T_size size=CHAR_BIT*sizeof(typename arrayInfo<T_i>::type))
						{
							AERRCHECKINTEGER(T_i,"ishftc","i");
							AERRCHECKINTEGER(T_shift,"ishftc","shift");
							AERRCHECKINTEGER(T_size,"ishftc","size");
							using type = typename arrayInfo<T_i>::type;
							using typeShift = typename arrayInfo<T_shift>::type;
							using typeSize = typename arrayInfo<T_size>::type;
							return makeLazy<ishftc_core<type,typeShift,typeSize>>
								(std::forward<T_i>(i),std::forward<T_shift>(shift),std::forward<T_size>(size));
						}

					/**
					 * Fortran IS_CONTIGUOUS
					 */
					template <typename T_a>
						bool is_contiguous(T_a &&a)
						{
							static_assert(arrayInfo<T_a>::value,"a must be an array in is_contiguous");
							return a.getContiguous();
						}

					/**
					 * Fortran KIND(https://fortranwiki.org/fortran/show/kind)
					 * This is not the same as the Fortran version because there is no equivalent
					 * Instead it is a using statement that gives you a type matching
					 * either the underlying type in an array or the type of the item called with
					 */
					template<typename T_i>
						using kind = arrayInfo<T_i>::type;

					/**
					 * Fortran LBOUND(https://fortranwiki.org/fortran/show/lbound)
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T_i=double>
						auto lbound(T_i &&array){
#ifdef RANK1_INQUIRY_IS_VALUE
							if constexpr(arrayInfo<T_i>::rank==1) {
								kind result = array.getLB(1);
								return result;
							}
							else {
#endif
								constexpr int rank = arrayInfo<T_i>::rank;
								Array<kind,1> result(rank);
								for (int i=1;i<=rank;++i){
									result(i) = array.getLB(i);
								}
								return result;
#ifdef RANK1_INQUIRY_IS_VALUE
							}
#endif
						}

					/**
					 * Fortran LBOUND(https://fortranwiki.org/fortran/show/lbound)
					 * with dim
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T_i=double>
						auto lbound(T_i &&array, int dim){
							kind result = array.getLB(dim);
							return result;
						}

					/**
					 * Fortran LEADZ(https://fortranwiki.org/fortran/show/leadz)
					 */
					template <typename T_i>
						auto leadz(T_i &&i)
						{
							AERRCHECKINTEGER(T_i,"leadz","i");
							auto l = [](typename arrayInfo<std::decay_t<T_i>>::type const &i)-> decltype(auto)
							{
								return std::countl_zero(std::make_unsigned_t<std::decay_t<decltype(i)>>(i));
							};
							return makeLazy<l>(std::forward<T_i>(i));
						}

					/**
					 * Fortran LEN(https://fortranwiki.org/fortran/show/len)
					 * This version is for fixed size character arrays
					 * and returns that maximum size of the string
					 */
					template <typename kind=FAR_UNSIGNED_INDEX_TYPE, typename T_string=char[1], 
						typename = std::enable_if_t<isFixedSizeCharArray<typename arrayInfo<std::remove_reference_t<T_string>>::unrefType>::value>>
							kind len(FAR_UNUSED T_string &&string)
							{
								return isFixedSizeCharArray<typename arrayInfo<std::remove_reference_t<T_string>>::unrefType>::size;
							}

					/**
					 * Fortran LEN(https://fortranwiki.org/fortran/show/len)
					 * This version is for std::string and returns an array of the same rank and shape
					 * as the passed array with the length of the string in each element
					 */
					template <typename kind = FAR_UNSIGNED_INDEX_TYPE, typename T_string = std::string, 
						std::enable_if_t<std::is_same_v<std::string, typename arrayInfo<T_string>::type>
							|| std::is_same_v<const std::string, typename arrayInfo<T_string>::type>
							, int> = 0>
							auto len(T_string &&string)
							{
								auto l = [](typename arrayInfo<std::decay_t<T_string>>::type const &string)-> kind
								{
									return string.length();
								};
								return makeLazy<l>(std::forward<T_string>(string));

							}

					/**
					 * Fortran LEN(https://fortranwiki.org/fortran/show/len)
					 * This version is for char* and returns an array of the same rank and shape
					 * as the passed array with the length of the string in each element
					 */
					template <typename kind = FAR_UNSIGNED_INDEX_TYPE, typename T_string = char*,  
						std::enable_if_t<
							!isFixedSizeCharArray<typename arrayInfo<std::remove_reference_t<T_string>>::unrefType>::value &&(
									std::is_same_v<typename arrayInfo<std::decay_t<T_string>>::type,char*>
									|| std::is_same_v<typename arrayInfo<std::decay_t<T_string>>::type,const char*>)
							, int> = 0>
							auto len(T_string &&string)
							{
								auto l = [](typename arrayInfo<std::decay_t<T_string>>::type const &string)-> kind
								{
									return strlen(string);
								};
								return makeLazy<l>(std::forward<T_string>(string));

							}


					/**
					 * Fortran LEN_TRIM(https://fortranwiki.org/fortran/show/len_trim)
					 * This version is for fixed length character strings
					 * In Fortran this returns a string with trailing spaces removed, but in C/C++
					 * trailing spaces are probably wanted since null termination is used
					 * This returns the length to the null terminator
					 */
					template <typename kind=FAR_UNSIGNED_INDEX_TYPE, typename T_string=char*,
						typename = std::enable_if_t<isFixedSizeCharArray<typename arrayInfo<std::remove_reference_t<T_string>>::type>::value>>
							auto len_trim(T_string &&string)
							{
								auto l = [](char* const string)-> kind
								{
									const char * src = string;
									kind last_non_space = 0, index = 0;
									while(*src){
										if (!std::isspace((int)*src))last_non_space=index;
										src++;
										index++;
									}
									return last_non_space+1;
								};
								return makeLazy<l>(std::forward<T_string>(string));
							}

					/**
					 * Fortran LEN_TRIM(https://fortranwiki.org/fortran/show/len_trim)
					 * This version is for character pointers
					 * In Fortran this returns a string with trailing spaces removed, but in C/C++
					 * trailing spaces are probably wanted since null termination is used
					 * This returns the length to the null terminator
					 */
					template <typename kind = FAR_UNSIGNED_INDEX_TYPE, typename T_string = char*,  
						std::enable_if_t<
							!isFixedSizeCharArray<typename arrayInfo<std::remove_reference_t<T_string>>::type>::value &&(
									std::is_same_v<typename arrayInfo<std::decay_t<T_string>>::type,char*>
									|| std::is_same_v<typename arrayInfo<std::decay_t<T_string>>::type,const char*>)
							, int> = 0>
							auto len_trim(T_string &&string)
							{
								auto l = [](typename arrayInfo<std::decay_t<T_string>>::type const &string)-> kind
								{ const char * src = string;
									kind last_non_space = 0, index = 0;
									while(*src){
										if (!std::isspace((int)*src))last_non_space=index;
										src++;
										index++;
									}
									return last_non_space+1;
								};
								return makeLazy<l>(std::forward<T_string>(string));
							}

					/**
					 * Fortran LEN_TRIM(https://fortranwiki.org/fortran/show/len_trim)
					 * This version is for std::string
					 * In Fortran this returns a string with trailing spaces removed, but in C/C++
					 * trailing spaces are probably wanted since null termination is used
					 * This returns the length returned by std::string::length
					 */
					template <typename kind=FAR_UNSIGNED_INDEX_TYPE, typename T_string=std::string,
						std::enable_if_t<std::is_same_v<typename arrayInfo<std::decay_t<T_string>>::type,std::string> 
							|| std::is_same_v<typename arrayInfo<std::decay_t<T_string>>::type,std::string>
							,int> =0>
							auto len_trim(T_string &&string)
							{
								auto l = [](typename arrayInfo<std::decay_t<T_string>>::type const &string)-> kind
								{
									const char * src = string.c_str();
									kind last_non_space = 0, index = 0;
									while(*src){
										if (!std::isspace((int)*src))last_non_space=index;
										src++;
										index++;
									}	
									return last_non_space+1;
								};
								return makeLazy<l>(std::forward<T_string>(string));
							}

					/**
					 * Fortran LGE(https://fortranwiki.org/fortran/show/lge)
					 */
					template <typename T_string1, typename T_string2>
						auto lge(T_string1 &&string1, T_string2 &&string2)
						{
							using T_inner1 = typename arrayInfo<T_string1>::type;
							using T_inner2 = typename arrayInfo<T_string2>::type;
							AERRCHECKSTRINGY(T_string1,"lge","string1");
							AERRCHECKSTRINGY(T_string2,"lge","string2");
							auto l = [](T_inner1 const &s1, 
									T_inner2 const &s2)-> decltype(auto)
							{
								const char *c1, *c2;
								FAR_SIGNED_INDEX_TYPE sz1, sz2;
								if constexpr(std::is_same_v<std::decay_t<decltype(s1)>,std::string>){
									c1=s1.c_str();
									sz1=s1.length();
								} else {
									c1=s1;
									sz1=strlen(c1);
								}
								if constexpr(std::is_same_v<std::decay_t<decltype(s2)>,std::string>){
									c2=s2.c_str();
									sz2=s2.length();
								} else {
									c2=s2;
									sz2=strlen(c2);
								}
								FAR_SIGNED_INDEX_TYPE mn = std::min(sz1,sz2);
								for (FAR_SIGNED_INDEX_TYPE i=0;i<mn;++i){
									//Convert from local collating sequence to ascii (likely does nothing, but that
									//guarantee is what separates this from <=
									char ac1 = char_to_ascii[c1[i]];
									char ac2 = char_to_ascii[c2[i]];
									if (ac1!=ac2) return ac1>=ac2;
								}
								return true;
							};
							return makeLazy<l>(std::forward<T_string1>(string1),std::forward<T_string2>(string2));
						}

					/**
					 * Fortran LGT(https://fortranwiki.org/fortran/show/lgt)
					 */
					template <typename T_string1, typename T_string2>
						auto lgt(T_string1 &&string1, T_string2 &&string2)
						{
							AERRCHECKSTRINGY(T_string1,"lgt","string1");
							AERRCHECKSTRINGY(T_string2,"lgt","string2");	
							auto l = [](typename arrayInfo<std::decay_t<T_string1>>::type const &s1,
									typename arrayInfo<std::decay_t<T_string2>>::type const &s2)-> decltype(auto)
							{
								const char *c1, *c2;
								FAR_SIGNED_INDEX_TYPE sz1, sz2;
								if constexpr(std::is_same_v<std::decay_t<decltype(s1)>,std::string>){
									c1=s1.c_str();
									sz1=s1.length();
								} else {
									c1=s1;
									sz1=strlen(c1);
								}
								if constexpr(std::is_same_v<std::decay_t<decltype(s2)>,std::string>){
									c2=s2.c_str();
									sz2=s2.length();
								} else {
									c2=s2;
									sz2=strlen(s2);
								}
								FAR_SIGNED_INDEX_TYPE mn = std::min(sz1,sz2);
								for (FAR_SIGNED_INDEX_TYPE i=0;i<mn;++i){
									//Convert from local collating sequence to ascii (likely does nothing, but that
									//guarantee is what separates this from <=
									char ac1 = char_to_ascii[c1[i]];
									char ac2 = char_to_ascii[c2[i]];
									if (ac1!=ac2) return ac1>ac2;
								}
								return false;
							};
							return makeLazy<l>(std::forward<T_string1>(string1),std::forward<T_string2>(string2));
						}

					/**
					 * Fortran LLE(https://fortranwiki.org/fortran/show/lle)
					 */
					template <typename T_string1, typename T_string2>
						auto lle(T_string1 &&string1, T_string2 &&string2)
						{
							AERRCHECKSTRINGY(T_string1,"lle","string1");
							AERRCHECKSTRINGY(T_string2,"lle","string2");
							auto l = [](typename arrayInfo<std::decay_t<T_string1>>::type const &s1,
									typename arrayInfo<std::decay_t<T_string2>>::type const &s2)-> decltype(auto)
							{
								const char *c1, *c2;
								FAR_SIGNED_INDEX_TYPE sz1, sz2;
								if constexpr(std::is_same_v<std::decay_t<decltype(s1)>,std::string>){
									c1=s1.c_str();
									sz1=s1.length();
								} else {
									c1=s1;
									sz1=strlen(c1);
								}
								if constexpr(std::is_same_v<std::decay_t<decltype(s2)>,std::string>){
									c2=s2.c_str();
									sz2=s2.length();
								} else {
									c2=s2;
									sz2=strlen(s2);
								}
								FAR_SIGNED_INDEX_TYPE mn = std::min(sz1,sz2);
								for (FAR_SIGNED_INDEX_TYPE i=0;i<mn;++i){
									//Convert from local collating sequence to ascii (likely does nothing, but that
									//guarantee is what separates this from <=
									char ac1 = char_to_ascii[c1[i]];
									char ac2 = char_to_ascii[c2[i]];
									if (ac1!=ac2) return ac1<=ac2;
								}
								return true;
							};
							return makeLazy<l>(std::forward<T_string1>(string1),std::forward<T_string2>(string2));
						}

					/**
					 * Fortran LLT(https://fortranwiki.org/fortran/show/llt)
					 */
					template <typename T_string1, typename T_string2>
						auto llt(T_string1 &&string1, T_string2 &&string2)
						{
							AERRCHECKSTRINGY(T_string1,"llt","string1");
							AERRCHECKSTRINGY(T_string2,"llt","string2");
							auto l = [](typename arrayInfo<std::decay_t<T_string1>>::type const &s1,
									typename arrayInfo<std::decay_t<T_string2>>::type const &s2)-> decltype(auto)
							{
								const char *c1, *c2;
								FAR_SIGNED_INDEX_TYPE sz1, sz2;
								if constexpr(std::is_same_v<std::decay_t<decltype(s1)>,std::string>){
									c1=s1.c_str();
									sz1=s1.length();
								} else {
									c1=s1;
									sz1=strlen(c1);
								}
								if constexpr(std::is_same_v<std::decay_t<decltype(s2)>,std::string>){
									c2=s2.c_str();
									sz2=s2.length();
								} else {
									c2=s2;
									sz2=strlen(s2);
								}
								FAR_SIGNED_INDEX_TYPE mn = std::min(sz1,sz2);
								for (FAR_SIGNED_INDEX_TYPE i=0;i<mn;++i){
									//Convert from local collating sequence to ascii (likely does nothing, but that
									//guarantee is what separates this from <=
									char ac1 = char_to_ascii[c1[i]];
									char ac2 = char_to_ascii[c2[i]];
									if (ac1!=ac2) return ac1<ac2;
								}
								return false;
							};
							return makeLazy<l>(std::forward<T_string1>(string1),std::forward<T_string2>(string2));
						}

					/**
					 * Fortran LOG(https://fortranwiki.org/fortran/show/log)
					 */
					template <typename T>
						auto log(T &&x)
						{
							AERRCHECKREAL(T,"log","x");
							auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
							{  return std::log(x); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran LOG_GAMMA(https://fortranwiki.org/fortran/show/log_gamma)
					 */
					template <typename T>
						auto log_gamma(T &&x)
						{
							AERRCHECKREAL(T,"log_gamma","x");
							auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
							{  return std::lgamma(x); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran LOG10(https://fortranwiki.org/fortran/show/log10)
					 */
					template <typename T>
						auto log10(T &&x)
						{
							AERRCHECKREAL(T,"log10","x");
							auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
							{  return std::log10(x); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran doesn't have LOG2, but C++ does so I have included it
					 */
					template <typename T>
						auto log2(T &&x)
						{
							AERRCHECKREAL(T,"log2","x");
							auto l = [](typename arrayInfo<std::decay_t<T>>::type const &x) -> decltype(auto)
							{  return std::log2(x); };
							return makeLazy<l>(std::forward<T>(x));
						}

					/**
					 * Fortran LOGICAL(https://fortranwiki.org/fortran/show/logical)
					 * C++ doesn't have different bool kinds, but it CAN use an integer
					 * as a logical so this function is designed for that.
					 */
					template <typename kind=bool, typename T=int>
						auto logical(T &&l)
						{
							auto f = [](typename arrayInfo<std::decay_t<T>>::type const &l) -> kind
							{  return static_cast<kind>(l); };
							return makeLazy<f>(std::forward<T>(l));
						}


					template<typename T>
						auto nullify(T &&l)
						{
							l.pointTo();
						}


					/**
					 * Fortran MASKL(https://fortranwiki.org/fortran/show/maskl)
					 */
					template <typename kind=int, typename T_i=int>
						auto maskl(T_i &&i)
						{
							AERRCHECKINTEGER(T_i,"maskl","i");
							using ukind = std::make_unsigned_t<kind>;
							using type = typename arrayInfo<T_i>::type;
							auto l = [](type const &i) -> decltype(auto)
							{
								return kind(i <= 0 ? 0 : (i >= type(CHAR_BIT * sizeof(ukind)) ? ~(ukind)0 : ~((ukind(1) << (CHAR_BIT*sizeof(ukind) - i)) - 1)));
							};
							return makeLazy<l>(std::forward<T_i>(i));
						}

					/**
					 * Fortran MASKR(https://fortranwiki.org/fortran/show/maskr)
					 */
					template <typename kind=int, typename T_i=int>
						auto maskr(T_i &&i)
						{
							AERRCHECKINTEGER(T_i,"maskr","i");
							using ukind = std::make_unsigned_t<kind>;
							using type = typename arrayInfo<T_i>::type;
							auto l = [](type const &i) -> decltype(auto)
							{
								return kind(i <= 0 ? 0 : (i >= type(CHAR_BIT * sizeof(ukind)) ? ~(ukind)0 : ((ukind(1) << i) - 1)));
							};
							return makeLazy<l>(std::forward<T_i>(i));
						}


					/**
					 * Fortran MATMUL(https://fortranwiki.org/fortran/show/matmul)
					 */
					template<typename T_matrix_a, typename T_matrix_b>
						auto matmul(T_matrix_a && matrix_a, T_matrix_b && matrix_b){
							using indexer = typename arrayInfo<std::decay_t<T_matrix_a>>::index;
							using core = typename arrayInfo<std::decay_t<T_matrix_a>>::type;
							if constexpr(std::is_same_v<core,bool>){
								return indexer::template matmul<[](const core& A, const core& B, core &C){C = C || (A&&B);}>
									(std::forward<T_matrix_a>(matrix_a), std::forward<T_matrix_b>(matrix_b));
							} else {
								return indexer::template matmul<[](const core& A, const core& B, core &C){C = C + (A*B);}>
									(std::forward<T_matrix_a>(matrix_a), std::forward<T_matrix_b>(matrix_b));
							}
						}


					/**
					 * Fortran MAX(https://fortranwiki.org/fortran/show/max)
					 */
					template<typename T_first, typename... T_others>
						auto max(T_first && first, T_others&&... others){
							return makeLazy<pmax<typename arrayInfo<T_first>::type, typename arrayInfo<T_others>::type...>>
								(std::forward<T_first>(first),std::forward<T_others>(others)...);
						}


					/**
					 * Fortran MAXEXPONENT(https://fortranwiki.org/fortran/show/maxexponent)
					 */
					template<typename T>
						auto maxexponent(FAR_UNUSED T && item){
							return std::numeric_limits<typename arrayInfo<T>::type>::max_exponent;
						}

					/**
					 * Fortran MAXLOC(https://fortranwiki.org/fortran/show/maxloc)
					 * Array, kind and back only
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T=void, bool natural=false>
						auto maxloc(T &&array, bool back=false){
							ERRCHECKREALINTEGERSTRING(T,"maxloc","array");
							using coreType = typename arrayInfo<T>::type;
							using locType = typename N_ary_tuple_type<kind,arrayInfo<T>::rank>::type;
							using rrType = indexInfo<arrayInfo<T>::rank>;
							struct rop{
								coreType accum;
								locType loc;
								bool active=false;
							};
							auto parallelCombine = [](rrType &rri, const T &src, rop &dest){
								auto val = getItem(src,rri);
								if (!dest.active || (val>dest.accum)){
									dest.active=true;
									dest.accum=val;
									if constexpr(natural){
										dest.loc = rri.SBindices;
									} else {
										dest.loc = rri.indices;
									}
								}
							};
							auto init = [](rop &dest){
								dest.active=false;
								assignTuple(dest.loc,defaultLB-1);
							};
							if (!back){
								auto serialCombine = [](const rop &local, rop &global){
									if (local.active &&((!global.active || indexTupleLT(local.loc,global.loc)) && (local.accum==global.accum ||
													local.accum>global.accum))){
										global.accum=local.accum;
										global.loc=local.loc;
										global.active=true;
									}
								};
								auto r = reduction<parallelCombine,serialCombine, init, nullFunction>(std::forward<T>(array));
								Array<kind,1> loc(arrayInfo<T>::rank);
								assignArrayTuple(loc,r.loc);
								return loc;
							} else {
								auto serialCombine = [](const rop &local, rop &global){
									if (local.active &&((!global.active || indexTupleGT(local.loc,global.loc)) && (local.accum==global.accum ||
													local.accum>global.accum))){
										global.accum=local.accum;
										global.loc=local.loc;
										global.active=true;
									}
								};
								auto r = reduction<parallelCombine,serialCombine, init, nullFunction,true>(std::forward<T>(array));
								Array<kind,1> loc(arrayInfo<T>::rank);
								assignArrayTuple(loc,r.loc);
								return loc;
							}
						}


					/**
					 * Fortran MAXLOC(https://fortranwiki.org/fortran/show/maxloc)
					 * Array, mask, kind and back
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T=void, typename T_mask=void, bool natural=false>
						auto maxloc_with_mask(T &&array, T_mask &&mask, bool back = false){
							ERRCHECKREALINTEGERSTRING(T,"maxloc_with_mask","array");
							AERRCHECKTRUTHY(T_mask,"maxloc_with_mask","mask");
							using coreType = typename arrayInfo<T>::type;
							using locType = typename N_ary_tuple_type<kind,arrayInfo<T>::rank>::type;
							using rrType = indexInfo<arrayInfo<T>::rank>;
							struct rop{
								coreType accum;
								locType loc;
								bool active;
							};
							auto parallelCombine = [](rrType &rri, const T &src, rop &dest, const T_mask &mask){
								auto val = getItem(src,rri);
								if (getItem(mask,rri) && (!dest.active || (val>dest.accum))){
									dest.active=true;
									dest.accum=val;
									if constexpr(natural){
										dest.loc = rri.SBindices;
									} else {
										dest.loc = rri.indices;
									}
								}
							};
							auto init = [](rop &dest){
								dest.active=false;
								assignTuple(dest.loc,defaultLB-1);
							};
							if (!back){
								auto serialCombine = [](const rop &local, rop &global){
									if (local.active &&((!global.active || indexTupleLT(local.loc,global.loc)) && (local.accum==global.accum ||
													local.accum>global.accum))){
										global.accum=local.accum;
										global.loc=local.loc;
										global.active=true;
									}
								};
								auto r = reduction<parallelCombine,serialCombine, init, nullFunction>(std::forward<T>(array),std::forward<T_mask>(mask));
								Array<kind,1> loc(arrayInfo<T>::rank);
								assignArrayTuple(loc,r.loc);
								return loc;
							} else {
								auto serialCombine = [](const rop &local, rop &global){
									if (local.active &&((!global.active || indexTupleGT(local.loc,global.loc)) && (local.accum==global.accum ||
													local.accum>global.accum))){
										global.accum=local.accum;
										global.loc=local.loc;
										global.active=true;
									}
								};
								auto r = reduction<parallelCombine,serialCombine, init, nullFunction,true>(std::forward<T>(array),std::forward<T_mask>(mask));
								Array<kind,1> loc(arrayInfo<T>::rank);
								assignArrayTuple(loc,r.loc);
								return loc;
							}
						}


					/**
					 * Fortran MAXLOC(https://fortranwiki.org/fortran/show/maxloc)
					 * Array, kind and back
					 * with direction
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T=void, bool natural=false>
						auto maxloc(T &&array, int direction, bool back=false){
							ERRCHECKREALINTEGERSTRING(T,"maxloc","array");
							constexpr int rank = arrayInfo<T>::rank;
							using coreType = typename arrayInfo<T>::type;
							using reduction = indexInfo<rank>;
							using rtype = Array<kind, rank-1>;
							struct meta{
								bool active;
								coreType accum;
							};
							using marray = Array<meta,rank-1>;

							auto parallelCombine = [](const reduction& rri, const T& array, rtype &result, marray &ma){
								auto val = getItem(array,rri);
								meta &m = getItem(ma,rri);
								if (!m.active || (val>m.accum)){
									m.active=true;
									m.accum=val;
									getItem(result,rri)=getTupleLevel(rri.SBindices,rri.direction);
								}
							};

							auto initResult = [](kind &result){
								result = defaultLB-1;
							};

							auto initMeta = [](meta &m){
								m.active=false;
							};

							if (!back){
								auto serialCombine = [](const kind &local, const meta &lm, kind &global, meta &gm){
									if ((!gm.active || lm.active) && ((lm.accum>gm.accum || lm.accum == gm.accum) && local<global)){
										gm.active=true;
										gm.accum=lm.accum;
										global=local;
									}
								};
								auto r = reductionWithDirectionAndMeta<parallelCombine, serialCombine, initResult, initMeta, nullFunction>(direction, std::forward<T>(array));
								return r;
							} else {
								auto serialCombine = [](const kind &local, const meta &lm, kind &global, meta &gm){
									if (!(gm.active || lm.active) && ((lm.accum>gm.accum || lm.accum == gm.accum) && local>global)){
										gm.active=true;
										gm.accum=lm.accum;
										global=local;
									}
								};
								auto r = reductionWithDirectionAndMeta<parallelCombine, serialCombine, initResult, initMeta, nullFunction,true>(direction, std::forward<T>(array));
								return r;
							}

						}

					/**
					 * Fortran MAXLOC(https://fortranwiki.org/fortran/show/maxloc)
					 * Array, mask, kind and back
					 * with direction
					 */

					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T=void, typename T_mask=void, bool natural=false>
						auto maxloc_with_mask(T &&array, int direction, T_mask &&mask, bool back=false){
							ERRCHECKREALINTEGERSTRING(T,"maxloc_with_mask","array");
							AERRCHECKTRUTHY(T_mask,"maxloc_with_mask","mask");
							constexpr int rank = arrayInfo<T>::rank;
							using coreType = typename arrayInfo<T>::type;
							using reduction = indexInfo<rank>;
							using rtype = Array<kind, rank-1>;
							struct meta{
								bool active;
								coreType accum;
							};
							using marray = Array<meta,rank-1>;

							auto parallelCombine = [](const reduction& rri, const T& array, rtype &result, marray &ma, const T_mask &mask){
								auto val = getItem(array,rri);
								meta &m = getItem(ma,rri);
								auto mval = getItem(mask,rri);
								if (mval && (!m.active || (val>m.accum))){
									m.active=true;
									m.accum=val;
									getItem(result,rri)=getTupleLevel(rri.SBindices,rri.direction);
								}
							};

							auto initResult = [](kind &result){
								result = defaultLB-1;
							};

							auto initMeta = [](meta &m){
								m.active=false;
							};

							if (!back){
								auto serialCombine = [](const kind &local, const meta &lm, kind &global, meta &gm){
									if ((!gm.active || lm.active) && ((lm.accum>gm.accum || lm.accum == gm.accum) && local<global)){
										gm.active=true;
										gm.accum=lm.accum;
										global=local;
									}
								};
								auto r = reductionWithDirectionAndMeta<parallelCombine, serialCombine, initResult, initMeta, nullFunction>(direction, std::forward<T>(array), std::forward<T_mask>(mask));
								return r;
							} else {
								auto serialCombine = [](const kind &local, const meta &lm, kind &global, meta &gm){
									if ((!gm.active || lm.active) && ((lm.accum>gm.accum || lm.accum == gm.accum) && local>global)){
										gm.active=true;
										gm.accum=lm.accum;
										global=local;
									}
								};
								auto r = reductionWithDirectionAndMeta<parallelCombine, serialCombine, initResult, initMeta, nullFunction,true>(direction, std::forward<T>(array), std::forward<T_mask>(mask));
								return r;
							}
						}

					/**
					 * Fortran MAXVAL(https://fortranwiki.org/fortran/show/maxval)
					 */
					template <typename T>
						auto maxval(T &&array)
						{
							ERRCHECKREALINTEGERSTRING(T,"maxval","array");
							using coreType = typename arrayInfo<T>::type;
							struct rop{
								coreType accum;
								bool active=false;
							};
							auto parallelCombine = [](const coreType &src, rop &dest){
								if (!dest.active || (src>dest.accum)){
									dest.active=true;
									dest.accum=src;
								}
							};
							auto serialCombine = [](const rop &local, rop &global){
								if (local.active && (local.accum>global.accum || !global.active)){
									global.accum=local.accum;
									global.active=true;
								}
							};
							auto initial = [](rop &dest){
								dest.active=false;
								//This is how Fortran does it, but in C++ can only do this for arithmetic types
								//For non-arithmetic types, the default constructor is used, so up to the user to
								//determine that a default constructed object is detected as "no maximum"
								if constexpr(std::is_arithmetic_v<coreType>){
									dest.accum = std::numeric_limits<coreType>::lowest();
								}
							};
							auto r = reduction<parallelCombine,serialCombine, initial, defaultFunction>(std::forward<T>(array));
							return r.accum;
						}

					/**
					 * Fortran MAXVAL(https://fortranwiki.org/fortran/show/maxval)
					 * with mask
					 */
					template <typename T, typename T_mask>
						auto maxval_with_mask(T &&array, T_mask &&mask)
						{
							ERRCHECKREALINTEGERSTRING(T,"maxval_with_mask","array");
							AERRCHECKTRUTHY(T_mask,"maxval_with_mask","mask");
							using coreType = typename arrayInfo<T>::type;
							using maskType = typename arrayInfo<T_mask>::type;
							struct rop{
								coreType accum;
								bool active=false;
							};
							auto parallelCombine = [](const coreType &src, rop &dest, const maskType &mask){
								if (mask &&(!dest.active || (src>dest.accum))){
									dest.active=true;
									dest.accum=src;
								}
							};
							auto serialCombine = [](const rop &local, rop &global){
								if (local.active && (local.accum>global.accum || !global.active)){
									global.accum=local.accum;
									global.active=true;
								}
							};
							auto initial = [](rop &dest){
								dest.active=false;
								//This is how Fortran does it, but in C++ can only do this for arithmetic types
								//For non-arithmetic types, the default constructor is used, so up to the user to
								//determine that a default constructed object is detected as "no maximum"
								if constexpr(std::is_arithmetic_v<coreType>){
									dest.accum = std::numeric_limits<coreType>::lowest();
								}
							};
							auto r = reduction<parallelCombine,serialCombine, initial, defaultFunction>(std::forward<T>(array),std::forward<T_mask>(mask));
							return r.accum;
						}

					/**
					 * Fortran MAXVAL(https://fortranwiki.org/fortran/show/maxval)
					 * With dim
					 */
					template <typename T>
						auto maxval(T &&array, int dim)
						{
							ERRCHECKREALINTEGERSTRING(T,"maxval","array");
							using coreType = typename arrayInfo<T>::type;
							struct meta{
								bool active;
							};

							auto parallelCombine = [](const coreType& val, coreType &result, meta &m){
								if (!m.active || (val>result)){
									m.active=true;
									result=val;
								}
							};

							auto initResult = [](coreType &result){
								if constexpr(std::is_arithmetic_v<coreType>){
									result = std::numeric_limits<coreType>::lowest();
								}

							};

							auto initMeta = [](meta &m){
								m.active=false;
							};

							auto serialCombine = [](const coreType &local, const meta &lm, coreType &global, meta &gm){
								if (lm.active && (local>global || !gm.active)){
									gm.active=true;
									global=local;
								}
							};
							auto r = reductionWithDirectionAndMeta<parallelCombine, serialCombine, initResult, initMeta>(dim, std::forward<T>(array));
							return r;
						}

					/**
					 * Fortran MAXVAL(https://fortranwiki.org/fortran/show/maxval)
					 * With dim
					 */
					template <typename T, typename T_mask>
						auto maxval_with_mask(T &&array, int dim, T_mask &&mask)
						{
							ERRCHECKREALINTEGERSTRING(T,"maxval","array");
							ERRCHECKTRUTHY(T_mask,"maxval","mask");
							using coreType = typename arrayInfo<T>::type;
							using maskType = typename arrayInfo<T_mask>::type;
							struct meta{
								bool active;
							};

							auto parallelCombine = [](const coreType& val, coreType &result, meta &m, const maskType &mask){
								if (mask && (!m.active || (val>result))){
									m.active=true;
									result=val;
								}
							};

							auto initResult = [](coreType &result){
								if constexpr(std::is_arithmetic_v<coreType>){
									result = std::numeric_limits<coreType>::lowest();
								}

							};

							auto initMeta = [](meta &m){
								m.active=false;
							};

							auto serialCombine = [](const coreType &local, const meta &lm, coreType &global, meta &gm){
								if (lm.active && (local>global || !gm.active)){
									gm.active=true;
									global=local;
								}
							};
							auto r = reductionWithDirectionAndMeta<parallelCombine, serialCombine, initResult, initMeta>(dim, std::forward<T>(array), std::forward<T_mask>(mask));
							return r;
						}

					/**
					 * Fortran MERGE(https://fortranwiki.org/fortran/show/merge)
					 */
					template<typename T_tsource, typename T_fsource, typename T_mask>
						auto merge(T_tsource &&tsource, T_fsource &&fsource, T_mask &&mask){
							AERRCHECKTRUTHY(T_mask,"merge","mask");
							using ixtype = const std::remove_pointer_t<decltype(getIndexer(tsource,fsource,mask))>;
							//using idtype = N_ary_tuple_type_t<FAR_SIGNED_INDEX_TYPE,maxRank<T_tsource,T_fsource,T_mask>>;
							constexpr int rank = maxRank<T_tsource,T_fsource,T_mask>;
							using rtype = std::conditional_t<!std::is_void_v<ixtype>,indexInfo<rank>,int>;
							if constexpr (std::is_void_v<ixtype>){
								return mask?tsource:fsource;
							} else {
								auto l = [](const rtype &rt, const T_tsource &tsource, const T_fsource &fsource, const T_mask &mask){
									return callDefaultBase(mask,rt.indices)?callDefaultBase(tsource,rt.indices):callDefaultBase(fsource,rt.indices);
								};
								//return Array<int,1>(2);
								return makeLazy<l>(std::forward<T_tsource>(tsource),std::forward<T_fsource>(fsource),std::forward<T_mask>(mask));
							}
						}

					/**
					 * Fortran MERGE_BITS(https://fortranwiki.org/fortran/show/merge_bits)
					 */
					template<typename T_i, typename T_j, typename T_mask>
						auto merge_bits(T_i &&i, T_j &&j, T_mask &&mask){
							AERRCHECKINTEGER(T_i,"merge_bits","i");
							AERRCHECKINTEGER(T_i,"merge_bits","i");
							auto l = [](const arrayInfo<T_i>::type&i, const arrayInfo<T_j>::type&j, const arrayInfo<T_mask>::type&mask){
								return (i & mask) | (j & ~mask);
							};
							return makeLazy<l>(std::forward<T_i>(i),std::forward<T_j>(j),std::forward<T_mask>(mask));
						}

					/**
					 * Fortran MIN(https://fortranwiki.org/fortran/show/min)
					 */
					template<typename T_first, typename... T_others>
						auto min(T_first && first, T_others&&... others){
							return makeLazy<pmin<typename arrayInfo<T_first>::type, typename arrayInfo<T_others>::type...>>
								(std::forward<T_first>(first),std::forward<T_others>(others)...);
						}

					/**
					 * Fortran MINEXPONENT(https://fortranwiki.org/fortran/show/minexponent)
					 */
					template<typename T>
						auto minexponent(FAR_UNUSED T && item){
							return std::numeric_limits<typename arrayInfo<T>::type>::min_exponent;
						}

					/**
					 * Fortran MINLOC(https://fortranwiki.org/fortran/show/minloc)
					 * Array, kind and back only
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T=void, bool natural=false>
						auto minloc(T &&array, bool back=false){
							ERRCHECKREALINTEGERSTRING(T,"minloc","array");
							using coreType = typename arrayInfo<T>::type;
							using locType = typename N_ary_tuple_type<kind,arrayInfo<T>::rank>::type;
							using rrType = indexInfo<arrayInfo<T>::rank>;
							struct rop{
								coreType accum;
								locType loc;
								bool active=false;
							};
							auto parallelCombine = [](rrType &rri, const T &src, rop &dest){
								auto val = getItem(src,rri);
								if (!dest.active || (val<dest.accum)){
									dest.active=true;
									dest.accum=val;
									if constexpr(natural){
										dest.loc = rri.SBindices;
									} else {
										dest.loc = rri.indices;
									}
								}
							};
							auto init = [](rop &dest){
								dest.active=false;
								assignTuple(dest.loc,defaultLB-1);
							};
							if (!back){
								auto serialCombine = [](const rop &local, rop &global){
									if (local.active &&((!global.active || indexTupleLT(local.loc,global.loc)) && (local.accum==global.accum ||
													local.accum<global.accum))){
										global.accum=local.accum;
										global.loc=local.loc;
										global.active=true;
									}
								};
								auto r = reduction<parallelCombine,serialCombine, init, nullFunction>(std::forward<T>(array));
								Array<kind,1> loc(arrayInfo<T>::rank);
								assignArrayTuple(loc,r.loc);
								return loc;
							} else {
								auto serialCombine = [](const rop &local, rop &global){
									if (local.active &&((!global.active || indexTupleGT(local.loc,global.loc)) && (local.accum==global.accum ||
													local.accum<global.accum))){
										global.accum=local.accum;
										global.loc=local.loc;
										global.active=true;
									}
								};
								auto r = reduction<parallelCombine,serialCombine, init, nullFunction,true>(std::forward<T>(array));
								Array<kind,1> loc(arrayInfo<T>::rank);
								assignArrayTuple(loc,r.loc);
								return loc;
							}
						}

					/**
					 * Fortran MINLOC(https://fortranwiki.org/fortran/show/minloc)
					 * Array, mask, kind and back
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T=void, typename T_mask=void, bool natural=false>
						auto minloc_with_mask(T &&array, T_mask &&mask, bool back = false){
							ERRCHECKREALINTEGERSTRING(T,"minloc_with_mask","array");
							AERRCHECKTRUTHY(T_mask,"minloc_with_mask","mask");
							using coreType = typename arrayInfo<T>::type;
							using locType = typename N_ary_tuple_type<kind,arrayInfo<T>::rank>::type;
							using rrType = indexInfo<arrayInfo<T>::rank>;
							struct rop{
								coreType accum;
								locType loc;
								bool active;
							};
							auto parallelCombine = [](rrType &rri, const T &src, rop &dest, const T_mask &mask){
								auto val = getItem(src,rri);
								if (getItem(mask,rri) && (!dest.active || (val<dest.accum))){
									dest.active=true;
									dest.accum=val;
									if constexpr(natural){
										dest.loc = rri.SBindices;
									} else {
										dest.loc = rri.indices;
									}
								}
							};
							auto init = [](rop &dest){
								dest.active=false;
								assignTuple(dest.loc,defaultLB-1);
							};
							if (!back){
								auto serialCombine = [](const rop &local, rop &global){
									if (local.active &&((!global.active || indexTupleLT(local.loc,global.loc)) && (local.accum==global.accum ||
													local.accum<global.accum))){
										global.accum=local.accum;
										global.loc=local.loc;
										global.active=true;
									}
								};
								auto r = reduction<parallelCombine,serialCombine, init, nullFunction>(std::forward<T>(array),std::forward<T_mask>(mask));
								Array<kind,1> loc(arrayInfo<T>::rank);
								assignArrayTuple(loc,r.loc);
								return loc;
							} else {
								auto serialCombine = [](const rop &local, rop &global){
									if (local.active &&((!global.active || indexTupleGT(local.loc,global.loc)) && (local.accum==global.accum ||
													local.accum<global.accum))){
										global.accum=local.accum;
										global.loc=local.loc;
										global.active=true;
									}
								};
								auto r = reduction<parallelCombine,serialCombine, init, nullFunction,true>(std::forward<T>(array),std::forward<T_mask>(mask));
								Array<kind,1> loc(arrayInfo<T>::rank);
								assignArrayTuple(loc,r.loc);
								return loc;
							}
						}

					/**
					 * Fortran MINLOC(https://fortranwiki.org/fortran/show/minloc)
					 * Array, kind and back
					 * with direction
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T=void, bool natural=false>
						auto minloc(T &&array, int direction, bool back=false){
							ERRCHECKREALINTEGERSTRING(T,"minloc","array");
							constexpr int rank = arrayInfo<T>::rank;
							using coreType = typename arrayInfo<T>::type;
							using reduction = indexInfo<rank>;
							using rtype = Array<kind, rank-1>;
							struct meta{
								bool active;
								coreType accum;
							};
							using marray = Array<meta,rank-1>;

							auto parallelCombine = [](const reduction& rri, const T& array, rtype &result, marray &ma){
								auto val = getItem(array,rri);
								meta &m = getItem(ma,rri);
								if (!m.active || (val<m.accum)){
									m.active=true;
									m.accum=val;
									getItem(result,rri)=getTupleLevel(rri.SBindices,rri.direction);
								}
							};

							auto initResult = [](kind &result){
								result = defaultLB-1;
							};

							auto initMeta = [](meta &m){
								m.active=false;
							};

							if (!back){
								auto serialCombine = [](const kind &local, const meta &lm, kind &global, meta &gm){
									if (lm.active && (!gm.active || (lm.accum<=gm.accum && local<global))){
										gm.active=true;
										gm.accum=lm.accum;
										global=local;
									}
								};
								auto r = reductionWithDirectionAndMeta<parallelCombine, serialCombine, initResult, initMeta, nullFunction>(direction, std::forward<T>(array));
								return r;
							} else {
								auto serialCombine = [](const kind &local, const meta &lm, kind &global, meta &gm){
									if (lm.active && (!gm.active || (lm.accum<=gm.accum && local>global))){
										gm.active=true;
										gm.accum=lm.accum;
										global=local;
									}
								};
								auto r = reductionWithDirectionAndMeta<parallelCombine, serialCombine, initResult, initMeta, nullFunction,true>(direction, std::forward<T>(array));
								return r;
							}
						}

					/**
					 * Fortran MINLOC(https://fortranwiki.org/fortran/show/minloc)
					 * Array, mask, kind and back
					 * with direction
					 */

					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T=void, typename T_mask=void, bool natural=false>
						auto minloc_with_mask(T &&array, int direction, T_mask&& mask, bool back=false){
							ERRCHECKREALINTEGERSTRING(T,"maxloc_with_mask","array");
							AERRCHECKTRUTHY(T_mask,"maxloc_with_mask","mask");
							constexpr int rank = arrayInfo<T>::rank;
							using coreType = typename arrayInfo<T>::type;
							using reduction = indexInfo<rank>;
							using rtype = Array<kind, rank-1>;
							struct meta{
								bool active;
								coreType accum;
							};
							using marray = Array<meta,rank-1>;

							auto parallelCombine = [](const reduction& rri, const T& array, rtype &result, marray &ma, const T_mask &mask){
								auto val = getItem(array,rri);
								meta &m = getItem(ma,rri);
								auto mval = getItem(mask,rri);
								if (mval && (!m.active || (val<m.accum))){
									m.active=true;
									m.accum=val;
									getItem(result,rri)=getTupleLevel(rri.SBindices,rri.direction);
								}
							};

							auto initResult = [](kind &result){
								result = defaultLB-1;
							};

							auto initMeta = [](meta &m){
								m.active=false;
							};

							if (!back){
								auto serialCombine = [](const kind &local, const meta &lm, kind &global, meta &gm){
									if (lm.active && (!gm.active || (lm.accum<=gm.accum && local<global))){
										gm.active=true;
										gm.accum=lm.accum;
										global=local;
									}
								};
								auto r = reductionWithDirectionAndMeta<parallelCombine, serialCombine, initResult, initMeta, nullFunction>(direction, std::forward<T>(array), std::forward<T_mask>(mask));
								return r;
							} else {
								auto serialCombine = [](const kind &local, const meta &lm, kind &global, meta &gm){
									if (lm.active && (!gm.active || (lm.accum<=gm.accum && local>global))){
										gm.active=true;
										gm.accum=lm.accum;
										global=local;
									}
								};
								auto r = reductionWithDirectionAndMeta<parallelCombine, serialCombine, initResult, initMeta, nullFunction,true>(direction, std::forward<T>(array), std::forward<T_mask>(mask));
								return r;
							}
						}

					/**
					 * Fortran MINVAL(https://fortranwiki.org/fortran/show/maxval)
					 */
					template <typename T>
						auto minval(T &&array)
						{
							ERRCHECKREALINTEGERSTRING(T,"minval","array");
							using coreType = typename arrayInfo<T>::type;
							struct rop{
								coreType accum;
								bool active=false;
							};
							auto parallelCombine = [](const coreType &src, rop &dest){
								if (!dest.active || (src<dest.accum)){
									dest.active=true;
									dest.accum=src;
								}
							};
							auto serialCombine = [](const rop &local, rop &global){
								if (local.active && ((local.accum<global.accum || !global.active))){
									global.accum=local.accum;
									global.active=true;
								}
							};
							auto initial = [](rop &dest){
								dest.active=false;
								//This is how Fortran does it, but in C++ can only do this for arithmetic types
								//For non-arithmetic types, the default constructor is used, so up to the user to
								//determine that a default constructed object is detected as "no maximum"
								if constexpr(std::is_arithmetic_v<coreType>){
									dest.accum = std::numeric_limits<coreType>::lowest();
								}
							};
							auto r = reduction<parallelCombine,serialCombine, initial, defaultFunction>(std::forward<T>(array));
							return r.accum;
						}

					/**
					 * Fortran MINVAL(https://fortranwiki.org/fortran/show/minval)
					 * with mask
					 */
					template <typename T, typename T_mask>
						auto minval_with_mask(T &&array, T_mask &&mask)
						{
							ERRCHECKREALINTEGERSTRING(T,"minval_with_mask","array");
							AERRCHECKTRUTHY(T_mask,"minval_with_mask","mask");
							using coreType = typename arrayInfo<T>::type;
							using maskType = typename arrayInfo<T_mask>::type;
							struct rop{
								coreType accum;
								bool active=false;
							};
							auto parallelCombine = [](const coreType &src, rop &dest, const maskType &mask){
								if (mask &&(!dest.active || (src<dest.accum))){
									dest.active=true;
									dest.accum=src;
								}
							};
							auto serialCombine = [](const rop &local, rop &global){
								if (local.active && (local.accum<global.accum || !global.active)){
									global.accum=local.accum;
									global.active=true;
								}
							};
							auto initial = [](rop &dest){
								dest.active=false;
								//This is how Fortran does it, but in C++ can only do this for arithmetic types
								//For non-arithmetic types, the default constructor is used, so up to the user to
								//determine that a default constructed object is detected as "no maximum"
								if constexpr(std::is_arithmetic_v<coreType>){
									dest.accum = std::numeric_limits<coreType>::lowest();
								}
							};
							auto r = reduction<parallelCombine,serialCombine, initial, defaultFunction>(std::forward<T>(array),std::forward<T_mask>(mask));
							return r.accum;
						}

					/**
					 * Fortran MINVAL(https://fortranwiki.org/fortran/show/minval)
					 * With dim
					 */
					template <typename T>
						auto minval(T &&array, int dim)
						{
							ERRCHECKREALINTEGERSTRING(T,"minval","array");
							using coreType = typename arrayInfo<T>::type;
							struct meta{
								bool active;
							};

							auto parallelCombine = [](const coreType& val, coreType &result, meta &m){
								if (!m.active || (val<result)){
									m.active=true;
									result=val;
								}
							};

							auto initResult = [](coreType &result){
								if constexpr(std::is_arithmetic_v<coreType>){
									result = std::numeric_limits<coreType>::lowest();
								}

							};

							auto initMeta = [](meta &m){
								m.active=false;
							};

							auto serialCombine = [](const coreType &local, const meta &lm, coreType &global, meta &gm){
								if (!gm.active || (lm.active && local>global)){
									gm.active=true;
									global=local;
								}
							};
							auto r = reductionWithDirectionAndMeta<parallelCombine, serialCombine, initResult, initMeta>(dim, std::forward<T>(array));
							return r;
						}

					/**
					 * Fortran MINVAL(https://fortranwiki.org/fortran/show/minval)
					 * With dim
					 */
					template <typename T, typename T_mask>
						auto minval_with_mask(T &&array, int dim, T_mask &&mask)
						{
							ERRCHECKREALINTEGERSTRING(T,"minval","array");
							ERRCHECKTRUTHY(T_mask,"minval","mask");
							using coreType = typename arrayInfo<T>::type;
							using maskType = typename arrayInfo<T_mask>::type;
							struct meta{
								bool active;
							};

							auto parallelCombine = [](const coreType& val, coreType &result, meta &m, const maskType &mask){
								if (mask && (!m.active || (val<result))){
									m.active=true;
									result=val;
								}
							};

							auto initResult = [](coreType &result){
								if constexpr(std::is_arithmetic_v<coreType>){
									result = std::numeric_limits<coreType>::lowest();
								}

							};

							auto initMeta = [](meta &m){
								m.active=false;
							};

							auto serialCombine = [](const coreType &local, const meta &lm, coreType &global, meta &gm){
								if (!gm.active || (lm.active && local>global)){
									gm.active=true;
									global=local;
								}
							};
							auto r = reductionWithDirectionAndMeta<parallelCombine, serialCombine, initResult, initMeta>(dim, std::forward<T>(array),std::forward<T_mask>(mask));
							return r;
						}

					/**
					 * Fortran MOD(https://fortranwiki.org/fortran/show/mod)
					 */
					template <typename T_a, typename T_p>
						auto mod(T_a &&a, T_p &&p)
						{
							AERRCHECKREALINTEGER(T_a,"mod","a");
							AERRCHECKREALINTEGER(T_p,"mod","p");
							auto l = [](typename arrayInfo<std::decay_t<T_a>>::type const &a,
									typename arrayInfo<std::decay_t<T_p>>::type const &p
									) -> decltype(auto)
							{ 
								return a - (std::trunc((double)a/(double)p))*p;
							};
							return makeLazy<l>(std::forward<T_a>(a),std::forward<T_p>(p));
						}

					/**
					 * Fortran MODULO(https://fortranwiki.org/fortran/show/modulo)
					 */
					template <typename T_a, typename T_p>
						auto modulo(T_a &&a, T_p &&p)
						{
							AERRCHECKREALINTEGER(T_a,"mod","a");
							AERRCHECKREALINTEGER(T_p,"mod","p");
							auto l = [](typename arrayInfo<std::decay_t<T_a>>::type const &a,
									typename arrayInfo<std::decay_t<T_p>>::type const &p
									)
							{
								return a - (std::floor((double)a/(double)p))*p;
							};
							return makeLazy<l>(std::forward<T_a>(a),std::forward<T_p>(p));
						}

					/**
					 * Fortran MOVE_ALLOC(https://fortranwiki.org/fortran/show/move_alloc)
					 */
					template<typename T>
						void move_alloc(T&&src, T&&dest){
							dest = std::move(src);
						}

					/**
					 * Fortran MVBITS(https://fortranwiki.org/fortran/show/mvbits)
					 */
					template<typename T_from, typename T_frompos, typename T_len, typename T_to, typename T_topos>
						void mvbits(T_from &&from, T_frompos &&frompos, T_len &&len, T_to &&to, T_topos &&topos){
							AERRCHECKINTEGER(T_from,"mvbits","from");
							AERRCHECKINTEGER(T_frompos,"mvbits","frompos");
							AERRCHECKINTEGER(T_len,"mvbits","len");
							AERRCHECKINTEGER(T_to,"mvbits","to");
							AERRCHECKINTEGER(T_topos,"mvbits","topos");
							using fromitem = arrayInfo<T_from>::type;
							using frompositem = arrayInfo<T_frompos>::type;
							using lenitem = arrayInfo<T_len>::type;
							using toitem = arrayInfo<T_to>::type;
							using topositem = arrayInfo<T_topos>::type;
							constexpr bool anyArray = arrayInfo<T_from>::value || arrayInfo<T_to>::value || arrayInfo<T_frompos>::value || arrayInfo<T_len>::value || arrayInfo<T_topos>::value;
							static_assert(std::is_same_v<fromitem,toitem>,"From and To must be of the same kind in mvbits");
							//To can be an array when from isn't, but From can't be an array if To isn't
							static_assert((arrayInfo<T_from>::value == arrayInfo<T_to>::value) || arrayInfo<T_to>::value, "To must be an array if From is an array in mvbits");

							auto l = [](const fromitem &from, const frompositem &frompos, const lenitem &len, toitem &to, const topositem &topos){
								to = ((to & (((~(((toitem(1)<<len)-1)<<topos))))) | (((from & ((toitem(1)<<len)-1)<<frompos)>>frompos))<<topos);
							};
							if constexpr(anyArray){
								forEach<l>(std::forward<T_from>(from),std::forward<T_frompos>(frompos),std::forward<T_len>(len),std::forward<T_to>(to),std::forward<T_topos>(topos));
							} else {
								l(std::forward<T_from>(from),std::forward<T_frompos>(frompos), std::forward<T_len>(len),
										std::forward<T_to>(to), std::forward<T_topos>(topos));
							}
						}

					/**
					 * Fortran NEAREST(https://fortranwiki.org/fortran/show/nearest)
					 */
					template <typename T_x, typename T_s>
						auto nearest(T_x &&x, T_s &&s)
						{
							AERRCHECKREAL(T_x,"nearest","x");
							AERRCHECKREAL(T_s,"nearest","s");
							using xitem = arrayInfo<T_x>::type;
							using sitem = arrayInfo<T_s>::type;
							static_assert(is_floating_point<T_x>, "x must be floating point in nearest");
							auto l = [](xitem const &x, sitem const &s) -> auto
							{
								auto dir = (s>=0?std::numeric_limits<xitem>::max():std::numeric_limits<xitem>::lowest());
								return std::nextafter(x,dir);
							};
							return makeLazy<l>(std::forward<T_x>(x),std::forward<T_s>(s));
						}

					/**
					 * Fortran NINT(https://fortranwiki.org/fortran/show/nint)
					 */
					template <typename kind = FAR_SIGNED_INDEX_TYPE, typename T_x=int>
						auto nint(T_x &&x)
						{
							AERRCHECKREALCOMPLEX(T_x,"nint","x");
							using xitem = arrayInfo<T_x>::type;
							auto l = [](xitem const &x) -> kind
							{
								if constexpr(is_complex<xitem>){
									return static_cast<kind>(std::llround(x.real()));
								} else {
									return static_cast<kind>(std::llround(x));
								}
							};
							return makeLazy<l>(std::forward<T_x>(x));
						}

					/**
					 * Fortran NORM2(https://fortranwiki.org/fortran/show/norm2)
					 */
					template <typename T_array>
						auto norm2(T_array &&array)
						{
							AERRCHECKREAL(T_array,"norm2","array");
							using type = typename arrayInfo<std::decay_t<T_array>>::type;
							return std::sqrt(reduction<[](const type &src, type &dest){dest+=src*src;},
									[](const type &src, type &dest){dest+=src*src;},
									[](type &dest){dest=0;}
									>(std::forward<T_array>(array)));
						}

					/**
					 * Fortran NORM2(https://fortranwiki.org/fortran/show/norm2)
					 * with dim
					 */
					template <typename T_array>
						auto norm2(T_array &&array, int dim)
						{
							AERRCHECKREAL(T_array,"norm2","array");
							using type = typename arrayInfo<std::decay_t<T_array>>::type;
							Array<type,arrayInfo<T_array>::rank-1> result;
							result = reductionWithDirection<[](const type &src, type &dest){dest+=src*src;},
										 [](const type &src, type &dest){dest+=src*src;},
										 [](type &dest){dest=0;}
							>(dim,std::forward<T_array>(array));

							forEach<[](type &dest){dest=std::sqrt(dest);}>(result);
							return result;
						}

					/**
					 * Fortran NOT(https://fortranwiki.org/fortran/show/not)
					 */
					template <typename T_i>
						auto f_not(T_i &&i)
						{
							AERRCHECKINTEGER(T_i,"not","i");
							using type = typename arrayInfo<T_i>::type;
							return makeLazy<[](type const &i) -> type
							{
								return ~i;
							}>(std::forward<T_i>(i));
						}

					/**
					 * Create a 1D array of items all given the value one
					 */
					template <typename T = double, bounds_check_state bounds = bc_default, typename... T_sizes>
						Array<T, sizeof...(T_sizes), bounds> ones(T_sizes... sizes)
						{ 
							Array<T, sizeof...(T_sizes), bounds> space(sizes...);
							space = 1;
							return space;
						}


					/**
					 * Fortran OUT_OF_RANGE
					 */
					template <typename mold_type, typename T_array>
						auto out_of_range(T_array &&array)
						{
							using type = typename arrayInfo<std::decay_t<T_array>>::type;
							using compare_type = long double;
							// Produce an array with the summed direction removed
							Array<bool,arrayInfo<T_array>::rank> dest;
							dest.mold(array);
							forEach<[]
								(bool &dest, const type &x) { 
									//If the range of the mold is greater than the range of the type, then it is always in range
									if constexpr(static_cast<compare_type>(std::numeric_limits<mold_type>::lowest()) <= static_cast<compare_type>(std::numeric_limits<type>::lowest()) && 
											static_cast<compare_type>(std::numeric_limits<mold_type>::max()) >= static_cast<compare_type>(std::numeric_limits<type>::max())){
										dest = false;
										return;
									}

									if constexpr(std::is_integral_v<mold_type>){
										type rounded;
										if constexpr(std::is_floating_point_v<type>){
										 	rounded = (x > 0) ? std::floor(x) : std::ceil(x);
										} else {rounded = x;}
										if (std::isinf(rounded) || std::isnan(rounded)) {
											dest = true;
										} else {
											if ((static_cast<compare_type>(rounded) < static_cast<compare_type>(std::numeric_limits<mold_type>::lowest()) || 
													static_cast<compare_type>(rounded) > static_cast<compare_type>(std::numeric_limits<mold_type>::max()))) 
											{
												dest = true;
											} else {
												dest = false;
											}
										}
									} else {
										if (std::isinf(x) || std::isnan(x)) {
											dest = true;
										} else {
											if ((static_cast<compare_type>(x) < static_cast<compare_type>(std::numeric_limits<mold_type>::lowest()) || 
													static_cast<compare_type>(x) > static_cast<compare_type>(std::numeric_limits<mold_type>::max()))) 
											{
												dest = true;
											} else {
												dest = false;
											}
										}
									}
								}>(dest, array);

							return dest;
						}

					/**
					 * Fortran OUT_OF_RANGE
					 */
					template <typename T_array, typename T_mold>
						auto out_of_range(T_array &&array, FAR_UNUSED T_mold &&mold)
						{
							return out_of_range<typename arrayInfo<T_mold>::type>(std::forward<T_array>(array));
						}

					/**
					 * Fortran OUT_OF_RANGE
					 */
					template <typename mold_type, typename T_array>
						auto out_of_range(T_array &&array, bool round)
						{
							using type = typename arrayInfo<std::decay_t<T_array>>::type;
							using compare_type = long double;
							// Produce an array with the summed direction removed
							Array<bool,arrayInfo<T_array>::rank> dest;
							dest.mold(array);
							forEach<[]
								(bool &dest, const type &x, bool round)
								{ 
									if constexpr(std::is_integral_v<mold_type>){
										if (round) {
											type rounded;
											if constexpr(std::is_floating_point_v<type>){
											 	rounded = std::round(x);
											} else {rounded = x;}
											if (std::isinf(rounded) || std::isnan(rounded)) {
												dest = true;
											} else if ((static_cast<compare_type>(rounded) < static_cast<compare_type>(std::numeric_limits<mold_type>::lowest()) || 
													static_cast<compare_type>(rounded) > static_cast<compare_type>(std::numeric_limits<mold_type>::max())))
											{
												dest = true;
											} else {
												dest = false;
											}
										} else {
											type rounded;
											if constexpr(std::is_floating_point_v<type>){
											 	rounded = (x > 0) ? std::floor(x) : std::ceil(x);
											} else {rounded = x;}
											if (std::isinf(rounded) || std::isnan(rounded)) {
												dest = true;
											} else {
												if ((static_cast<compare_type>(rounded) < static_cast<compare_type>(std::numeric_limits<mold_type>::lowest()) || 
														static_cast<compare_type>(rounded) > static_cast<compare_type>(std::numeric_limits<mold_type>::max())))
												{
													dest = true;
												} else {
													dest = false;
												}
											}
										}
									} else {
										if (std::isinf(x) || std::isnan(x)) {
											dest = true;
										} else {
											type rounded = std::round(x);
											if ((static_cast<compare_type>(rounded) < static_cast<compare_type>(std::numeric_limits<mold_type>::lowest()) || 
													static_cast<compare_type>(rounded) > static_cast<compare_type>(std::numeric_limits<mold_type>::max())))
											{
												dest = true;
											} else {
												dest = false;
											}
										}
									}
								}>(dest, array, round);

							return dest;
						}


					/**
					 * Fortran OUT_OF_RANGE
					 */
					template <typename T_array, typename T_mold>
						auto out_of_range(T_array &&array, FAR_UNUSED T_mold &&mold, bool round)
						{
							return out_of_range<typename arrayInfo<T_mold>::type>(std::forward<T_array>(array), round);
						}


					/**
					 * Fortran PACK(https://fortranwiki.org/fortran/show/pack)
					 * This version isn't in the Fortran standard and is equivalent to pack(array,.TRUE.)
					 * It is included here because it feels like a more complete interface
					 */
					template <typename T_array>
						auto pack(T_array &&array)
						{
							using type = arrayInfo<T_array>::type;
							constexpr int rank = arrayInfo<T_array>::rank;
							auto packed = array.template buildVariant<type,1>();
							packed.allocate(array.getSize());
							FAR_SIGNED_INDEX_TYPE index=defaultLB;

							auto mapper = [] (const indexInfo<rank> &info, const T_array &array, decltype(packed) &dest, FAR_SIGNED_INDEX_TYPE &index){
								return std::make_tuple(std::ref(getItem(array,info)), std::ref(dest), std::ref(index));};

							auto packer = [] (const type&src, decltype(packed) &dest, FAR_SIGNED_INDEX_TYPE &index){
								dest(index) = src;
								index++;
							};

							typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
							array.getIndexer().template forEachCore<packer, mapper>(t1,t2,rdNone,array, packed, index);
							return packed;
						}


					/**
					 * Fortran PACK(https://fortranwiki.org/fortran/show/pack)
					 * With mask parameter. Returns as many elements as match mask
					 */
					template <typename T_array, typename T_mask>
						auto pack(T_array &&array, T_mask &&mask)
						{
							using type = arrayInfo<T_array>::type;
							using mask_type = arrayInfo<T_mask>::type;
							constexpr int rank = arrayInfo<T_array>::rank;
							auto packed = array.template buildVariant<type,1>();
							//If mask is an array then count true elements
							//otherwise use the size of the array
							FAR_SIGNED_INDEX_TYPE items;
							if constexpr(arrayInfo<T_mask>::value){
								items = count(mask);
							} else {
								items = array.getSize();
							}
							packed.allocate(items);
							if (items==0) return packed;
							FAR_SIGNED_INDEX_TYPE index=packed.getLB(1);

							auto mapper = [] (const indexInfo<rank> &info, const T_array &array, const T_mask &mask, decltype(packed) &dest, FAR_SIGNED_INDEX_TYPE &index){
								return std::make_tuple(getItem(array,info), getItem(mask,info), std::ref(dest), std::ref(index));};

							auto packer = [] (const type&src, const mask_type&mask, decltype(packed) &dest, FAR_SIGNED_INDEX_TYPE &index){
								if (mask){
									dest(index) = src;
									index++;
								}
							};

							typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
							array.getIndexer().template forEachCore<packer, mapper>(t1, t2, rdNone, array, mask, packed, index);
							return packed;
						}

					/**
					 * Fortran PACK(https://fortranwiki.org/fortran/show/pack)
					 * With mask parameter and vector parameter
					 */
					template <typename T_array, typename T_mask, typename T_vector>
						auto pack(T_array &&array, T_mask &&mask, T_vector &&vector)
						{
							using type = arrayInfo<T_array>::type;
							using mask_type = arrayInfo<T_mask>::type;
							constexpr int rank = arrayInfo<T_array>::rank;
							auto packed = array.template buildVariant<type,1>();
							//If mask is an array then count true elements
							//otherwise use the size of the array
							FAR_SIGNED_INDEX_TYPE items = vector.getSize();
							packed.allocate(items);
							if (items==0) return packed;
							FAR_SIGNED_INDEX_TYPE index=packed.getLB(1);

							auto mapper = [] (const indexInfo<rank> &info, const T_array &array, const T_mask &mask, decltype(packed) &dest, FAR_SIGNED_INDEX_TYPE &index){
								return std::make_tuple(getItem(array,info), getItem(mask,info), std::ref(dest), std::ref(index));};
							auto packer = [] (const type&src, const mask_type&mask, decltype(packed) &dest, FAR_SIGNED_INDEX_TYPE &index){
								if (mask){
									dest(index) = src;
									index++;
								}
							};

							typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
							array.getIndexer().template forEachCore<packer, mapper>(t1, t2, rdNone, array, mask, packed, index);
							//Have now packed everything from array, now pack from mask
							for (;index<=vector.getUB(1);++index){
								packed(index)=vector(index);
							}
							return packed;
						}

					/**
					 * Fortran PARITY(https://fortranwiki.org/fortran/show/parity)
					 */
					template <typename T>
						auto parity(T &&source)
						{
							AERRCHECKTRUTHY(T,"parity","mask")
								using type = typename arrayInfo<std::decay_t<T>>::type;
							return reduction<[](const type &src, type &dest){dest=src?!dest:dest;},
										 defaultFunction,
										 [](type &dest){dest=type(false);}
							>(std::forward<T>(source));
						}


					/**
					 * Fortran PARITY(https://fortranwiki.org/fortran/show/parity)
					 * with dir
					 */
					template <typename itype=FAR_UNSIGNED_INDEX_TYPE, typename T=double>
						auto parity(T &&source, int dim)
						{
							AERRCHECKTRUTHY(T,"parity","mask")
								using type = typename arrayInfo<std::decay_t<T>>::type;
							return reductionWithDirection<[](const type &src, type &dest){dest=src?!dest:dest;},
										 defaultFunction,
										 [](type &dest){dest=type(false);}
							>(dim, std::forward<T>(source));
						}

					/**
					 * Fortran POPCNT(https://fortranwiki.org/fortran/show/popcnt)
					 */
					template <typename T_i>
						auto popcnt(T_i &&i)
						{
							AERRCHECKINTEGER(T_i,"popcnt","i");
							auto l = [](typename arrayInfo<T_i>::type const &i) -> decltype(auto)
							{
								return std::popcount(i);
							};
							return makeLazy<l>(std::forward<T_i>(i));
						}

					/**
					 * Fortran POPPAR(https://fortranwiki.org/fortran/show/poppar)
					 */
					template <typename T_i>
						auto poppar(T_i &&i)
						{
							AERRCHECKINTEGER(T_i,"poppar","i");
							auto l = [](typename arrayInfo<T_i>::type const &i) -> decltype(auto)
							{ 
								return std::popcount(i) % 2;
							};
							return makeLazy<l>(std::forward<T_i>(i));
						}


					/**
					 * Fortran PRECISION(https://fortranwiki.org/fortran/show/precision)
					 */
					template <typename T_x>
						auto precision(FAR_UNUSED T_x &&x)
						{
							AERRCHECKREALCOMPLEX(T_x,"precision","x");
							if constexpr(!is_complex<T_x>){
								return std::numeric_limits<typename arrayInfo<T_x>::type>::digits10;
							} else {
								return std::numeric_limits<typename arrayInfo<T_x>::type::value_type>::digits10;
							}
						}

					/**
					 * Fortran PRODUCT(https://fortranwiki.org/fortran/show/product)
					 */
					template <typename T>
						auto product(T &&source)
						{
							ERRCHECKREALINTEGERCOMPLEX(T,"product","array");
							using type = typename arrayInfo<std::decay_t<T>>::type;
							return reduction<[](const type &src, type &dest){dest*=src;},
										 defaultFunction,
										 [](type &dest){dest=type(1);}
							>(std::forward<T>(source));
						}

					/**
					 * Fortran PRODUCT(https://fortranwiki.org/fortran/show/product)
					 * with mask
					 */
					template <typename T, typename T_mask>
						auto product_with_mask(T &&source, T_mask &&mask)
						{
							ERRCHECKREALINTEGERCOMPLEX(T,"product","array");
							AERRCHECKTRUTHY(T_mask,"product","mask");
							using type = typename arrayInfo<std::decay_t<T>>::type;
							using masktype = typename arrayInfo<std::decay_t<T_mask>>::type;
							return reduction<[](const type &src, type &dest, const masktype &mask){if (mask) dest*=src;},
										 defaultFunction,
										 [](type &dest){dest=type(1);}
							>(std::forward<T>(source),std::forward<T_mask>(mask));
						}

					/**
					 * Fortran PRODUCT(https://fortranwiki.org/fortran/show/product)
					 * with dim
					 */
					template <typename T>
						auto product(T &&source,const int dim)
						{
							ERRCHECKREALINTEGERCOMPLEX(T,"product","array");
							using type = typename arrayInfo<std::decay_t<T>>::type;
							return reductionWithDirection<[](const type &src, type &dest){dest*=src;},
										 defaultFunction,
										 [](type &dest){dest=type(1);}
							>(dim, std::forward<T>(source));

						}

					/**
					 * Fortran PRODUCT(https://fortranwiki.org/fortran/show/product)
					 * with dim and mask
					 */
					template <typename T, typename T_mask>
						auto product_with_mask(T &&source, const int dim, T_mask &&mask)
						{
							ERRCHECKREALINTEGERCOMPLEX(T,"product","array");
							AERRCHECKTRUTHY(T_mask,"product","mask");
							using type = typename arrayInfo<std::decay_t<T>>::type;
							using masktype = typename arrayInfo<std::decay_t<T_mask>>::type;
							return reductionWithDirection<[](const type &src, type &dest, const masktype &mask){if(mask)dest*=src;},
										 defaultFunction,
										 [](type &dest){dest=type(1);}
							>(dim, std::forward<T>(source),std::forward<T_mask>(mask));
						}

					/**
					 * POW. In Fortran raising to the power of is an operator, but it isn't in C++
					 * This implements POW to implement it in C++
					 */
					template <typename T_x, typename T_power>
						auto pow(T_x &&x, T_power &&power)
						{
							AERRCHECKREALINTEGERCOMPLEX(T_x,"pow","x");
							AERRCHECKREALINTEGER(T_power,"pow","power");
							using type = typename arrayInfo<std::decay_t<T_x>>::type;
							using powertype = typename arrayInfo<std::decay_t<T_power>>::type;
							auto l = [](type const &x, 
									powertype const &power) -> type
							{
								return std::pow(x,power);
							};
							return makeLazy<l>(std::forward<T_x>(x), std::forward<T_power>(power));
						}


					/**
					 * Fortran RADIX(https://fortranwiki.org/fortran/show/radix)
					 */
					template <typename T>
						constexpr auto radix(FAR_UNUSED T &&source)
						{
							return std::numeric_limits<typename arrayInfo<T>::type>::radix;
						}

					

					template<typename T_Engine = FAR_DEFAULT_GLOBAL_RNG>
						struct RNGHandle{
							private:
								// Helper struct to check for existence of ::state_size
								template <typename Engine, typename = void>
									struct has_state_size : std::false_type {};

								template <typename Engine>
									struct has_state_size<Engine, std::void_t<decltype(Engine::state_size)>> : std::true_type {};

								// Primary template for non-Mersenne Twister engines
								template <typename Engine, bool = has_state_size<Engine>::value>
									struct SeedSize {
										static constexpr FAR_UNSIGNED_INDEX_TYPE value = sizeof(typename Engine::result_type) / sizeof(uint32_t);
									};

								// Specialization for engines with ::state_size
								template <typename Engine>
									struct SeedSize<Engine, true> {
										static constexpr FAR_UNSIGNED_INDEX_TYPE value = Engine::state_size;
									};

							public:
								//The actual random engine
								T_Engine engine{std::random_device{}()};

								static constexpr auto randomSeedSize(){
									return SeedSize<T_Engine>::value;
								}
								void randomSeedGet(std::ostream &s){
									s << engine;
								}

								void randomSeedPut(std::istream &s){
									s >> engine;
								}

								void randomSeedPut(const Array<uint32_t,1> &seed){
									std::seed_seq seq(seed.begin(), seed.end());
									engine.seed(seq);
								}

						};

					//Concept for an RNGHandle however templated
					template<typename T>
						concept RNGHandleConcept = requires(T a){
							{a.engine};
							{a.randomSeedGet(std::declval<std::ostream&>())};
							{a.randomSeedPut(std::declval<std::istream&>())};
							{a.randomSeedPut(std::declval<Array<uint32_t,1>&>())};
							{a.randomSeedSize()};
						};

					namespace{
						inline static RNGHandle globalRNGHandle;
					};

					template<typename T=double, typename... T_indices>
						auto random_array(T_indices... sizes){
							Array<T,sizeof...(T_indices)> result;
							result.allocate(sizes...);
							random_number(result);
							return result;
						}

					template<typename T=double, typename T_handle, typename... T_indices>
						auto random_array(T_handle &handle, T_indices... sizes){
							Array<T,sizeof...(T_indices)> result;
							result.allocate(sizes...);
							random_number(handle,result);
							return result;
						}

					/**
					 * RANDOM_DISTRIBUTION. An extension to RANDOM_NUMBER for generating random numbers from a
					 * provided distribution. Uses a provided RNGHandle object to allow multiple RNGs
					 */
					template <typename T_handle, typename T_harvest, typename T_dist>
						auto random_distribution(T_handle &handle, T_harvest &harvest, T_dist &&dist)
						{
							AERRCHECKREAL(T_harvest,"random_distribution","harvest");
							using typehandle = typename arrayInfo<T_handle>::type;
							using type = typename arrayInfo<T_harvest>::type;
							using typedist = typename arrayInfo<T_dist>::type;
							if constexpr(arrayInfo<T_harvest>::value){
								forEach<[](type &harvest, typedist &dist, typehandle &handle) {
									harvest = dist(handle.engine);
								}>(std::forward<T_harvest>(harvest),std::forward<T_dist>(dist), handle);
							} else {
								harvest = dist(handle.engine);
							}
						}

					/**
					 * RANDOM_DISTRIBUTION. An extension to RANDOM_NUMBER for generating random numbers from a
					 * provided distribution. Uses global RNG
					 */
					template <typename T_harvest, typename T_dist>
						auto random_distribution(T_harvest &harvest, T_dist &&dist)
						{
							random_distribution(globalRNGHandle,harvest,dist);
						}

					/**
					 * RANDOM_NORMAL. An extension to RANDOM_NUMBER for generating normally distributed
					 * numbers with a specified mean and standard deviation
					 */
					template <typename T_handle, typename T_harvest, typename T_mean, typename T_stdev>
						auto random_normal(T_handle &handle, T_harvest &harvest, const T_mean &&mean, const T_stdev &&stdev)
						{
							AERRCHECKREAL(T_harvest,"random_normal","harvest");
							AERRCHECKREAL(T_mean,"random_normal","mean");
							AERRCHECKREAL(T_stdev,"random_notmal","stdev");
							using typehandle = typename arrayInfo<T_handle>::type;
							using type = typename arrayInfo<T_harvest>::type;
							using typemean = typename arrayInfo<T_mean>::type;
							using typestdev = typename arrayInfo<T_stdev>::type;
							if constexpr(arrayInfo<T_harvest>::value){
								forEach<[](typehandle &handle, type &harvest, const typemean &mean, const typestdev &stdev) {
									std::normal_distribution<type> distribution{mean,stdev};
									harvest = distribution(handle.engine);
								}>(std::forward<T_handle>(handle), std::forward<T_harvest>(harvest),std::forward<const T_mean>(mean),std::forward<const T_stdev>(stdev));
							} else {
								std::normal_distribution<type> distribution{mean,stdev};
								harvest = distribution(handle.engine);
							}
						}

					/**
					 * RANDOM_NORMAL. An extension to RANDOM_NUMBER for generating normally distributed
					 * numbers with a specified mean and standard deviation
					 */
					template <typename T_harvest, typename T_mean, typename T_stdev>
						void random_normal(T_harvest &harvest, const T_mean &&mean, const T_stdev &&stdev){
							random_normal(globalRNGHandle,harvest,std::forward<const T_mean>(mean),std::forward<const T_stdev>(stdev));
						}


					/**
					 * Fortran RANDOM_NUMBER(https://fortranwiki.org/fortran/show/random_number)
					 * Uses a provided RNGHandle object to allow multiple RNGs
					 */
					template <typename T_handle, typename T_harvest>
						void random_number(T_handle &handle, T_harvest &harvest)
						{
							AERRCHECKREAL(T_harvest,"random_number","harvest");
							using typehandle = typename arrayInfo<T_handle>::type;
							using type = typename arrayInfo<T_harvest>::type;
							if constexpr(arrayInfo<T_harvest>::value){
								forEach<[](type &harvest, typehandle &handle) {
									static std::uniform_real_distribution<type> distribution(0.0, 1.0);
									harvest = distribution(handle.engine);
								}>(harvest,handle);
							} else {
								static std::uniform_real_distribution<type> distribution(0.0, 1.0);
								harvest = distribution(handle.engine);
							}
						}


					/**
					 * Fortran RANDOM_NUMBER(https://fortranwiki.org/fortran/show/random_number)
					 */
					template <typename T_harvest>
						void random_number(T_harvest &harvest)
						{
							random_number(globalRNGHandle, harvest);
						}

					/**
					 * Fortran RANDOM_SEED(https://fortranwiki.org/fortran/show/random_seed)
					 * Implements the "size" version of RANDOM_SEED
					 */
					template<RNGHandleConcept T_handle>
						void random_seed_size(FAR_UNUSED T_handle &handle, FAR_SIGNED_INDEX_TYPE &size){
							size = arrayInfo<T_handle>::type::randomSeedSize();
						}

					/**
					 * Fortran RANDOM_SEED(https://fortranwiki.org/fortran/show/random_seed)
					 * Implements the "size" version of RANDOM_SEED
					 */
					inline void random_seed_size(FAR_SIGNED_INDEX_TYPE &size){
						size = arrayInfo<decltype(globalRNGHandle)>::type::randomSeedSize();
					}

					/**
					 * Fortran RANDOM_SEED(https://fortranwiki.org/fortran/show/random_seed)
					 * Implements the "size" version of RANDOM_SEED. Since I am having to change
					 * the interface anyway, I'll make this one have a C++ style return
					 */
					template<RNGHandleConcept T_handle>
						constexpr FAR_SIGNED_INDEX_TYPE random_seed_size(FAR_UNUSED T_handle &handle){
							return arrayInfo<T_handle>::type::randomSeedSize();
						}

					/**
					 * Fortran RANDOM_SEED(https://fortranwiki.org/fortran/show/random_seed)
					 * Implements the "size" version of RANDOM_SEED. Since I am having to change
					 * the interface anyway, I'll make this one have a C++ style return
					 */
					constexpr FAR_SIGNED_INDEX_TYPE random_seed_size(){
						return random_seed_size(globalRNGHandle);
					}

					/**
					 * Fortran RANDOM_SEED(https://fortranwiki.org/fortran/show/random_seed)
					 * Implements the "put" version of RANDOM_SEED for an array
					 */
					template<typename T_handle>
						void random_seed_put(T_handle &handle, Array<uint32_t,1> &data){
							using handletype = typename arrayInfo<T_handle>::type;
							if constexpr(arrayInfo<T_handle>::value){
								constexpr int rank = arrayInfo<T_handle>::rank;
								typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
								handle.getIndexer().template forEachCore<[](handletype &handle, Array<uint32_t,1> &data){
									handle.randomSeedPut(data);
								},
									[](const indexInfo<rank> &rii, std::decay_t<T_handle>& handle, Array<uint32_t,1> &data){
										return std::make_tuple(std::ref(getItem(handle,rii)),std::ref(data));
									}
								>(t1,t2,rdNone,handle,data);
							} else {
								handle.randomSeedPut(data);
							}
						}

					/**
					 * Fortran RANDOM_SEED(https://fortranwiki.org/fortran/show/random_seed)
					 * Implements the "put" version of RANDOM_SEED for an array
					 */
					template<typename T_handle, typename std::enable_if_t<arrayInfo<T_handle>::value, int>* = nullptr>
						void random_seed_put(T_handle &handle, Array<uint32_t,arrayInfo<T_handle>::rank+1> &data){
							using handletype = typename arrayInfo<T_handle>::type;
							constexpr int rank = arrayInfo<T_handle>::rank;
							if constexpr(arrayInfo<T_handle>::value){
								typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
								handle.getIndexer().template forEachCore<
									[](const indexInfo<rank> &rii, handletype &handle, Array<uint32_t,rank+1> &data){
										auto index = std::tuple_cat(rii.SBindices,std::make_tuple(Range()));
										handle.randomSeedPut(std::apply([&data](auto&...params){return data(params...);},index));
									},
									[](const indexInfo<rank> &rii, std::decay_t<T_handle>& handle, Array<uint32_t,rank+1> &data){
										return std::make_tuple(rii, std::ref(getItem(handle,rii)), std::ref(data));
									}
								>(t1,t2,rdNone,handle,data);
							} else {
								handle.randomSeedPut(data);
							}
						}

					/**
					 * Fortran RANDOM_SEED(https://fortranwiki.org/fortran/show/random_seed)
					 * Implements the "put" version of RANDOM_SEED for an array
					 */
					inline void random_seed_put(Array<uint32_t,1> &data){
						random_seed_put(globalRNGHandle,data);
					}


					/**
					 * Fortran RANDOM_SEED(https://fortranwiki.org/fortran/show/random_seed)
					 * Implements the "get" version of RANDOM_SEED to a stringstream
					 */
					template<typename T_handle>
						void random_seed_get(T_handle &handle, std::ostream &s){
							using handletype = typename arrayInfo<T_handle>::type;
							if constexpr(arrayInfo<T_handle>::value){
								forEach<[](handletype &handle, std::ostream &s) {
									handle.randomSeedGet(s);
								}>(handle,s);
							} else {
								handle.randomSeedGet(s);
							}
						}

					/**
					 * Fortran RANDOM_SEED(https://fortranwiki.org/fortran/show/random_seed)
					 * Implements the "get" version of RANDOM_SEED to a stringstream
					 */
					inline void random_seed_get(std::ostream &s){
						random_seed_get(globalRNGHandle,s);
					}


					/**
					 * Fortran RANDOM_SEED(https://fortranwiki.org/fortran/show/random_seed)
					 * Implements the "put" version of RANDOM_SEED to a stringstream
					 */
					template<typename T_handle>
						void random_seed_put(T_handle &handle, std::istream &s){
							using handletype = typename arrayInfo<T_handle>::type;
							if constexpr(arrayInfo<T_handle>::value){
								forEach<[](handletype &handle, std::istream &s) {
									handle.randomSeedPut(s);
								}>(handle,s);
							} else {
								handle.randomSeedPut(s);
							}
						}


					/**
					 * Fortran RANDOM_SEED(https://fortranwiki.org/fortran/show/random_seed)
					 * Implements the "put" version of RANDOM_SEED to a stringstream
					 */
					inline void random_seed_put(std::istream &s){
						random_seed_put(globalRNGHandle,s);
					}

					/**
					 * Fortran RANGE(https://fortranwiki.org/fortran/show/range)
					 */
					template <typename T_x>
						constexpr auto range(FAR_UNUSED T_x &&x)
						{
							if constexpr(!is_complex<T_x>){
								return std::numeric_limits<typename arrayInfo<T_x>::type>::max_exponent10-1;
							} else {
								return std::numeric_limits<typename arrayInfo<T_x>::type::value_type>::max_exponent10-1;
							}
						}

					/**
					 * Fortran RANK(https://fortranwiki.org/fortran/show/rank)
					 */
					template <typename T_x>
						constexpr auto rank(FAR_UNUSED T_x &&x)
						{
							return arrayInfo<T_x>::rank;
						}


					/**
					 * Fortran REAL(https://fortranwiki.org/fortran/show/real)
					 */
					template <typename kind=double, typename T_x=int>
						auto real(T_x &&x)
						{
							AERRCHECKREALINTEGERCOMPLEX(T_x,"real","x");
							auto l = [](typename arrayInfo<T_x>::type const &i) -> decltype(auto)
							{
								return static_cast<kind>(i);
							};
							return makeLazy<l>(std::forward<T_x>(x));
						}


					/**
					 * REALLOCATE function. Equivalant to testing for allocated, deallocating if true and then allocating
					 */
					template<typename... T>
						void reallocate(T&& ...items){
							alloc_middle<true>(items...);
						}


					/**
					 * Fortran REDUCE
					 */

					template<typename T, typename T_op>
						auto reduce(T&& array, T_op &&operation, FAR_UNUSED bool ordered=true){
							using atype = arrayInfo<T>::type;
							atype result={};
							bool active = false;

							forEachSerial<[](atype &dest, const atype &src, T_op &op, bool& active){
								if (!active){
									dest = src;
									active = true;
								} else {
									dest = op(dest,src);
								}
							}>(result,array,operation,active);
							if (!active) throw std::out_of_range("Initial sequence is empty");
							return result;
						}

					/**
					 * Fortran REDUCE with mask
					 */
					template<typename T, typename T_op, typename T_mask>
						auto reduce_with_mask(T&& array, T_op &&operation, T_mask &&mask, FAR_UNUSED bool ordered=true){
							using atype = arrayInfo<T>::type;
							using masktype = arrayInfo<T_mask>::type;
							atype result={};
							bool active = false;

							forEachSerial<[](atype &dest, const atype &src, const masktype &mask, T_op &op, bool& active){
								if (mask){
									if (!active){
										dest = src;
										active = true;
									} else {
										dest = op(dest,src);
									}
								}
							}>(result,array,mask,operation,active);
							if (!active) throw std::out_of_range("Initial sequence is empty");
							return result;
						}

					/**
					 * Fortran REDUCE
					 * with identity
					 */
					template<typename T, typename T_op, typename T_identity>
						auto reduce_with_identity(T&& array, T_op &&operation, T_identity &&identity, FAR_UNUSED bool ordered=true){
							using atype = arrayInfo<T>::type;
							static_assert(!arrayInfo<T_identity>::value,"identity must be scalar in reduce_with_identity");
							static_assert(std::is_assignable_v<atype&,std::decay_t<T_identity>>,"identity must be assignable to array in reduce_with_identity");
							atype result={};
							bool active = false;

							forEachSerial<[](atype &dest, const atype &src, T_op &op, bool& active){
								if (!active){
									dest = src;
									active = true;
								} else {
									dest = op(dest,src);
								}
							}>(result,array,operation,active);
							if (!active) result = identity;
							return result;
						}

					/**
					 * Fortran REDUCE
					 * with identity and mask
					 */
					template<typename T, typename T_op, typename T_mask, typename T_identity>
						auto reduce_with_mask_and_identity(T&& array, T_op &&operation, T_mask &&mask, T_identity &&identity, FAR_UNUSED bool ordered=true){
							using atype = arrayInfo<T>::type;
							using mtype = arrayInfo<T_mask>::type;
							static_assert(!arrayInfo<T_identity>::value,"identity must be scalar in reduce_with_mask_and_identity");
							static_assert(std::is_assignable_v<atype&,std::decay_t<T_identity>>,"identity must be assignable to array in reduce_with_mask_and_identity");
							atype result={};
							bool active = false;

							forEachSerial<[](atype &dest, const atype &src, const mtype &mask, T_op &op, bool& active){
								if (mask) {
									if (!active){
										dest = src;
										active = true;
									} else {
										dest = op(dest,src);
									}
								}
							}>(result,array,mask,operation,active);
							if (!active) result = identity;
							return result;
						}


					/**
					 * Fortran REDUCE
					 * with direction
					 */

					template<typename T, typename T_op>
						auto reduce(T&& array, T_op &&operation, int dim, FAR_UNUSED bool ordered=true){
							using atype = arrayInfo<T>::type;
							auto result = reduce_rank<atype>(array, dim);
							auto active = reduce_rank<bool>(array,dim);
							if (array.getSize()==0) throw std::out_of_range("Initial sequence is empty");
							active = false;
							forEachReduceSerial<[](atype &dest, const atype &src, T_op &op, bool& active){
								if (!active){
									dest = src;
									active = true;
								} else {
									dest = op(dest,src);
								}
							}>(dim,result,array,operation,active);
							if (!all(active)) throw std::out_of_range("Initial sequence is empty");
							return result;
						}

					/**
					 * Fortran REDUCE
					 * with direction and identity
					 */

					template<typename T, typename T_op, typename T_identity>
						auto reduce_with_identity(T&& array, T_op &&operation, int dim, T_identity &&identity, FAR_UNUSED bool ordered=true){
							using atype = arrayInfo<T>::type;
							static_assert(!arrayInfo<T_identity>::value,"identity must be scalar in reduce_with_identity");
							static_assert(std::is_assignable_v<atype&,std::decay_t<T_identity>>,"identity must be assignable to array in reduce_with_identity");
							auto result = reduce_rank<atype>(array, dim);
							auto active = reduce_rank<bool>(array,dim);
							active = false;
							forEachReduceSerial<[](atype &dest, const atype &src, T_op &op, bool& active){
								if (!active){
									dest = src;
									active = true;
								} else {
									dest = op(dest,src);
								}
							}>(dim,result,array,operation,active);
							forEachSerial<[](atype &dest, const bool &active, const T_identity &identity){
								if (!active) dest = identity;
							}> (result,active,identity);
							return result;
						}

					/**
					 * Fortran REDUCE
					 * with mask and direction
					 */

					template<typename T, typename T_op, typename T_mask>
						auto reduce_with_mask(T&& array, T_op &&operation, int dim, T_mask &&mask, FAR_UNUSED bool ordered=true){
							using atype = arrayInfo<T>::type;
							using mask_type = arrayInfo<T_mask>::type;
							auto result = reduce_rank<atype>(array, dim);
							auto active = reduce_rank<bool>(array,dim);
							active = false;
							forEachReduceSerial<[](atype &dest, const atype &src, const mask_type &mask, T_op &op, bool& active){
								if (mask) {
									if (!active){
										dest = src;
										active = true;
									} else {
										dest = op(dest,src);
									}
								}
							}>(dim,result,array,mask,operation,active);
							if (!all(active)) throw std::out_of_range("Initial sequence is empty");
							return result;
						}

					/**
					 * Fortran REDUCE
					 * with mask, direction and identity
					 */
					template<typename T, typename T_op, typename T_mask, typename T_identity>
						auto reduce_with_mask_and_identity(T&& array, T_op &&operation, int dim, T_mask&&mask, T_identity &&identity, FAR_UNUSED bool ordered=true){
							using atype = arrayInfo<T>::type;
							using mask_type = arrayInfo<T_mask>::type;
							static_assert(!arrayInfo<T_identity>::value,"identity must be scalar in reduce_with_identity");
							static_assert(std::is_assignable_v<atype&,std::decay_t<T_identity>>,"identity must be assignable to array in reduce_with_identity");
							auto result = reduce_rank<atype>(array, dim);
							auto active = reduce_rank<bool>(array,dim);
							active = false;
							forEachReduceSerial<[](atype &dest, const atype &src, const mask_type &mask, T_op &op, bool& active){
								if (mask){
									if (!active){
										dest = src;
										active = true;
									} else {
										dest = op(dest,src);
									}
								}
							}>(dim,result,array,mask,operation,active);
							forEach<[](atype &dest, const bool &active, const T_identity &identity){
								if (!active) dest = identity;
							}> (result,active,identity);
							return result;
						}

					/**
					 * Take reference to data in an array and return an array of pointers to the data
					 */
					template<typename T>
						auto reference(T &ptr){
							static_assert(!std::is_base_of_v<array_operator,std::decay_t<T>>
								&& !(std::is_pointer_v<typename arrayInfo<T>::type> || std::is_reference_v<typename arrayInfo<T>::type>)
								,"Cannot take references to lazy arrays unless they return references or pointers");
							auto l = [](typename arrayInfo<T>::type &ptr){
								return &ptr;
							};
							return makeLazy<l>(ptr);
						}


					/**
					 * Fortran REPEAT(https://fortranwiki.org/fortran/show/repeat)
					 * Source is std::string
					 */
					template<typename T>
						std::basic_string<T> repeat(const std::basic_string<T> &source, FAR_UNSIGNED_INDEX_TYPE size){
							std::string result;
							result.reserve(source.size()*size);
							for (FAR_UNSIGNED_INDEX_TYPE i=0;i<size;++i){
								result+=source;
							}
							return result;
						}

					/**
					 * Fortran REPEAT(https://fortranwiki.org/fortran/show/repeat)
					 * Source is char*
					 */
					template<typename T>
						std::basic_string<T> repeat(const T *source, FAR_UNSIGNED_INDEX_TYPE size){
							std::string result;
							FAR_UNSIGNED_INDEX_TYPE len = strlen(source);
							result.reserve(len*size);
							for (FAR_UNSIGNED_INDEX_TYPE i=0;i<size;++i){
								result+=source;
							}
							return result;
						}

					/**
					 * Fortran RESHAPE(https://fortranwiki.org/fortran/show/reshape)
					 * using parameters for sizes
					 */
					template<typename T_src, typename... T_sizes>
						auto reshape(T_src &&source, T_sizes... sizes)
						{
							using atype = std::invoke_result_t<decltype(std::decay_t<T_src>::template buildVariant<typename arrayInfo<T_src>::type,sizeof...(T_sizes)>)>;
							static_assert(sizeof...(T_sizes)>0,"reshape must have at least one size parameter");
							atype result(sizes...);
							if (result.getSize() > source.getSize())
								throw std::out_of_range("Unable to reshape array - number of elements after reshape larger than number of elements before reshape");
							if (result.getSize() == source.getSize() && std::is_rvalue_reference_v<decltype(source)> && source.getContiguous() && !source.getAssociated()){
								result.free_core();
								result.shouldDelete=true;
								result.rdata=source.rdata;
								result.bdata=source.bdata;
								source.rdata=nullptr;
								source.rdata=nullptr;
							} else {
								auto its = source.begin();
								auto itr = result.begin();
								auto itse = source.end();
								auto itre = result.end();
								while (its!=itse && itr!=itre){
									(*itr)=(*its);
									itr++;
									its++;
								}
							}
							return result;
						}

					/**
					 * Fortran RESHAPE(https://fortranwiki.org/fortran/show/reshape)
					 * using array of sizes
					 */
					template<typename T_src, typename T_kind, FAR_UNSIGNED_INDEX_TYPE size>
						auto reshape(T_src &&source,std::array<T_kind,size> sizes)
						{
							using atype = std::invoke_result_t<decltype(std::decay_t<T_src>::template buildVariant<typename arrayInfo<T_src>::type,size>)>;
							atype result;
							result.allocate(sizes);
							if (result.getSize() > source.getSize())
								throw std::out_of_range("Unable to reshape array - number of elements after reshape larger than number of elements before reshape");
							if (result.getSize() == source.getSize() && std::is_rvalue_reference_v<decltype(source)> && source.getContiguous() && !source.getAssociated()){
								result.free_core();
								result.shouldDelete = true;
								result.rdata=source.rdata;
								result.bdata=source.bdata;
								source.rdata=nullptr;
								source.rdata=nullptr;
							} else {
								auto its = source.begin();
								forEach<[](arrayInfo<atype>::type &dest, std::decay_t<decltype(its)> &iterator){
									dest = *iterator++;
								}> (result,its);
							}
							return result;
						}

					/**
					 * Fortran RESHAPE(https://fortranwiki.org/fortran/show/reshape)
					 * with pad
					 */
					template<typename T_src, typename T_kind, FAR_UNSIGNED_INDEX_TYPE size, typename T_pad, std::enable_if_t<arrayInfo<T_pad>::value || isStdArray<std::decay_t<T_pad>>::value,int>* =nullptr>
						auto reshape_with_pad(T_src &&source,std::array<T_kind,size> sizes, T_pad &&pad)
						{
							using atype = std::invoke_result_t<decltype(std::decay_t<T_src>::template buildVariant<typename arrayInfo<T_src>::type,size>)>;
							atype result;
							result.allocate(sizes);
							if (result.getSize() == source.getSize() && std::is_rvalue_reference_v<decltype(source)> && source.getContiguous() && !source.getAssociated()){
								result.free_core();
								result.shouldDelete = true;
								result.rdata=source.rdata;
								result.bdata=source.bdata;
								source.rdata=nullptr;
								source.rdata=nullptr;
							} else {
								auto its = source.begin();
								auto ite = source.end();
								auto pts = pad.begin();
								auto ptc = pad.begin();
								auto pte = pad.end();
								forEach<[](arrayInfo<atype>::type &dest, std::decay_t<decltype(its)> &iterator, std::decay_t<decltype(its)> &end, 
										std::decay_t<decltype(pts)> &padstart, std::decay_t<decltype(ptc)> &padit, std::decay_t<decltype(pte)> &padend){
									if (iterator<end){
										dest = *iterator++;
									} else {
										dest = *padit++;
									}
									if (padit==padend) padit=padstart;
								}> (result,its,ite,pts,ptc,pte);
							}
							return result;
						}

					/**
					 * Fortran RESHAPE(https://fortranwiki.org/fortran/show/reshape)
					 * with pad (scalar pad, extension to Fortran)
					 */
					template<typename T_src, typename T_kind, FAR_UNSIGNED_INDEX_TYPE size, typename T_pad, std::enable_if_t<!arrayInfo<T_pad>::value & !isStdArray<std::decay_t<T_pad>>::value,int>* =nullptr>
						auto reshape_with_pad(T_src &&source,std::array<T_kind,size> sizes, T_pad &&pad)
						{
							using atype = std::invoke_result_t<decltype(std::decay_t<T_src>::template buildVariant<typename arrayInfo<T_src>::type,size>)>;
							atype result;
							result.allocate(sizes);
							if (result.getSize() == source.getSize() && std::is_rvalue_reference_v<decltype(source)> && source.getContiguous() && !source.getAssociated()){
								result.free_core();
								result.shouldDelete = true;
								result.rdata=source.rdata;
								result.bdata=source.bdata;
								source.rdata=nullptr;
								source.rdata=nullptr;
							} else {
								auto its = source.begin();
								auto ite = source.end();
								forEach<[](arrayInfo<atype>::type &dest, std::decay_t<decltype(its)> &iterator, std::decay_t<decltype(its)> &end, std::decay_t<decltype(pad)> &pad){
									if (iterator<end){
										dest = *iterator++;
									} else {
										dest = pad;
									}
								}> (result,its,ite,pad);
							}
							return result;
						}

					/**
					 * Fortran RESHAPE(https://fortranwiki.org/fortran/show/reshape)
					 * with order
					 */
					template<typename T_src, typename T_kind, FAR_UNSIGNED_INDEX_TYPE size, typename T_order>
						auto reshape_with_order(T_src &&source,std::array<T_kind,size> sizes, T_order &&order)
						{
							constexpr static int rank = size;
							using atype = std::invoke_result_t<decltype(std::decay_t<T_src>::template buildVariant<typename arrayInfo<T_src>::type,size>)>;
							FAR_UNSIGNED_INDEX_TYPE orderedSizes[rank];
							auto getOrder = [&order](int i){
								if constexpr(arrayInfo<T_order>::value){
									return order.getItem(i)-1;
								} else {
									return order[i]-1;
								}
							};
							for (int i=0;i<rank;++i){
								orderedSizes[i] = sizes[getOrder(i)];
							}
							typename atype::index_type it(orderedSizes); // Create an indexer with the indices in incrementing order so that I can get and advance the indices
							typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type otuple, dtuple;
							bool finished=false;
							atype result;
							result.allocate(sizes);
							if (result.getSize() > source.getSize())
								throw std::out_of_range("Unable to reshape array - number of elements after reshape larger than number of elements before reshape");
							for (int i=0;i<rank;++i) getTupleLevel(otuple,i) = defaultLB;
							auto sit = source.begin();
							while(!finished){
								//Descramble the index tuple back into the real array
								for (int i=0;i<rank;++i){
									getTupleLevel(dtuple,getOrder(i)) = getTupleLevel(otuple,i);
								}
								//Copy the data across
								callDefaultBase(result,dtuple) = *sit++;
								//Move the simple incrementing tuple
								finished = it.advanceTuple(otuple);
							}
							return result;
						}

					/**
					 * Fortran RESHAPE(https://fortranwiki.org/fortran/show/reshape)
					 * with pad and order
					 */
					template<typename T_src, typename T_kind, FAR_UNSIGNED_INDEX_TYPE size, typename T_order, typename T_pad, std::enable_if_t<arrayInfo<T_pad>::value || isStdArray<std::decay_t<T_pad>>::value,int>* =nullptr>
						auto reshape_with_order_and_pad(T_src &&source,std::array<T_kind,size> sizes, T_order &&order, T_pad &&pad)
						{
							constexpr static int rank = size;
							using atype = std::invoke_result_t<decltype(std::decay_t<T_src>::template buildVariant<typename arrayInfo<T_src>::type,size>)>;
							FAR_UNSIGNED_INDEX_TYPE orderedSizes[rank];
							auto getOrder = [&order](int i){
								if constexpr(arrayInfo<T_order>::value){
									return order.getItem(i)-1;
								} else {
									return order[i]-1;
								}
							};
							for (int i=0;i<rank;++i){
								orderedSizes[i] = sizes[getOrder(i)];
							}
							typename atype::index_type it(orderedSizes);
							typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type otuple, dtuple;
							bool finished=false;
							atype result;
							result.allocate(sizes);
							for (int i=0;i<rank;++i) getTupleLevel(otuple,i) = defaultLB;
							auto sit = source.begin();
							auto site = source.end();
							auto pst = pad.begin();
							auto pit = pad.begin();
							auto pte = pad.end();
							while(!finished){
								for (int i=0;i<rank;++i){
									getTupleLevel(dtuple,getOrder(i)) = getTupleLevel(otuple,i);
								}
								if (sit<site){
									callDefaultBase(result,dtuple) = *sit++;
								} else {
									callDefaultBase(result,dtuple) = *pit++;
								}
								if (pit==pte) pit=pst;
								finished = it.advanceTuple(otuple);
							}
							return result;
						}

					/**
					 * Fortran RESHAPE(https://fortranwiki.org/fortran/show/reshape)
					 * with pad and order. Scalar pad (extension to Fortran)
					 */
					template<typename T_src, typename T_kind, FAR_UNSIGNED_INDEX_TYPE size, typename T_order, typename T_pad, std::enable_if_t<!arrayInfo<T_pad>::value & !isStdArray<std::decay_t<T_pad>>::value,int>* =nullptr>
						auto reshape_with_order_and_pad(T_src &&source,std::array<T_kind,size> sizes, T_order &&order, T_pad &&pad)
						{
							constexpr static int rank = size;
							using atype = std::invoke_result_t<decltype(std::decay_t<T_src>::template buildVariant<typename arrayInfo<T_src>::type,size>)>;
							FAR_UNSIGNED_INDEX_TYPE orderedSizes[rank];
							auto getOrder = [&order](int i){
								if constexpr(arrayInfo<T_order>::value){
									return order.getItem(i)-1;
								} else {
									return order[i]-1;
								}
							};
							for (int i=0;i<rank;++i){
								orderedSizes[i] = sizes[getOrder(i)];
							}
							typename atype::index_type it(orderedSizes);
							typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type otuple, dtuple;
							bool finished=false;
							atype result;
							result.allocate(sizes);
							for (int i=0;i<rank;++i) getTupleLevel(otuple,i) = defaultLB;
							auto sit = source.begin();
							auto site = source.end();
							while(!finished){
								for (int i=0;i<rank;++i){
									getTupleLevel(dtuple,getOrder(i)) = getTupleLevel(otuple,i);
								}
								if (sit<site){
									callDefaultBase(result,dtuple) = *sit++;
								} else {
									callDefaultBase(result,dtuple) = pad;
								}
								finished = it.advanceTuple(otuple);
							}
							return result;
						}

					/**
					 * Fortran RRSPACING(https://fortranwiki.org/fortran/show/rrspacing)
					 */
					template <typename T_x>
						auto rrspacing(T_x &&x)
						{
							AERRCHECKREAL(T_x,"rrspacing","x")
								using type = arrayInfo<T_x>::type;
							auto l = [](type const &x) -> decltype(auto)
							{ 
								int exponent;
								auto radix = std::numeric_limits<type>::radix;
								auto epsilon = std::numeric_limits<type>::epsilon();
								type fraction = std::frexp(x, &exponent);
								return std::abs(fraction) * radix/epsilon;
							};
							return makeLazy<l>(std::forward<T_x>(x));
						}

					/**
					 * Fortran SAME_TYPE_AS(https://fortranwiki.org/fortran/show/same_type_as)
					 */
					template<typename T1, typename T2>
						bool constexpr same_type_as(FAR_UNUSED T1 &&a, FAR_UNUSED T2 &&b){
							return std::is_same_v<std::decay_t<T1>,std::decay_t<T2>>;
						}

					/**
					 * Fortran SCALE(https://fortranwiki.org/fortran/show/scale)
					 */
					template <typename T_x, typename T_i>
						auto scale(T_x &&x, T_i &&i)
						{ 
							AERRCHECKREAL(T_x,"scale","x");
							ERRCHECKINTEGER(T_i,"scale","i");
							auto l = [](typename arrayInfo<T_x>::type const &x,
									typename arrayInfo<T_i>::type const &i) -> decltype(auto)
							{
								return std::scalbn(x,i);
							};
							return makeLazy<l>(std::forward<T_x>(x), std::forward<T_i>(i));
						}

					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T_str=int, typename T_set=int, typename T_back=int>
						auto scan(T_str &&string, T_set &&set, T_back &&back){
							AERRCHECKSTRINGY(T_str,"scan","string");
							AERRCHECKSTRINGY(T_set,"scan","set");
							AERRCHECKTRUTHY(T_back,"scan","back");
							using type = typename arrayInfo<T_str>::type;
							using setType = typename arrayInfo<T_set>::type;
							using backtype = typename arrayInfo<T_back>::type;
							auto l = [](const type &string, const setType &set, const backtype &back){
								kind result=0;
								using strtype = std::conditional_t<is_stdstring_v<type>,type,std::basic_string<std::remove_cv_t<std::remove_pointer_t<type>>>>;
								using settype = std::conditional_t<is_stdstring_v<setType>,setType,std::basic_string<std::remove_cv_t<std::remove_pointer_t<setType>>>>;
								strtype sstr = string;
								settype sset = set;
								if (sset.size()==0) return result;
								if (back){
									auto pos = sstr.find_last_of(sset);
									if (pos!=std::string::npos){
										result = pos+1;
									}
								} else {
									auto pos = sstr.find_first_of(sset);
									if (pos!=std::string::npos){
										result = pos+1;
									}
								}
								return result;
							};

							return makeLazy<l>(std::forward<T_str>(string), std::forward<T_set>(set), std::forward<T_back>(back));
						}

					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T_str=int, typename T_set=int>
						auto scan(T_str &&string, T_set &&set){
							AERRCHECKSTRINGY(T_str,"scan","string");
							AERRCHECKSTRINGY(T_set,"scan","set");
							return scan<kind>(std::forward<T_str>(string),std::forward<T_set>(set),false);
						}


					namespace{
						// Helper struct to map the number of digits to the corresponding integer type
						template<int Digits, typename T = void>
							struct int_kind_helper;

						template<int Digits>
							struct int_kind_helper<Digits, typename std::enable_if<(Digits <= std::numeric_limits<int8_t>::digits10)>::type> {
#ifdef FAR_EXACT_INT_KIND
								using type = int8_t;
#else
								using type = int_fast8_t;
#endif
							};

						template<int Digits>
							struct int_kind_helper<Digits, typename std::enable_if<(Digits > std::numeric_limits<int8_t>::digits10) && (Digits <= std::numeric_limits<int16_t>::digits10)>::type> {
#ifdef FAR_EXACT_INT_KIND
								using type = int16_t;
#else
								using type = int_fast16_t;
#endif

							};

						template<int Digits>
							struct int_kind_helper<Digits, typename std::enable_if<(Digits > std::numeric_limits<int16_t>::digits10) && (Digits <= std::numeric_limits<int32_t>::digits10)>::type> {
#ifdef FAR_EXACT_INT_KIND
								using type = int32_t;
#else
								using type = int_fast32_t;
#endif

							};

						template<int Digits>
							struct int_kind_helper<Digits, typename std::enable_if<(Digits > std::numeric_limits<int32_t>::digits10) && (Digits <= std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::digits10)>::type> {
#ifdef FAR_EXACT_INT_KIND
								using type = FAR_SIGNED_INDEX_TYPE;
#else
								using type = int_fast64_t;
#endif

							};

#ifdef FAR_INT128
						template<int Digits>
							struct int_kind_helper<Digits, typename std::enable_if<(Digits > std::numeric_limits<FAR_SIGNED_INDEX_TYPE>::digits10) && (Digits <= std::numeric_limits<FAR_INT128>::digits10)>::type> {
								using type = FAR_INT128;
							};	
#endif

						template<int P=0, int R=0, typename T=void>
							struct real_kind_helper;

						template<int P, int R>
							struct real_kind_helper<P, R, typename std::enable_if<(P <= std::numeric_limits<float>::digits10) && (R <= std::numeric_limits<float>::max_exponent10-1)>::type> {
								using type = float;
							};

						template<int P, int R>
							struct real_kind_helper<P, R, typename std::enable_if<((P <= std::numeric_limits<double>::digits10) && (R <= std::numeric_limits<double>::max_exponent10-1)) &&
							(P > std::numeric_limits<float>::digits10 || R> std::numeric_limits<float>::max_exponent10-1)>::type> {
								using type = double;
							};

						template<int P, int R>
							struct real_kind_helper<P, R, typename std::enable_if<
							((std::numeric_limits<double>::digits10 < (std::numeric_limits<long double>::digits10) || std::numeric_limits<double>::max_exponent10 < std::numeric_limits<long double>::max_exponent10) &&
							 (P <= std::numeric_limits<long double>::digits10) && (R <= std::numeric_limits<long double>::max_exponent10-1)) &&
							(P > std::numeric_limits<double>::digits10 || R> std::numeric_limits<double>::max_exponent10-1)>::type> {
								using type = long double;
							};

					};

					template<int R>
						using selected_int_kind = typename int_kind_helper<R>::type;

					template<int R>
						using selected_uint_kind = std::make_unsigned_t<typename int_kind_helper<R>::type>;

					template<int P, int R>
						using selected_real_kind = typename real_kind_helper<P,R>::type;

					template<int P>
						using selected_real_kind_p = typename real_kind_helper<P,0>::type;

					template<int R>
						using selected_real_kind_r = typename real_kind_helper<0,R>::type;

					template<int P, int R>
						using selected_complex_kind = std::complex<selected_real_kind<P,R>>;

					template<int P>
						using selected_complex_kind_p = std::complex<selected_real_kind_p<P>>;

					template<int R>
						using selected_complex_kind_r = std::complex<selected_real_kind_r<R>>;

					/**
					 * Fortran SET_EXPONENT(https://fortranwiki.org/fortran/show/set_exponent)
					 */
					template <typename T_x, typename T_i>
						auto set_exponent(T_x &&x, T_i &&i)
						{
							AERRCHECKREAL(T_x,"set_exponent","x");
							AERRCHECKINTEGER(T_i,"set_exponent","i");
							using xtype = arrayInfo<T_x>::type;
							using itype = arrayInfo<T_i>::type;
							auto l = [](xtype const &x,
									itype const &i) -> decltype(auto)
							{
								auto f = fraction(x);
								return std::scalbn(f,i);
							};
							return makeLazy<l>(std::forward<T_x>(x), std::forward<T_i>(i));
						}

					/**
					 * Fortran SHAPE(https://fortranwiki.org/fortran/show/shape)
					 */
					template <typename kind = FAR_SIGNED_INDEX_TYPE, typename T_source=int>
						auto shape(T_source &&source)
						{
							Array<kind,1> result;
							if constexpr(!arrayInfo<T_source>::value){
								result.allocate(0);
								return result;
							} else {
								constexpr int rank = arrayInfo<T_source>::rank;
								result.allocate(rank);
								for (int i=1;i<=rank;++i){
									result(i)=source.getRankSize(i);
								}
								return result;
							}
						}

					/**
					 * Fortran SHIFTA(https://fortranwiki.org/fortran/show/shifta)
					 */
					template <typename T_i, typename T_shift>
						auto shifta(T_i &&i, T_shift &&shift)
						{
							AERRCHECKINTEGER(T_i,"shifta","i");
							AERRCHECKINTEGER(T_shift,"shifta","shift");
							using type = typename arrayInfo<T_i>::type;
							using shiftType = typename arrayInfo<T_shift>::type;
							//Algebraic shift for both signed and unsigned integers
							auto l =[](type const &i, shiftType const &shift) -> decltype(auto)
							{
								if constexpr(std::is_signed_v<type>){
									//Signed integers already have algebraic shift
									return i>>shift;
								} else {
									//This could be done correctly, but at the moment, just bit
									//cast to signed, do the shift and bit cast back
									std::make_signed_t<type> si = std::bit_cast<std::make_signed_t<type>>(i);
									si >>= shift;
									return std::bit_cast<type>(si);
								}
							};
							return makeLazy<l>(std::forward<T_i>(i),std::forward<T_shift>(shift));
						}

					/**
					 * Fortran SHIFTL(https://fortranwiki.org/fortran/show/shiftl)
					 */
					template <typename T_i, typename T_shift>
						auto shiftl(T_i &&i, T_shift &&shift)
						{
							AERRCHECKINTEGER(T_i,"shiftl","i");
							AERRCHECKINTEGER(T_shift,"shiftl","shift");
							using type = typename arrayInfo<T_i>::type;
							using shiftType = typename arrayInfo<T_shift>::type;
							return makeLazy<shiftl_core<type,shiftType>>(i,shift);
						}

					/**
					 * Fortran SHIFTR(https://fortranwiki.org/fortran/show/shiftr)
					 */
					template <typename T_i, typename T_shift>
						auto shiftr(T_i &&i, T_shift &&shift)
						{
							AERRCHECKINTEGER(T_i,"shiftr","i");
							AERRCHECKINTEGER(T_shift,"shiftr","shift");
							using type = typename arrayInfo<T_i>::type;
							using shiftType = typename arrayInfo<T_shift>::type;
							return makeLazy<shiftr_core<type,shiftType>>(i,shift);
						}

					/**
					 * Fortran SIGN(https://fortranwiki.org/fortran/show/sign)
					 */
					template <typename T_a, typename T_b>
						auto sign(T_a &&a, T_b &&b)
						{
							AERRCHECKREALINTEGERCOMPLEX(T_a,"sign","a");
							AERRCHECKREALINTEGERCOMPLEX(T_b,"sign","b");
							using type_a = typename arrayInfo<T_a>::type;
							using type_b = typename arrayInfo<T_b>::type;
							if constexpr(!std::is_floating_point_v<type_b>){
								return makeLazy<[](const type_a &a, const type_b &b){return b>=0?std::abs(a):-std::abs(a);}>(a,b);
							} else {
								//For floats/doubles must use signbit
								//otherwise -0.0 will be compared as zero
								return makeLazy<[](const type_a &a, const type_b &b){return std::signbit(b)?-std::abs(a):std::abs(a);}>(a,b);
							}
						}

					/**
					 * Fortran SIN(https://fortranwiki.org/fortran/show/sin)
					 */
					template <typename T_x>
						auto sin(T_x &&x)
						{
							AERRCHECKREALCOMPLEX(T_x,"sin","x");
							using type = typename arrayInfo<T_x>::type;
							return makeLazy<[](const type&x){return std::sin(x);}>(std::forward<T_x>(x));
						}

					/**
					 * Fortran SIND(sin in degrees)
					 */
					template <typename T_x>
						auto sind(T_x &&x)
						{
							AERRCHECKREALCOMPLEX(T_x,"sind","x");
							using type = typename arrayInfo<T_x>::type;
							return makeLazy<[](const type&x){return std::sin(x/180.0*M_PI);}>(std::forward<T_x>(x));
						}

					/**
					 * Fortran SINH(https://fortranwiki.org/fortran/show/sinh)
					 */
					template <typename T_x>
						auto sinh(T_x &&x)
						{
							AERRCHECKREALCOMPLEX(T_x,"sinh","x")
								using type = typename arrayInfo<T_x>::type;
							return makeLazy<[](const type&x){return std::sinh(x);}>(std::forward<T_x>(x));
						}

					/**
					 * Fortran SINPI(https://fortranwiki.org/fortran/show/sinpi)
					 */
					template <typename T_x>
						auto sinpi(T_x &&x)
						{
							AERRCHECKREALCOMPLEX(T_x,"sinpi","x");
							using type = typename arrayInfo<T_x>::type;
							return makeLazy<[](const type&x){return std::sin(x*M_PI);}>(std::forward<T_x>(x));
						}

					/**
					 * Fortran SIZE(https://fortranwiki.org/fortran/show/size)
					 */
					template <typename kind=FAR_SIGNED_INDEX_TYPE, typename T_array>
						kind size(T_array &&array)
						{
							static_assert(arrayInfo<T_array>::value,"size must be called with an array");
							return array.getSize();
						}

					/**
					 * Fortran SIZE(https://fortranwiki.org/fortran/show/size)
					 * with dim
					 */
					template <typename kind=FAR_SIGNED_INDEX_TYPE, typename T_array>
						kind size(T_array &&array, int dim)
						{
							static_assert(arrayInfo<T_array>::value,"size must be called with an array");
							return array.getRankSize(dim);
						}

					/**
					 * Fortran SNGL(https://fortranwiki.org/fortran/show/sngl)
					 * with dim
					 */
					template <typename T_a>
						auto sngl(T_a &&a)
						{
							AERRCHECKREALINTEGER(T_a,"sngl","a");
							using type = typename arrayInfo<T_a>::type;
							return makeLazy<[](const type& a){return float(a);}>(a);
						}

					/**
					 * Fortran SPACING(https://fortranwiki.org/fortran/show/spacing)
					 */
					template <typename T_x>
						auto spacing(T_x &&x)
						{
							AERRCHECKREAL(T_x,"spacing","x");
							using type = typename arrayInfo<T_x>::type;
							auto l = [](type const &x) -> decltype(auto)
							{
								if (x==std::numeric_limits<type>::infinity() || x==-std::numeric_limits<type>::infinity()){
									return std::numeric_limits<type>::quiet_NaN();
								}
								if (std::isnan(x)){
									return std::numeric_limits<type>::quiet_NaN();
								}
								if (x==0.0){
									return std::numeric_limits<type>::min();
								}
								if (x>0) {
									return std::nextafter(x, std::numeric_limits<type>::max()) - x;
								} else {
									return x - std::nextafter(x, std::numeric_limits<type>::lowest());
								}
							};
							return makeLazy<l>(std::forward<T_x>(x));
						}

					/**
					 * Fortran SPREAD(https://fortranwiki.org/fortran/show/spread)
					 */
					template<typename T_source>
						auto spread(T_source&& source, int dim, FAR_SIGNED_INDEX_TYPE copies)
						{
							constexpr int rank = arrayInfo<T_source>::rank;
							using coreType = arrayInfo<T_source>::type;
							if constexpr(rank==0){
								contiguousIndexBuilder<1> fib(copies);
								return lazyArray<[](const coreType &t){return t;},contiguousIndexBuilder<1>, coreType>(fib,source);
							} else {
								using indexType = contiguousIndexBuilder<rank+1>;
								using indexInfoType = indexInfo<rank+1>;
								int originIndex=1;
								FAR_UNSIGNED_INDEX_TYPE sizes[rank+1];
								for (int i=0;i<rank+1;++i){
									sizes[i] = ((i+1)==dim)?copies:source.getRankSize(originIndex++);
								}
								indexType ib(sizes);
								auto l = [](indexInfoType info,
										const std::decay_t<T_source> &source,
										int dim) -> auto
								{
									typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type caller;
									reduceTuple(info.indices, caller, dim);
									return std::apply([&source](auto&...params){return source.defaultBase(params...);},caller);
								};
								auto ret = lazyArray<l, indexType, std::decay_t<T_source>, int>(ib,source,dim-1);
								return ret;
							}
						}

					/**
					 * Fortran SQRT(https://fortranwiki.org/fortran/show/sqrt)
					 */
					template <typename T_x>
						auto sqrt(T_x &&x)
						{
							AERRCHECKREALCOMPLEX(T_x,"sqrt","x");
							using type = typename arrayInfo<T_x>::type;
							return makeLazy<[](const type& x){return std::sqrt(x);}>(std::forward<T_x>(x));
						}

					/**
					 * Fortran STORAGE_SIZE(https://fortranwiki.org/fortran/show/storage_size)
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T>
						kind storage_size(FAR_UNUSED T &&t){
							using type = typename arrayInfo<T>::type;
							return sizeof(type) * CHAR_BIT;
						}

					/**
					 * Fortran SUM(https://fortranwiki.org/fortran/show/sum)
					 */
					template <typename T_array>
						auto sum(T_array &&array)
						{
							ERRCHECKREALINTEGERCOMPLEX(T_array,"sum","array");
							using type = typename arrayInfo<T_array>::type;
							return reduction<[](const type &src, type &dest){dest+=src;},
										 [](const type &src, type &dest){dest+=src;},
										 [](type &dest){dest=0.0;}
							>(array);
						}

					/**
					 * Fortran SUM(https://fortranwiki.org/fortran/show/sum)
					 * with mask
					 */
					template <typename T_array, typename T_mask>
						auto sum_with_mask(T_array &&array, T_mask &&mask)
						{
							ERRCHECKREALINTEGERCOMPLEX(T_array,"sum","array");
							AERRCHECKTRUTHY(T_mask,"sum","mask");
							using type = typename arrayInfo<std::decay_t<T_array>>::type;
							using typeMask = typename arrayInfo<std::decay_t<T_mask>>::type;
							return reduction<[](const type &src, type &dest, typeMask &mask){dest+=mask?src:0;},
										 [](const type &src, type&dest){dest+=src;},
										 [](type &dest){dest=0.0;}
							>(array,mask);
						}

					/**
					 * Fortran SUM(https://fortranwiki.org/fortran/show/sum)
					 * with dim
					 */
					template <typename T_array>
						auto sum(T_array &&array, int dim)
						{
							ERRCHECKREALINTEGERCOMPLEX(T_array,"sum","array");
							using type = typename arrayInfo<std::decay_t<T_array>>::type;
							return reductionWithDirection<[](const type &src, type &dest){dest+=src;}>(dim,array);
						}

					template <typename T_array>
						auto sumSerial(T_array &&array, int dim)
						{
							ERRCHECKREALINTEGERCOMPLEX(T_array,"sum","array");
							using type = typename arrayInfo<std::decay_t<T_array>>::type;
							auto result = reduce_rank(array, dim);
							forEachSerial<[](type &src){src={};}>(result);
							forEachReduceSerial<[](type &src, type &dest){dest+=src;}>(dim,array,result);
							return result;
						}

					/**
					 * Fortran SUM(https://fortranwiki.org/fortran/show/sum)
					 * with dim and mask
					 */
					template <typename T_array, typename T_mask>
						auto sum_with_mask(T_array &&array, int dim, T_mask &&mask)
						{
							ERRCHECKREALINTEGERCOMPLEX(T_array,"sum","array");
							AERRCHECKTRUTHY(T_mask,"sum","mask");
							using type = typename arrayInfo<std::decay_t<T_array>>::type;
							using typeMask = typename arrayInfo<std::decay_t<T_mask>>::type;
							// Produce an array with the summed direction removed
							auto result = reduce_rank(array, dim);
							forEach<[](type &src){src={};}>(result);
							forEachReduce<[](const type &src, type &dest, const typeMask&mask ){dest+=mask?src:0;}>(dim,array,result,mask);
							return result;
						}

					/**
					 * Fortran TAN(https://fortranwiki.org/fortran/show/tan)
					 */
					template <typename T_x>
						auto tan(T_x &&x)
						{
							AERRCHECKREALCOMPLEX(T_x,"tan","x");
							using type = typename arrayInfo<T_x>::type;
							return makeLazy<[](const type& x){return std::tan(x);}>(std::forward<T_x>(x));
						}

					/**
					 * Fortran TAND(https://fortranwiki.org/fortran/show/tand)
					 */
					template <typename T_x>
						auto tand(T_x &&x)
						{
							AERRCHECKREALCOMPLEX(T_x,"tand","x");
							using type = typename arrayInfo<T_x>::type;
							return makeLazy<[](const type& x){return std::tan(x/180.0*M_PI);}>(std::forward<T_x>(x));
						}

					/**
					 * Fortran TANH(https://fortranwiki.org/fortran/show/tanh)
					 */
					template <typename T_x>
						auto tanh(T_x &&x)
						{
							AERRCHECKREALCOMPLEX(T_x,"tanh","x");
							using type = typename arrayInfo<T_x>::type;
							return makeLazy<[](const type& x){return std::tanh(x);}>(std::forward<T_x>(x));
						}

					/**
					 * Fortran TANPI
					 */
					template <typename T_x>
						auto tanpi(T_x &&x)
						{
							AERRCHECKREALCOMPLEX(T_x,"tanpi","x");
							using type = typename arrayInfo<T_x>::type;
							return makeLazy<[](const type& x){return std::tan(x*M_PI);}>(std::forward<T_x>(x));
						}

					/**
					 * Fortran TINY(https://fortranwiki.org/fortran/show/tiny)
					 */
					template <typename T_x>
						auto tiny(FAR_UNUSED T_x &&x)
						{
							return std::numeric_limits<typename arrayInfo<T_x>::type>::min();
						}

					/**
					 * Fortran TRAILZ(https://fortranwiki.org/fortran/show/trailz)
					 */
					template <typename T_i>
						auto trailz(T_i &&i)
						{ AERRCHECKINTEGER(T_i,"trailz","i");
							using type = typename arrayInfo<T_i>::type;
							return makeLazy<[](type x){
								if (x==0) return (int)(CHAR_BIT * sizeof(x));
								int count = 0;
								while (!(x & 1) && count < static_cast<int>(CHAR_BIT * sizeof(x))){
									x>>=1;
									count++;
								}
								return count;
							}>(std::forward<T_i>(i));
						}

					/**
					 * Fortran TRANSFER(https://fortranwiki.org/fortran/show/transfer)
					 */
					template <typename T_source, typename T_mold>
						auto transfer(T_source &&source, T_mold &&mold)
						{
							using sourceType = arrayInfo<T_source>::type;
							using moldType = arrayInfo<T_mold>::type;

							if constexpr(arrayInfo<T_mold>::value){
								//Mold is an array, so result is an array large enough to store all of source
								if constexpr(arrayInfo<T_source>::value){
									//Source is also an array
									Array<moldType,1> result(source.getSize() * FAR_SIGNED_INDEX_TYPE(std::ceil(float(sizeof(sourceType))/float(sizeof(moldType)))));
									if (source.getContiguous()){
										//Source is contiguous, so can just copy directly 
										std::memcpy((char*)result.data(),(char*)source.data(),source.getSize() * sizeof(sourceType));
										return result;
									} else {
										//Source is non contiguous, so have to move it element by element (could in theory be done in contiguous
										//blocks, but I'm not doing that right now)
										char* dest = (char*)(result.data());
										forEach<[](const sourceType &src, char**dest){
											std::memcpy(*dest,(char*)&src,sizeof(sourceType));
											*dest+=sizeof(sourceType);
										}>(source,&dest);
										return result;
									}
								} else {
									//Source is scalar
									Array<moldType,1> result(FAR_SIGNED_INDEX_TYPE(std::ceil(float(sizeof(sourceType))/float(sizeof(moldType)))));
									std::memcpy((char*)result.data(),(char*)&source,sizeof(source));
									return result;
								}
							} else {
								//Mold is not an array
								if constexpr(arrayInfo<T_source>::value){
									//Source is an array
									FAR_UNSIGNED_INDEX_TYPE copySize = std::min(sizeof(mold),sizeof(sourceType)*source.getSize());
									if (source.getContiguous()){
										moldType result{};
										std::memcpy((char*)(&result),(char*)(source.data()),copySize);
										return result;
									} else {
										moldType result{};
										char *dest = (char*)(&result);
										forEach<[](const sourceType &src, char**dest, FAR_UNSIGNED_INDEX_TYPE &remaining){
											if (((remaining) > 0)){
												std::memcpy(*dest,(char*)&src,sizeof(sourceType));
												*dest+=sizeof(sourceType);
												remaining -= sizeof(sourceType);
												return iaContinue;
											} else {
												return iaTerminate;
											}
										}>(source, &dest, copySize);
										return result;
									}
								} else {
									//Neither source nor mold are arrays
									//Simply copy being careful not to buffer overrun on either side
									moldType result{};
									FAR_UNSIGNED_INDEX_TYPE copySize = std::min(sizeof(moldType),sizeof(sourceType));
									std::memcpy((char*)&result,(char*)&source,copySize);
									return result;
								}
							}
						}


					/**
					 * Fortran TRANSFER(https://fortranwiki.org/fortran/show/transfer)
					 */
					template <typename T_source, typename T_mold>
						auto transfer(T_source &&source, FAR_UNUSED T_mold &&mold, FAR_SIGNED_INDEX_TYPE size)
						{
							using sourceType = arrayInfo<T_source>::type;
							using moldType = arrayInfo<T_mold>::type;

							if constexpr(arrayInfo<T_source>::value){
								//Source is also an array
								FAR_SIGNED_INDEX_TYPE sourceSize = source.getSize() * FAR_SIGNED_INDEX_TYPE(std::ceil(float(sizeof(sourceType))/float(sizeof(moldType))));
								FAR_SIGNED_INDEX_TYPE copySize = std::min(FAR_SIGNED_INDEX_TYPE(sourceSize*sizeof(sourceType)),FAR_SIGNED_INDEX_TYPE(size*sizeof(moldType)));
								Array<moldType,1> result(size);
								if (source.getContiguous()){
									//Source is contiguous, so can just copy directly
									std::memcpy((char*)result.data(),(char*)source.data(),copySize);
									return result;
								} else {
									//Source is non contiguous, so have to move it element by element (could in theory be done in contiguous
									//blocks, but I'm not doing that right now)
									char* dest = (char*)(result.data());
									forEach<[](const sourceType &src, char**dest, FAR_SIGNED_INDEX_TYPE &remaining){
										if (((remaining) > 0)){
											std::memcpy(*dest,(char*)&src,sizeof(sourceType));
											*dest+=sizeof(sourceType);
											remaining -= sizeof(sourceType);
											return iaContinue;
										} else {
											return iaTerminate;
										}
									}>(source,&dest,copySize);
									return result;
								}
							} else {
								//Source is scalar
								Array<moldType,1> result(FAR_SIGNED_INDEX_TYPE(std::ceil(float(sizeof(sourceType))/float(sizeof(moldType)))));
								std::memcpy((char*)result.data(),(char*)&source,sizeof(source));
								return result;
							}
						}

					/**
					 * Fortran TRANSPOSE(https://fortranwiki.org/fortran/show/transpose)
					 */
					template<typename T_matrix>
						auto transpose(T_matrix&& source)
						{
							constexpr int rank = arrayInfo<T_matrix>::rank;
							static_assert(rank==2, "Rank must be 2 in transpose");
							using indexInfoType = indexInfo<rank>;

							auto l = [](indexInfoType info,
									const std::decay_t<T_matrix> &source) -> auto
							{
								return source.defaultBase(std::get<1>(info.indices),std::get<0>(info.indices));
							};
							typename arrayInfo<T_matrix>::index ib(source.getRankSize(2),source.getRankSize(1));
							return lazyArray<l,typename arrayInfo<T_matrix>::index,T_matrix>(ib,std::forward<T_matrix>(source));
						}


					/**
					 * Fortran UBOUND(https://fortranwiki.org/fortran/show/ubound)
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T_i=double>
						auto ubound(T_i &&array){
#ifdef RANK1_INQUIRY_IS_VALUE
							if constexpr(arrayInfo<T_i>::rank==1) {
								kind result = array.getUxB(1);
								return result;
							}
							else {
#endif
								constexpr int rank = arrayInfo<T_i>::rank;
								Array<kind,1> result(rank);
								for (int i=1;i<=rank;++i){
									result(i) = array.getUB(i);
								}
								return result;
#ifdef RANK1_INQUIRY_IS_VALUE
							}
#endif
						}

					/**
					 * Fortran UBOUND(https://fortranwiki.org/fortran/show/ubound)
					 * with dim
					 */
					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T_i=double>
						auto ubound(T_i &&array, int dim){
							kind result = array.getUB(dim);
							return result;
						}

					/**
					 * Fortran UNPACK(https://fortranwiki.org/fortran/show/unpack)
					 * With mask parameter and vector parameter
					 */
					template <typename T_vector, typename T_mask, typename T_field>
						auto unpack(T_vector &&vector, T_mask &&mask, T_field &&field)
						{
							using type = arrayInfo<T_vector>::type;
							using mask_type = arrayInfo<T_mask>::type;
							using field_type = arrayInfo<T_field>::type;
							constexpr int rank = arrayInfo<T_mask>::rank;
							static_assert(arrayInfo<T_vector>::rank==1,"vector must be rank 1 in unpack");
							auto unpacked = vector.template buildVariant<type,rank>();
							//If mask is an array then count true elements
							//otherwise use the size of the array
							FAR_SIGNED_INDEX_TYPE items = vector.getSize();
							unpacked.mold(mask);
							if (items==0){
								unpacked = field;
								return unpacked;
							}
							FAR_SIGNED_INDEX_TYPE index=vector.getLB(1);

							auto mapper = [] (const indexInfo<rank> &info, const T_vector &vector, const T_mask &mask, const T_field &field, 
									decltype(unpacked) &dest, FAR_SIGNED_INDEX_TYPE &index){
								return std::make_tuple(std::ref(vector), getItem(mask,info), getItem(field,info), std::ref(getItem(dest,info)), std::ref(index));
							};
							auto unpacker = [] (std::decay_t<T_vector> &vector, const mask_type&mask, const field_type&field, type &dest, FAR_SIGNED_INDEX_TYPE &index){
								if (mask){
									dest = vector(index);
									index++;
								} else {
									dest = field;
								}
							};
							typename N_ary_tuple_type<FAR_SIGNED_INDEX_TYPE,rank>::type t1,t2;
							mask.getIndexer().template forEachCore<unpacker, mapper>(t1, t2, rdNone, vector, mask, field, unpacked, index);
							return unpacked;
						}


					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T_str=int, typename T_set=int, typename T_back=int>
						auto verify(T_str &&string, T_set &&set, T_back &&back){
							AERRCHECKSTRINGY(T_str,"verify","string");
							AERRCHECKSTRINGY(T_set,"verify","set");
							AERRCHECKTRUTHY(T_back,"verify","back");
							using type = typename arrayInfo<T_str>::type;
							using setType = typename arrayInfo<T_set>::type;
							using backtype = typename arrayInfo<T_back>::type;
							auto l = [](const type &string, const setType &set, const backtype &back){
								kind result=0;
								using strtype = std::conditional_t<is_stdstring_v<type>,type,std::basic_string<std::remove_cv_t<std::remove_pointer_t<type>>>>;
								using settype = std::conditional_t<is_stdstring_v<setType>,setType,std::basic_string<std::remove_cv_t<std::remove_pointer_t<setType>>>>;
								strtype sstr = string;
								settype sset = set;
								if (sset.size()==0) return result;
								if (back){
									auto pos = sstr.find_last_not_of(sset);
									if (pos!=std::string::npos){
										result = pos+1;
									}
								} else {
									auto pos = sstr.find_first_not_of(sset);
									if (pos!=std::string::npos){
										result = pos+1;
									}
								}
								return result;
							};

							return makeLazy<l>(std::forward<T_str>(string), std::forward<T_set>(set), std::forward<T_back>(back));
						}

					template<typename kind=FAR_SIGNED_INDEX_TYPE, typename T_str=int, typename T_set=int>
						auto verify(T_str &&string, T_set &&set){
							AERRCHECKSTRINGY(T_str,"verify","string");
							AERRCHECKSTRINGY(T_set,"verify","set");
							return verify<kind>(std::forward<T_str>(string),std::forward<T_set>(set),false);
						}

					/**
					 * An attempt at replicating Fortran's WHERE construct
					 */
					template <typename T_condition, typename T_where, typename... T_params>
						void where(T_condition &&condition, T_where &&where, T_params&&... params)
						{
							forEach<([](T_condition &condition, T_where &where, arrayInfo<T_params>::type&... params)
									{if (condition(params...)){where(params...);} })>(condition,where,params...);
						}


					/**
					 * An attempt at replicating Fortran's WHERE/ELSEWHERE construct
					 */
					template <typename T_condition, typename T_where, typename T_elsewhere, typename... T_params>
						void where_elsewhere(T_condition &&condition, T_where &&where, T_elsewhere &&elsewhere, T_params&&... params)
						{
							forEach<([](T_condition &condition, T_where &where, T_elsewhere &elsewhere, arrayInfo<T_params>::type&... params)
									{if (condition(params...)){where(params...);} else {elsewhere(params...);} })>(condition,where,elsewhere,params...);
						}

					/**
					 * Create a 1D array of items all given the value zero
					 */
					template <typename T = double, bounds_check_state bounds = bc_default, typename... T_sizes>
						Array<T, sizeof...(T_sizes), bounds> zeros(T_sizes... sizes)
						{
							Array<T, sizeof...(T_sizes), bounds> space(sizes...);
							space = 0;
							return space;
						}

					//************************* OPERATORS ********************************
					//


					/**
					 * Addition operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator+(T_lhs &&lhs, T_rhs &&rhs)
						{
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							static_assert(arrayInfo<lhs_base>::rank==arrayInfo<rhs_base>::rank || arrayInfo<rhs_base>::rank==0 || arrayInfo<lhs_base>::rank==0,"Addition of arrays with different ranks is not allowed");
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1  +  op2; };
							return makeLazy<l>(std::forward<T_lhs>(lhs),std::forward<T_rhs>(rhs));
						}

					/**
					 * Unary plus for array
					 */
					template <typename T>
						FAR_INLINE auto operator+(T &&lhs)
						{
							using lhs_base = std::decay_t<T>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op) -> decltype(auto)
							{  return -op; };
							return makeLazy<l>(std::forward<T>(lhs));
						}

					/**
					 * Unit increment operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename = std::enable_if_t<arrayInfo<T_lhs>::value>>
						FAR_INLINE auto& operator++(T_lhs &&lhs)
						{
							using lhs_base = std::decay_t<T_lhs>;
							auto l = [](typename arrayInfo<lhs_base>::type &op)
							{  op++; };
							forEach<l>(lhs);
							return lhs;
						}

					/**
					 * Unit increment operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename = std::enable_if_t<arrayInfo<T_lhs>::value>>
						FAR_INLINE auto operator++(T_lhs &&lhs,int)
						{
							using lhs_base = std::decay_t<T_lhs>;
							lhs_base val = lhs;
							auto l = [](typename arrayInfo<lhs_base>::type &op)
							{  op++; };
							forEach<l>(lhs);
							return val;
						}

					/**
					 * N increment operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<arrayInfo<T_lhs>::value>>
						FAR_INLINE auto& operator+=(T_lhs &&lhs, T_rhs &&rhs)
						{
							lhs=lhs+rhs;
							return lhs;
						}

					/**
					 * Subtraction operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator-(T_lhs &&lhs, T_rhs &&rhs)
						{
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							using index = firstIndex<T_lhs, T_rhs>;
							static_assert(arrayInfo<lhs_base>::rank==arrayInfo<rhs_base>::rank || arrayInfo<rhs_base>::rank==0 || arrayInfo<lhs_base>::rank==0,"Subtraction of arrays with different ranks is not allowed");
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1  -  op2; };
							auto a = lazyArray<l, index, decltype(lhs), decltype(rhs)>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return a;
						}

					/**
					 * Unary minus for array
					 */
					template <typename T>
						FAR_INLINE auto operator-(T &&lhs)
						{
							using lhs_base = std::decay_t<T>;
							using index = arrayInfo<T>::index;
							auto l = [](typename arrayInfo<lhs_base>::type const &op) -> decltype(auto)
							{  return -op; };
							auto a = lazyArray<l, index, decltype(lhs)>(std::forward<T>(lhs));
							return a;
						}


					/**
					 * Unit decrement operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename = std::enable_if_t<arrayInfo<T_lhs>::value>>
						FAR_INLINE auto& operator--(T_lhs &&lhs)
						{
							using lhs_base = std::decay_t<T_lhs>;
							auto l = [](typename arrayInfo<lhs_base>::type &op)
							{  op--; };
							forEach<l>(lhs);
							return lhs;
						}

					/**
					 * Unit decrement operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename = std::enable_if_t<arrayInfo<T_lhs>::value>>
						FAR_INLINE auto operator--(T_lhs &&lhs,int)
						{
							using lhs_base = std::decay_t<T_lhs>;
							lhs_base val = lhs;
							auto l = [](typename arrayInfo<lhs_base>::type &op)
							{  op--; };
							forEach<l>(lhs);
							return val;
						}

					/**
					 * N decrement operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<arrayInfo<T_lhs>::value>>
						FAR_INLINE auto& operator-=(T_lhs &&lhs, const T_rhs &&rhs)
						{
							lhs=lhs-rhs;
							return lhs;
						}

					/**
					 * Multiplication operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator*(T_lhs &&lhs, T_rhs &&rhs)
						{
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							using index = firstIndex<T_lhs, T_rhs>;
							static_assert(arrayInfo<lhs_base>::rank==arrayInfo<rhs_base>::rank || arrayInfo<rhs_base>::rank==0 || arrayInfo<lhs_base>::rank==0,"Multiplication of arrays with different ranks is not allowed");
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1  *  op2; };
							auto a = lazyArray<l, index, decltype(lhs), decltype(rhs)>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return a;
						}

					/**
					 * N multiply operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<arrayInfo<T_lhs>::value>>
						FAR_INLINE auto& operator*=(T_lhs &&lhs, T_rhs &&rhs)
						{
							lhs=lhs*rhs;
							return lhs;
						}

					/**
					 * Division operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator/(T_lhs &&lhs, T_rhs &&rhs)
						{
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							using index = firstIndex<T_lhs, T_rhs>;
							static_assert(arrayInfo<lhs_base>::rank==arrayInfo<rhs_base>::rank || arrayInfo<rhs_base>::rank==0 || arrayInfo<lhs_base>::rank==0,"Division of arrays with different ranks is not allowed");
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1  /  op2; };
							auto a = lazyArray<l, index, decltype(lhs), decltype(rhs)>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return a;
						}

					/**
					 * N divide operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<arrayInfo<T_lhs>::value>>
						FAR_INLINE auto& operator/=(T_lhs &&lhs, T_rhs &&rhs)
						{
							lhs=lhs/rhs;
							return lhs;
						}

					/**
					 * Modulo division operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator%(T_lhs &&lhs, T_rhs &&rhs)
						{
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							static_assert(arrayInfo<lhs_base>::rank==arrayInfo<rhs_base>::rank || arrayInfo<rhs_base>::rank==0 || arrayInfo<lhs_base>::rank==0,"Modulo of arrays with different ranks is not allowed");
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1  %  op2; };
							auto a = makeLazy<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return a;
						}

					/**
					 * Modulo. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<arrayInfo<T_lhs>::value>>
						FAR_INLINE auto& operator%=(T_lhs &&lhs, T_rhs &&rhs)
						{
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							static_assert(arrayInfo<lhs_base>::rank==arrayInfo<rhs_base>::rank || arrayInfo<rhs_base>::rank==0 || arrayInfo<lhs_base>::rank==0,"Modulo of arrays with different ranks is not allowed");
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{ op1  %=  op2; };
							forEach<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return lhs;

						}


					/**
					 * <operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator<(T_lhs &&lhs, T_rhs &&rhs)
						{
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							using index = firstIndex<T_lhs, T_rhs>;
							static_assert(arrayInfo<lhs_base>::rank==arrayInfo<rhs_base>::rank || arrayInfo<rhs_base>::rank==0 || arrayInfo<lhs_base>::rank==0,"Comparison of arrays with different ranks is not allowed");
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1 < op2; };
							auto a = lazyArray<l, index, decltype(lhs), decltype(rhs)>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return a;
						}

					/**
					 * <=operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator<=(T_lhs &&lhs, T_rhs &&rhs)
						{
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							using index = firstIndex<T_lhs, T_rhs>;
							static_assert(arrayInfo<lhs_base>::rank==arrayInfo<rhs_base>::rank || arrayInfo<rhs_base>::rank==0 || arrayInfo<lhs_base>::rank==0,"Comparison of arrays with different ranks is not allowed");
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1 <= op2; };
							auto a = lazyArray<l, index, decltype(lhs), decltype(rhs)>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return a;
						}

					/**
					 * >operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator>(T_lhs &&lhs, T_rhs &&rhs)
						{
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							using index = firstIndex<T_lhs, T_rhs>;
							static_assert(arrayInfo<lhs_base>::rank==arrayInfo<rhs_base>::rank || arrayInfo<rhs_base>::rank==0 || arrayInfo<lhs_base>::rank==0,"Comparison of arrays with different ranks is not allowed");
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{return op1 > op2; };
							auto a = lazyArray<l, index, decltype(lhs), decltype(rhs)>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return a;
						}

					/**
					 * >=operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator>=(T_lhs &&lhs, T_rhs &&rhs)
						{
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							using index = firstIndex<T_lhs, T_rhs>;
							static_assert(arrayInfo<lhs_base>::rank==arrayInfo<rhs_base>::rank || arrayInfo<rhs_base>::rank==0 || arrayInfo<lhs_base>::rank==0,"Comparison of arrays with different ranks is not allowed");
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1 >= op2; };
							auto a = lazyArray<l, index, decltype(lhs), decltype(rhs)>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return a;
						}

					/**
					 * ==operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator==(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator==");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1 == op2; };
							return makeLazy<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
						}

					/**
					 * !=operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator!=(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator!=");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1 != op2; };
							return makeLazy<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
						}

					/**
					 * || operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator||(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator||");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1 || op2; };
							return makeLazy<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
						}

					/**
					 * | operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator|(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator|");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1 | op2; };
							return makeLazy<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
						}

					/**
					 * |= operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs,T_rhs>>>
						FAR_INLINE T_lhs& operator|=(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator|");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type &op1, typename arrayInfo<rhs_base>::type const &op2)
							{  op1 |= op2; };
							forEach<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return lhs;
						}

					/**
					 * && operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator&&(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator&&");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1 && op2; };
							return makeLazy<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
						}

					/**
					 * & operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator&(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator&");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1 & op2; };
							return makeLazy<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
						}

					/**
					 * &= operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs,T_rhs>>>
						FAR_INLINE void operator&=(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator|");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type &op1, typename arrayInfo<rhs_base>::type const &op2)
							{  op1 &= op2; };
							forEach<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
						}

					/**
					 * ^ operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator^(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator&");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1 ^ op2; };
							return makeLazy<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
						}

					/**
					 * ^= operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto& operator^=(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator&");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  op1 ^= op2; };
							return forEach<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return lhs;
						}

					/**
					 * ! operator. Only triggers if is an Array or lazyArray
					 */
					template <typename T, typename = std::enable_if_t<std::is_base_of_v<array_base,clean<T>>>>
						FAR_INLINE auto operator!(T &&op)
						{
							using base = arrayInfo<std::decay_t<T>>::type;
							auto l = [](base const &op) -> decltype(auto)
							{  return !op; };
							return makeLazy<l>(std::forward<T>(op));
						}

					/**
					 * ~ operator. Only triggers if is an Array or lazyArray
					 */
					template <typename T, typename = std::enable_if_t<std::is_base_of_v<array_base,clean<T>>>>
						FAR_INLINE auto operator~(T &&op)
						{
							using base = arrayInfo<std::decay_t<T>>::type;
							auto l = [](base const &op) -> base
							{  return ~op; };
							return makeLazy<l>(std::forward<T>(op));
						}

					/**
					 * << operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator<<(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator&");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1 << op2; };
							return makeLazy<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
						}

					/**
					 * <<= operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto& operator<<=(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator&");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  op1 <<= op2; };
							forEach<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return lhs;
						}


					/**
					 * << operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto operator>>(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator&");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  return op1 >> op2; };
							return makeLazy<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
						}

					/**
					 * << operator. Only triggers if at least one side is an Array or lazyArray
					 */
					template <typename T_lhs, typename T_rhs, typename = std::enable_if_t<oneSideArray<T_lhs, T_rhs>>>
						FAR_INLINE auto& operator>>=(T_lhs &&lhs, T_rhs &&rhs)
						{
							constexpr int rank_lhs = arrayInfo<T_lhs>::rank;
							constexpr int rank_rhs = arrayInfo<T_rhs>::rank;
							static_assert(rank_lhs==rank_rhs || rank_rhs==0,"Rank of lhs and rhs must be the same or rhs must be scalar in operator&");
							using lhs_base = std::decay_t<T_lhs>;
							using rhs_base = std::decay_t<T_rhs>;
							auto l = [](typename arrayInfo<lhs_base>::type const &op1, typename arrayInfo<rhs_base>::type const &op2) -> decltype(auto)
							{  op1 >>= op2; };
							forEach<l>(std::forward<T_lhs>(lhs), std::forward<T_rhs>(rhs));
							return lhs;
						}


					// Non Fortran helpers
					template <typename T>
						auto forceArray(T &&source)
						{
							auto A = std::decay_t<T>:: template buildVariant<typename arrayInfo<std::decay_t<T>>::type, arrayInfo<std::decay_t<T>>::rank,
									 arrayInfo<std::decay_t<T>>::bounds_check>();
							A=source;
							return A;
						}

					template <typename T>
						auto forcePointer(T &&source)
						{
							auto A = std::decay_t<T>:: template buildVariant<typename arrayInfo<std::decay_t<T>>::type, arrayInfo<std::decay_t<T>>::rank,
									 arrayInfo<std::decay_t<T>>::bounds_check>();
							A.pointTo(source);
							A.pointerCopy=true;
							return A;
						}

					/**
					 * Create an array from a source lazy array
					 */
					template <typename T, typename std::enable_if_t<std::is_base_of_v<array_operator, std::decay_t<T>>,int> = 0>
						auto toArray(T &&source)
						{
							return forceArray(std::forward<T>(source));
						}

					/**
					 * Pass through an array if it is not lazy and an lvalue
					 */
					template <typename T, typename std::enable_if_t<std::is_base_of_v<array_actual, std::decay_t<T>>,int> = 0>
						auto& toArray(T &source)
						{
							return std::forward<T&>(source);
						}

					template <typename T, typename std::enable_if_t<std::is_base_of_v<array_actual, std::decay_t<T>>,int> = 0>
						auto&& toArray(T &&source)
						{
							return std::forward<T&&>(source);
						}

						/**
						 * Convert a C++ iterable to a far::Array
						 */
					template <isIterable T>
						auto toArray(const T &source)
						{
							using coreType = std::decay_t<T>::value_type;
							Array<coreType,1> result(source.size());
							FAR_SIGNED_INDEX_TYPE index = result.getLB(1);
							for (auto &item : source)
							{
								//If the item has an assignment operator, use it
								if constexpr (std::is_assignable_v<coreType,decltype(item)>)
								{
									result(index) = item;
								}
								else
								{
									//Otherwise, use the constructor
									new(&result(index)) coreType(item);
								}
								index++;
							}
							return result;
						}

					/**
					 * Create a contiguous array from a source
					 */
					template <typename T=int, typename std::enable_if_t<std::is_base_of_v<array_base, std::decay_t<T>>,int> = 0>
						auto toContiguousArray(T &&source)
						{
							auto A = std::decay_t<T>:: template buildVariant<typename arrayInfo<std::decay_t<T>>::type, arrayInfo<std::decay_t<T>>::rank,
									 arrayInfo<std::decay_t<T>>::bounds_check,true>();
							A=source;
							return A;
						}

						/**
						 * Convert a C++ iterable to a contiguous far::Array
						 */
					template <isIterable T>
						auto toContiguousArray(const T &source)
						{
							using coreType = std::decay_t<T>::value_type;
							contiguousArray<coreType,1> result(source.size());
							FAR_SIGNED_INDEX_TYPE index = result.getLB(1);
							for (auto &item : source)
							{
								//If the item has an assignment operator, use it
								if constexpr (std::is_assignable_v<coreType,decltype(item)>)
								{
									result(index) = item;
								}
								else
								{
									//Delete the existing item
									result(index).~coreType();
									//Otherwise, use the constructor
									new(&result(index)) coreType(item);
								}
								index++;
							}
							return result;
						}

					/**
					 * Create an array with a specific indexer from a source
					 */
					template<template<int rank, bounds_check_state bc> class index_type, bool contiguous = false, typename T,  typename std::enable_if_t<std::is_base_of_v<array_base, std::decay_t<T>>,int> = 0>
						decltype(auto) toIndexedArray(T&& src)
						{
							Array<typename arrayInfo<T>::type, arrayInfo<T>::rank, arrayInfo<T>::bounds_check,index_type,contiguous> result;
							result = src;
							return result;
						}

					/** 
					 * Create a Fortran indexed array from a source
					 */
					template <typename T, bool contiguous =false, typename std::enable_if_t<std::is_base_of_v<array_base, std::decay_t<T>>,int> = 0>
						auto toFortranArray(T &&source)
						{
							return toIndexedArray<FortranIndexBuilder, contiguous>(std::forward<T>(source));
						}

					/** 
					 * Create a C indexed array from a source
					 */
					template <typename T, bool contiguous =false, typename std::enable_if_t<std::is_base_of_v<array_base, std::decay_t<T>>,int> = 0>
						auto toCArray(T &&source)
						{
							return toIndexedArray<CIndexBuilder, contiguous>(std::forward<T>(source));
						}

						/**
						 * Convert a C++ iterable to a far::Array with a specific indexer
						 */
					template <template<int rank, bounds_check_state bc> class index_type, isIterable T, bool contiguous=false>
						auto toIndexedArray(const T &source)
						{
							using coreType = std::decay_t<T>::value_type;
							Array<coreType,1,bc_default,index_type,contiguous> result(source.size());
							FAR_SIGNED_INDEX_TYPE index = result.getLB(1);
							for (auto &item : source)
							{
								result(index) = item;
								index++;
							}
							return result;
						}

						/**
						 * Convert a C++ iterable to a far::Array with a Fortran Indexer
						 */
					template <isIterable T, bool contiguous=false>
						auto toFortranArray(const T &source)
						{
							return toIndexedArray<FortranIndexBuilder>(source);
						}

						/**
						 * Convert a C++ iterable to a far::Array with a C Indexer
						 */
					template <isIterable T, bool contiguous=false>
						auto toCArray(const T &source)
						{
							return toIndexedArray<CIndexBuilder>(source);
						}



					template<typename T, bool copy_in = true, bool copy_out = true>
						class contiguous{
							using core_type = arrayInfo<T>::type;
							constexpr static int rank =  arrayInfo<T>::rank;
							constexpr static bounds_check_state bc = arrayInfo<T>::bounds_check;
							using stripped_type = std::invoke_result_t<decltype(std::decay_t<T>:: template buildVariant<core_type,rank,bc,false>)>;
							using contiguous_type = std::invoke_result_t<decltype(std::decay_t<stripped_type>:: template buildVariant<core_type, rank, bc, true>)>;
							stripped_type *bound=nullptr;
							stripped_type *forward=nullptr;
							contiguous_type *contig=nullptr;
							template<bool errorOnFail>
							void copyBackInner(){
								if constexpr(copy_out){
									if (bound && (forward!=bound)) (*bound)=(*forward);
								} else if constexpr(errorOnFail){
									static_assert(!errorOnFail,"Requesting manual copy_out but copy_out not specified");
								}
							}
							public:
							contiguous(T& src){
								bind(src);
							}
							contiguous(const contiguous& src){
								bind(*src.bound);
							}
							contiguous(contiguous&& src){
								bound = src.bound;
								forward = src.forward;
								contig = src.contig;
								src.bound = nullptr;
								src.forward = nullptr;
								src.contig = nullptr;
							}
							void copyBack(){
								copyBackInner<true>();
							}
							void unbind() {
								copyBackInner<false>();
								if (forward&&(!bound || forward!=bound)) delete forward;
								if (contig) delete contig;
								bound=nullptr;
								forward=nullptr;
								contig=nullptr;
							}
							void bind(T& src){
								unbind();
								if constexpr(!std::is_base_of_v<array_operator,std::decay_t<T>>){
									bound = &src;
								}
								if (src.getContiguous()){
									//If the array is contiguous and not associated then
									//it must be an lvalue, so can just point to it
									if (!src.getAssociated()){
										forward = bound;
									} else {
										//Otherwise we need to copy it
										if constexpr(std::is_base_of_v<array_actual, std::decay_t<T>>){
											forward = new stripped_type();
											forward->pointTo(src);
										} else {
											forward = new stripped_type();
											forward->mold(src);
											if constexpr(copy_in){
												forEachCopyTight(*forward,src);
											}
										}
									}
								} else {
									forward = new stripped_type();
									forward->mold(src);
									if constexpr(copy_in){
										forEachCopyTight(*forward,src);
									}
								}
								contig = new contiguous_type();
								contig->pointTo(*forward);
							}
							~contiguous(){
								unbind();
							}
							core_type* data() const {
								return forward->data();
							}
							stripped_type & toArray() const {
								return *forward;
							}
							contiguous_type & toContiguous() const {
								return *contig;
							}
							core_type* toPointer() const {
								return forward->data();
							}

							operator stripped_type&() const {return toArray();}
							template<typename T_inner = stripped_type, std::enable_if_t<!arrayInfo<T_inner>::forceContiguous, int> = 0>
								operator contiguous_type&() const {return toContiguous();}
							operator core_type*() const {return toPointer();}
							auto& operator[](FAR_SIGNED_INDEX_TYPE index) const {return toPointer()[index];}
						};

					template<typename T>
						using contiguousCopyIn = contiguous<T,true,false>;

					template<typename T>
						using contiguousCopyOut = contiguous<T,false,true>;

					template<typename T>
						using contiguousCopyInOut = contiguous<T,true,true>;

					template<typename T, bool copyIn = true, bool copyOut = true>
						auto makeContiguous(T&& src){
							if constexpr (!std::is_base_of_v<array_operator,std::decay_t<T>>){
								return contiguous<T, copyIn, copyOut>(src);
							} else {
								static_assert(!copyOut,"Cannot copy out to array expression");
								return contiguous<T,copyIn,false>(src);
							}
						}

					template<typename T>
						auto makeCopyIn(T&& src){
							return makeContiguous<T,true,false>(std::forward<T>(src));
						}

					template<typename T>
						auto makeCopyOut(T&& src){
							return makeContiguous<T, false, true>(std::forward<T>(src));
						}

					template<typename T>
						auto makeCopyInOut(T&& src){
							return makeContiguous<T, true, true>(std::forward<T>(src));
						}

					/**
					 * Create a 1D array of items spanning a range
					 */
					template <typename T = double, bounds_check_state bounds = bc_default>
						Array<T, 1, bounds> linspace(T min, T max, FAR_SIGNED_INDEX_TYPE N)
						{
							Array<T, 1, bounds> space(N);
							double delta = (max - min) / double(N - 1);
							for (FAR_SIGNED_INDEX_TYPE i = defaultLB; i < defaultLB+N; ++i)
							{
								space(i) = delta * double(i-defaultLB) + min;
							}
							return space;
						}

					/**
					 * Create a log spaced array of items spanning a range
					 */
					template <typename T = double, bounds_check_state bounds = bc_default>
						Array<T, 1, bounds> logspace(T min, T max, FAR_SIGNED_INDEX_TYPE N)
						{
							Array<T, 1, bounds> space(N);
							double delta = (std::log10(max) - std::log10(min)) / double(N - 1);
							for (FAR_SIGNED_INDEX_TYPE i = defaultLB; i < defaultLB+N; ++i)
							{
								space(i) = std::pow(10.0, delta * double(i-defaultLB) + std::log10(min));
							}
							return space;
						}

					inline auto& fortFile(int n) {
						static std::map<int, std::fstream> streamMap;
						std::fstream *stream;
						FAR_OMP_CRITICAL
						{
							stream = &streamMap[n];
							if (!stream->is_open()) {
								std::string filename = "fort." + std::to_string(n);
								stream->open(filename, std::ios::in | std::ios::out | std::ios::trunc);
							}
						}
						return *stream;
					}

					//Produce an array the same size as a given array with one dimension removed
					template <typename T_result, typename T>
						auto reduce_rank(T &&source, int dir)
						{
							if constexpr(arrayInfo<T>::rank==1){
								using src_type = std::conditional_t<std::is_void_v<T_result>,typename arrayInfo<T>::type, T_result>;
								return src_type();
							} else {
								using src_type = std::conditional_t<std::is_void_v<T_result>,typename arrayInfo<T>::type, T_result>;
								constexpr bounds_check_state bounds = arrayInfo<T>::bounds_check;
								constexpr int rank = arrayInfo<T>::rank;
								using dest_type = std::invoke_result_t<decltype(std::decay_t<T>:: template buildVariant<src_type, rank-1, bounds>)>;
								dest_type r;
								FAR_UNSIGNED_INDEX_TYPE SZ[rank-1];
								int idest = 0;
								// Allocate the array to the right size
								for (int src = 0; src < rank; ++src)
								{
									if (src != dir-1)
									{
										SZ[idest] = source.getRankSize(src+1);
										idest++; 
									}
								}
								r.allocate(SZ);
								return r;
							}
						}

					// Define the << operator for arrays
					template<typename T_os, typename T, typename = std::enable_if_t<arrayInfo<T>::value || arrayInfo<T_os>::value>>
						inline T_os& operator<<(T_os &os, T&& obj) {
							using core = typename arrayInfo<T>::type;
							forEach<[](core& item, std::ostream& os){
								if constexpr(std::is_same_v<core,bool>){
#ifndef C_BOOL_PRINT
									os << (item?"T":"F") << " ";
#else
									os << item << " ";
#endif
								} else {
									os << item << " ";
								}
							}>(obj, os);
							return os;
						}

					// Define the >> operator for arrays
					template<typename T_os, typename T, typename = std::enable_if_t<arrayInfo<T>::value || arrayInfo<T_os>::value>>
						inline T_os& operator>>(std::istream&& is, T&& obj) {
							using core = typename arrayInfo<T>::type;
							forEach<[](core& item, std::istream& is){
								if constexpr(!std::is_same_v<core,bool>){
									(is >> item);
								} else {
									//Fortran outputs T and F rather than 1 and 0 for true and false
									//This makes FAR++ able to read that output
									char c;
									(is >> c);
									if (c=='1' or c=='t' or c=='T') {
										item = true;
									} else {
										item = false;
									}
								}
							}>(obj, is);
							return is;
						}

					/**
					 * Object to print a rank 2 array as a grid of rows and columns
					 * only works for small arrays!
					 */
					class gridPrint{
						const void *data;
						void(*actualPrint)(const void* data, std::ostream &out);
						template<typename T>
							static void corePrint(const void* data, std::ostream &os)
							{
								const T* objp = static_cast<const T*>(data);
								const T& obj = *objp;
								std::string s;
								FAR_UNSIGNED_INDEX_TYPE mx=0;
								forEach<[](const typename arrayInfo<T>::type &item, FAR_UNSIGNED_INDEX_TYPE &mx){
									std::string s = std::to_string(item);
									mx=std::max(mx,s.size()+1);
								}>(obj,mx);
								++mx; //Increment mx so there is a space
								FAR_SIGNED_INDEX_TYPE LB1 = obj.getLB(1), UB1 = obj.getUB(1);
								FAR_SIGNED_INDEX_TYPE LB2 = obj.getLB(2), UB2 = obj.getUB(2);
								for (FAR_SIGNED_INDEX_TYPE i=LB2-1;i<=UB2;++i){
									s=std::to_string(i);
									if (i==LB2-1) s="";
									if (s.size() < mx) s.insert(s.begin(), (mx-1) - s.size(), ' ');
									s=s + " |";
									os << s;
								}
								os << "\n";
								for (FAR_UNSIGNED_INDEX_TYPE i=0;i<mx;++i) os << " ";
								for (FAR_UNSIGNED_INDEX_TYPE i=0;i<mx*(UB2-LB2+2);++i) os << "-";
								os << "\n";
								for (FAR_SIGNED_INDEX_TYPE j=LB1;j<=UB1;++j){
									for (FAR_SIGNED_INDEX_TYPE i=LB2;i<=UB2;++i){
										if (i==LB2){
											s=std::to_string(j);
											if (s.size() < mx) s.insert(s.begin(), (mx-1) - s.size(), ' ');
											s=s+" |";
											os << s;
										}
										s=std::to_string(obj(j,i));
										if (s.size() < mx) s.insert(s.begin(), mx - s.size(), ' ');
										s+="|";
										os << s;
									}
									if (j!=UB1) os << "\n";
								}
							}

						template<typename T>
							static void corePrintNN(const void* data, std::ostream &os)
							{
								const T* objp = static_cast<const T*>(data);
								const T& obj = *objp;
								FAR_UNSIGNED_INDEX_TYPE mx=0;
								forEach<[](typename arrayInfo<T>::type &item, FAR_UNSIGNED_INDEX_TYPE &mx){
									std::stringstream ss;
									ss << item;
									mx=std::max(mx,ss.str().size()+1);
								}>(obj,mx);
								FAR_SIGNED_INDEX_TYPE LB1 = obj.getLB(1), UB1 = obj.getUB(1);
								FAR_SIGNED_INDEX_TYPE LB2 = obj.getLB(2), UB2 = obj.getUB(2);
								for (FAR_SIGNED_INDEX_TYPE i=LB2-1;i<=UB2;++i){
									std::string s=std::to_string(i);
									if (i==LB2-1) s="";
									if (s.size() < mx) s.insert(s.begin(), (mx-1) - s.size(), ' ');
									s=s + " |";
									os << s;
								}
								os << "\n";
								for (FAR_UNSIGNED_INDEX_TYPE i=0;i<mx;++i) os << " ";
								for (FAR_UNSIGNED_INDEX_TYPE i=0;i<mx*(UB2-LB2+2);++i) os << "-";
								os << "\n";

								for (FAR_SIGNED_INDEX_TYPE j=LB1;j<=UB1;++j){
									for (FAR_SIGNED_INDEX_TYPE i=LB2;i<=UB2;++i){
										std::stringstream ss;
										std::string s;
										if (i==LB2){
											ss << j;
											s = ss.str();
											if (s.size() < mx) s.insert(s.begin(), (mx-1) - s.size(), ' ');
											s=s+" |";
											os << s;
										}
										ss.str("");
										ss << obj(j,i);
										s=ss.str();
										if (s.size() < mx) s.insert(s.begin(), mx - s.size(), ' ');
										s+="|";
										os << s;

									}
									if (j!=UB2) os << "\n";
								}
							}

						public:
						template<typename T, typename = std::enable_if_t<arrayInfo<T>::value 
							&& arrayInfo<T>::rank==2>>
							gridPrint(const T& toPrint)
							{this->data = &toPrint;
								if constexpr(std::is_arithmetic_v<typename arrayInfo<T>::type>){
									this->actualPrint = corePrint<std::decay_t<T>>;
								} else {
									this->actualPrint = corePrintNN<std::decay_t<T>>;}
							}
						void callPrint(std::ostream &os) const {this->actualPrint(this->data,os);}
					};

					/**
					 * insertion operator for gridPrint object
					 */
					inline std::ostream& operator << (std::ostream & os, const gridPrint &pp){
						pp.callPrint(os);
						return os;
					}

					template<typename T_array, typename T_item, typename... T_others>
						void fillArray(T_array& array, FAR_SIGNED_INDEX_TYPE toFill, T_item&& item, T_others&&... others){
							if constexpr(std::is_base_of_v<array_base,std::decay_t<T_item>>){
								for (auto &el:item){
									array(toFill++) = el;
								}
							} else {
								array(toFill++) = item;
							}
							if constexpr(sizeof...(others)>0) fillArray(array, toFill, others...);
						}

					template<typename T_item, typename ...T_others>
						FAR_SIGNED_INDEX_TYPE sizeList(T_item&& item, T_others&&... others){
							FAR_SIGNED_INDEX_TYPE sz=0;
							if constexpr(std::is_base_of_v<array_base,std::decay_t<T_item>>){
								sz+= item.getSize();
							} else {
								++sz;
							}
							if constexpr(sizeof...(others)>0) {
								sz+=sizeList(others...);
							}
							return sz;
						}

					//Construct a literal array from a comma separated list of values
					template<typename T_item, typename... T_others>
						auto literal(T_item &&item, T_others&&... others){
							Array<typename arrayInfo<T_item>::type,1> array(sizeList(item,others...));
							fillArray(array, defaultLB, item, others...);
							return array;
						}

					template<typename T_first, typename... T_values>
						auto constLiteral(T_first &&first, T_values&&... values){
							std::array<std::decay_t<T_first>, sizeof...(values)+1> a = {first,values...};
							return a;
						}

					template<typename T, typename T_src>
						auto lazy_cast(T_src&& src){
							auto l = [](arrayInfo<T_src>::type const &op) -> decltype(auto)
							{  return static_cast<T>(op); };
							return makeLazy<l>(std::forward<T_src>(src));
						}

					template<typename T, typename T_src>
						auto array_cast(T_src&& src){
							auto r = src.template buildVariant<T>();
							r=src;
							return r;
						}
			};
#endif
