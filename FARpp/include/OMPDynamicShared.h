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
#ifndef OMPDYNAMICSHARED_H
#define OMPDYNAMICSHARED_H

#include <iostream>
#include <map>
#include <array>
#ifdef _OPENMP
#include <omp.h>
#endif
/**
 * @namespace ompdynamicshared
 * @brief Namespace containing classes for dynamic shared memory management in OpenMP.
 */

/**
 * @class dynamicShared
 * @brief Manages dynamic shared memory for OpenMP parallel regions.
 * 
 * @tparam T The type of the shared variable.
 * @tparam perThread If true, each thread gets its own instance of T; otherwise, all threads share a single instance.
 * @tparam slots The number of slots available for shared variables.
 * 
 * This class provides a mechanism to dynamically allocate and manage shared memory in OpenMP parallel regions.
 * It supports both per-thread and shared instances of the variable type T.
 * 
 * @note This class requires OpenMP support. If OpenMP is not available, it falls back to a single-threaded implementation.
 * 
 * @fn dynamicShared::dynamicShared(Ts&&... params)
 * @brief Constructor that initializes the shared variable.
 * 
 * @param params Parameters to initialize the shared variable.
 * 
 * @fn dynamicShared::dynamicShared(dynamicShared &&other) noexcept
 * @brief Move constructor.
 * 
 * @param other The dynamicShared object to move from.
 * 
 * @fn dynamicShared::dynamicShared(const dynamicShared &other)
 * @brief Copy constructor.
 * 
 * @param other The dynamicShared object to copy from.
 * 
 * @fn T* dynamicShared::getAllData()
 * @brief Returns a pointer to the shared data.
 * 
 * @return Pointer to the shared data.
 * 
 * @fn T& dynamicShared::getData()
 * @brief Returns a reference to the shared data for the current thread.
 * 
 * @return Reference to the shared data.
 * 
 * @fn dynamicShared::operator T&()
 * @brief Conversion operator to T&.
 * 
 * @return Reference to the shared data.
 * 
 * @fn void dynamicShared::operator= (const T& other)
 * @brief Assignment operator to assign a value to the shared data.
 * 
 * @param other The value to assign.
 * 
 * @fn dynamicShared::~dynamicShared()
 * @brief Destructor that releases the shared variable.
 */

/**
 * @typedef perThreadVariable
 * @brief Alias for dynamicShared with perThread set to true.
 * 
 * @tparam T The type of the shared variable.
 * @tparam slots The number of slots available for shared variables.
 */

/**
 * @typedef sharedVariable
 * @brief Alias for dynamicShared with perThread set to false.
 * 
 * @tparam T The type of the shared variable.
 * @tparam slots The number of slots available for shared variables.
 */
namespace ompdynamicshared{
	template<typename T, bool perThread, int slots>
		class dynamicShared{
#ifdef _OPENMP
			struct meta{
				int count = 0;
				int uid =  0;
				T* data=nullptr;
			};
			static inline std::map<int,std::array<bool,slots>> uids;
			static inline std::map<int,int> activeUID;
			static inline std::map<int,std::array<meta,slots>> data;
			meta *myMeta=nullptr;

			static int getID(){
				for (int i=0;i<slots;++i){
					if (!uids[omp_get_team_num()][i]) {
						uids[omp_get_team_num()][i] = true;
						return i+1;
					}
				}
				throw std::out_of_range("");
			}

			static void releaseID(int id){
				uids[omp_get_team_num()][id-1] = false;
			}

			template <typename... T_params>
				static meta* getVar(T_params&&... params){
					int team = omp_get_team_num();
					std::cout << team << std::endl;
					meta* res;
#pragma omp critical
					{
						if (activeUID[team]==0) activeUID[team]=getID();
						auto& m = data[team][activeUID[team]];
						m.uid = activeUID[team];
						if (!m.data) {
							if constexpr(perThread) {
								m.data = (T*)std::malloc(sizeof(T)*omp_get_num_threads());
							} else {
								m.data = new T(params...);
							}
						}
						if constexpr(perThread) new(&m.data[omp_get_thread_num()])T(params...);
						if(++m.count==omp_get_num_threads()) {
							activeUID[team]=0;
						}
						res = &m;
					}
#pragma omp barrier
					return res;
				}

			static void releaseVar(int uid){
				int team = omp_get_team_num();
				auto& m = data[team][uid];
				if constexpr(perThread) m.data[omp_get_thread_num()].~T();
				if(--m.count == 0){
					if (perThread){
						std::free(m.data);
					} else {
						delete m.data;
					}
					m.data = nullptr;
					releaseID(uid);
				}
#pragma omp barrier
			}
#else
			T*data;
#endif
			public:
			template<typename... Ts>
				dynamicShared(Ts&&... params){
#ifdef _OPENMP
					this->myMeta = getVar(params...);
#else
					this->data = new T(params...);
#endif
				}

				dynamicShared(dynamicShared &&other) noexcept{
#ifdef _OPENMP
          this->myMeta = other.myMeta;
#else
          this->data = other.data;
#endif
				}

        dynamicShared(const dynamicShared &other){
#ifdef _OPENMP
          this->myMeta = other.myMeta;
					this->myMeta->count++;
#else
          this->data = other.data;
#endif
        }

			T* getAllData(){
#ifdef _OPENMP
				return myMeta->data;
#else
				return data;
#endif
			}

			T& getData() {
#ifdef _OPENMP
				if constexpr(perThread){
					return myMeta->data[omp_get_thread_num()];
				} else {
					return *myMeta->data;
				}
#else
				return *data;
#endif
			}

			operator T&(){
				return getData();
			}

			void operator= (const T& other){
				getData() = other;
			}

			~dynamicShared(){
#ifdef _OPENMP
				releaseVar(myMeta->uid);
#else
				delete data;
#endif
			}
		};

	template<typename T, int slots=2048>
		using perThreadVariable=dynamicShared<T,true,slots>;

	template<typename T, int slots=2048>
		using sharedVariable=dynamicShared<T,false,slots>;
}

#endif
