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
#ifndef FAR_DEBUG_H
#define FAR_DEBUG_H

#include <iostream>

  template<typename T_index, int level=0>
    inline void printTuple(T_index &&index){
      constexpr int rank = std::tuple_size_v<std::decay_t<T_index>>;
      if constexpr(level<rank-1){
        std::cout << std::get<level>(index) << " ";
        printTuple<T_index,level+1>(std::forward<T_index>(index));
      } else {
        std::cout << std::get<level>(index);
      }
    }

  template<typename T, typename... T_others>
    inline void printPackTypes(){
      if constexpr(std::is_const_v<T>){
        std::cout << "const ";
      }
      std::cout << demangle(typeid(T).name());
      if constexpr(std::is_rvalue_reference_v<T>){
        std::cout << "&&";
      } else {
        if constexpr(std::is_reference_v<T>){
          std::cout << "&";
        }
      }
      if constexpr(sizeof...(T_others)>0){
        std::cout << " || ";
        printPackTypes<T_others...>();
      } else {
        std::cout << "\n";
      }
    }

  template<typename T, typename... T_others>
    inline void printPack(T&& c, T_others... o){
			//std::cout << c;
      if constexpr(sizeof...(T_others)>0){
        std::cout << " || ";
        printPack<T_others...>(std::forward<T_others>(o)...);
      } else {
				std::cout << c;
        std::cout << "\n";
      }
    }

  template<int level, typename T1, typename... Ts>
    void printRV(T1&& current, Ts&&... params){
      std::cout << level << " is rvalue " << std::is_rvalue_reference_v<decltype(current)> << "\n";
      if constexpr(sizeof...(params)>0) printRV<level+1>(std::forward<Ts>(params)...);
    }

#endif
