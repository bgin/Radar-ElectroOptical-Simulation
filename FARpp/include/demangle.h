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
#ifndef DEMANGLE_H
#define DEMANGLE_H
#include <string>
#include <stdlib.h>
#if __has_include(<cxxabi.h>)
	#include <cxxabi.h>
	#define HAS_DEMANGLE 1
#else
	#define HAS_DEMANGLE 0
#endif


inline std::string demangle(std::string mangled){
	int status;
#ifdef HAS_DEMANGLE
	char *dem = abi::__cxa_demangle(mangled.c_str(),nullptr,nullptr,&status);
	std::string s = dem;
	free(dem);
	return s;
#else
	std::string s="";
#ifndef DML_NO_WARNING
	static bool warning=false;
	if (!warning){
		s="***WARNING*** Unable to demangle names\n"
		warning=true;
	}
#endif
	s+= mangled;
#endif
}

template<typename T>
inline std::string demangle(){
	std::string s=demangle(typeid(T).name());
	if (std::is_reference_v<T>) {
		if (std::is_rvalue_reference_v<T>) s+="&&";
		else s+="&";
	}
	if (std::is_const_v<T>) s="const " + s;
	return s;
}

#endif
