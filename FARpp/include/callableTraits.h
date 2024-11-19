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
/**
 * @file callableTraits.h
 * @brief Provides utilities to extract type information from callable objects.
 *
 * This header defines the `callableTraits` structure template, which can be used to
 * determine the return type and parameter types of various callable objects, including
 * function pointers, member function pointers, functors, and lambdas.
 */

namespace far {

/**
 * @brief Base template for callableTraits.
 *
 * This template is used as a fallback for types that are not recognized as callables.
 * It defines a static constexpr boolean `value` set to false and a `params` type alias
 * to a tuple containing a single `int`.
 *
 * @tparam T The type to be inspected.
 * @tparam void A SFINAE parameter to enable partial specialization.
 */
template <typename T, typename = void>
struct callableTraits {
    static constexpr bool value = false; ///< Indicates whether the type is a callable.
    using params = std::tuple<int>; ///< Placeholder tuple type for non-callables.
};

/**
 * @brief Partial specialization for callable types.
 *
 * This specialization is used for types that have an `operator()` member function.
 * It inherits from `callableTraits` instantiated with the type of the `operator()`.
 *
 * @tparam T The type to be inspected.
 */
template <typename T>
struct callableTraits<T, std::void_t<decltype(&T::operator())>> : callableTraits<decltype(&T::operator())> {};

/**
 * @brief Specialization for const member function pointers.
 *
 * This specialization extracts the return type and parameter types from const member
 * function pointers.
 *
 * @tparam C The class type.
 * @tparam R The return type of the member function.
 * @tparam Args The parameter types of the member function.
 */
template <typename C, typename R, typename... Args>
struct callableTraits<R(C::*)(Args...) const> {
    using type = R; ///< The return type of the callable.
    using params = std::tuple<Args...>; ///< A tuple of the parameter types.
    static constexpr bool value = true; ///< Indicates that the type is a callable.
};

/**
 * @brief Specialization for non-const member function pointers.
 *
 * This specialization extracts the return type and parameter types from non-const member
 * function pointers.
 *
 * @tparam C The class type.
 * @tparam R The return type of the member function.
 * @tparam Args The parameter types of the member function.
 */
template <typename C, typename R, typename... Args>
struct callableTraits<R(C::*)(Args...)> {
    using type = R; ///< The return type of the callable.
    using params = std::tuple<Args...>; ///< A tuple of the parameter types.
    static constexpr bool value = true; ///< Indicates that the type is a callable.
};

/**
 * @brief Specialization for function pointers.
 *
 * This specialization extracts the return type and parameter types from function pointers.
 *
 * @tparam R The return type of the function.
 * @tparam Args The parameter types of the function.
 */
template <typename R, typename... Args>
struct callableTraits<R(*)(Args...)> {
    using type = R; ///< The return type of the callable.
    using params = std::tuple<Args...>; ///< A tuple of the parameter types.
    static constexpr bool value = true; ///< Indicates that the type is a callable.
};

/**
 * @brief Specialization for std::function types.
 *
 * This specialization extracts the return type and parameter types from std::function types.
 *
 * @tparam R The return type of the function.
 * @tparam Args The parameter types of the function.
 */
template <typename R, typename... Args>
struct callableTraits<R(Args...)> {
    using type = R; ///< The return type of the callable.
    using params = std::tuple<Args...>; ///< A tuple of the parameter types.
    static constexpr bool value = true; ///< Indicates that the type is a callable.
};

} // namespace far
