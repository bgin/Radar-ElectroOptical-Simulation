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
#ifndef LBREFWRAPPER
#define LBREFWRAPPER

namespace far{

    /**
     * @brief A wrapper class to allow for late binding of references
     * @details This class is used to allow for late binding of references. 
     * This is useful when you want to pass a reference or have an array of references.
     * The first time you assign a value to the reference, the reference is bound to the value
     * and all subsequent operations are done on the bound reference.
     * Includes casting operators to allow for easy conversion to the reference type.
     * 
     * @tparam T The type of the reference to be wrapped.
     */
    template<typename T>
    class LBrefWrapper{
        private:
        T* ptr=nullptr; ///< Pointer to the reference being wrapped.

        public:
        using type = T; ///< Type alias for the wrapped reference type.

        /**
         * @brief Default constructor. Initializes the pointer to nullptr.
         */
        LBrefWrapper():ptr(nullptr){}

        /**
         * @brief Constructor that initializes the wrapper with a pointer.
         * @param ptr Pointer to the reference to be wrapped.
         */
        LBrefWrapper(T* ptr):ptr(ptr){}

        /**
         * @brief Constructor that initializes the wrapper with a reference.
         * @param ref Reference to be wrapped.
         */
        LBrefWrapper(T& ref):ptr(&ref){}

        /**
         * @brief Constructor that initializes the wrapper with a const reference.
         * @param ref Const reference to be wrapped.
         */
        LBrefWrapper(const T& ref):ptr(const_cast<T*>(&ref)){}

        /**
         * @brief Copy constructor.
         * @param ref Another LBrefWrapper object to copy from.
         */
        LBrefWrapper(const LBrefWrapper& ref):ptr(ref.ptr){}

        /**
         * @brief Move constructor.
         * @param ref Another LBrefWrapper object to move from.
         */
        LBrefWrapper(LBrefWrapper&& ref):ptr(ref.ptr){
            ref.ptr=nullptr;
        }

        /**
         * @brief Binds the wrapper to a reference.
         * @param ref Reference to bind to.
         */
        void bind(T& ref){
            ptr = &ref;
        }

        /**
         * @brief Binds the wrapper to a const reference.
         * @param ref Const reference to bind to.
         */
        void bind(const T& ref){
            ptr = const_cast<T*>(&ref);
        }

        /**
         * @brief Binds the wrapper to a pointer.
         * @param ref Pointer to bind to.
         */
        void bind(T* ref){
            ptr = ref;
        }

        /**
         * @brief Assignment operator for non-const lvalue reference.
         * @param value Reference to assign.
         * @return The assigned reference.
         */
        T& operator=(T& value){
            if(ptr!=nullptr){
                *ptr = value;
            } else {
                ptr=&value;
            }
            return value;
        }

        /**
         * @brief Assignment operator for const lvalue reference.
         * @param value Const reference to assign.
         * @return The assigned reference.
         */
        T& operator=(const T& value){
            if(ptr!=nullptr){
                *ptr = value;
            } else {
                ptr=const_cast<T*>(&value);
            }
            return *ptr;
        }

        /**
         * @brief Assignment operator for rvalue reference.
         * @param value Rvalue reference to assign.
         * @return The assigned reference.
         * @throws std::runtime_error if the wrapper is not bound to a reference.
         */
        T& operator=(T&& value){
            if(ptr!=nullptr){
                *ptr = value;
            } else {
                throw std::runtime_error("Cannot have late bound reference to rvalue");
            }
            return *ptr;
        }

        /**
         * @brief Assignment operator for const rvalue reference.
         * @param value Const rvalue reference to assign.
         * @return The assigned reference.
         * @throws std::runtime_error if the wrapper is not bound to a reference.
         */
        T& operator=(const T&& value){
            if(ptr!=nullptr){
                *ptr = value;
            } else {
                throw std::runtime_error("Cannot have late bound reference to rvalue");
            }
            return *ptr;
        }

        /**
         * @brief Assignment operator for another LBrefWrapper object.
         * @param ref Another LBrefWrapper object to assign from.
         * @return The assigned LBrefWrapper object.
         */
        LBrefWrapper& operator=(const LBrefWrapper& ref){
            ptr = ref.ptr;
            return *this;
        }

        /**
         * @brief Conversion operator to the reference type.
         * @return The bound reference.
         * @throws std::runtime_error if the wrapper is not bound to a reference.
         */
        operator T&(){
            if (!ptr) throw std::runtime_error("Single use reference not set");
            return *ptr;
        }

        /**
         * @brief Conversion operator to the const reference type.
         * @return The bound const reference.
         * @throws std::runtime_error if the wrapper is not bound to a reference.
         */
        operator const T&() const{
            if (!ptr) throw std::runtime_error("Single use reference not set");
            return *ptr;
        }

        /**
         * @brief Gets the bound reference.
         * @return The bound reference.
         * @throws std::runtime_error if the wrapper is not bound to a reference.
         */
        T& get(){
            if (!ptr) throw std::runtime_error("Single use reference not set");
            return *ptr;
        }
    };

    /**
     * @brief Trait class to check if a type is an LBrefWrapper.
     * @tparam T The type to check.
     */
    template<typename T>
    class isLBrefWrapper{
        public:
        static constexpr bool value = false; ///< Indicates if the type is an LBrefWrapper.
        using type = T; ///< Type alias for the checked type.
    };

    /**
     * @brief Specialization of isLBrefWrapper for LBrefWrapper types.
     * @tparam T The wrapped reference type.
     */
    template<typename T>
    class isLBrefWrapper<LBrefWrapper<T>>{
        public:
        static constexpr bool value = true; ///< Indicates if the type is an LBrefWrapper.
        using type = T&; ///< Type alias for the wrapped reference type.
    };
}

#endif
