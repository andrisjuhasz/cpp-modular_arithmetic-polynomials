# Modular Arithmetic & Polynomial Library in C++

A C++ project built to deepen my understanding of modern C++, algebraic abstractions, and the mathematical foundations used in cryptography and proof systems.

This repository implements:

- custom `std::vector` utilities
- finite field arithmetic (`ℤ/pℤ`)
- univariate polynomial arithmetic
- polynomial representation conversion:
  - monomial basis
  - Lagrange basis
- caching/precomputation optimizations
- demo examples

The long-term goal of this project is to better understand the algebra behind:

- zero-knowledge proofs
- polynomial commitment schemes
- finite field cryptography

---

## Features

### STL Vector Utilities

File: `vector_utilities.hpp`

Adds utility functionality on top of `std::vector`, including:

- coefficient-vector multiplication (convolution-style multiplication)
- scalar multiplication
- formatted printing

Example:

```cpp
std::vector<int> a = {1,2};
std::vector<int> b = {3,4};

auto result = a * b;
// result = {3,10,8}
```

### Finite Field Arithmetic (`IntMod<P>`)

Template-based modular arithmetic over finite fields.

Supports:

- addition/subtraction
- multiplication/division
- modular inversion
- fast exponentiation
- operator overloading
- normalization of negative/out-of-range inputs
- exception handling for invalid operations

Example:

```cpp
IntMod<13> a(15);   // -> 2
IntMod<13> b(-3);   // -> 10

std::cout << a * b << std::endl;
```

### Univariate Polynomial Class (`UnivPoly<T>`)

A generic polynomial class that works with arbitrary coefficient types (most commonly `IntMod<P>`).

Supports:

- polynomial addition
- polynomial multiplication
- scalar multiplication
- evaluation
- degree computation
- normalization

Example:

```cpp
UnivPoly<IntMod<13>> p({...});
```

#### Lagrange Representation

`UnivPoly<T>` also supports conversion between:

- monomial representation
- Lagrange representation

with caching/precomputation for repeated interpolation.

This functionality is relevant for polynomial commitments, Reed-Solomon encoding, and zero-knowledge proof systems.

---

## Project Structure

```
.
├── vector_utilities.hpp
├── modular_arithmetic.hpp
├── modular_arithmetic.tpp
├── polynomials.hpp
├── polynomials.tpp
├── demo.cpp
├── CMakeLists.txt
```

---

## Build & Run

```bash
mkdir build
cd build
cmake ..
make
./demo
```

---

## Why I Built This

I created this project to:

- improve my C++ skills
- practice templates and operator overloading
- build stronger portfolio projects for software engineering roles

---

## Future Work

- Reed-Solomon encoding
- zk-proof primitives
- unit tests
