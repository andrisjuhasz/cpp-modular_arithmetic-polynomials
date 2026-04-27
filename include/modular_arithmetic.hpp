#pragma once
#include <ostream>

// IntMod<P> represents integers modulo P.
// 
// IMPORTANT:
// This implementation assumes that P is a prime number.
// 
// Using a non-prime modulus may lead to incorrect behavior for division
// and inverse operations.

template<int P>
class IntMod{
    static_assert(P > 1, "Modulus must be greater than 1");

    public:
        // --- Constructors ---
        IntMod();
        explicit IntMod(int n);

        // --- Comparison operators ---
        bool operator==(const IntMod& other) const;
        bool operator!=(const IntMod<P>& other) const;
            
        // --- Arithmetic operators --- 
        IntMod& operator+=(const IntMod& other);
        IntMod& operator++();
        IntMod operator++(int);
        IntMod operator+(const IntMod& other) const;

        IntMod& operator-=(const IntMod& other);
        IntMod operator-(const IntMod& other) const;
        IntMod operator-() const;

        IntMod& operator*=(const IntMod& other);
        IntMod operator*(const IntMod& other) const;

        IntMod Inverse() const;
        IntMod operator/(const IntMod& other) const;

        IntMod FastExponentiation(int exp) const;

        // --- Output stream operator ---
        friend std::ostream& operator<< (std::ostream& os, const IntMod<P>& obj)
        {
            return os << obj.nMod_;
        }

        // --- Access value ---
        int Value() const;

    private:
        int nMod_;
};

#include <modular_arithmetic.tpp>
