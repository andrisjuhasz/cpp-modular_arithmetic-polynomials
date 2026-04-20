#pragma once
#include <ostream>

template<int P>
class IntMod{
    private:
        int nMod_;
    public:
        //Constructors
        IntMod();
        explicit IntMod(int n);

        //Operators
        bool operator==(const IntMod& other) const;
        IntMod& operator+=(const IntMod& other);
        IntMod& operator++();
        IntMod operator+(const IntMod& other) const;
        IntMod& operator-=(const IntMod& other);
        IntMod operator-(const IntMod& other) const;
        IntMod operator-() const;
        IntMod& operator*=(const IntMod& other);
        IntMod operator*(const IntMod& other) const;
        IntMod operator*(int) const; 
        IntMod operator/(const IntMod& other) const;
        friend std::ostream& operator<< (std::ostream& os, const IntMod<P>& obj)
        {
            return os << obj.nMod_;
        } //friend keyword -> "free" function

        //Miscellaneous
        IntMod FastExponentiation(int exp) const;
        IntMod Inverse() const;
};

#include <modular_arithmetic.tpp>
