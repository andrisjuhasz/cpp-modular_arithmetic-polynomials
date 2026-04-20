#pragma once
#include <vector>
#include <initializer_list>
#include <unordered_map>

//only for types T that are constructible from 1 (see e.g. Evaluate)
template <class T=int>
class UnivPoly
{
    private:
        std::vector<T> monomCoeffVect_; //always initialized
        mutable std::unordered_map<size_t, std::vector<T>> lagrangeCoeffVects_; //(logical constness)

        //Store coefficient vectors in a "normalize" form, in particular nonempty:
        static std::vector<T> MonomNormalize_(std::vector<T> list);
        static std::vector<T> MonomNormalize_(std::initializer_list<T> list);//returns argument minus trailing zeros and {T{}} if it would be empty
        static std::vector<T> LagrangeNormalize_(std::vector<T> list);//returns {T{}} if argument is empty
        static std::vector<T> LagrangeNormalize_(std::initializer_list<T> list); 

        //Monom2Lagrange:
        //Calculate Lagrange coefficients for a given degree (from monomial coefficients)
        std::vector<T> CalculateLagrangeCoeffs_(size_t deg) const; //Used in SetLagrangeCoeffs

        //Lagrange2Monom:
        //Calculate, update and store the (a priori) monomial coefficients of Lagrange bases polynomials:
        static std::vector<bool> lagrangeBasesAprioriMonomCoeffsInitialized_;
        static std::vector<std::vector<std::vector<T>>> lagrangeBasesAprioriMonomCoeffs_;

        static std::vector<T> PuncturedInterval_(size_t n, size_t i); //used in CalculateLagrangeBasisAprioriMonomCoeffs_
        static std::vector<std::vector<T>> CalculateLagrangeBasisAprioriMonomCoeffs_(size_t lagrangecoeffvectsize); //for the  corresponding degree(=lagrangecoeffvectsize-1)
        static void UpdateLagrangeBasesAprioriMonomCoeffs_(size_t lagrangecoeffvectsize); //updates lagrangeBasesAprioriMonomCoeffsInitialized_ and lagrangeBasesAprioriMonomCoeffs_ if necessary

        //Calculate (a priori) monomial coefficients (from Lagrange coefficients)
        static std::vector<T> Lagrange2AprioriMonomVect_(const std::vector<T>& lagrangecoeffvect);

    public:
        //Constructors: In monomial basis
        UnivPoly();
        UnivPoly(std::initializer_list<T> list); //not standard: MonomNormalize_
        UnivPoly(std::vector<T> coefficientList); //not standars: MonomNormalize_
        //??? UnivPoly(const std::vector<T>& coefficientList) ~ FromLagrangeCoeffs;
        //Create polynomial from Lagrange coefficients (and calculate monomial coefficients)
        static UnivPoly FromLagrangeCoeffs(std::initializer_list<T> apriorilagrangeinitlist);
        static UnivPoly FromLagrangeCoeffs(const std::vector<T>& apriorilagrangevect);

        //Calculate and store Lagrange coefficients for a given degree
        void SetLagrangeCoeffs(size_t degree) const;

        //Operators
        UnivPoly& operator+=(const UnivPoly& other);
        UnivPoly operator+(const UnivPoly& other) const;
        UnivPoly& operator*=(const UnivPoly& other);
        UnivPoly& operator*=(const T& num);
        UnivPoly operator*(const UnivPoly& other) const;
        UnivPoly operator*(const T& num) const;

        //Miscellaneous
        size_t Degree() const; //0 for the constant 0 polynomial
        T Evaluate(const T& arg) const; //only for types that are constructible from 1 
        T Evaluate(int arg) const;
        void PrintMonomCoeffs() const;
        void PrintLagrangeCoeffVect(size_t deg) const;
        void PrintCalculatedLagrangeCoeffVects() const;
};

#include "polynomials.tpp"
