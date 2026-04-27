#pragma once
#include <vector>
#include <initializer_list>
#include <unordered_map>

// Intended to represent univariate polynomials
// - for coefficient types T that are constructible from integers, and
// - support field operations: +, -, *, / (or equivalently multiplicative inverse for nonzero elements).
//
// Not thread-safe: lazy initialization of shared Lagrange basis caches
// assumes single-threaded use or external synchronization.

template <class T=int>
class UnivPoly
{
    public:
        // --- Constructors ---
        // --- From monomial coefficients ---
        UnivPoly(); //zero polynomial
        explicit UnivPoly(std::initializer_list<T> list); //store monomial coefficients in normalized form
        explicit UnivPoly(std::vector<T> coeffvect); //store monomial coefficients in normalized form
        // --- From Lagrange coefficients ---
        static UnivPoly FromLagrangeCoeffs(std::initializer_list<T> apriorilagrangeinitlist); //store Lagrange coefficients (nonempty) and calculate monomial coefficients
        static UnivPoly FromLagrangeCoeffs(const std::vector<T>& apriorilagrangevect);
        static UnivPoly FromLagrangeCoeffs(std::vector<T>&& apriorilagrangevect);

        // --- Lagrange coefficients: Calculate and store for a given degree ---
        void SetLagrangeCoeffs(size_t degree) const;

        // --- Coefficient access ---
        const std::vector<T>& MonomialCoeffs() const noexcept;
        const std::vector<T>* TryGetLagrangeCoeffs(size_t deg) const noexcept;

        // --- Arithmetic operators --- 
        UnivPoly& operator+=(const UnivPoly& other);
        UnivPoly operator+(const UnivPoly& other) const;
        UnivPoly& operator*=(const UnivPoly& other);
        UnivPoly& operator*=(const T& num);
        UnivPoly operator*(const UnivPoly& other) const;
        UnivPoly operator*(const T& num) const;

        // --- Evaluate ---
        T Evaluate(const T& arg) const; //only for types that are constructible from 1 
        T Evaluate(int arg) const;

        // --- Basic properties ---
        size_t Degree() const noexcept; //0 for the constant 0 polynomial

    private:
        // --- Internal representation: Coefficients in monomial and Lagrange bases ---
        std::vector<T> monomCoeffs_; //always initialized
        mutable std::unordered_map<size_t, std::vector<T>> lagrangeCoeffsCache_; //(logical constness)

        // ---  Monomial Normalization ---        
        static std::vector<T> MonomNormalize_(std::vector<T> list); //Takes ownership of the vector, removes trailing zeros and makes it nonempty 
        static std::vector<T> MonomNormalize_(std::initializer_list<T> list);
        // ---  Lagrange Normalization ---        
        static std::vector<T> LagrangeNormalize_(std::vector<T> list); // Takes ownership of the vector and makes it nonempty
        static std::vector<T> LagrangeNormalize_(std::initializer_list<T> list); 

        // --- Monomial to Lagrange basis ---
        std::vector<T> CalculateLagrangeCoeffs_(size_t deg) const; //Calculate Lagrange coefficients for a given degree 

        // --- Lagrange Bases: Monomial coefficients ---
        static std::vector<bool> lagrangeBasesAprioriMonomCoeffsInitialized_; //Initialization of Lagrange bases
        static std::vector<std::vector<std::vector<T>>> lagrangeBasesAprioriMonomCoeffs_; //(a priori) monomial coefficients in Lagrange basis
        
        static std::vector<std::vector<T>> CalculateLagrangeBasisAprioriMonomCoeffs_(size_t lagrangecoeffvectsize); //Calculate the (a priori) monomial coefficients in Lagrange bases for the  corresponding degree
        static void UpdateLagrangeBasesAprioriMonomCoeffs_(size_t lagrangecoeffvectsize); //updates lagrangeBasesAprioriMonomCoeffsInitialized_ and lagrangeBasesAprioriMonomCoeffs_ if necessary

        // --- Lagrange to Monomial basis ---
        static std::vector<T> Lagrange2AprioriMonomCoeffs_(const std::vector<T>& lagrangecoeffvect); //Calculate (a priori) monomial coefficients from Lagrange coefficients

        // --- Internal helpers ---
        static std::vector<size_t> PuncturedInterval_(size_t n, size_t i); 
};

#include "polynomials.tpp"
