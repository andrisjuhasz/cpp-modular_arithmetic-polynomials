#include <iostream>
#include <vector>
#include "stdvector_extension.hpp"
#include "modular_arithmetic.hpp"
#include "polynomials.hpp"


template<class T>
void PrintMonomCoeffs(const UnivPoly<T>& poly) {
    std::cout << "Monomial coefficients: ";
    vectorext::operator<<(std::cout,poly.MonomialCoeffs());
    std::cout << '\n';
}

template<class T>
void PrintCachedLagrangeCoeffs(const UnivPoly<T>& poly, size_t deg) {
    if (const auto* coeffs = poly.TryGetLagrangeCoeffs(deg)) {
        std::cout << "Cached Lagrange coefficients for degree "
                  << deg << ": ";
        vectorext::operator<<(std::cout,*coeffs);
        std::cout << '\n';
    }
    else {
        std::cout << "No cached Lagrange coefficients for degree "
                  << deg << '\n';
    }
}

int main()
{
    // --- MODULAR ARITHMETICS ---
    std::cout << "MODULAR ARITHMETICS DEMO: \n";

    IntMod<7> a(10);
    IntMod<7> b(-4);
    IntMod<7> c(5);

    std::cout << "Constructor normalization:\n";
    std::cout << "a = " << a << '\n';
    std::cout << "b = " << b << '\n';
    std::cout << "c = " << c << "\n\n";

    std::cout << "Arithmetic operators:\n";
    std::cout << "a + c = " << a + c << '\n';
    std::cout << "a - c = " << a - c << '\n';
    std::cout << "a * c = " << a * c << "\n\n";

    std::cout << "Division and inverse:\n";
    std::cout << "c inverse = " << c.Inverse() << '\n';
    std::cout << "a / c = " << a / c << "\n\n";

    std::cout << "Increment operators:\n";
    std::cout << "a before ++a: " << a << '\n';
    std::cout << "++a = " << ++a << '\n';

    std::cout << "a before a++: " << a << '\n';
    std::cout << "a++ = " << a++ << '\n';
    std::cout << "a after a++: " << a << "\n\n";

    std::cout << "Comparison operators:\n";
    if (a == b)
        std::cout << "a == b\n";
    else
        std::cout << "a != b\n";

    std::cout << "\nException example:\n";
    try {
        IntMod<7> zero(0);
        std::cout << zero.Inverse() << '\n';
    }
    catch (const std::runtime_error& e) {
        std::cout << "Caught exception: " << e.what() << '\n';
    }

    std::cout << std::string(80, '-') << std::endl;

    // --- UNIVARIATE POLYNOMIALS ---
    std::cout << "UNIVARIATE POLYNOMIALS DEMO: \n";
    // --- poly1 ---
    std::cout <<"poly1: constructed from monomial coefficients (default): \n";
    UnivPoly<IntMod<13>> poly1;
    PrintMonomCoeffs(poly1);
    try{
        poly1.SetLagrangeCoeffs(3);
    } catch(const std::exception& e) {
        std::cout << e.what() << '\n';
    }
    PrintCachedLagrangeCoeffs(poly1,0);
    PrintCachedLagrangeCoeffs(poly1,3);

    std::cout << std::endl;

    // --- poly2 ---
    std::cout <<"poly2: constructed from monomial coefficients (initializer_list): \n";
    UnivPoly<IntMod<13>> poly2{IntMod<13>(15),IntMod<13>(-2),IntMod<13>(7)};
    PrintMonomCoeffs(poly2);
    try{
        poly2.SetLagrangeCoeffs(1);
    } catch(const std::exception& e) {
        std::cout << e.what() << '\n';
    }

    try{
        poly2.SetLagrangeCoeffs(4);
    } catch(const std::exception& e) {
        std::cout << e.what() << '\n';
    }
    PrintCachedLagrangeCoeffs(poly2,1);
    PrintCachedLagrangeCoeffs(poly2,4);

    std::cout << std::endl;

    // --- poly3 ---
    std::vector<IntMod<13>> vect3{IntMod<13>(4),IntMod<13>(18)};
    std::cout <<"poly3: constructed from monomial coefficients (vector): ";
    vectorext::operator<<(std::cout,vect3);
    std::cout << std::endl;
    UnivPoly<IntMod<13>> poly3(vect3);
    PrintMonomCoeffs(poly3);

    std::cout << std::endl;

    // --- poly4 ---
    std::initializer_list<IntMod<13>> list4{IntMod<13>(1),IntMod<13>(2),IntMod<13>(2),IntMod<13>(1)};
    std::cout<<"poly4: constructed from Lagrange coefficients (initializer_list): ";
    for (const auto& item : list4) {
        std::cout << item << ", ";
    }
    std::cout << std::endl;
    UnivPoly<IntMod<13>> poly4=UnivPoly<IntMod<13>>::FromLagrangeCoeffs(list4);
    // UnivPoly<IntMod<13>> poly4=UnivPoly<IntMod<13>>::FromLagrangeCoeffs({IntMod<13>(1),IntMod<13>(2),IntMod<13>(2),IntMod<13>(1)});
    PrintMonomCoeffs(poly4);
    PrintCachedLagrangeCoeffs(poly4,2);
    PrintCachedLagrangeCoeffs(poly4,3);

    std::cout << std::endl;

    // --- poly5 ---
    std::vector<IntMod<13>> vect5{IntMod<13>(-4),IntMod<13>(14),IntMod<13>(5)};
    std::cout<<"poly5: constructed from Lagrange coefficients (vector): ";
    vectorext::operator<<(std::cout,vect5);
    std::cout << std::endl;
    UnivPoly<IntMod<13>> poly5=UnivPoly<IntMod<13>>::FromLagrangeCoeffs(vect5);
    std::cout << std::endl;

    // --- Polynomial arithmetics ---
    std::cout << "Some operations on polynomials: " << std::endl;
    std::cout << "For poly3*pol4: ";
    PrintMonomCoeffs(poly3*poly4);
}

