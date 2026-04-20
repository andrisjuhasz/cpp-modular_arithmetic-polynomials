#include <iostream>
#include <vector>
#include "stdvector_extension.hpp"
#include "modular_arithmetic.hpp"
#include "polynomials.hpp"

int main()
{
    std::cout <<"poly1: constructed from monomial coefficients (default): " << std::endl;
    UnivPoly<IntMod<13>> poly1;
    poly1.PrintMonomCoeffs();
    poly1.SetLagrangeCoeffs(3);
    poly1.PrintCalculatedLagrangeCoeffVects();
    poly1.PrintLagrangeCoeffVect(0);
    poly1.PrintLagrangeCoeffVect(3);
    std::cout << std::endl;

    std::cout <<"poly2: constructed from monomial coefficients (initializer_list): " << std::endl;
    UnivPoly<IntMod<13>> poly2{IntMod<13>(15),IntMod<13>(-2),IntMod<13>(7)};
    poly2.PrintMonomCoeffs();
    poly2.SetLagrangeCoeffs(1);
    std::cout << std::endl;

    std::cout <<"poly3: constructed from monomial coefficients (vector): ";
    std::vector<IntMod<13>> vect3{IntMod<13>(4),IntMod<13>(18)};
    vectorext::operator<<(std::cout,vect3);
    std::cout<<std::endl;
    UnivPoly<IntMod<13>> poly3(vect3);
    poly3.PrintMonomCoeffs();
    poly3.SetLagrangeCoeffs(1);
    poly3.SetLagrangeCoeffs(4);
    poly3.PrintCalculatedLagrangeCoeffVects();
    std::cout << std::endl;

    std::cout<<"poly4: constructed from Lagrange coefficients (initializer_list): "<<std::endl;
    UnivPoly<IntMod<13>> poly4=UnivPoly<IntMod<13>>::FromLagrangeCoeffs({IntMod<13>(1),IntMod<13>(2),IntMod<13>(2),IntMod<13>(1)});
    poly4.PrintCalculatedLagrangeCoeffVects();
    std::cout << poly4.Degree() << std::endl;
    poly4.PrintMonomCoeffs();
    std::cout << std::endl;

    std::cout<<"poly5: constructed from Lagrange coefficients (vector): ";
    std::vector<IntMod<13>> vect5{IntMod<13>(-4),IntMod<13>(14),IntMod<13>(5)};
    vectorext::operator<<(std::cout,vect5);
    std::cout << std::endl;
    UnivPoly<IntMod<13>> poly5=UnivPoly<IntMod<13>>::FromLagrangeCoeffs(vect5);
    poly5.PrintCalculatedLagrangeCoeffVects();
    poly5.PrintMonomCoeffs();
    std::cout << std::endl;

     std::cout << "Some operations on polynomials: " << std::endl;
     std::cout << "Monomial coefficients of poly3*pol4: ";
     (poly3*poly4).PrintMonomCoeffs();
     std::cout << std::endl;
}

