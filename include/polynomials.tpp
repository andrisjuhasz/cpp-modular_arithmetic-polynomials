#include <iostream>
#include "stdvector_extension.hpp"

template<class T>
std::vector<T> UnivPoly<T>::MonomNormalize_(std::initializer_list<T> list){
    if(list.size()==0) { return std::vector<T>{T{}}; }
    else
    {
        const T* begin=list.begin();
        const T* end=list.end();
        while(begin!=end && *(end-1)==T{})
        {
            --end;
        }
        return begin==end ? std::vector<T>{T{}} : std::vector<T>(begin,end);
    }
}

template <class T>
std::vector<T> UnivPoly<T>::MonomNormalize_(std::vector<T> arg){
    if(arg.size()==0) { return std::vector<T>{T{}}; }
    else
    {
        while(arg.size()>1 && arg.back()==T{})
        {arg.pop_back();}
        return arg.size()==0? std::vector<T>{T{}} : arg;
    }
}

template<class T>
std::vector<T> UnivPoly<T>::LagrangeNormalize_(std::initializer_list<T> list){
    std::vector<T> result(list);
    if (result.empty()) {
        result.emplace_back();
    }
    return result;
}

template <class T>
std::vector<T> UnivPoly<T>::LagrangeNormalize_(std::vector<T> arg){
    if (arg.empty()) {
        arg.emplace_back();
    }
    return arg;
}

template<class T>
std::vector<bool> UnivPoly<T>::lagrangeBasesAprioriMonomCoeffsInitialized_;

template<class T>
std::vector<std::vector<std::vector<T>>> UnivPoly<T>::lagrangeBasesAprioriMonomCoeffs_;

template <class T>
UnivPoly<T>::UnivPoly(): monomCoeffVect_{T{}} {}

template <class T>
UnivPoly<T>::UnivPoly(std::initializer_list<T> coefflist):
    monomCoeffVect_(MonomNormalize_(coefflist)) {} //initializer_list not movable

template <class T>
UnivPoly<T>::UnivPoly(std::vector<T> coeffvect):
    monomCoeffVect_(MonomNormalize_(std::move(coeffvect))) {}

template <class T>
UnivPoly<T> UnivPoly<T>::FromLagrangeCoeffs(std::initializer_list<T> apriorilagrangeinitlist){
    std::vector<T> lagrangecoeffvect(LagrangeNormalize_(apriorilagrangeinitlist));
    UnivPoly<T> resultpoly;
    resultpoly.lagrangeCoeffVects_.try_emplace(lagrangecoeffvect.size()-1, lagrangecoeffvect);
    resultpoly.monomCoeffVect_=MonomNormalize_(Lagrange2AprioriMonomVect_(lagrangecoeffvect));
    return resultpoly;
}

template <class T>
UnivPoly<T> UnivPoly<T>::FromLagrangeCoeffs(const std::vector<T>& apriorilagrangevect){
    std::vector<T> lagrangecoeffvect(LagrangeNormalize_(std::move(apriorilagrangevect)));
    UnivPoly<T> resultpoly;
    resultpoly.lagrangeCoeffVects_.try_emplace(lagrangecoeffvect.size()-1, lagrangecoeffvect);
    resultpoly.monomCoeffVect_=MonomNormalize_(Lagrange2AprioriMonomVect_(lagrangecoeffvect));
    return resultpoly;
}

template <class T>
size_t UnivPoly<T>::Degree() const {
    return monomCoeffVect_.size()-1;
}

template <class T>
std::vector<T> UnivPoly<T>::PuncturedInterval_(size_t n, size_t i){
    std::vector<T> vec(n);
    for (size_t k = 0; k < i; ++k)
        vec[k] = T(k);
    for (size_t k = i; k < n; ++k)
        vec[k] = T(k + 1);
    return vec;
}

template<class T>
std::vector<std::vector<T>> UnivPoly<T>::CalculateLagrangeBasisAprioriMonomCoeffs_(size_t lagrangecoeffvectsize){
    std::vector<std::vector<T>> resultVectOfLagrBasApriMonomCoeffs;
    T tempDenominator;
    std::vector<T> tempLagrPolyMonomCoeffsVect;
    for(size_t i=0;i<lagrangecoeffvectsize;i++){
        tempDenominator=T(1);
        tempLagrPolyMonomCoeffsVect=std::vector{T(1)};
        for(auto element:PuncturedInterval_(lagrangecoeffvectsize-1,i)){
            tempDenominator*=(T(static_cast<int>(i))-element);
            tempLagrPolyMonomCoeffsVect=vectorext::operator*(tempLagrPolyMonomCoeffsVect,std::vector{-element,T(1)});
        }
        vectorext::MultiliedByScalar(tempLagrPolyMonomCoeffsVect,tempDenominator.Inverse());
        resultVectOfLagrBasApriMonomCoeffs.push_back(tempLagrPolyMonomCoeffsVect);
    }
    return resultVectOfLagrBasApriMonomCoeffs;
} 

//only call in Lagrange2AprioriMonomVect_ ensures that argument > 0
template<class T>
void UnivPoly<T>::UpdateLagrangeBasesAprioriMonomCoeffs_(size_t lagrangecoeffvectsize){
    size_t apriorideg=lagrangecoeffvectsize-1;
    if (lagrangeBasesAprioriMonomCoeffsInitialized_.size() <= apriorideg){
        lagrangeBasesAprioriMonomCoeffsInitialized_.resize(apriorideg + 1, false);
        lagrangeBasesAprioriMonomCoeffs_.resize(apriorideg + 1);
    }
    if (!lagrangeBasesAprioriMonomCoeffsInitialized_[apriorideg]) {
        lagrangeBasesAprioriMonomCoeffs_[apriorideg]=CalculateLagrangeBasisAprioriMonomCoeffs_(lagrangecoeffvectsize);
        lagrangeBasesAprioriMonomCoeffsInitialized_[apriorideg] = true;
    }
}

//it will only be called by SetLagrangeCoeffs
template<class T>
std::vector<T> UnivPoly<T>::CalculateLagrangeCoeffs_(size_t deg) const{
    std::vector<T> resultVect(deg+1);
    for(int arg(0);arg<=static_cast<int>(deg);++arg)
        resultVect[arg] = Evaluate(arg);
    return resultVect;
}

template<class T>
void UnivPoly<T>::SetLagrangeCoeffs(size_t deg) const{
    if (deg < this->Degree()) {
        std::cerr << "Warning: SetLagrangeCoeffs argument ("<< deg<< ") is below allowed threshold (degree)" << std::endl;
        return;
    }
    lagrangeCoeffVects_.try_emplace(deg, CalculateLagrangeCoeffs_(deg));
}

//only called with proper (LagrangeNormalize_'d) argument 
template<class T>
std::vector<T> UnivPoly<T>::Lagrange2AprioriMonomVect_(const std::vector<T>& lagrangecoeffvect){
    size_t lagrangecoeffvectsize(lagrangecoeffvect.size());
    UpdateLagrangeBasesAprioriMonomCoeffs_(lagrangecoeffvectsize);
    std::vector<std::vector<T>> lagrangebasismonomcoeffs(lagrangeBasesAprioriMonomCoeffs_[lagrangecoeffvectsize-1]);
    std::vector<T> resultmonomcoeffs(lagrangecoeffvectsize,T{});
    for(size_t i=0;i<lagrangecoeffvectsize;i++){
        for(size_t j=0;j<lagrangecoeffvectsize;j++){
            resultmonomcoeffs[j]+=lagrangecoeffvect[i]*lagrangebasismonomcoeffs[i][j];
        }
    }
    return resultmonomcoeffs;
}

template<class T>
UnivPoly<T> UnivPoly<T>::operator+(const UnivPoly<T>& other) const{
     const std::vector<T>& longer =
        (monomCoeffVect_.size() > other.monomCoeffVect_.size()) ?
        monomCoeffVect_ : other.monomCoeffVect_;
     const std::vector<T>& shorter =
        (monomCoeffVect_.size() > other.monomCoeffVect_.size()) ?
        other.monomCoeffVect_ : monomCoeffVect_;
    std::vector<T> result = longer;  // copy longer
    for (std::size_t i = 0; i < shorter.size(); ++i) result[i] += shorter[i];
    return UnivPoly<T>(result);
}

template<class T>
UnivPoly<T>& UnivPoly<T>::operator+=(const UnivPoly<T>& other){
     *this = *this + other;
    return *this;
}

template<class T>
UnivPoly<T> UnivPoly<T>::operator*(const UnivPoly<T>& other) const{
    std::vector<T> resultCoefficientList(monomCoeffVect_.size()+other.monomCoeffVect_.size()-1);
    for (std::size_t i=0; i<monomCoeffVect_.size();++i)
        for(std::size_t j=0; j<other.monomCoeffVect_.size();++j)
            resultCoefficientList[i+j]+=monomCoeffVect_[i]*other.monomCoeffVect_[j];
    return UnivPoly<T>(resultCoefficientList);
}

template<class T>
UnivPoly<T>& UnivPoly<T>::operator*=(const UnivPoly<T>& other){
    *this = *this * other;
    return *this;
}

template<class T>
UnivPoly<T> UnivPoly<T>::operator*(const T& num) const{
    UnivPoly<T> temp=*this;
    temp*=num;
    return temp;
}

template<class T>
UnivPoly<T>& UnivPoly<T>::operator*=(const T& num){
    for (auto& element : monomCoeffVect_) element *= num;
    return *this;
}

template<class T>
T UnivPoly<T>::Evaluate(int arg) const{
    int power(1);
    T result{};
    for (const auto& coefficient : monomCoeffVect_) {
        result+=coefficient*power; 
        power*=arg;
    }
    return result;
}

template<class T>
T UnivPoly<T>::Evaluate(const T& arg) const{
    T power(1); 
    T result{};
    for (const auto& coefficient : monomCoeffVect_) {
        result+=coefficient*power; 
        power*=arg;
    }
   return result;
}

template<class T>
void UnivPoly<T>::PrintMonomCoeffs() const{
    vectorext::operator<<(std::cout,monomCoeffVect_);
    std::cout << std::endl;
}

template<class T>
void UnivPoly<T>::PrintLagrangeCoeffVect(size_t deg) const{
    if (lagrangeCoeffVects_.contains(deg)) {
        vectorext::operator<<(std::cout,lagrangeCoeffVects_[deg]);
    }
    else std::cout << "Lagrange Coefficient List Uninitialized";
    std::cout << std::endl;
}

template<class T>
void UnivPoly<T>::PrintCalculatedLagrangeCoeffVects() const{
    for (const auto& [deg, lagrangecoeffvect] : lagrangeCoeffVects_) {
        vectorext::operator<<(std::cout, lagrangecoeffvect);
        std::cout << std::endl;
    }
}

