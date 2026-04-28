#include <string>
#include <stdexcept>
#include "vector_utilities.hpp"

template<class T>
std::vector<bool> UnivPoly<T>::lagrangeBasesAprioriMonomCoeffsInitialized_;

template<class T>
std::vector<std::vector<std::vector<T>>> UnivPoly<T>::lagrangeBasesAprioriMonomCoeffs_;

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
std::vector<T> UnivPoly<T>::MonomNormalize_(std::vector<T> coeffs){
    if(coeffs.size()==0) { return std::vector<T>{T{}}; }
    else
    {
        while(coeffs.size()>1 && coeffs.back()==T{})
        {coeffs.pop_back();}
        return coeffs.size()==0? std::vector<T>{T{}} : coeffs;
    }
}

template<class T>
std::vector<T> UnivPoly<T>::LagrangeNormalize_(std::initializer_list<T> coefflist){
    std::vector<T> result(coefflist);
    if (result.empty()) {
        result.emplace_back();
    }
    return result;
}

template <class T>
std::vector<T> UnivPoly<T>::LagrangeNormalize_(std::vector<T> coeffs){
    if (coeffs.empty()) {
        coeffs.emplace_back();
    }
    return coeffs;
}

template <class T>
UnivPoly<T>::UnivPoly(): monomCoeffs_{T{}} {}

template <class T>
UnivPoly<T>::UnivPoly(std::initializer_list<T> coefflist):
    monomCoeffs_(MonomNormalize_(coefflist)) {} //initializer_list not movable

template <class T>
UnivPoly<T>::UnivPoly(std::vector<T> coeffvect):
    monomCoeffs_(MonomNormalize_(std::move(coeffvect))) {}

template <class T>
UnivPoly<T> UnivPoly<T>::FromLagrangeCoeffs(std::initializer_list<T> apriorilagrangeinitlist){
    std::vector<T> lagrangecoeffvect(LagrangeNormalize_(apriorilagrangeinitlist));
    UnivPoly<T> resultpoly;
    resultpoly.lagrangeCoeffsCache_.try_emplace(lagrangecoeffvect.size()-1, lagrangecoeffvect);
    resultpoly.monomCoeffs_=MonomNormalize_(Lagrange2AprioriMonomCoeffs_(lagrangecoeffvect));
    return resultpoly;
}

template<class T>
UnivPoly<T> UnivPoly<T>::FromLagrangeCoeffs(const std::vector<T>& apriorilagrangevect) {
    std::vector<T> normalized = LagrangeNormalize_(apriorilagrangevect);
    UnivPoly<T> resultpoly;
    resultpoly.monomCoeffs_=MonomNormalize_(Lagrange2AprioriMonomCoeffs_(normalized));
    resultpoly.lagrangeCoeffsCache_.try_emplace(normalized.size()-1,std::move(normalized));
    return resultpoly;
}

template<class T>
UnivPoly<T> UnivPoly<T>::FromLagrangeCoeffs(std::vector<T>&& apriorilagrangevect) {
    std::vector<T> normalized = LagrangeNormalize_(std::move(apriorilagrangevect));
    UnivPoly<T> resultpoly;
    resultpoly.monomCoeffs_ =
        MonomNormalize_(Lagrange2AprioriMonomCoeffs_(normalized));
    resultpoly.lagrangeCoeffsCache_.try_emplace(
        normalized.size() - 1,
        std::move(normalized)
    );
    return resultpoly;
}

template <class T>
size_t UnivPoly<T>::Degree() const noexcept {
    return monomCoeffs_.size()-1;
}

template <class T>
std::vector<size_t> UnivPoly<T>::PuncturedInterval_(size_t n, size_t i){
    std::vector<size_t> vec(n);
    for (size_t k = 0; k < i; ++k)
        vec[k] = k;
    for (size_t k = i; k < n; ++k)
        vec[k] = k + 1;
    return vec;
}

template<class T>
std::vector<std::vector<T>> UnivPoly<T>::CalculateLagrangeBasisAprioriMonomCoeffs_(size_t lagrangecoeffvectsize){
    std::vector<std::vector<T>> resultVectOfLagrBasApriMonomCoeffs;
    resultVectOfLagrBasApriMonomCoeffs.reserve(lagrangecoeffvectsize);
    T tempDenominator;
    std::vector<T> tempLagrPolyMonomCoeffsVect;
    for(size_t i=0;i<lagrangecoeffvectsize;++i){
        tempDenominator=T(1);
        tempLagrPolyMonomCoeffsVect=std::vector{T(1)};
        for(auto element:PuncturedInterval_(lagrangecoeffvectsize-1,i)){
            tempDenominator*=T(static_cast<int>(i)-static_cast<int>(element));
            tempLagrPolyMonomCoeffsVect=vectorext::operator*(tempLagrPolyMonomCoeffsVect,
                   std::vector{T(-static_cast<int>(element)),T(1)});
        }
        vectorext::MultiliedByScalar(tempLagrPolyMonomCoeffsVect,tempDenominator.Inverse());
        resultVectOfLagrBasApriMonomCoeffs.push_back(tempLagrPolyMonomCoeffsVect);
    }
    return resultVectOfLagrBasApriMonomCoeffs;
} 

//only call in Lagrange2AprioriMonomCoeffs_ ensures that argument > 0
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
    if (deg < Degree()) {
        throw std::invalid_argument("SetLagrangeCoeffs: deg (" +
                std::to_string(deg) +
                ") is smaller than polynomial degree (" +
                std::to_string(Degree()) +
                ")"
                );
    }
    lagrangeCoeffsCache_.try_emplace(deg, CalculateLagrangeCoeffs_(deg));
}

//only called with proper (LagrangeNormalize_'d) argument 
template<class T>
std::vector<T> UnivPoly<T>::Lagrange2AprioriMonomCoeffs_(const std::vector<T>& lagrangecoeffvect){
    size_t lagrangecoeffvectsize(lagrangecoeffvect.size());
    UpdateLagrangeBasesAprioriMonomCoeffs_(lagrangecoeffvectsize);
    const auto& lagrangebasismonomcoeffs=lagrangeBasesAprioriMonomCoeffs_[lagrangecoeffvectsize-1];
    std::vector<T> resultmonomcoeffs(lagrangecoeffvectsize,T{});
    for(size_t i=0;i<lagrangecoeffvectsize;++i){
        for(size_t j=0;j<lagrangecoeffvectsize;++j){
            resultmonomcoeffs[j]+=lagrangecoeffvect[i]*lagrangebasismonomcoeffs[i][j];
        }
    }
    return resultmonomcoeffs;
}

template<class T>
UnivPoly<T> UnivPoly<T>::operator+(const UnivPoly<T>& other) const{
     const std::vector<T>& longer =
        (monomCoeffs_.size() > other.monomCoeffs_.size()) ?
        monomCoeffs_ : other.monomCoeffs_;
     const std::vector<T>& shorter =
        (monomCoeffs_.size() > other.monomCoeffs_.size()) ?
        other.monomCoeffs_ : monomCoeffs_;
    std::vector<T> result = longer;  // copy longer
    for (std::size_t i = 0; i < shorter.size(); ++i) result[i] += shorter[i];
    return UnivPoly<T>(result);
}

template<class T>
UnivPoly<T>& UnivPoly<T>::operator+=(const UnivPoly<T>& other) {
    if (other.Degree() > Degree()) {
        monomCoeffs_.resize(other.Degree() + 1, T{});
    }
    for (size_t i = 0; i <= other.Degree(); ++i) {
        monomCoeffs_[i] += other.monomCoeffs_[i];
    }
    monomCoeffs_ = MonomNormalize_(std::move(monomCoeffs_));
    lagrangeCoeffsCache_.clear();
    return *this;
}

template<class T>
UnivPoly<T> UnivPoly<T>::operator*(const UnivPoly<T>& other) const{
    std::vector<T> resultCoefficientList(monomCoeffs_.size()+other.monomCoeffs_.size()-1);
    for (std::size_t i=0; i<monomCoeffs_.size();++i)
        for(std::size_t j=0; j<other.monomCoeffs_.size();++j)
            resultCoefficientList[i+j]+=monomCoeffs_[i]*other.monomCoeffs_[j];
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
UnivPoly<T>& UnivPoly<T>::operator*=(const T& num) {
    for (auto& c : monomCoeffs_) {
        c *= num;
    }
    monomCoeffs_ = MonomNormalize_(std::move(monomCoeffs_));
    lagrangeCoeffsCache_.clear();
    return *this;
}

template<class T>
T UnivPoly<T>::Evaluate(const T& arg) const{
    T power(1); 
    T result{};
    for (const auto& coefficient : monomCoeffs_) {
        result+=coefficient*power; 
        power*=arg;
    }
   return result;
}

template<class T>
T UnivPoly<T>::Evaluate(int arg) const{
    return Evaluate(T(arg));
}

template<class T>
const std::vector<T>& UnivPoly<T>::MonomialCoeffs() const noexcept {
    return monomCoeffs_;
}

template<class T>
const std::vector<T>* UnivPoly<T>::TryGetLagrangeCoeffs(size_t deg) const noexcept {
    auto it = lagrangeCoeffsCache_.find(deg);
    if (it == lagrangeCoeffsCache_.end()) {
        return nullptr;
    }
    return &it->second;
}
