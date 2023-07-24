cmake_minimum_required(VERSION 3.25)
project(LearningProject)

set(CMAKE_CXX_STANDARD 14)

include_directories(BathTemperature)
include_directories(BrownianMotion)
include_directories(Convergence)
include_directories(DomainSize)
include_directories(FirstTry)
include_directories(GPUTries)
include_directories(HarmonicOscillator)
include_directories(Header)
include_directories(LatticeModel)
include_directories(Odeint)
include_directories(SimpleInteraction)
include_directories(Testing)
include_directories(TwoStepMethod)

add_executable(LearningProject
        BathTemperature/Bath.cpp
        BathTemperature/Bath.h
        BrownianMotion/main.cpp
        BrownianMotion/main.h
        Convergence/ConvergenceCheck.cpp
        Convergence/ConvergenceCheckDiffusion.cpp
        Convergence/ConvergenceCheckDiffusion.h
        DomainSize/CalcCorrelationFunction.cpp
        DomainSize/CalcDomainSize.cpp
        DomainSize/CalcDomainSize.h
        DomainSize/FourierLatticeTrafo.cpp
        DomainSize/FourierLatticeTrafo.h
        DomainSize/matplotlibcpp.h
        DomainSize/matplotlibtest.cpp
        FirstTry/first_try.cpp
        FirstTry/kzm_silicon_2d.cpp
        FirstTry/Lattice.cpp
        FirstTry/main.cpp
        FirstTry/matrices.cpp
        FirstTry/matrices.h
        FirstTry/point_type.hpp
        FirstTry/rng.cpp
        FirstTry/rng.h
        FirstTry/SolarSystem.cpp
        FirstTry/Testclass.cpp
        FirstTry/Testclass.h
        FullInteraction/main.cpp
        GPUTries/GPUTry.cpp
        GPUTries/GPUTry.cu
        GPUTries/GPUTry.h
        HarmonicOscillator/main.cpp
        HarmonicOscillator/main.h
        "Header/Helpfunctions and Classes.h"
        Header/Solvers.h
        Header/Systems.h
        LatticeModel/main.cpp
        LatticeModel/main.h
        Odeint/OdeintTest.cpp
        Odeint/OdeintTest.h
        SimpleInteraction/main.cpp
        SimpleInteraction/main.h
        SimpleInteraction/RightAndDown.cpp
        SimpleInteraction/RightAndDown.h
        Testing/EigenCurveFitting.cpp
        Testing/EigenCurveFitting.h
        Testing/fftwMinimal.cpp
        Testing/fftwMinimal.h
        Testing/Polymorphism.cpp
        Testing/Polymorphism.h
        TwoStepMethod/main.cpp
        TwoStepMethod/main.h CurveFitting/Curvefitting.cpp CurveFitting/Curvefitting.h DomainSize/QuenchCorrLength.cpp DomainSize/QuenchCorrLength.h DomainSize/LatticeFFT.cpp DomainSize/LatticeFFT.h DomainSize/LatticeFFTError.cpp)
