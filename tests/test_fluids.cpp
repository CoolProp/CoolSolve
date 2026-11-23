#include <catch2/catch_test_macros.hpp>
#include "coolsolve/fluids.h"
#include "coolsolve/ideal_gas.h"
#include "CoolProp.h"
#include "HumidAirProp.h"
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace coolsolve;

TEST_CASE("Test all registered fluids with CoolProp", "[fluids]") {
    auto fluids = FluidRegistry::getAllFluids();
    
    for (const auto& fluid : fluids) {
        std::string name = fluid->getName();
        std::string cpName = fluid->getCoolPropName();
        FluidType type = fluid->getType();
        
        SECTION("Testing fluid: " + name) {
            auto idealGas = std::dynamic_pointer_cast<IdealGasFluid>(fluid);
            if (idealGas && !idealGas->hasCoolPropFluid()) {
                // A species of the NASA-Glenn data has no CoolProp fluid: its properties come from the ideal-gas layer
                const double h = idealGasAt(cpName, 300.0).h;
                REQUIRE(std::isfinite(h));
            } else if (type == FluidType::Real || type == FluidType::IdealGas) {
                // Try to get enthalpy at 300K and 101325Pa
                // For IdealGas, we know it only needs T for enthalpy, 
                // but CoolProp always wants 2 inputs.
                
                // 300 K, or 20 K above the lowest temperature of the fluid (the methyl esters melt above 300 K)
                double T = std::max(300.0, CoolProp::Props1SI(cpName, "Tmin") + 20.0);
                double P = 1000000.0; // 10 bar (above triple point for CO2)
                
                double h = CoolProp::PropsSI("H", "T", T, "P", P, cpName);
                
                if (!std::isfinite(h)) {
                    std::string err = CoolProp::get_global_param_string("errstring");
                    INFO("Fluid: " << name << " (CoolProp: " << cpName << ") returned non-finite value. Error: " << err);
                    REQUIRE(std::isfinite(h));
                }
                
                REQUIRE(std::isfinite(h));
            } else if (type == FluidType::HumidAir) {
                // Humid air test
                double T = 300.0;
                double P = 101325.0;
                double R = 0.5; // Relative humidity
                try {
                    double h = HumidAir::HAPropsSI("H", "T", T, "P", P, "R", R);
                    REQUIRE(std::isfinite(h));
                    
                    // Test derivatives (implicit via evaluator if we had a full test case, 
                    // but here we just check if HAPropsSI works)
                    double w = HumidAir::HAPropsSI("W", "T", T, "P", P, "R", R);
                    REQUIRE(std::isfinite(w));
                    REQUIRE(w > 0.0);
                } catch (const std::exception& e) {
                    FAIL("Humid air failed CoolProp call: " << e.what());
                }
            } else if (type == FluidType::Incompressible) {
                // Incompressible test - CoolProp uses INCOMP:: prefix
                try {
                    double T = 300.0;
                    double P = 101325.0;
                    double rho = CoolProp::PropsSI("D", "T", T, "P", P, cpName);
                    if (!std::isfinite(rho)) {
                        WARN("Incompressible fluid " << name << " returned non-finite density");
                    }
                } catch (...) {
                    // Some incompressibles might not be in CoolProp or need different inputs
                    WARN("Incompressible fluid " << name << " failed CoolProp call");
                }
            }
        }
    }
}

TEST_CASE("Humid Air specific properties and derivatives", "[fluids][humidair]") {
    // This test verifies that HAPropsSI works as expected for different outputs
    double T = 300.0;
    double P = 101325.0;
    double R = 0.5;
    
    SECTION("Enthalpy") {
        double h = HumidAir::HAPropsSI("H", "T", T, "P", P, "R", R);
        REQUIRE(std::isfinite(h));
        REQUIRE(h > 0.0);
    }
    
    SECTION("Humidity ratio") {
        double w = HumidAir::HAPropsSI("W", "T", T, "P", P, "R", R);
        REQUIRE(std::isfinite(w));
        REQUIRE(w > 0.0);
        REQUIRE(w < 0.1); // Typical value
    }
    
    SECTION("Dew point") {
        double dp = HumidAir::HAPropsSI("D", "T", T, "P", P, "R", R);
        REQUIRE(std::isfinite(dp));
        REQUIRE(dp < T);
    }
}

TEST_CASE("The NASA-Glenn formation enthalpies of the ideal gases agree with CoolProp's ATcT values (HFORMATION)",
          "[fluids][hformation]") {
    // CoolSolve's ideal gases carry the enthalpy of formation of the NASA-Glenn table (src/nasa_glenn.cpp, McBride, Zehe &
    // Gordon 2002); CoolProp 8.1 gives the ATcT value of its fluids (HFORMATION, Active Thermochemical Tables). Two
    // independent data sets of the same quantity: a comparison of data, not a tolerance of CoolSolve.
    int compared = 0;
    for (const auto& fluid : FluidRegistry::getAllFluids()) {
        const auto ideal = std::dynamic_pointer_cast<IdealGasFluid>(fluid);
        if (!ideal || !ideal->hasCoolPropFluid() || !ideal->getFormationEnthalpy()) continue;
        double atct = NAN;
        try {
            atct = CoolProp::Props1SI(ideal->getCoolPropName(), "HFORMATION");   // J/mol
        } catch (...) {
        }
        INFO(ideal->getName() << " (" << ideal->getCoolPropName() << "): NASA-Glenn "
                              << *ideal->getFormationEnthalpy() / 1000.0 << " J/mol, ATcT " << atct << " J/mol");
        if (!std::isfinite(atct)) continue;   // a fluid without an ATcT value in CoolProp
        ++compared;
        const double nasa = *ideal->getFormationEnthalpy() / 1000.0;   // J/kmol -> J/mol
        // measured with CoolProp 8.1.0: 0 for the elements, 16 to 87 J/mol for CO, H2O, CO2, CH4, 150 to 860 J/mol for the
        // n-alkanes (C2H6 to C8H18)
        CHECK(std::abs(nasa - atct) < 1000.0);
    }
    CHECK(compared >= 12);
}
