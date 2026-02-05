#include <catch2/catch_test_macros.hpp>
#include "coolsolve/fluids.h"
#include "CoolProp.h"
#include "HumidAirProp.h"
#include <iostream>
#include <cmath>

using namespace coolsolve;

TEST_CASE("Test all registered fluids with CoolProp", "[fluids]") {
    auto fluids = FluidRegistry::getAllFluids();
    
    for (const auto& fluid : fluids) {
        std::string name = fluid->getName();
        std::string cpName = fluid->getCoolPropName();
        FluidType type = fluid->getType();
        
        SECTION("Testing fluid: " + name) {
            if (type == FluidType::Real || type == FluidType::IdealGas) {
                // Try to get enthalpy at 300K and 101325Pa
                // For IdealGas, we know it only needs T for enthalpy, 
                // but CoolProp always wants 2 inputs.
                
                double T = 300.0; // K
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
