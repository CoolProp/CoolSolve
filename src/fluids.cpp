#include "coolsolve/fluids.h"
#include <algorithm>
#include <stdexcept>
#include <set>

namespace coolsolve {

std::map<std::string, std::shared_ptr<Fluid>> FluidRegistry::registry_;
bool FluidRegistry::initialized_ = false;

void FluidRegistry::initialize() {
    if (initialized_) return;
    
    // --- Real Fluids ---
    auto addReal = [](const std::string& name, const std::string& cpName) {
        auto fluid = std::make_shared<RealFluid>(name, cpName);
        registry_[name] = fluid;
        // Add lowercase version for case-insensitive lookup
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower != name) registry_[lower] = fluid;
    };
    
    addReal("Water", "Water");
    addReal("Steam", "Water");
    addReal("Steam_IAPWS", "Water");
    addReal("R718", "Water");
    addReal("R134a", "R134a");
    addReal("R245fa", "R245fa");
    addReal("Ammonia", "Ammonia");
    addReal("R717", "Ammonia");
    addReal("Propane", "Propane");
    addReal("R290", "Propane");
    addReal("Ethane", "Ethane");
    addReal("Methane", "Methane");
    addReal("Nitrogen", "Nitrogen");
    addReal("Oxygen", "Oxygen");
    addReal("Argon", "Argon");
    addReal("CarbonDioxide", "CarbonDioxide");
    addReal("R744", "CarbonDioxide");
    addReal("Toluene", "Toluene");
    addReal("R22", "R22");
    addReal("R32", "R32");
    addReal("R404A", "R404A");
    addReal("R407C", "R407C");
    addReal("R410A", "R410A");
    addReal("R507A", "R507A");
    addReal("R600a", "IsoButane");
    addReal("IsoButane", "IsoButane");
    addReal("n-Butane", "n-Butane");
    addReal("R600", "n-Butane");
    addReal("n-Pentane", "n-Pentane");
    addReal("IsoPentane", "Isopentane");
    addReal("n-Hexane", "n-Hexane");
    addReal("IsoHexane", "Isohexane");
    addReal("n-Heptane", "n-Heptane");
    addReal("n-Octane", "n-Octane");
    addReal("IsoOctane", "Isooctane");
    addReal("CycloHexane", "CycloHexane");
    addReal("Acetone", "Acetone");
    addReal("Benzene", "Benzene");
    addReal("Butene", "1-Butene");
    addReal("cis-2-Butene", "cis-2-Butene");
    addReal("trans-2-Butene", "trans-2-Butene");
    addReal("IsoButene", "IsoButene");
    addReal("CycloPentane", "CycloPentane");
    addReal("D4", "D4");
    addReal("D5", "D5");
    addReal("D6", "D6");
    addReal("Deuterium", "Deuterium");
    addReal("HeavyWater", "HeavyWater");
    addReal("DiethylEther", "DiethylEther");
    addReal("DimethylEther", "DimethylEther");
    addReal("EthylBenzene", "EthylBenzene");
    addReal("Ethylene", "Ethylene");
    addReal("Fluorine", "Fluorine");
    addReal("Helium", "Helium");
    addReal("H2", "Hydrogen");
    addReal("Hydrogen", "Hydrogen");
    addReal("HydrogenChloride", "HydrogenChloride");
    addReal("HydrogenSulfide", "HydrogenSulfide");
    addReal("Krypton", "Krypton");
    addReal("MDM", "MDM");
    addReal("MD2M", "MD2M");
    addReal("MD3M", "MD3M");
    addReal("MD4M", "MD4M");
    addReal("MM", "MM");
    addReal("Neon", "Neon");
    addReal("Neopentane", "Neopentane");
    addReal("NitrousOxide", "NitrousOxide");
    addReal("Novec649", "Novec649");
    addReal("OrthoHydrogen", "OrthoHydrogen");
    addReal("ParaHydrogen", "ParaHydrogen");
    addReal("Propylene", "Propylene");
    addReal("R11", "R11");
    addReal("R12", "R12");
    addReal("R13", "R13");
    addReal("R14", "R14");
    addReal("R23", "R23");
    addReal("R32", "R32");
    addReal("R41", "R41");
    addReal("R113", "R113");
    addReal("R114", "R114");
    addReal("R115", "R115");
    addReal("R116", "R116");
    addReal("R123", "R123");
    addReal("R124", "R124");
    addReal("R125", "R125");
    addReal("R134a", "R134a");
    addReal("R141b", "R141b");
    addReal("R142b", "R142b");
    addReal("R143a", "R143a");
    addReal("R152a", "R152a");
    addReal("R161", "R161");
    addReal("R218", "R218");
    addReal("R227ea", "R227ea");
    addReal("R236ea", "R236ea");
    addReal("R236fa", "R236fa");
    addReal("R245ca", "R245ca");
    addReal("R245fa", "R245fa");
    addReal("R365mfc", "R365mfc");
    addReal("R1233zd(E)", "R1233zd(E)");
    addReal("R1234yf", "R1234yf");
    addReal("R1234ze(E)", "R1234ze(E)");
    addReal("R1234ze(Z)", "R1234ze(Z)");
    addReal("R1243zf", "R1243zf");
    addReal("RC318", "RC318");
    addReal("SES36", "SES36");
    addReal("SulfurDioxide", "SulfurDioxide");
    addReal("SulfurHexafluoride", "SulfurHexafluoride");
    addReal("Xenon", "Xenon");
    
    // --- Unsupported Real Fluids (EES has them, CoolProp might not or they are special) ---
    auto addUnsupported = [](const std::string& name, const std::string& reason) {
        auto fluid = std::make_shared<UnsupportedFluid>(name, reason);
        registry_[name] = fluid;
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower != name) registry_[lower] = fluid;
    };
    
    addUnsupported("Acetylene", "Not in CoolProp standard fluids");
    addUnsupported("Chlorine", "Not in CoolProp standard fluids");
    addUnsupported("Chloroethene", "Not in CoolProp standard fluids");
    addUnsupported("EthyleneOxide", "Not in CoolProp standard fluids");
    addUnsupported("Helium3", "Not in CoolProp standard fluids");
    addUnsupported("Potassium", "Liquid metal, not in CoolProp");
    addUnsupported("Sodium", "Liquid metal, not in CoolProp");
    addUnsupported("IsoOctane", "Isooctane is not in CoolProp standard fluids");
    addUnsupported("IsoHexane", "Isohexane is not in CoolProp standard fluids");
    
    // --- Ideal Gases ---
    auto addIdeal = [](const std::string& name, const std::string& cpName) {
        auto fluid = std::make_shared<IdealGasFluid>(name, cpName);
        registry_[name] = fluid;
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower != name) registry_[lower] = fluid;
    };
    
    addIdeal("Air", "Air");
    addIdeal("N2", "Nitrogen");
    addIdeal("O2", "Oxygen");
    addIdeal("H2", "Hydrogen");
    addIdeal("He", "Helium");
    addIdeal("Ar", "Argon");
    addIdeal("CO2", "CarbonDioxide");
    addIdeal("CO", "CarbonMonoxide");
    addIdeal("CH4", "Methane");
    addIdeal("C2H6", "Ethane");
    addIdeal("C3H8", "Propane");
    
    // --- Real Fluid Air_ha (EES real gas) ---
    addReal("Air_ha", "Air");

    // --- Humid Air (EES humid air) ---
    registry_["airh2o"] = std::make_shared<HumidAirFluid>();
    
    // --- Incompressible Fluids (Examples) ---
    auto addIncomp = [](const std::string& name) {
        auto fluid = std::make_shared<IncompressibleFluid>(name);
        registry_[name] = fluid;
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower != name) registry_[lower] = fluid;
    };
    addIncomp("Aluminum");
    addIncomp("Copper");
    
    // --- Mixtures (Examples) ---
    auto addMix = [](const std::string& name) {
        auto fluid = std::make_shared<MixtureFluid>(name);
        registry_[name] = fluid;
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower != name) registry_[lower] = fluid;
    };
    addMix("NH3H2O");
    addMix("LiBrH2O");
    
    initialized_ = true;
}

std::shared_ptr<Fluid> FluidRegistry::getFluid(const std::string& name) {
    if (!initialized_) initialize();
    
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    
    auto it = registry_.find(lower);
    if (it != registry_.end()) {
        return it->second;
    }
    
    // If not found, we could try to see if it's a valid CoolProp fluid directly
    // but for now, we only support the ones in the registry.
    return nullptr;
}

void FluidRegistry::addFluid(std::shared_ptr<Fluid> fluid) {
    if (!initialized_) initialize();
    registry_[fluid->getName()] = fluid;
    std::string lower = fluid->getName();
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    registry_[lower] = fluid;
}

std::vector<std::shared_ptr<Fluid>> FluidRegistry::getAllFluids() {
    if (!initialized_) initialize();
    std::vector<std::shared_ptr<Fluid>> fluids;
    std::set<std::shared_ptr<Fluid>> uniqueFluids;
    for (const auto& [name, fluid] : registry_) {
        uniqueFluids.insert(fluid);
    }
    for (const auto& fluid : uniqueFluids) {
        fluids.push_back(fluid);
    }
    return fluids;
}

} // namespace coolsolve
