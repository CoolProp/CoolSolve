#include "coolsolve/fluids.h"
#include "AbstractState.h"
#include "DataStructures.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <stdexcept>
#include <set>

namespace coolsolve {

std::map<std::string, std::shared_ptr<Fluid>> FluidRegistry::registry_;
bool FluidRegistry::initialized_ = false;
std::map<std::string, std::shared_ptr<Fluid>> FluidRegistry::mixtures_;
std::shared_mutex FluidRegistry::mixturesMutex_;

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
    // Standard enthalpies of formation of the gas at 298.15 K [J/kmol] (EES/JANAF convention: the elements
    // in their reference state, N2 O2 H2 He Ar, are 0) and standard molar entropies of the ideal gas at
    // 298.15 K and 1 bar [J/kmol-K] (third law). Both come from one source, the NASA-Glenn thermodynamic
    // data (McBride, Zehe, Gordon, NASA/TP-2002-211556; github.com/nasa/cea, data/thermo.inp): the heat of
    // formation is the field of the species header (J/mol, x 1000), the entropy the 9-coefficient polynomial
    // evaluated at 298.15 K with R = 8.314462618 J/mol-K (the polynomials' own H(298.15 K) differs from the
    // header by at most 2.2 J/mol, the R of the fit). It is the data family of the EES ideal gases, which
    // use it with R = 8.314 (EES help: CO2 -393 486, H2O -241 811, CH4 -74 595 kJ/kmol). The rounded JANAF
    // values of Cengel & Boles (table A-26) used before differed by 5 to 830 kJ/kmol (CO2 -393 520, CO
    // -110 530, CH4 -74 850, C2H6 -84 680, C3H8 -103 850, H2O -241 820); NIST-JANAF entropies differ from the
    // NASA ones by at most 120 J/kmol-K (CH4). The lower heating value of CO is 282 974.8 kJ/kmol.
    // `Air` has neither: EES's Air is not a JANAF species, it keeps CoolProp's enthalpy and entropy references.
    auto addIdeal = [](const std::string& name, const std::string& cpName,
                       std::optional<double> formationEnthalpy = std::nullopt,
                       std::optional<double> standardEntropy = std::nullopt) {
        auto fluid = std::make_shared<IdealGasFluid>(name, cpName, 101325.0, formationEnthalpy, standardEntropy);
        registry_[name] = fluid;
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower != name) registry_[lower] = fluid;
    };
    
    addIdeal("Air", "Air");
    addIdeal("N2", "Nitrogen", 0.0, 191608.6);
    addIdeal("O2", "Oxygen", 0.0, 205148.3);
    addIdeal("H2", "Hydrogen", 0.0, 130680.3);
    addIdeal("He", "Helium", 0.0, 126153.0);
    addIdeal("Ar", "Argon", 0.0, 154845.8);
    addIdeal("CO2", "CarbonDioxide", -393510e3, 213786.2);
    addIdeal("CO", "CarbonMonoxide", -110535.196e3, 197658.7);
    addIdeal("CH4", "Methane", -74600e3, 186370.0);
    addIdeal("C2H6", "Ethane", -83851.544e3, 229219.8);
    addIdeal("C3H8", "Propane", -104680e3, 270313.3);
    
    // H2O as ideal gas (water vapor for combustion/chemistry): formation enthalpy and entropy of the vapor.
    // Use a low dummy pressure (100 Pa) so CoolProp stays in the gas phase
    // even at low temperatures (saturation T at 100 Pa ≈ −23 °C).
    {
        auto fluid = std::make_shared<IdealGasFluid>("H2O", "Water", 100.0, -241826e3, 188828.0);
        registry_["h2o"] = fluid;
    }
    
    // --- Real Fluid Air_ha (EES real gas) ---
    addReal("Air_ha", "Air");

    // --- Humid Air (EES humid air) ---
    registry_["airh2o"] = std::make_shared<HumidAirFluid>();
    
    // --- Incompressible Fluids ---
    // Pure incompressibles (no concentration): solids, pure liquids
    auto addIncomp = [](const std::string& name, const std::string& cpKey = "") {
        auto fluid = std::make_shared<IncompressibleFluid>(name, /*isSolution=*/false, cpKey);
        registry_[name] = fluid;
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower != name) registry_[lower] = fluid;
    };
    // Solutions: glycols, brines, secondary refrigerants (need concentration)
    auto addIncompSolution = [](const std::string& name, const std::string& cpKey) {
        auto fluid = std::make_shared<IncompressibleFluid>(name, /*isSolution=*/true, cpKey);
        registry_[name] = fluid;
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower != name) registry_[lower] = fluid;
    };
    
    addIncomp("Aluminum");
    addIncomp("Copper");
    
    // Glycols
    addIncompSolution("EthyleneGlycol", "MEG");
    addIncompSolution("EG", "MEG");
    addIncompSolution("MEG", "MEG");
    addIncompSolution("PropyleneGlycol", "MPG");
    addIncompSolution("PG", "MPG");
    addIncompSolution("MPG", "MPG");
    
    // Brines
    addIncompSolution("CaCl2", "MCA");
    addIncompSolution("MCA", "MCA");
    addIncompSolution("NaCl", "MNA");
    addIncompSolution("MNA", "MNA");
    addIncompSolution("LiCl", "MLI");
    addIncompSolution("MLI", "MLI");
    
    // Secondary refrigerants / solvents (solutions in water)
    addIncompSolution("Methanol", "MMA");
    addIncompSolution("MMA", "MMA");
    addIncompSolution("EthanolSolution", "MEA");
    addIncompSolution("MEA", "MEA");
    addIncompSolution("Glycerol", "MGL");
    addIncompSolution("MGL", "MGL");
    
    // Pure incompressible
    addIncomp("LiBr", "LiBr");   // pure LiBr correlations
    addIncomp("SeaWater", "MITSW");
    
    // --- Mixtures (unsupported in CoolProp with EES parameterization) ---
    addUnsupported("NH3H2O",
        "Ammonia-water mixture with EES correlation is not available in CoolProp. "
        "Use REFPROP via \"REFPROP::Ammonia&Water\" backend or provide an external correlation.");
    addUnsupported("LiBrH2O",
        "Lithium-bromide-water solution with EES correlation is not available in CoolProp. "
        "Use an external correlation or the incompressible 'LiBr' fluid for pure LiBr data.");
    
    initialized_ = true;
}

std::shared_ptr<Fluid> FluidRegistry::getFluid(const std::string& name, std::string* error) {
    if (!initialized_) initialize();
    
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    
    auto it = registry_.find(lower);
    if (it != registry_.end()) {
        return it->second;
    }
    
    // Not in the registry: a CoolProp mixture `A[x]&B[y]` is built from the string. Any other
    // CoolProp fluid is only reachable through the registry.
    if (name.find('&') != std::string::npos) return getMixture(name, error);
    return nullptr;
}

namespace {

std::string trimmed(const std::string& s) {
    const auto first = s.find_first_not_of(" \t");
    if (first == std::string::npos) return {};
    return s.substr(first, s.find_last_not_of(" \t") - first + 1);
}

bool equalsNoCase(const std::string& a, const std::string& b) {
    return a.size() == b.size() &&
           std::equal(a.begin(), a.end(), b.begin(), [](char x, char y) { return std::tolower(x) == std::tolower(y); });
}

std::string numberText(double v, int digits) {
    std::ostringstream oss;
    oss.imbue(std::locale::classic());
    oss << std::setprecision(digits) << v;
    return oss.str();
}

/// Mole fraction written in brackets: a plain number, nothing else.
bool parseMoleFraction(const std::string& text, double& value) {
    std::istringstream iss(trimmed(text));
    iss.imbue(std::locale::classic());
    iss >> value;
    return !iss.fail() && iss.eof() && std::isfinite(value);
}

}  // namespace

std::string MixtureFluid::getComponentList() const {
    std::string list;
    for (const auto& c : components_) list += (list.empty() ? "" : "&") + c;
    return list;
}

std::shared_ptr<Fluid> FluidRegistry::getMixture(const std::string& name, std::string* error) {
    {
        std::shared_lock<std::shared_mutex> lock(mixturesMutex_);
        auto it = mixtures_.find(name);
        if (it != mixtures_.end()) return it->second;
    }
    auto fail = [&](const std::string& reason) -> std::shared_ptr<Fluid> {
        if (error) *error = "Invalid mixture '" + name + "': " + reason;
        return nullptr;
    };

    // Optional backend prefix: only the default backend, HEOS, is accepted
    std::string text = trimmed(name);
    const auto colons = text.find("::");
    if (colons != std::string::npos) {
        const std::string prefix = trimmed(text.substr(0, colons));
        if (!equalsNoCase(prefix, "HEOS"))
            return fail("the backend prefix '" + prefix + "::' is not supported (use 'HEOS::' or none)");
        text = trimmed(text.substr(colons + 2));
    }

    // Components with their mole fractions: NAME[x]
    std::vector<std::string> written, components;
    std::vector<double> fractions;
    for (size_t start = 0;;) {
        const auto amp = text.find('&', start);
        const std::string part = trimmed(text.substr(start, amp == std::string::npos ? amp : amp - start));
        if (part.empty()) return fail("empty component (stray '&')");
        const auto open = part.find('[');
        if (open == std::string::npos)
            return fail("component '" + part + "' has no mole fraction (write 'A[x]&B[y]', e.g. 'R245fa[0.5]&R134a[0.5]')");
        const std::string compName = trimmed(part.substr(0, open));
        if (compName.empty()) return fail("component name missing before '" + part.substr(open) + "'");
        if (part.back() != ']')
            return fail("the mole fraction of component '" + compName + "' is not closed by ']'");
        double x = 0.0;
        const std::string xText = part.substr(open + 1, part.size() - open - 2);
        if (!parseMoleFraction(xText, x))
            return fail("the mole fraction '" + xText + "' of component '" + compName + "' is not a number");
        if (x < 0.0 || x > 1.0)
            return fail("the mole fraction " + numberText(x, 10) + " of component '" + compName +
                        "' must be between 0 and 1");
        written.push_back(compName);
        fractions.push_back(x);
        if (amp == std::string::npos) break;
        start = amp + 1;
    }
    double sum = 0.0;
    for (double x : fractions) sum += x;
    if (std::abs(sum - 1.0) > 1e-6)
        return fail("the mole fractions sum to " + numberText(sum, 10) + ", not 1");
    for (double& x : fractions) x /= sum;

    // Components: CoolProp names (the EES names of the registry, case-insensitive, are accepted too)
    for (const auto& compName : written) {
        std::string cpName = compName;
        auto known = getFluid(compName);
        if (known && known->getType() == FluidType::Real) cpName = known->getCoolPropName();
        try {
            std::unique_ptr<CoolProp::AbstractState> probe(CoolProp::AbstractState::factory("HEOS", cpName));
        } catch (const std::exception&) {
            return fail("unknown component '" + compName + "', not a CoolProp fluid name");
        }
        if (std::find(components.begin(), components.end(), cpName) != components.end())
            return fail("component '" + compName + "' appears twice");
        components.push_back(cpName);
    }

    // Like CoolProp, drop the components with a zero mole fraction: a pure composition (MM = 1 in
    // 'A[1]&B[0]') is the pure fluid, with exactly its properties and its functions
    for (size_t i = components.size(); i-- > 0;) {
        if (fractions[i] > 10 * std::numeric_limits<double>::epsilon()) continue;
        components.erase(components.begin() + i);
        fractions.erase(fractions.begin() + i);
    }
    if (components.size() == 1) {
        auto pure = std::make_shared<RealFluid>(components[0], components[0]);
        std::unique_lock<std::shared_mutex> lock(mixturesMutex_);
        return mixtures_.emplace(name, pure).first->second;
    }

    std::string canonical;
    for (size_t i = 0; i < components.size(); ++i)
        canonical += (i ? "&" : "") + components[i] + "[" + numberText(fractions[i], 15) + "]";
    auto mixture = std::make_shared<MixtureFluid>(canonical, components, fractions);
    try {   // binary interaction parameters of every pair, among others
        std::unique_ptr<CoolProp::AbstractState> probe(
            CoolProp::AbstractState::factory("HEOS", mixture->getComponentList()));
        probe->set_mole_fractions(fractions);
    } catch (const std::exception& e) {
        return fail(std::string("CoolProp cannot build it: ") + e.what());
    }

    std::unique_lock<std::shared_mutex> lock(mixturesMutex_);
    mixtures_.emplace(canonical, mixture);
    return mixtures_.emplace(name, mixtures_.at(canonical)).first->second;
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

ReferenceState IdealGasFluid::getReferenceState(const std::string& backend) const {
    if (!formationEnthalpy_ && !standardEntropy_) return {0.0, 0.0};
    {
        std::shared_lock<std::shared_mutex> lock(offsetMutex_);
        auto it = referenceByBackend_.find(backend);
        if (it != referenceByBackend_.end()) return it->second;
    }
    // The formation enthalpy and the standard entropy apply at 25 °C (1 atm, resp. 1 bar).
    constexpr double T_REFERENCE_K = 298.15;
    constexpr double P_STANDARD_PA = 1.0e5;          // standard-state pressure of the entropy tables
    constexpr double P_IDEAL_LIMIT_PA = 100.0;       // CoolProp's equation of state is an ideal gas here
    ReferenceState reference;
    try {
        std::unique_ptr<CoolProp::AbstractState> state(
            CoolProp::AbstractState::factory(backend, cpName_));
        if (formationEnthalpy_) {
            // temperature-only calls use the dummy pressure
            state->update(CoolProp::PT_INPUTS, dummyPressure_, T_REFERENCE_K);
            const double formationPerKg = *formationEnthalpy_ / (1000.0 * state->molar_mass());
            reference.h_offset = formationPerKg - state->hmass();
        }
        if (standardEntropy_) {
            // Third-law entropy at (T, p) of the ideal gas: s0 - R ln(p/1 bar). Compared with CoolProp at a
            // pressure low enough for the residual entropy to vanish, so that the offset is a pure reference
            // constant; at any pressure, CoolProp then adds its real-gas residual to the ideal-gas entropy.
            // A backend that refuses 100 Pa (tabular tables) is called at the dummy pressure instead.
            double p = P_IDEAL_LIMIT_PA;
            try {
                state->update(CoolProp::PT_INPUTS, p, T_REFERENCE_K);
            } catch (const std::exception&) {
                p = dummyPressure_;
                state->update(CoolProp::PT_INPUTS, p, T_REFERENCE_K);
            }
            const double molarMass = state->molar_mass();                              // kg/mol
            const double s0PerMolK = *standardEntropy_ / 1000.0;                       // J/(mol K)
            const double idealPerMolK = s0PerMolK - state->gas_constant() * std::log(p / P_STANDARD_PA);
            reference.s_offset = idealPerMolK / molarMass - state->smass();
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("Cannot compute the reference state (enthalpy of formation, standard entropy) of '" +
                                 name_ + "' with backend '" + backend + "': " + e.what());
    }
    std::unique_lock<std::shared_mutex> lock(offsetMutex_);
    referenceByBackend_.emplace(backend, reference);
    return reference;
}

} // namespace coolsolve
