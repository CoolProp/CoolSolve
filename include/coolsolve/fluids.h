#pragma once

#include <string>
#include <sstream>
#include <iomanip>
#include <locale>
#include <vector>
#include <memory>
#include <map>
#include <mutex>
#include <optional>
#include <shared_mutex>
#include "coolsolve/units.h"

namespace coolsolve {

enum class FluidType {
    Real,
    IdealGas,
    Incompressible,
    HumidAir,
    Mixture,
    AmmoniaWater,   ///< the EES fluid NH3H2O: the IAPWS G4-01 formulation (nh3h2o.h), not CoolProp
    Unknown
};

struct ReferenceState {
    double h_offset = 0.0; // J/kg
    double s_offset = 0.0; // J/kg-K
};

class Fluid {
public:
    virtual ~Fluid() = default;
    
    virtual std::string getName() const = 0;
    virtual std::string getCoolPropName() const = 0;
    virtual FluidType getType() const = 0;
    virtual int getMinInputs() const = 0;
    
    // Get reference state offsets (EES - CoolProp) for the given CoolProp backend
    virtual ReferenceState getReferenceState(const std::string& backend = "HEOS") const {
        (void)backend;
        return {0.0, 0.0};
    }
};

class RealFluid : public Fluid {
public:
    RealFluid(const std::string& name, const std::string& cpName) 
        : name_(name), cpName_(cpName) {}
        
    std::string getName() const override { return name_; }
    std::string getCoolPropName() const override { return cpName_; }
    FluidType getType() const override { return FluidType::Real; }
    int getMinInputs() const override { return 2; }

private:
    std::string name_;
    std::string cpName_;
};

/**
 * @brief Ideal-gas species of EES (N2, CO2, H2O, ...), evaluated from the ideal-gas part (alpha0) of the CoolProp
 * equation of state of the substance (src/ideal_gas.cpp): pv = RT exactly, h and u functions of T only.
 *
 * Optionally carries the standard enthalpy of formation and the standard (third-law) entropy, which put the
 * enthalpy and the entropy of the species on the EES/JANAF basis: h(25 °C) equals the enthalpy of formation
 * (elements: 0), s(25 °C, 1 bar) the absolute entropy of the tables, so that heats and entropies of reaction,
 * Gibbs energies and exergies of combustion can be built from species properties (see getReferenceState()).
 */
class IdealGasFluid : public Fluid {
public:
    /**
     * @param formationEnthalpy  standard enthalpy of formation at 25 °C [J/kmol];
     *                           unset: CoolProp's own enthalpy reference is kept
     * @param standardEntropy    standard molar entropy at 298.15 K and 1 bar, ideal gas [J/kmol-K];
     *                           unset: CoolProp's own entropy reference is kept
     */
    IdealGasFluid(const std::string& name, const std::string& cpName,
                  std::optional<double> formationEnthalpy = std::nullopt,
                  std::optional<double> standardEntropy = std::nullopt)
        : name_(name), cpName_(cpName),
          formationEnthalpy_(formationEnthalpy), standardEntropy_(standardEntropy) {}

    /// Standard enthalpy of formation at 25 °C [J/kmol], if the species has one.
    std::optional<double> getFormationEnthalpy() const { return formationEnthalpy_; }

    /// Standard molar entropy at 298.15 K and 1 bar (ideal gas) [J/kmol-K], if the species has one.
    std::optional<double> getStandardEntropy() const { return standardEntropy_; }

    /**
     * @brief Reference-state offsets (EES - CoolProp) of the enthalpy and of the entropy.
     *
     * Enthalpy: `dHf/M - h0(25 °C)`, so that h(25 °C) = enthalpy of formation. Entropy: `s°/M - s0(25 °C, 1 bar)`.
     * h0 and s0 are the ideal-gas part of the CoolProp equation of state of the substance: the offsets are pure
     * reference constants, without any real-gas residual. They are computed once (thread-safe) and added to every
     * enthalpy / internal energy / entropy returned for this fluid; enthalpy, internal-energy and entropy inputs
     * are corrected the other way round. A species without a formation enthalpy (standard entropy) keeps CoolProp's
     * enthalpy (entropy) reference: no offset, no cost. The CoolProp backend does not matter (ideal gases always
     * use the HEOS equation of state): the argument is that of the base class.
     */
    ReferenceState getReferenceState(const std::string& backend = "HEOS") const override;

    /// The real fluid of the registry that has the same CoolProp equation of state (`Nitrogen` for `N2`, `Water` for
    /// `H2O`, `Air_ha` for `Air`), or an empty string when there is none (`CO`).
    std::string getRealFluidName() const;

    std::string getName() const override { return name_; }
    std::string getCoolPropName() const override { return cpName_; }
    FluidType getType() const override { return FluidType::IdealGas; }
    int getMinInputs() const override { return 1; }

private:
    std::string name_;
    std::string cpName_;
    std::optional<double> formationEnthalpy_;
    std::optional<double> standardEntropy_;
    mutable std::shared_mutex offsetMutex_;
    mutable std::optional<ReferenceState> reference_;
};

class HumidAirFluid : public Fluid {
public:
    std::string getName() const override { return "Air_ha"; }
    std::string getCoolPropName() const override { return "HumidAir"; }
    FluidType getType() const override { return FluidType::HumidAir; }
    int getMinInputs() const override { return 3; } // T, P, and one of (W, R, D, ...)
};

class IncompressibleFluid : public Fluid {
public:
    // `isSolution` indicates whether the fluid requires a mass-fraction
    // concentration (e.g. MEG, MPG, brines, LiBr solution, sea water).  A pure incompressible
    // does not need one; none is registered at present (the solids Aluminum and Copper are
    // unsupported fluids: CoolProp has none).
    // `coolPropKey` optionally overrides the CoolProp key used for PropsSI.
    IncompressibleFluid(const std::string& name,
                        bool isSolution = false,
                        const std::string& coolPropKey = "")
        : name_(name),
          coolPropKey_(coolPropKey.empty() ? name : coolPropKey),
          isSolution_(isSolution) {}
    std::string getName() const override { return name_; }
    std::string getCoolPropName() const override { return "INCOMP::" + coolPropKey_; }
    std::string getCoolPropKey() const { return coolPropKey_; }
    bool isSolution() const { return isSolution_; }
    FluidType getType() const override { return FluidType::Incompressible; }
    int getMinInputs() const override { return 2; }

    /// Pressure [Pa] of a property call whose only state input is the temperature (EES: `cp('EG',T=T,C=30)`).
    /// CoolProp's density, specific heat, viscosity and conductivity of an incompressible do not depend on the
    /// pressure; its enthalpy, entropy and internal energy do (a p/rho term), so they are those at this pressure.
    static constexpr double kInternalPressureSI = 101325.0;
    
    /// Build a CoolProp fluid name for a given mass fraction (0..1), the form that PropsSI reads
    /// (`INCOMP::MEG[0.3]`). Only meaningful when isSolution() is true.
    /// The fraction is written with every digit it has: the shortest text of 15 to 17 significant digits that CoolProp
    /// reads back as exactly `massFraction` (the default stream precision, 6 digits, rounded it: CS-BUG-INCOMP-CONC-AD).
    /// The evaluator asks CoolProp's cached state for the base name getCoolPropName() and sets the fraction on it;
    /// this name is for PropsSI and for the messages.
    std::string getCoolPropNameWithFraction(double massFraction) const {
        std::string text;
        for (int digits = 15; digits <= 17; ++digits) {
            std::ostringstream number;
            number.imbue(std::locale::classic());
            number << std::setprecision(digits) << massFraction;
            text = number.str();
            if (digits == 17) break;   // 17 digits always read back exactly
            std::istringstream back(text);
            back.imbue(std::locale::classic());
            double value = 0.0;
            if ((back >> value) && value == massFraction) break;
        }
        return "INCOMP::" + coolPropKey_ + "[" + text + "]";
    }

private:
    std::string name_;
    std::string coolPropKey_;
    bool isSolution_;
};

/**
 * @brief A CoolProp mixture written `A[x]&B[y]` (mole fractions in brackets, two or more
 * components, optional `HEOS::` prefix), e.g. `'R245fa[0.5]&R134a[0.5]'`.
 *
 * Not in the fixed registry: created and checked by FluidRegistry::getFluid() from the string.
 * The name (also the CoolProp name) is the canonical form `A[x]&B[y]` without prefix, spaces
 * and EES aliases; the properties come from CoolProp's mixture model. Components with a zero
 * mole fraction are dropped, as CoolProp does: a pure composition is a RealFluid, not a
 * MixtureFluid, so a MixtureFluid always has two or more components.
 */
class MixtureFluid : public Fluid {
public:
    MixtureFluid(std::string name, std::vector<std::string> components, std::vector<double> moleFractions)
        : name_(std::move(name)), components_(std::move(components)),
          moleFractions_(std::move(moleFractions)) {}
    std::string getName() const override { return name_; }
    std::string getCoolPropName() const override { return name_; }
    FluidType getType() const override { return FluidType::Mixture; }
    int getMinInputs() const override { return 2; }

    /// CoolProp names of the components, in the order of the string.
    const std::vector<std::string>& getComponents() const { return components_; }
    /// Mole fractions (sum 1), in the order of the components.
    const std::vector<double>& getMoleFractions() const { return moleFractions_; }
    /// Component names joined by '&', the fluid name of CoolProp's AbstractState factory.
    std::string getComponentList() const;

private:
    std::string name_;
    std::vector<std::string> components_;
    std::vector<double> moleFractions_;
};

/**
 * @brief The ammonia-water mixture of EES (`NH3H2O`), evaluated by the module `nh3h2o` (IAPWS Guideline G4-01), not by
 * CoolProp. The state is fixed by three of T, P, X (mass fraction of ammonia), Q (vapour mass fraction), H, S, U, V.
 */
class AmmoniaWaterFluid : public Fluid {
public:
    std::string getName() const override { return "NH3H2O"; }
    std::string getCoolPropName() const override { return ""; }
    FluidType getType() const override { return FluidType::AmmoniaWater; }
    int getMinInputs() const override { return 3; }
};

class UnsupportedFluid : public Fluid {
public:
    UnsupportedFluid(const std::string& name, const std::string& reason) 
        : name_(name), reason_(reason) {}
    std::string getName() const override { return name_; }
    std::string getCoolPropName() const override { return ""; }
    FluidType getType() const override { return FluidType::Unknown; }
    int getMinInputs() const override { return 0; }
    std::string getReason() const { return reason_; }

private:
    std::string name_;
    std::string reason_;
};

/**
 * @brief Standard atomic weight of a chemical element, the argument of MOLARMASS(C), MOLARMASS(H), ...
 *
 * The symbol is read without regard to case, as every EES name. These are not fluids: a name that is a
 * registered fluid (N2, He, CO, ...) is looked up in FluidRegistry first, and only the other names are tried here.
 *
 * @param symbol  element symbol such as "C", "H", "Fe"
 * @return the atomic weight in kg/kmol (IUPAC 2005, Wieser, Pure Appl. Chem. 78 (2006) 2051), or nothing when
 *         the symbol is not an element with a standard atomic weight (also the elements without stable nuclide)
 */
std::optional<double> elementAtomicWeight(const std::string& symbol);

class FluidRegistry {
public:
    /**
     * @brief Look a fluid up by name (case-insensitive); a name containing '&' is a CoolProp mixture.
     *
     * Mixtures are parsed, checked against CoolProp (components, binary parameters) and cached by
     * string, once per process; pure fluids pay nothing for it.
     *
     * @param error  when given, receives the reason why a mixture string was rejected (fractions
     *               not summing to 1, unknown component, ...); left empty for other unknown names
     * @return the fluid, or nullptr when the name is unknown or the mixture is invalid
     */
    static std::shared_ptr<Fluid> getFluid(const std::string& name, std::string* error = nullptr);
    static void addFluid(std::shared_ptr<Fluid> fluid);
    static std::vector<std::shared_ptr<Fluid>> getAllFluids();
    
private:
    static std::shared_ptr<Fluid> getMixture(const std::string& name, std::string* error);

    static std::map<std::string, std::shared_ptr<Fluid>> registry_;
    static bool initialized_;
    static std::map<std::string, std::shared_ptr<Fluid>> mixtures_;   ///< by string as written and by canonical name
    static std::shared_mutex mixturesMutex_;
    static void initialize();
};

} // namespace coolsolve
