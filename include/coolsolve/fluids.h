#pragma once

#include <string>
#include <vector>
#include <memory>
#include <map>
#include "coolsolve/units.h"

namespace coolsolve {

enum class FluidType {
    Real,
    IdealGas,
    Incompressible,
    HumidAir,
    Mixture,
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
    
    // Get reference state offsets (EES - CoolProp)
    virtual ReferenceState getReferenceState() const { return {0.0, 0.0}; }
    
    // Check if a property depends on pressure for this fluid
    virtual bool propertyDependsOnPressure(const std::string& prop) const {
        return true; 
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

class IdealGasFluid : public Fluid {
public:
    IdealGasFluid(const std::string& name, const std::string& cpName) 
        : name_(name), cpName_(cpName) {}
        
    std::string getName() const override { return name_; }
    std::string getCoolPropName() const override { return cpName_; }
    FluidType getType() const override { return FluidType::IdealGas; }
    int getMinInputs() const override { return 1; }
    
    bool propertyDependsOnPressure(const std::string& prop) const override {
        std::string p = prop;
        for (auto& c : p) c = std::tolower(c);
        if (p == "h" || p == "enthalpy" || p == "u" || p == "internalenergy" || p == "cp" || p == "cv") {
            return false;
        }
        return true;
    }

private:
    std::string name_;
    std::string cpName_;
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
    IncompressibleFluid(const std::string& name) : name_(name) {}
    std::string getName() const override { return name_; }
    std::string getCoolPropName() const override { return "INCOMP::" + name_; }
    FluidType getType() const override { return FluidType::Incompressible; }
    int getMinInputs() const override { return 2; }

private:
    std::string name_;
};

class MixtureFluid : public Fluid {
public:
    MixtureFluid(const std::string& name) : name_(name) {}
    std::string getName() const override { return name_; }
    std::string getCoolPropName() const override { return name_; }
    FluidType getType() const override { return FluidType::Mixture; }
    int getMinInputs() const override { return 2; }

private:
    std::string name_;
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

class FluidRegistry {
public:
    static std::shared_ptr<Fluid> getFluid(const std::string& name);
    static void addFluid(std::shared_ptr<Fluid> fluid);
    static std::vector<std::shared_ptr<Fluid>> getAllFluids();
    
private:
    static std::map<std::string, std::shared_ptr<Fluid>> registry_;
    static bool initialized_;
    static void initialize();
};

} // namespace coolsolve
