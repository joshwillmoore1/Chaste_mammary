#ifndef MAMDELTANOTCHODESYSTEM_HPP_
#define MAMDELTANOTCHODESYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"

/**
 * Represents the notch model basal stem cell for the mammary organoid
 * Joshua W. Moore - 2020
 */
class MamDeltaNotchOdeSystem : public AbstractOdeSystem
{
private:

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

    ///\todo extract model parameters as member variables

public:

    /**
     * Default constructor.
     *
     * @param stateVariables optional initial conditions for state variables (only used in archiving)
     */
    MamDeltaNotchOdeSystem(std::vector<double> stateVariables=std::vector<double>());

    /**
     * Destructor.
     */
    ~MamDeltaNotchOdeSystem();

    /**
     * Compute the RHS of the Moore et al. system of ODEs.
     *
     * Returns a vector representing the RHS of the ODEs at each time step, y' = [y1' ... yn'].
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time used to evaluate the RHS.
     * @param rY value of the solution vector used to evaluate the RHS.
     * @param rDY filled in with the resulting derivatives (using  Moores et al. system of equations).
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MamDeltaNotchOdeSystem)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a MamDeltaNotchOdeSystem.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const MamDeltaNotchOdeSystem * t, const unsigned int file_version)
{
    const std::vector<double>& state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/**
 * De-serialize constructor parameters and initialise a MamDeltaNotchOdeSystem.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, MamDeltaNotchOdeSystem * t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new(t)MamDeltaNotchOdeSystem(state_variables);
}
}
} // namespace ...

#endif /*MAMDELTANOTCHODESYSTEM_HPP_*/