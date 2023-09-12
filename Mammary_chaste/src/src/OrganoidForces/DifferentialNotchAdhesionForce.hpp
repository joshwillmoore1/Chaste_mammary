#ifndef DIFFERENTIALNOTCHADHESIONFORCE_HPP_
#define DIFFERENTIALNOTCHADHESIONFORCE_HPP_

#include "GeneralisedLinearSpringForce.hpp"

/**
 * A class for a simple two-body differential adhesion force law between
 * labelled and unlabelled cells (as defined by the CellLabel cell
 * property).
 *
 * Designed for use in node and mesh-based simulations.
 *
 * \todo #2266 - throw exceptions if using other cell population objects?
 * \todo #2266 - override CalculateForceBetweenNodes() to use a default rest length of 1.0 for all springs?
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class DifferentialNotchAdhesionForce : public GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>
{
private :

    /**
     * A scalar value which determines the maximum spring constant
     * multiplier.
     *
     * Defaults to 1.0 in the constructor.
     *
     * returned by the method VariableSpringConstantMultiplicationFactor()
     * in the parent class GeneralisedLinearSpringForce.
     */
    double mMaximumSpringConstantMultiplier;

    /**
     * A scalar determining the spring constajt multiplier decay used
     * in the overridden method VariableSpringConstantMultiplicationFactor().
     *
     * Defaults to 1.0 in the constructor.
     */
    double mAdhesionDecayConstant;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mMaximumSpringConstantMultiplier;
        archive & mAdhesionDecayConstant;
    }

public :

    /**
     * Constructor.
     */
    DifferentialNotchAdhesionForce();

    /**
     * Overridden VariableSpringConstantMultiplicationFactor() method.
     *
     * This method takes account of the distinct spring constants present
     * for homotypic (labelled-labelled and unlabelled-unlabelled) and
     * heterotypic (labelled-unlabelled) interactions between neighbouring
     * cells.
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @param isCloserThanRestLength whether the neighbouring nodes lie closer than the rest length of their connecting spring
     *
     * @return the multiplication factor.
     */
    double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                      unsigned nodeBGlobalIndex,
                                                      AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                      bool isCloserThanRestLength);

    /**
     * @return #mHomotypicLabelledSpringConstantMultiplier.
     */
    double GetMaximumSpringConstantMultiplier();

    /**
     * Set mHomotypicLabelledSpringConstantMultiplier.
     *
     * @param TypeSpringMaximumSpringConstantMultiplier the new value of mHomotypicLabelledSpringConstantMultiplier
     */
    void SetMaximumSpringConstantMultiplier(double TypeSpringMaximumSpringConstantMultiplier);

    /**
     * @return #mHeterotypicSpringConstantMultiplier.
     */
    double GetAdhesionDecayConstant();

    /**
     * Set mHeterotypicSpringConstantMultiplier.
     *
     * @param AdhesionDecayConstant the new value of mHeterotypicSpringConstantMultiplier
     */
    void SetAdhesionDecayConstant(double AdhesionDecayConstant);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(DifferentialNotchAdhesionForce)

#endif /*DIFFERENTIALNOTCHADHESIONFORCE_HPP_*/
