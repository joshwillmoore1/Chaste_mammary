#ifndef ORGANOIDSPRINGFORCE_HPP_
#define ORGANOIDSPRINGFORCE_HPP_

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
class OrganoidSpringForce : public GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>
{
private :

    /**
     * A scalar determining the relative spring constant for homotypic
     * interactions between neighbouring labelled cells, used in the
     * overridden method VariableSpringConstantMultiplicationFactor().
     *
     * Defaults to 1.0 in the constructor.
     *
     * Note that for homotypic interactions between neighbouring
     * unlabelled cells, we use the multiplier value 1.0 that is
     * returned by the method VariableSpringConstantMultiplicationFactor()
     * in the parent class GeneralisedLinearSpringForce.
     */
    double mHomotypicTypeSpringConstantMultiplier;

    /**
     * A scalar determining the relative spring constant for heterotypic
     * (labelled-unlabelled) interactions between neighbouring cells, used
     * in the overridden method VariableSpringConstantMultiplicationFactor().
     *
     * Defaults to 1.0 in the constructor.
     */
    double mHeterotypicSpringConstantMultiplier;

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
        archive & mHomotypicTypeSpringConstantMultiplier;
        archive & mHeterotypicSpringConstantMultiplier;
    }

public :

    /**
     * Constructor.
     */
    OrganoidSpringForce();

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
    double GetHomotypicTypeSpringConstantMultiplier();

    /**
     * Set mHomotypicLabelledSpringConstantMultiplier.
     *
     * @param TypeSpringConstantMultiplier the new value of mHomotypicLabelledSpringConstantMultiplier
     */
    void SetHomotypicTypeSpringConstantMultiplier(double TypeSpringConstantMultiplier);

    /**
     * @return #mHeterotypicSpringConstantMultiplier.
     */
    double GetHeterotypicSpringConstantMultiplier();

    /**
     * Set mHeterotypicSpringConstantMultiplier.
     *
     * @param heterotypicSpringConstantMultiplier the new value of mHeterotypicSpringConstantMultiplier
     */
    void SetHeterotypicSpringConstantMultiplier(double heterotypicSpringConstantMultiplier);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OrganoidSpringForce)

#endif /*ORGANOIDSPRINGFORCE_HPP_*/
