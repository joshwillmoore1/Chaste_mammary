#ifndef COLLIERFIXEDTRACKINGMODIFIER_HPP_
#define COLLIERFIXEDTRACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <complex> 
#include "AbstractCellBasedSimulationModifier.hpp"
#include "MyoEpiCellProliferativeType.hpp"
#include <AbstractTwoBodyInteractionForce.hpp>

/**
 * A modifier class in which the mean levels of Delta in neighbouring cells
 * are computed and stored in CellData. To be used in conjunction with Delta
 * Notch cell cycle models.
 */
template<unsigned DIM>
class CollierFixedTrackingModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    /** Needed for serialization. */
    double mSameCellW1;


    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    CollierFixedTrackingModifier();

    /**
     * Destructor.
     */
    virtual ~CollierFixedTrackingModifier();

     /**
55	     * Set the polarisation gradient bound constant for the cells.
56	     *
57	     * @param samecellw1 the movement parameter to use
58	     */
	    void SetSameCellW1(double samecellw1);
	
	    /**
62	     * Get the polarisation coefficient.
63	     *
64	     * @return mMovementParameter
65	     */
	    double GetSameCellW1();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to compute the mean level of Delta in each cell's neighbours and store these in the CellData.
     *
     * Note: If using a CaBasedCellPopulation, we assume a Moore neighbourhood and unit carrying capacity.
     * If a cell has no neighbours (such as an isolated cell in a CaBasedCellPopulation), we store the
     * value -1 in the CellData.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CollierFixedTrackingModifier)

#endif /*COLLIERFIXEDTRACKINGMODIFIER_HPP_*/
