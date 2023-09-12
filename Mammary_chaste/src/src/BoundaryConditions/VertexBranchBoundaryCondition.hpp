#ifndef VERTEXBRANCHBOUNDARYCONDITION_HPP_
#define VERTEXBRANCHBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"


#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * Update this
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class VertexBranchBoundaryCondition : public AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>
{
private:

    /**
     * Centre of the boundary.
     */
    c_vector<double, SPACE_DIM> mCentre;

    /**
     * Size of circular boundary.
     */
    double mBoundaryRadius;

    /**
     * Oscilating radius.
     */
    double mOscRadius;

    /**
     * Whether to jiggle the cells on the bottom surface, initialised to false
     * in the constructor.
     */
    bool mUseJiggledNodesOnPlane;

    /** Needed for serialization. */
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
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM> >(*this);
        //archive & mUseJiggledNodesOnPlane;
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param centre a point on the boundary plane
     * @param BoundaryRadius the outward-facing unit normal vector to the boundary plane
     * @param oscRadius the outward-facing unit normal vector to the boundary plane
     */
    VertexBranchBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
                           c_vector<double, SPACE_DIM> centre,
                           double BoundaryRadius, double oscRadius);


    /**
     * @return #mCentre.
     */
    const c_vector<double, SPACE_DIM>& rGetCentre() const;

    /**
     * @return #mBoundaryRadius.
     */
    double rGetBoundaryRadius() const;


    /**
     * @return #mOscRadius.
     */
    double rGetOscRadius() const;


    /**
     * Set method for mUseJiggledNodesOnPlane
     *
     * @param useJiggledNodesOnPlane whether to jiggle the nodes on the surface of the plane, can help stop overcrowding on plane.
     */
    void SetUseJiggledNodesOnPlane(bool useJiggledNodesOnPlane);

    /** @return #mUseJiggledNodesOnPlane. */
    bool GetUseJiggledNodesOnPlane();

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexBranchBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a PlaneBoundaryCondition.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const VertexBranchBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;

    // Archive c_vectors one component at a time
    
    c_vector<double, SPACE_DIM> point = t->rGetCentre();
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        ar << point[i];
    }
    
    double radius  = t ->rGetBoundaryRadius();
    ar << radius; 

}

/**
 * De-serialize constructor parameters and initialize a PlaneBoundaryCondition.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, VertexBranchBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population;
    ar >> p_cell_population;

    // Archive c_vectors one component at a time
    
    c_vector<double, SPACE_DIM> point;
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        ar >> point[i];
    }
     
   double radius = t->rGetBoundaryRadius();
   double oscRadius = t->rGetOscRadius();
    

    // Invoke inplace constructor to initialise instance
    ::new(t)VertexBranchBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(p_cell_population, point, radius, oscRadius);
}
}
} // namespace ...
#endif /*VERTEXBRANCHBOUNDARYCONDITION_HPP_*/
