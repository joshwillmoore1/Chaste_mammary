#ifndef NORMALCENTREDIVISIONPLANERULE_HPP_
#define NORMALCENTREDIVISIONPLANERULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractVertexBasedDivisionRule.hpp"
#include "VertexBasedCellPopulation.hpp"

// Forward declaration prevents circular include chain
template<unsigned SPACE_DIM> class VertexBasedCellPopulation;
template<unsigned SPACE_DIM> class AbstractVertexBasedDivisionRule;

/**
 * A class to generate the short axis of a cell for vertex-based cell
 * populations, to be used in cell division. This is the default rule that
 * is used in most of the vertex-based simulations.
 *
 * The short axis is the eigenvector associated with the largest eigenvalue
 * of the moment of inertia of the cell's polygon.
 */
template <unsigned SPACE_DIM>
class NormalCentreDivisionPlaneRule : public AbstractVertexBasedDivisionRule<SPACE_DIM>
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
        archive & boost::serialization::base_object<AbstractVertexBasedDivisionRule<SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    NormalCentreDivisionPlaneRule()
    {
    }

    /**
     * Empty destructor.
     */
    virtual ~NormalCentreDivisionPlaneRule()
    {
    }

    /**
     * Overridden CalculateCellDivisionVector() method.
     *
     * Return the short axis of the existing cell, which will be used to
     * form the boundary between the daughter cells.
     *
     * @param pParentCell  The cell to divide
     * @param rCellPopulation  The vertex-based cell population
     * @return the division vector.
     */
    virtual c_vector<double, SPACE_DIM> CalculateCellDivisionVector(CellPtr pParentCell,
        VertexBasedCellPopulation<SPACE_DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NormalCentreDivisionPlaneRule)

#endif // NORMALCENTREDIVISIONPLANERULE_HPP_
