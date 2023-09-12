#include "VertexBranchBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexBranchBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::VertexBranchBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> *pCellPopulation,
                                                                                     c_vector<double, SPACE_DIM> centre,
                                                                                     double BoundaryRadius, double oscRadius)
    : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(pCellPopulation),
      mCentre(centre),
      mUseJiggledNodesOnPlane(false)
{
    assert(BoundaryRadius > 0.0);
    mBoundaryRadius = BoundaryRadius;

    assert(oscRadius > 0.0);
    mOscRadius = oscRadius;


}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM> &VertexBranchBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetCentre() const
{
    return mCentre;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexBranchBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetBoundaryRadius() const
{
    return mBoundaryRadius;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexBranchBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetOscRadius() const
{
    return mOscRadius;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexBranchBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::SetUseJiggledNodesOnPlane(bool useJiggledNodesOnPlane)
{
    mUseJiggledNodesOnPlane = useJiggledNodesOnPlane;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexBranchBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::GetUseJiggledNodesOnPlane()
{
    return mUseJiggledNodesOnPlane;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexBranchBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::ImposeBoundaryCondition(const std::map<Node<SPACE_DIM> *, c_vector<double, SPACE_DIM>> &rOldLocations)
{

    ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM> *>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("VertexBranchBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }

    assert((dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM> *>(this->mpCellPopulation)) || (SPACE_DIM == ELEMENT_DIM && (dynamic_cast<VertexBasedCellPopulation<SPACE_DIM> *>(this->mpCellPopulation))));

    // This is a magic number
    //double max_jiggle = 1e-4;

    if (SPACE_DIM != 1)
    {
        assert(SPACE_DIM == ELEMENT_DIM); // LCOV_EXCL_LINE
        assert(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM> *>(this->mpCellPopulation));

        c_vector<double, SPACE_DIM> ReferenceVector;
        ReferenceVector[0] = 1;
        ReferenceVector[1] = 0;

        //boundaryLinspace
        
        c_vector<double, 1000> BoundaryPointX;
        c_vector<double, 1000> BoundaryPointY;
        for (double i = 0; i < 1000; i++)
        {
            double tempT = 2 * 3.14159265359 * (i / 1000);
            double TempR = mBoundaryRadius + mOscRadius * sin(4 * tempT);
            int Index = (int)round(i);

            BoundaryPointX[Index] = TempR * cos(tempT) + mCentre[0];
            BoundaryPointY[Index] = TempR * sin(tempT) + mCentre[1];
        }

        // Iterate over all nodes and update their positions according to the boundary conditions
        unsigned num_nodes = this->mpCellPopulation->GetNumNodes();
        for (unsigned node_index = 0; node_index < num_nodes; node_index++)
        {
            Node<SPACE_DIM> *p_node = this->mpCellPopulation->GetNode(node_index);
            c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();

            c_vector<double, SPACE_DIM> node_location_trans = node_location - mCentre;

            double node_length_from_centre = norm_2(node_location - mCentre);

            if (node_length_from_centre > 1e-4)
            {

                double theta = atan2(node_location_trans[1], node_location_trans[0]);
                double OscilatingRad = mBoundaryRadius + mOscRadius * sin(4 * theta);

                double signed_dist_from_boundary = OscilatingRad - node_length_from_centre;

                if (signed_dist_from_boundary < 0.0)
                {
                    c_vector<double, SPACE_DIM> BoundaryPoint;
                    c_vector<double, SPACE_DIM> nearest_point;
                    double SmallestDist = 1e10;

                    for (double i = 0; i < 1000; i++)
                    {
                        BoundaryPoint[0] = BoundaryPointX[i];
                        BoundaryPoint[1] = BoundaryPointY[i];

                        double TempDist = norm_2(BoundaryPoint - node_location);

                        if (TempDist < SmallestDist)
                        {
                            SmallestDist = TempDist;
                            nearest_point = BoundaryPoint;
                        }
                    }

                    //check for node movement breaking the sim - cast downstream of abstract pop class to vertex
                    c_vector<double, SPACE_DIM> displacementVector = nearest_point - node_location;
                    dynamic_cast<VertexBasedCellPopulation<SPACE_DIM> *>(this->mpCellPopulation)->CheckForStepSizeException(node_index, displacementVector, 0.001);
                    p_node->rGetModifiableLocation() = nearest_point;
                }
            }
        }
    }

    else
    {
        // DIM == 1
        NEVER_REACHED;
        //DomeBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexBranchBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;
    c_vector<double, SPACE_DIM> ReferenceVector;
    ReferenceVector[0] = 1;
    ReferenceVector[1] = 0;
    unsigned num_nodes = this->mpCellPopulation->GetNumNodes();
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        Node<SPACE_DIM> *p_node = this->mpCellPopulation->GetNode(node_index);
        c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();
        c_vector<double, SPACE_DIM> node_location_trans = node_location - mCentre;

        double node_length_from_centre = norm_2(node_location - mCentre);

        if (node_length_from_centre > 1e-4)
        {

            double theta = atan2(node_location_trans[1], node_location_trans[0]);

            double OscilatingRad = mBoundaryRadius + mOscRadius * sin(4 * theta);

            double signed_dist_from_boundary = (OscilatingRad + 1e-2) - node_length_from_centre;
            if (signed_dist_from_boundary < 0)
            {
                std::cout << "EXCEPTION IN VERIFICATION" << std::endl;
                condition_satisfied = false;
                break;
            }
        }
    }
    return condition_satisfied;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexBranchBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream &rParamsFile)
{

    *rParamsFile << "\t\t\t<Centre>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mCentre[index] << ",";
    }
    *rParamsFile << mCentre[SPACE_DIM - 1] << "</Centre>\n";

    *rParamsFile << "\t\t\t<OscRadius>" << mOscRadius << "</OscRadius>\n";

    *rParamsFile << "\t\t\t<UseJiggledNodesOnPlane>" << mUseJiggledNodesOnPlane << "</UseJiggledNodesOnPlane>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class VertexBranchBoundaryCondition<1, 1>;
template class VertexBranchBoundaryCondition<1, 2>;
template class VertexBranchBoundaryCondition<2, 2>;
template class VertexBranchBoundaryCondition<1, 3>;
template class VertexBranchBoundaryCondition<2, 3>;
template class VertexBranchBoundaryCondition<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexBranchBoundaryCondition)
