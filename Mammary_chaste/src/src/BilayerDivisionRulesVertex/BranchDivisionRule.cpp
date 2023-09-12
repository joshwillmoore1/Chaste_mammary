#include "BranchDivisionRule.hpp"
#include "LumenERPositiveCellProliferativeType.hpp"
#include "LumenERNegativeCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "GhostCellProliferativeType.hpp"
#include "MyoEpiCellProliferativeType.hpp"
#include "BasalStemCellProliferativeType.hpp"

template <unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> BranchDivisionRule<SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    VertexBasedCellPopulation<SPACE_DIM> &rCellPopulation)
{

    std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(pParentCell);

    //loop through the neighbours of the cell - find the centroids of the neighbouring cells of different types
    //take the average of the centroids and return the division vector from the centre of the main cell and the averaged point
    if (!neighbour_indices.empty())
    {

        std::vector<double> NeighbourCentreX;
        std::vector<double> NeighbourCentreY;

        for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
             iter != neighbour_indices.end();
             ++iter)
        {

            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);

            if (pParentCell->template HasCellProperty<BasalStemCellProliferativeType>() || pParentCell->template HasCellProperty<MyoEpiCellProliferativeType>())
            {

                if (p_cell->template HasCellProperty<LumenERNegativeCellProliferativeType>() || p_cell->template HasCellProperty<LumenERPositiveCellProliferativeType>())
                {

                    // then get the centroid coordinates
                    c_vector<double, SPACE_DIM> neighbour_location = rCellPopulation.GetLocationOfCellCentre(p_cell);
                    NeighbourCentreX.push_back(neighbour_location[0]);
                    NeighbourCentreY.push_back(neighbour_location[1]);
                }
            }
            else
            {
                if (p_cell->template HasCellProperty<BasalStemCellProliferativeType>() || p_cell->template HasCellProperty<MyoEpiCellProliferativeType>())
                {
                    c_vector<double, SPACE_DIM> neighbour_location = rCellPopulation.GetLocationOfCellCentre(p_cell);
                    NeighbourCentreX.push_back(neighbour_location[0]);
                    NeighbourCentreY.push_back(neighbour_location[1]);
                }
            }
        }

        // if there are heterotypic connections pass the perpendicular vector to AB
        if (!NeighbourCentreX.empty())
        {
            double MeanX = std::accumulate(NeighbourCentreX.begin(), NeighbourCentreX.end(), 0.0) / NeighbourCentreX.size();
            double MeanY = std::accumulate(NeighbourCentreY.begin(), NeighbourCentreY.end(), 0.0) / NeighbourCentreX.size();
            c_vector<double, SPACE_DIM> MeanVec;
            MeanVec[0] = MeanX;
            MeanVec[1] = MeanY; 

            c_vector<double, SPACE_DIM> pCellCentre = rCellPopulation.GetLocationOfCellCentre(pParentCell);
            c_vector<double, SPACE_DIM> VectAB = pCellCentre - MeanVec;

            c_vector<double, SPACE_DIM> VectABNorm;

            if (pParentCell->template HasCellProperty<BasalStemCellProliferativeType>() || pParentCell->template HasCellProperty<MyoEpiCellProliferativeType>())
            {
                VectABNorm[0] = VectAB[0]/norm_2(VectAB);
                VectABNorm[1] = VectAB[1]/norm_2(VectAB);

            }
            else {
                // add noise

                if (RandomNumberGenerator::Instance()->ranf() > 0.8)
                {
                    VectABNorm[0] = VectAB[1]/norm_2(VectAB);
                    VectABNorm[1] = -VectAB[0]/norm_2(VectAB);

                }
                else
                {
                VectABNorm[0] = VectAB[0]/norm_2(VectAB);
                VectABNorm[1] = VectAB[1]/norm_2(VectAB);
                }

            }

            return VectABNorm;

        } else 
        {
            unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(pParentCell);
            return rCellPopulation.rGetMesh().GetShortAxisOfElement(elem_index);
        }


    }
    else
    {

        unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(pParentCell);
        return rCellPopulation.rGetMesh().GetShortAxisOfElement(elem_index);
    }

   

}

// Explicit instantiation
template class BranchDivisionRule<1>;
template class BranchDivisionRule<2>;
template class BranchDivisionRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BranchDivisionRule)
