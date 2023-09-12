#include "SmoothVertexBoundaries.hpp"
#include <AbstractTwoBodyInteractionForce.hpp>
#include "Cell.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "AbstractCellPopulation.hpp"
#include "MutableVertexMesh.hpp"
#include "VertexElement.hpp"
#include "Node.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
SmoothVertexBoundaries<DIM>::SmoothVertexBoundaries()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mSimulationEndTime(100),
      mSimulationDt(0.001)
{
}

template <unsigned DIM>
SmoothVertexBoundaries<DIM>::~SmoothVertexBoundaries()
{
}

template <unsigned DIM>
void SmoothVertexBoundaries<DIM>::SetSimulationEndTime(double simulationEndTime)
{
    assert(simulationEndTime > 0.0);
    mSimulationEndTime = simulationEndTime;
}

template <unsigned DIM>
double SmoothVertexBoundaries<DIM>::GetSimulationEndTime()
{
    return mSimulationEndTime;
}

template <unsigned DIM>
void SmoothVertexBoundaries<DIM>::SetSimulationDt(double simulationDt)
{
    assert(simulationDt> 0.0);
    mSimulationDt = simulationDt;
}

template <unsigned DIM>
double SmoothVertexBoundaries<DIM>::GetSimulationDt()
{
    return mSimulationDt;
}

template <unsigned DIM>
void SmoothVertexBoundaries<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM> &rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template <unsigned DIM>
void SmoothVertexBoundaries<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM> &rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template <unsigned DIM>
void SmoothVertexBoundaries<DIM>::UpdateCellData(AbstractCellPopulation<DIM, DIM> &rCellPopulation)
{

    //this should only occur at the final step
    double currentSimTim = SimulationTime::Instance()->GetTime();
    //double NumOfSteps = mSimulationEndTime/mSimulationDt;
    if (currentSimTim == (mSimulationEndTime - 2*mSimulationDt))
    {

        rCellPopulation.Update();
        //make sure this only operates on a vertex population
        assert(dynamic_cast<VertexBasedCellPopulation<DIM> *>(&rCellPopulation));
        // Helper variable that is a static cast of the cell population
        VertexBasedCellPopulation<DIM> *p_cell_population = static_cast<VertexBasedCellPopulation<DIM> *>(&rCellPopulation);
        MutableVertexMesh<DIM, DIM> &CurrentMesh = static_cast<VertexBasedCellPopulation<DIM> *>(&rCellPopulation)->rGetMesh();

        std::vector<unsigned> InterpolatedNodeA;
        std::vector<unsigned> InterpolatedNodeB;

        unsigned NewIndex = p_cell_population->GetNumNodes();
        //unsigned NewNodeIndexCounter = p_cell_population->rGetMesh().GetNodeIteratorEnd();
        for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter_A = p_cell_population->rGetMesh().GetNodeIteratorBegin();
             node_iter_A != p_cell_population->rGetMesh().GetNodeIteratorEnd();
             ++node_iter_A)
        {
            unsigned NodeAindex = node_iter_A->GetIndex();
            Node<DIM> *ThisNodeA = p_cell_population->GetNode(NodeAindex);
            // ITERATE ONLY OVER THE NEIGHBOURS OF A

            // Get set of neighbours.
            //std::vector<unsigned>& near_nodes_ind = node_iter_A->rGetNeighbouringNodeIndices();

            std::set<unsigned> near_nodes_ind = p_cell_population->GetNeighbouringNodeIndices(NodeAindex);

            for (std::set<unsigned>::iterator NodeBindex_iter = near_nodes_ind.begin();
                 NodeBindex_iter != near_nodes_ind.end(); ++NodeBindex_iter)
            {
                Node<DIM> *ThisNodeB  = p_cell_population->GetNode(*NodeBindex_iter);
                unsigned NodeBindex = ThisNodeB ->GetIndex();
                Node<DIM> *node_iter_B  = p_cell_population->GetNode(*NodeBindex_iter);
                if (NodeAindex != NodeBindex)
                {

                    std::set<unsigned> elements_containing_nodeA = node_iter_A->rGetContainingElementIndices();
                    std::set<unsigned> elements_containing_nodeB = node_iter_B->rGetContainingElementIndices();

                    // Find common elements
                    std::set<unsigned> shared_elements;
                    std::set_intersection(elements_containing_nodeA.begin(),
                                          elements_containing_nodeA.end(),
                                          elements_containing_nodeB.begin(),
                                          elements_containing_nodeB.end(),
                                          std::inserter(shared_elements, shared_elements.begin()));

                    // THEN IS THE BOUNDARY ELEMENT - SEE MUTABLE VERTEX CPP for adding in new nodes
                    if (shared_elements.size() == 1)
                    {   

                        
                        
                        //KEEP TRACK OF THE NODES ITERATED OVER
                        std::vector<unsigned>::iterator Finding_NodeB_in_SetA = std::find(InterpolatedNodeA.begin(), InterpolatedNodeA.end(), NodeBindex);
                        std::vector<unsigned>::iterator Finding_NodeA_in_SetB = std::find(InterpolatedNodeB.begin(), InterpolatedNodeB.end(), NodeAindex);

                        if (Finding_NodeB_in_SetA != InterpolatedNodeA.end() && Finding_NodeA_in_SetB != InterpolatedNodeB.end())
                        {
                            if (Finding_NodeB_in_SetA - InterpolatedNodeA.begin() != Finding_NodeA_in_SetB - InterpolatedNodeB.begin())
                            {   

                                CurrentMesh.DivideEdge(ThisNodeA , ThisNodeB);
                                
                            }
                        }
                        else
                        {
                            //pair has not met yet
                           // CurrentMesh.DivideEdge(ThisNodeA , ThisNodeB );
                        }

                        InterpolatedNodeA.push_back(NodeAindex);
                        InterpolatedNodeB.push_back(NodeBindex);

                        //end of just one common element...
                    }
                }
            }
        }
    }

    rCellPopulation.Update();
}

template <unsigned DIM>
void SmoothVertexBoundaries<DIM>::OutputSimulationModifierParameters(out_stream &rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class SmoothVertexBoundaries<1>;
template class SmoothVertexBoundaries<2>;
template class SmoothVertexBoundaries<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SmoothVertexBoundaries)