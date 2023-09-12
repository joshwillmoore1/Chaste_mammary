#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "CollierSrnModel.hpp"
#include "CollierTrackingModifier.hpp"
#include "CollierFixedTrackingModifier.hpp"
#include "LumenKiller.hpp"
#include "BasalStemCellProliferativeType.hpp"
#include "LumenERPositiveCellProliferativeType.hpp"
#include "LumenERNegativeCellProliferativeType.hpp"
#include "MyoEpiCellProliferativeType.hpp"
#include "GhostCellProliferativeType.hpp"
#include "OrganoidG1FixedCellCycleModel.hpp"
#include "CellData.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "WildTypeCellMutationState.hpp"
#include "RandomMotionForce.hpp"
#include "CellLabel.hpp"
#include "CellAgesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "NagaiHondaDiffAdTypeDependent.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "CellsGenerator.hpp"

//notch dependent adhesion
//#include "DifferentialNotchAdhesionForce.hpp"

//Normal division vector
#include "NormalCentreDivisionPlaneRule.hpp"
#include "ApicalBasalDivisionRule.hpp"
#include "LumenPressureForce.hpp"
#include "LumenApop.hpp"
#include "BoundaryCellTypeMod.hpp"

#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"

//boundary conditions
#include "VertexBranchBoundaryCondition.hpp"
#include "SmoothVertexBoundaries.hpp"
#include "BoundaryCurvatureCalculation.hpp"

class TestGrowingBilayerNotch : public AbstractCellBasedTestSuite
{
    /*
* This test is used to generate the growing bilayer with ghost cells
* for lumen formation. The polaried collier model is also applied.
*
*/

public:
    const double spring_const = 25;
    const double Random_pert = 0.003;

    const double Notch_hss = 0.1936;
    const double Delta_hss = 0.049;

    const double basal_Delta_IC = 0.2;
    const double luminal_Delta_IC = 0.1;

    const double basal_Notch_IC = 0.1;
    const double luminal_Notch_IC = 0.2;

    const int End_time = 100; //100;
    const double time_step = 0.01;
    const double sample_step = 50;

    const double lum_radius = 1;
    const int cells_across = 21;
    const int cells_up = 21;

    const double preBirthTime = 5.0;

    //CELL FORCE PARAMETERS
    //this is for target area
    const double areaEnergy = 30.0;

    //this is for rounded shapes
    const double SurfaceEnergy = 50;

    //differential adhesion
    const double BasalBasalEnergy = 1;
    const double LuminalLuminalEnergy = 1;
    const double GhostGhostEnergy = 1;

    const double BasalLuminalEnergy = 1;
    const double GhostLuminalEnergy = 1;
    const double GhostBasalEnergy = 0;

    const double GhostboundaryEnergy = 0;
    const double BasalboundaryEnergy = 3;
    const double LuminalboundaryEnergy = 1;

    // LUMINAL FORCE PARAS
    const double force_decay = 0.05; //0.3
    const double force_scale = 0.2;  //0.25

    const double SeedValue = 15; //6

    //Values for boundary
    const double circleRadius = 5 * lum_radius;
    const double oscRadius = 0;
    const int NumberOfBoundarySteps = 11;

        //see ipad notes for seed ordering

    void TestVertexBoundaryBilayersPolarised_1()
    {

        std::string TestoutputName = "Test_adaptive_polarity_growning_" + std::to_string(1);
        RandomNumberGenerator::Instance()->Reseed(4);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11 * 0.95);
        simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        /*MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            */

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

     void TestVertexBoundaryBilayersPolarised_2()
    {

        std::string TestoutputName = "Test_adaptive_polarity_growning_" + std::to_string(2);
        RandomNumberGenerator::Instance()->Reseed(2);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11 * 0.95);
        simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        /*MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            */

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

     void TestVertexBoundaryBilayersPolarised_3()
    {

        std::string TestoutputName = "Test_adaptive_polarity_growning_" + std::to_string(3);
        RandomNumberGenerator::Instance()->Reseed(3);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11 * 0.95);
        simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        /*MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            */

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

     void TestVertexBoundaryBilayersPolarised_4()
    {

        std::string TestoutputName = "Test_adaptive_polarity_growning_" + std::to_string(4);
        RandomNumberGenerator::Instance()->Reseed(7);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11 * 0.95);
        simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        /*MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            */

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

 void TestVertexBoundaryBilayersPolarised_5()
    {

        std::string TestoutputName = "Test_adaptive_polarity_growning_" + std::to_string(5);
        RandomNumberGenerator::Instance()->Reseed(8);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11 * 0.95);
        simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        /*MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            */

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

    

     void TestVertexBoundaryBilayersPolarised_6()
    {

        std::string TestoutputName = "Test_adaptive_polarity_growning_" + std::to_string(6);
        RandomNumberGenerator::Instance()->Reseed(10);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11 * 0.95);
        simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        /*MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            */

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

    
     void TestVertexBoundaryBilayersPolarised_7()
    {

        std::string TestoutputName = "Test_adaptive_polarity_growning_" + std::to_string(7);
        RandomNumberGenerator::Instance()->Reseed(11);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11 * 0.95);
        simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        /*MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            */

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

    

     void TestVertexBoundaryBilayersPolarised_8()
    {

        std::string TestoutputName = "Test_adaptive_polarity_growning_" + std::to_string(8);
        RandomNumberGenerator::Instance()->Reseed(12);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11 * 0.95);
        simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        /*MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            */

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

    
 void TestVertexBoundaryBilayersPolarised_9()
    {

        std::string TestoutputName = "Test_adaptive_polarity_growning_" + std::to_string(9);
        RandomNumberGenerator::Instance()->Reseed(14);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11 * 0.95);
        simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        /*MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            */

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

     void TestVertexBoundaryBilayersPolarised_10()
    {

        std::string TestoutputName = "Test_adaptive_polarity_growning_" + std::to_string(10);
        RandomNumberGenerator::Instance()->Reseed(15);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11 * 0.95);
        simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        /*MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            */

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

    

    /////////////////////// FIXED SIMULATIONS ///////////////////////////////////

   void TestVertexBoundaryBilayersFixed_1()
    {

        std::string TestoutputName = "Test_fixed_polarity_growning_" + std::to_string(1);
        RandomNumberGenerator::Instance()->Reseed(4);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

       // MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        //p_modifier->SetPolarisationParameter(0.11 * 0.95);
        //simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

     void TestVertexBoundaryBilayersFixed_2()
    {

        std::string TestoutputName = "Test_fixed_polarity_growning_" + std::to_string(2);
        RandomNumberGenerator::Instance()->Reseed(2);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        //MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        //p_modifier->SetPolarisationParameter(0.11 * 0.95);
        //simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

     void TestVertexBoundaryBilayersFixed_3()
    {

        std::string TestoutputName = "Test_fixed_polarity_growning_" + std::to_string(3);
        RandomNumberGenerator::Instance()->Reseed(3);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        //MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        //p_modifier->SetPolarisationParameter(0.11 * 0.95);
        //simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

     void TestVertexBoundaryBilayersFixed_4()
    {

        std::string TestoutputName = "Test_fixed_polarity_growning_" + std::to_string(4);
        RandomNumberGenerator::Instance()->Reseed(7);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        //MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        //p_modifier->SetPolarisationParameter(0.11 * 0.95);
        //simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

 void TestVertexBoundaryBilayersFixed_5()
    {

        std::string TestoutputName = "Test_fixed_polarity_growning_" + std::to_string(5);
        RandomNumberGenerator::Instance()->Reseed(8);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

       // MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
       // p_modifier->SetPolarisationParameter(0.11 * 0.95);
       // simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

    

     void TestVertexBoundaryBilayersFixed_6()
    {

        std::string TestoutputName = "Test_fixed_polarity_growning_" + std::to_string(6);
        RandomNumberGenerator::Instance()->Reseed(10);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

       // MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
       // p_modifier->SetPolarisationParameter(0.11 * 0.95);
       // simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
        

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

    
     void TestVertexBoundaryBilayersFixed_7()
    {

        std::string TestoutputName = "Test_adaptive_fixed_growning_" + std::to_string(7);
        RandomNumberGenerator::Instance()->Reseed(11);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        //MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        //p_modifier->SetPolarisationParameter(0.11 * 0.95);
        //simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

    

     void TestVertexBoundaryBilayersFixed_8()
    {

        std::string TestoutputName = "Test_fixed_polarity_growning_" + std::to_string(8);
        RandomNumberGenerator::Instance()->Reseed(12);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        //MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        //p_modifier->SetPolarisationParameter(0.11 * 0.95);
        //simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

    
 void TestVertexBoundaryBilayersFixed9()
    {

        std::string TestoutputName = "Test_fixed_polarity_growning_" + std::to_string(9);
        RandomNumberGenerator::Instance()->Reseed(14);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

       // MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
       // p_modifier->SetPolarisationParameter(0.11 * 0.95);
       // simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
            

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

     void TestVertexBoundaryBilayersFixed_10()
    {

        std::string TestoutputName = "Test_fixed_polarity_growning_" + std::to_string(10);
        RandomNumberGenerator::Instance()->Reseed(15);

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(GhostCellProliferativeType, p_ghost_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());

        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //AB polarity - use this for growing bilayer simulations
        MAKE_PTR(ApicalBasalDivisionRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double, 2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location - centre);

            if (dist_to_centre < Centre_tol)
            {
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;
            }
        }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        c_vector<double, 2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Get distance from centre of cell population
            c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double r = norm_2(location - centre_cell);
            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
            {
                //slight biase towards bilayer
                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                cell_iter->SetCellProliferativeType(p_BSC_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                //initial_conditions.push_back(basal_Notch_IC);
                //initial_conditions.push_back(basal_Delta_IC);
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 0.85);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }
            if (r <= lum_radius + 1e-5 && r > lum_radius - 0.5)
            {

                double pert1 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 0.8);
                double pert2 = 0.01 * (2 * RandomNumberGenerator::Instance()->ranf() - 1.2);
                //cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_Lneg_type);
                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions;
                initial_conditions.push_back(Notch_hss + pert1);
                initial_conditions.push_back(Delta_hss + pert2);
                p_srn_model->SetInitialConditions(initial_conditions);
                cell_iter->SetSrnModel(p_srn_model);

                double birth_time = RandomNumberGenerator::Instance()->ranf() * -preBirthTime;
                cell_iter->SetBirthTime(birth_time);
            }

            if (r > lum_radius + 1.2)
            {
                cell_iter->Kill();
            }

            if (r < lum_radius - 0.5)
            {

                CollierSrnModel *p_srn_model = new CollierSrnModel();
                std::vector<double> initial_conditions_override;
                cell_iter->SetCellProliferativeType(p_ghost_type);
                initial_conditions_override.push_back(0);
                initial_conditions_override.push_back(0);
                p_srn_model->SetInitialConditions(initial_conditions_override);
                cell_iter->SetSrnModel(p_srn_model);
                cell_iter->GetCellData()->SetItem("target area", 1);
            }
        }
        cell_population.RemoveDeadCells();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(TestoutputName);

        //ADAPTIVE POLARITY

        //MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        //p_modifier->SetPolarisationParameter(0.11 * 0.95);
        //simulator.AddSimulationModifier(p_modifier);

        //FIXED POLARITY
        MAKE_PTR(CollierFixedTrackingModifier<2>, p_modifier);
            p_modifier->SetSameCellW1(0.11 * 0.95);
            simulator.AddSimulationModifier(p_modifier);
        

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        MAKE_PTR(NagaiHondaDiffAdTypeDependent<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(areaEnergy);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);

        p_force->SetNagaiHondaBasalBasalAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaGhostGhostAdhesionEnergyParameter(GhostGhostEnergy);

        p_force->SetNagaiHondaBasalLuminalAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaBasalGhostAdhesionEnergyParameter(GhostBasalEnergy);
        p_force->SetNagaiHondaLuminalGhostAdhesionEnergyParameter(GhostLuminalEnergy);

        p_force->SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(BasalboundaryEnergy);
        p_force->SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(LuminalboundaryEnergy);
        p_force->SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(GhostboundaryEnergy);

        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(BoundaryCellTypeMod<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //Luminal pressure
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

        // Run simulation
        simulator.Solve();
    }

};