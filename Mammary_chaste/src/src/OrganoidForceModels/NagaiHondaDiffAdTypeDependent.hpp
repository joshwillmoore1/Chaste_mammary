#ifndef NAGAIHONDADIFFADTYPEDEPENDENT_HPP_
#define NAGAIHONDADIFFADTYPEDEPENDENT_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "NagaiHondaForce.hpp"

#include <iostream>

/**
 * A force class for use in vertex-based simulations, based on a model
 * model proposed by T. Nagai and H. Honda ("A dynamic cell model for
 * the formation of epithelial tissues", Philosophical Magazine Part B
 * 81:699-719) to include differential adhesion between normal and
 * labelled cells. To include differential adhesion we override the
 * GetAdhesionParameter() method.
 *
 * Each of the model parameter member variables are rescaled such that
 * mDampingConstantNormal takes the default value 1, whereas Nagai and
 * Honda (who denote the parameter by nu) take the value 0.01.
 */
template <unsigned DIM>
class NagaiHondaDiffAdTypeDependent : public NagaiHondaForce<DIM>
{
private:
    /**
     * Cell-cell adhesion energy parameter for two labelled cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaBasalLuminalAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for labelled and non-labelled cells.
     * Has has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaLuminalLuminalAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter for labelled cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaBasalBasalAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter for labelled cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaGhostGhostAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter for labelled cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaBasalBoundaryAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter for labelled cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaLuminalBoundaryAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter for labelled cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaGhostBoundaryAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter for labelled cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaLuminalGhostAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter for labelled cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaBasalGhostAdhesionEnergyParameter;

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive &boost::serialization::base_object<NagaiHondaForce<DIM>>(*this);
    }

public:
    /**
     * Constructor.
     */
    NagaiHondaDiffAdTypeDependent();

    /**
     * Destructor.
     */
    virtual ~NagaiHondaDiffAdTypeDependent();

    /**
     * Overridden GetAdhesionParameter() method.
     *
     * Get the adhesion parameter for the edge between two given nodes. Depends
     * on the type of cells attached to the elements.
     *
     * @param pNodeA one node
     * @param pNodeB the other node
     * @param rVertexCellPopulation reference to the cell population
     *
     * @return the adhesion parameter for this edge.
     */
    virtual double GetAdhesionParameter(Node<DIM> *pNodeA, Node<DIM> *pNodeB, VertexBasedCellPopulation<DIM> &rVertexCellPopulation);

    /**
     * @return mNagaiHondaLabelledCellCellAdhesionEnergyParameter
     */
    double GetNagaiHondaBasalLuminalAdhesionEnergyParameter();

    /**
     * @return mNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter
     */
    double GetNagaiHondaLuminalLuminalAdhesionEnergyParameter();

    /**
     * @return mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter
     */
    double GetNagaiHondaBasalBasalAdhesionEnergyParameter();

    /**
     * @return mNagaiHondaLabelledCellCellAdhesionEnergyParameter
     */
    double GetNagaiHondaGhostGhostAdhesionEnergyParameter();

    /**
     * @return mNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter
     */
    double GetNagaiHondaBasalBoundaryAdhesionEnergyParameter();

    /**
     * @return mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter
     */
    double GetNagaiHondaLuminalBoundaryAdhesionEnergyParameter();

    /**
     * @return mNagaiHondaLabelledCellCellAdhesionEnergyParameter
     */
    double GetNagaiHondaGhostBoundaryAdhesionEnergyParameter();

    /**
     * @return mNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter
     */
    double GetNagaiHondaLuminalGhostAdhesionEnergyParameter();

    /**
     * @return mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter
     */
    double GetNagaiHondaBasalGhostAdhesionEnergyParameter();
    //Setting functions

    /**
     * Set mNagaiHondaLabelledCellCellAdhesionEnergyParameter.
     *
     * @param NagaiHondaBasalLuminalAdhesionEnergyParameter the new value of mNagaiHondaLabelledCellCellAdhesionEnergyParameter
     */
    void SetNagaiHondaBasalLuminalAdhesionEnergyParameter(double NagaiHondaBasalLuminalAdhesionEnergyParameter);

    /**
     * Set mNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter.
     *
     * @param NagaiHondaLuminalLuminalAdhesionEnergyParameter the new value of mNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter
     */
    void SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(double NagaiHondaLuminalLuminalAdhesionEnergyParameter);

    /**
     * Set mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter.
     *
     * @param NagaiHondaBasalBasalAdhesionEnergyParameter the new value of mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter
     */
    void SetNagaiHondaBasalBasalAdhesionEnergyParameter(double NagaiHondaBasalBasalAdhesionEnergyParameter);

    /**
     * Set mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter.
     *
     * @param NagaiHondaGhostGhostAdhesionEnergyParameter the new value of mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter
     */
    void SetNagaiHondaGhostGhostAdhesionEnergyParameter(double NagaiHondaGhostGhostAdhesionEnergyParameter);

    /**
     * Set mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter.
     *
     * @param NagaiHondaBasalBoundaryAdhesionEnergyParameter the new value of mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter
     */
    void SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(double NagaiHondaBasalBoundaryAdhesionEnergyParameter);

    /**
     * Set mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter.
     *
     * @param NagaiHondaLuminalBoundaryAdhesionEnergyParameter the new value of mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter
     */
    void SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(double NagaiHondaLuminalBoundaryAdhesionEnergyParameter);

    /**
     * Set mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter.
     *
     * @param NagaiHondaGhostBoundaryAdhesionEnergyParameter the new value of mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter
     */
    void SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(double NagaiHondaGhostBoundaryAdhesionEnergyParameter);

    /**
     * Set mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter.
     *
     * @param NagaiHondaLuminalGhostAdhesionEnergyParameter the new value of mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter
     */
    void SetNagaiHondaLuminalGhostAdhesionEnergyParameter(double NagaiHondaLuminalGhostAdhesionEnergyParameter);

    /**
     * Set mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter.
     *
     * @param NagaiHondaBasalGhostAdhesionEnergyParameter the new value of mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter
     */
    void SetNagaiHondaBasalGhostAdhesionEnergyParameter(double NagaiHondaBasalGhostAdhesionEnergyParameter);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NagaiHondaDiffAdTypeDependent)

#endif /*NAGAIHONDADIFFADTYPEDEPENDENT_HPP_*/
