#ifndef LUMENERNEGATIVECELLPROLIFERATIVETYPE_HPP_
#define LUMENERNEGATIVECELLPROLIFERATIVETYPE_HPP_

#include "AbstractCellProliferativeType.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Subclass of AbstractCellProliferativeType defining a 'basal stem' cell.
 * This is the generic stem cell within the mammary organoid
 */
class LumenERNegativeCellProliferativeType : public AbstractCellProliferativeType
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell proliferative type.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProliferativeType>(*this);
    }

public:
    /**
     * Constructor.
     */
    LumenERNegativeCellProliferativeType();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(LumenERNegativeCellProliferativeType)

#endif /*LUMENERNEGATIVECELLPROLIFERATIVETYPE_HPP_*/
