#include "MyoEpiCellProliferativeType.hpp"

MyoEpiCellProliferativeType::MyoEpiCellProliferativeType()
    : AbstractCellProliferativeType(40)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(MyoEpiCellProliferativeType)
