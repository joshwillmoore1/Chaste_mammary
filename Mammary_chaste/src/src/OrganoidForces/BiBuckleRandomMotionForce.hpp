#ifndef BIBUCKLERANDOMMOTIONFORCE_HPP_
	#define BIBUCKLERANDOMMOTIONFORCE_HPP_
	
	#include "ChasteSerialization.hpp"
	#include <boost/serialization/base_object.hpp>
	
	#include "AbstractForce.hpp"
	#include "AbstractOffLatticeCellPopulation.hpp"
	#include "RandomNumberGenerator.hpp"
	
	/**
12	 * A force class to model random cell movement.
13	 */
	template<unsigned DIM>
	class BiBuckleRandomMotionForce : public AbstractForce<DIM>
	{
	private :
	
	    /**
20	     * Random Movement Parameter.
21	     */
	    double mMovementParameter;
		    /**
25	     * Archiving.
26	     */
	    friend class boost::serialization::access;
	    /**
29	     * Boost Serialization method for archiving/checkpointing.
30	     * Archives the object and its member variables.
31	     *
32	     * @param archive  The boost archive.
33	     * @param version  The current version of this class.
34	     */
	    template<class Archive>
	    void serialize(Archive & archive, const unsigned int version)
	    {
	        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
	        archive & mMovementParameter;
	    }
	
	public :
	
	    /**
45	     * Constructor.
46	     */
	    BiBuckleRandomMotionForce();
	
	    /**
50	     * Destructor.
51	     */
	    ~BiBuckleRandomMotionForce();
	
	    /**
55	     * Set the diffusion constant for the cells.
56	     *
57	     * @param movementParameter the movement parameter to use
58	     */
	    void SetMovementParameter(double movementParameter);
	
	    /**
62	     * Get the random motion coefficient.
63	     *
64	     * @return mMovementParameter
65	     */
	    double GetMovementParameter();
	
	    /**
69	     * Overridden AddForceContribution() method.
70	     *
71	     * @param rCellPopulation reference to the tissue
72	     *
73	     */
	    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);
	
	    /**
77	     * Overridden OutputForceParameters() method.
78	     *
79	     * @param rParamsFile the file stream to which the parameters are output
80	     */
	    void OutputForceParameters(out_stream& rParamsFile);
	};
	
	#include "SerializationExportWrapper.hpp"
	EXPORT_TEMPLATE_CLASS_SAME_DIMS(BiBuckleRandomMotionForce)
	
#endif /*BIBUCKLERANDOMMOTIONFORCE_HPP_*/