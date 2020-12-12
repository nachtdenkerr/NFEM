package fem_abstract;	

import java.util.List;

import fem_1D.GaussPoint;
import fem_1D.Node;
import shape_functions.ShapeFunctions;
import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVectorRO;

public abstract class Element{
	public abstract void enumerateDOFs();

	public abstract int[] getDOFNumbers();
	
	public abstract IMatrixRO computeStiffnessMatrix();
	
	public abstract double getEModulus();
	
	public abstract double getNuy();

	public abstract IMatrixRO getC();
	
	public abstract List<Node> getLON();
	
	public abstract int getStart();

	public abstract IMatrixRO computeStiffnessMatrix(ShapeFunctions sf);

	public abstract IVectorRO computeInternalForce(ShapeFunctions sf);

	public abstract void setDisplacement() ;

	public abstract void computeBMatrix(ShapeFunctions sf) ;

	public abstract List<GaussPoint> getLOGP() ;

	public abstract IVectorRO getDisplacement() ;

	public abstract double getForce() ;

	public abstract IVectorRO computeInternalForce() ;

	public abstract void assignPlasticGaussPoints(int i);
}