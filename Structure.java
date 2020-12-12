package fem_abstract;

import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVectorRO;

import java.util.List;

import fem_1D.LineLoad;
import fem_1D.Node;
import shape_functions.ShapeFunctions;

public abstract class Structure {
	
	//public List<Node> listOfNodes = new ArrayList<Node>();
	//public List<Element> listOfElements = new ArrayList<Element>();
	
	//private IMatrixRO kGlobal = new Array2DMatrix(this.enumerateDOFs(),this.enumerateDOFs());
	//private IVectorRO rGlobal = new Array1DVector(this.enumerateDOFs());
	//private IVectorRO uGlobal = new Array1DVector(this.getNumberOfNodes() * 3);
	
	public abstract Node addNode(double x1, double x2, double x3);
	
	public abstract List<Node> getLON();
	
	public abstract int getNumberOfNodes();
	
	public abstract Node getNode(int id);
	
	public abstract void addElement(Element e) ;
	
	public abstract Element getElement(int id);
	
	public abstract List<Element> getLOE();
	
	public abstract int getNumberOfElements();
	
	public abstract void setSF(ShapeFunctions sf);
	
	public abstract void setLineLoad (double q, Node n1, Node n2, String a);
	
	public abstract List<LineLoad> getLOLL() ;
	
	public abstract void addSpring (double k, Node n, IVectorRO dir);
			
	public abstract IMatrixRO getC();
		
	public abstract void assignStartIndex();
	
	public abstract int enumerateDOFs();
		
	public abstract IMatrixRO assembleStiffnessMatrix();
	
	public abstract IVectorRO solve();
	
	public abstract void setDisplacement();

	public abstract void selectPlasticDisplacements(IVectorRO delta_u) ;

	public abstract void selectDisplacements() ;

	public abstract IVectorRO assembleExternalForces() ;

	public abstract IVectorRO assembleInternalForces() ;

	public abstract IMatrixRO assemblePlasticStiffnessMatrix() ;

}