package fem_3D;

import iceb.jnumerics.IMatrix;
import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVector;
import iceb.jnumerics.IVectorRO;
import iceb.jnumerics.MatrixFormat;
import iceb.jnumerics.Vector3D;

import java.util.ArrayList;
import java.util.List;

import fem_1D.Force;
import fem_1D.LineLoad;
import fem_1D.Node;
import fem_1D.Spring;
import shape_functions.ShapeFunctions;
import solver.LUDecomposition;
import linalg.Array1DVector;
import linalg.Array2DMatrix;
import fem_abstract.Element;
import fem_abstract.Structure;

public class BrickStructure extends Structure{
	
	public List<Node> listOfNodes = new ArrayList<Node>();
	public List<Element> listOfElements = new ArrayList<Element>();
	public List<LineLoad> listOfLLoads = new ArrayList<LineLoad>();
	public List<Spring> listOfSprings = new ArrayList<Spring>();

	private IMatrixRO kGlobal;
	private IVectorRO rGlobal;
	private IVectorRO uGlobal;
	
	private double eModulus, nuy;
	private IMatrixRO C_matrix;
	private ShapeFunctions sf;
	
	public Node addNode(double x1, double x2, double x3) {
		Node n = new Node(x1,x2,x3);
		listOfNodes.add(n);
		return n;
	}
	
	public void addElement(Element e){
		this.listOfElements.add(e);
	}
	
	public void addParameters(double E, double n){
		this.eModulus = E;
		this.nuy = n;
		
		IMatrix C = new Array2DMatrix(6,6);
		double c = eModulus/((1+this.nuy)*(1-2*this.nuy));
		for (int i = 0; i < 3; i++){
			C.set(i, i, 1- this.nuy);
			C.set(i+3, i+3, (1 - 2*this.nuy)/2);
			for (int j = 0; j < 3; j++){
				if (i != j){
					C.set(i, j, this.nuy);
				}
			}
		}
		C_matrix = (IMatrix) C.multiply(c);
	}
	
	public IMatrixRO getC(){
		return this.C_matrix;
	}
	
	public void setSF(ShapeFunctions sf){
		this.sf = sf;
	}
	
	public ShapeFunctions getSF(){
		return this.sf;
	}

	public void setLineLoad(double q, Node n1, Node n2, String a) {
		LineLoad ll = new LineLoad (q,n1,n2,a);
		listOfLLoads.add(ll);
		double[] load = ll.computeLineLoad();
		Force f = new Force (load[0],load[1],load[2]);
		n1.setForce(f);
		n2.setForce(f);
		n1.setLineLoad(true);
		n2.setLineLoad(true);
		listOfLLoads.add(ll);
	}
	
	public List<LineLoad> getLOLL() {
		return listOfLLoads;
	}

	
	public void addSpring(double k, Node n, IVector dir) {
		Spring spr = new Spring (k,n,dir);
		listOfSprings.add(spr);
		
		IMatrix kSpring = new Array2DMatrix(3,3);
		IVectorRO c = new Array1DVector(dir.get(0)/dir.normTwo(),
				dir.get(1)/dir.normTwo(),dir.get(2)/dir.normTwo());
		kSpring = (IMatrix) c.dyadicProduct(c);
		
		int[] ndof = n.getDOFNumbers();
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				if (ndof[i] != -1 && ndof[j] != -1){
					((Array2DMatrix)kGlobal).add(ndof[i], ndof[j], kSpring.get(i, j));
				}
			}
		}
	}

	
	public List<Node> getLON() {
		return listOfNodes;
	}

	
	public int getNumberOfNodes() {
		return listOfNodes.size();
	}

	
	public Node getNode(int id) {
		return listOfNodes.get(id);
	}

	
	public List<Element> getLOE() {
		return listOfElements;
	}

	
	public int getNumberOfElements() {
		return listOfElements.size();
	}

	
	public Element getElement(int id) {
		return listOfElements.get(id);
	}

	
	public void assignStartIndex() {
		if (this.getNumberOfNodes() > 1){
			int max_index = 0;
			listOfNodes.get(0).setCounter(max_index);
			for (int i = 0; i < listOfNodes.size()-1; i++){
				int start_next = listOfNodes.get(i).enumerateDOFs(max_index);
				if (start_next >= max_index) {
					max_index = start_next;
					listOfNodes.get(i+1).setCounter(max_index);
				}
			}
		}
	}

	
	public int enumerateDOFs() {
		this.assignStartIndex();
		int start = 0;
		int nDOF = 0;
		for (int i = 0; i < this.getLON().size(); i++){
			int a = this.getLON().get(i).enumerateDOFs(start);
			nDOF = nDOF + a;
		} 
		return nDOF;
	}

		
	public IMatrixRO assembleStiffnessMatrix() {
		IMatrixRO kGlobal = new Array2DMatrix (this.enumerateDOFs(), this.enumerateDOFs()); 

		for (Element e : listOfElements){
			IMatrixRO kElement = e.computeStiffnessMatrix(this.sf);
			int[] eleDOF = e.getDOFNumbers();
			for (int i = 0; i < eleDOF.length; i++){
				for (int j = 0; j< eleDOF.length; j++){
					if (eleDOF[i] != -1 && eleDOF[j] != -1){
						((Array2DMatrix)kGlobal).add(eleDOF[i], eleDOF[j], kElement.get(i, j));
					}
				}
			}
		//this.kGlobal = kGlobal;
		}
		return kGlobal;
	}
	
	public IVectorRO assembleInternalForces () {
		IVectorRO riGlobal = new Array1DVector(this.enumerateDOFs());
		//from Node select DOF which is non -1, get the index, plus to the R_global with the same DoF
		
		for (Element e : listOfElements){
			IVectorRO riElement = e.computeInternalForce(this.sf);
			int[] eleDOF = e.getDOFNumbers();
			for (int i = 0; i < eleDOF.length; i++){
				if (eleDOF[i] != -1){
					((Array1DVector) riGlobal).add(eleDOF[i], riElement.get(i));
				}
			}
//			System.out.println("r_int " + MatrixFormat.format(riElement));
		}

		return riGlobal;
	}

	
	public IVectorRO solve() {
		IMatrixRO K = this.assembleStiffnessMatrix();
		IVectorRO r = this.assembleExternalForces();
		
		LUDecomposition lud = new LUDecomposition(K);
		IVectorRO u = lud.solveFor(r);
		uGlobal = u;
		return u;
	}
	
	public IVectorRO solvePlastic() {
		
		IMatrixRO K = this.assembleStiffnessMatrix();
		IVectorRO r = this.assembleExternalForces();
		IVectorRO ri = this.assembleInternalForces();

		LUDecomposition lud = new LUDecomposition(K);
		IVectorRO delta_u = lud.solveFor(r.subtract(ri));
		
		System.out.println("r_ext " + MatrixFormat.format(r));
		System.out.println("r_int " + MatrixFormat.format(ri));
		//System.out.println(MatrixFormat.format(K));

		return delta_u;
	}
	
	public void selectDisplacements() {
		IVectorRO solved_u = this.solve();
		for (Element e : this.getLOE()){
			for (Node n : e.getLON()){
				IVectorRO dis = new Array1DVector(3);
				int[] dof = n.getDOFNumbers();
				for (int i = 0; i < 3; i++){
					if (dof[i] == -1){
						((Array1DVector) dis).set(i,0);
					}
					else ((Array1DVector) dis).set(i,solved_u.get(dof[i]));
				}
				Vector3D displacement = new Vector3D(dis);
				n.setDisplacement(displacement);
			}
		}
	}
	
	public void selectPlasticDisplacements(IVectorRO delta_u) {
		for (Node n : this.getLON()){
			
			Vector3D old_displacement = n.getDisplacement();
			if (old_displacement == null){
				old_displacement = new Vector3D (0.0, 0.0, 0.0);
			}
			IVectorRO new_displacement = new Array1DVector(3);
			int[] node_DOF = n.getDOFNumbers();
			
			for (int i = 0; i < 3; i++){
				if (node_DOF[i] == -1){
					((Array1DVector) new_displacement).set(i,0);
				}
				else {
					((Array1DVector) new_displacement).set(i,old_displacement.get(i) + delta_u.get(node_DOF[i]));
				}
			}
			Vector3D displacement = new Vector3D(new_displacement);
			n.setDisplacement(displacement);
			
		}
	}
	
	/*public void computeStrain(){
		for (Element e : this.listOfElements){
			((Brick) e).computeStrain();
		}
	}*/

	public void setDisplacement(){
		for (Element e : this.listOfElements){
			e.setDisplacement();
		}
	}
	
	public void computeBMatrix(ShapeFunctions sf){
		for (Element e : this.listOfElements){
			e.computeBMatrix(sf);
		}
	}
	
	public void printResults(){
		System.out.println(" Solving K u = r");
		System.out.println(" Matrix K");
		System.out.println(MatrixFormat.format(this.kGlobal));
		System.out.println(" Vector r");
		System.out.println(MatrixFormat.format(this.rGlobal));
		
		// print result
		System.out.println(" Solution u");
		System.out.println(MatrixFormat.format(uGlobal));
	}

	@Override
	public void addSpring(double k, Node n, IVectorRO dir) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public IVectorRO assembleExternalForces() {
		IVectorRO rExternal = new Array1DVector(this.enumerateDOFs());
		//from Node select DOF which is non -1, get the index, plus to the R_global with the same DoF
		for (Node n : listOfNodes){
			int[] nDOF = n.getDOFNumbers();
			for(int i = 0; i < nDOF.length; i++){
				if (nDOF[i] != -1 && n.getForce() != null){
					((Array1DVector) rExternal).add(nDOF[i], n.getForce().getComponent(i));
				}
			}
		}
		return rExternal;
	}

	@Override
	public IMatrixRO assemblePlasticStiffnessMatrix() {
		IMatrixRO kGlobal = new Array2DMatrix (this.enumerateDOFs(), this.enumerateDOFs()); 

		for (Element e : listOfElements){
			IMatrixRO kElement = e.computeStiffnessMatrix(this.sf);
			int[] eleDOF = e.getDOFNumbers();
			for (int i = 0; i < eleDOF.length; i++){
				for (int j = 0; j< eleDOF.length; j++){
					if (eleDOF[i] != -1 && eleDOF[j] != -1){
						((Array2DMatrix)kGlobal).add(eleDOF[i], eleDOF[j], kElement.get(i, j));
					}
				}
			}
		//this.kGlobal = kGlobal;
		}
		return kGlobal;
	}

	
}
