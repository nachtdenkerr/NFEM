package fem_2D;

import iceb.jnumerics.IMatrixRO;
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

public class QuadStructure extends Structure{
	public List<Node> listOfNodes = new ArrayList<Node>();
	public List<Element> listOfElements = new ArrayList<Element>();
	public List<LineLoad> listOfLLoads = new ArrayList<LineLoad>();
	public List<Spring> listOfSprings = new ArrayList<Spring>();

	private IMatrixRO kGlobal = new Array2DMatrix(this.enumerateDOFs(),this.enumerateDOFs());
	private IVectorRO rGlobal = new Array1DVector(this.enumerateDOFs());
	private IVectorRO uGlobal = new Array1DVector(this.getNumberOfNodes() * 3);
	
	private ShapeFunctions sf;
	
	public Node addNode(double x1, double x2, double x3){
		Node n = new Node(x1,x2,x3);
		listOfNodes.add(n);
		return n;
	}
	
	public Element addElement(double e, double n, String c, double h, int n1, int n2, int n3, int n4){
		Node node1 = listOfNodes.get(n1);
		Node node2 = listOfNodes.get(n2);
		Node node3 = listOfNodes.get(n3);
		Node node4 = listOfNodes.get(n4);
		Quad ele = new Quad(e,n,c,h,node1,node2,node3,node4);
		listOfElements.add(ele);
		return ele;
	}
	
	public void addElement(Quad q){
		listOfElements.add(q);
		Node node1 = q.getLON().get(0);
		Node node2 = q.getLON().get(1);
		Node node3 = q.getLON().get(2);
		Node node4 = q.getLON().get(3);
		listOfNodes.add(node1);
		listOfNodes.add(node2);
		listOfNodes.add(node3);
		listOfNodes.add(node4);
	}
	
	public void setLineLoad (double q, Node n1, Node n2, String a){
		LineLoad ll = new LineLoad (q,n1,n2,a);
		listOfLLoads.add(ll);
		double[] load = ll.computeLineLoad();
		Force f = new Force (load[0],load[1],load[2]);
		n1.setForce(f);
		n2.setForce(f);
	}
	
	@Override
	public void addSpring(double k, Node n, IVectorRO dir) {
		Spring spr = new Spring (k,n,dir);
		listOfSprings.add(spr);
		
		IMatrixRO kSpring = new Array2DMatrix(3,3);
		IVectorRO c = new Array1DVector(dir.get(0)/dir.normTwo(),
				dir.get(1)/dir.normTwo(),dir.get(2)/dir.normTwo());
		kSpring = c.dyadicProduct(c);
		
		int[] ndof = n.getDOFNumbers();
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				if (ndof[i] != -1 && ndof[j] != -1){
					((Array2DMatrix)kGlobal).add(ndof[i], ndof[j], kSpring.get(i, j));
				}
			}
		}
	}
	
	public void computeB_Plastic_Matrix(){
		for (Element q : this.listOfElements){
			((Quad) q).computeB_Plastic_Matrix();
		}
	}
	
	@Override
	public IMatrixRO getC() {
		return this.getElement(0).getC();
	}
	
	public List<Node> getLON(){
		return listOfNodes;
	}
	
	public int getNumberOfNodes(){
		return listOfNodes.size();
	}
	
	public Node getNode(int id){
		return listOfNodes.get(id);
	}
	
	public List<Element> getLOE(){
		return listOfElements;
	}
	
	public int getNumberOfElements(){
		return listOfElements.size();
	}
	
	public Element getElement(int id){
		return listOfElements.get(id);
	}
	
	public void assignStartIndex(){
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
	
	public int enumerateDOFs(){
		this.assignStartIndex();
		int start = 0;
		int nDOF = 0;
		for (int i = 0; i < this.getLON().size(); i++){
			int a = this.getLON().get(i).enumerateDOFs(start);
			nDOF = nDOF + a;
		} 
		return nDOF;
	}
	
	public IMatrixRO assembleStiffnessMatrix(){
		int n = this.enumerateDOFs();
		IMatrixRO kGlobal = new Array2DMatrix (n,n); 
		//from element select DoF which is non -1, get the index, plus to the K_global with the same DoF
		
		for (Element e : listOfElements){
			IMatrixRO kElement = e.computeStiffnessMatrix();
			int[] eleDOF = e.getDOFNumbers();

			for (int i = 0; i < eleDOF.length; i++){
				for (int j = 0; j< eleDOF.length; j++){
					if (eleDOF[i] != -1 && eleDOF[j] != -1){
						((Array2DMatrix)kGlobal).add(eleDOF[i], eleDOF[j], kElement.get(i, j));
					}
				}
			}
		}
		this.kGlobal = kGlobal;
		return this.kGlobal;
	}
	
	public IMatrixRO assemblePlasticStiffnessMatrix(){
		int n = this.enumerateDOFs();
		IMatrixRO kGlobal = new Array2DMatrix (n,n); 
		//from element select DoF which is non -1, get the index, plus to the K_global with the same DoF
		
		for (Element e : listOfElements){
			IMatrixRO kElement = ((Quad) e).computePlasticStiffnessMatrix();
			int[] eleDOF = e.getDOFNumbers();

			for (int i = 0; i < eleDOF.length; i++){
				for (int j = 0; j< eleDOF.length; j++){
					if (eleDOF[i] != -1 && eleDOF[j] != -1){
						((Array2DMatrix)kGlobal).add(eleDOF[i], eleDOF[j], kElement.get(i, j));
					}
				}
			}
		}
		this.kGlobal = kGlobal;
		return this.kGlobal;
	}
	
	public IVectorRO assembleInternalForces () {
		IVectorRO riGlobal = new Array1DVector(this.enumerateDOFs());
		//from Node select DOF which is non -1, get the index, plus to the R_global with the same DoF
		
		for (Element e : listOfElements){
			IVectorRO riElement = e.computeInternalForce();
			int[] eleDOF = e.getDOFNumbers();
			for (int i = 0; i < eleDOF.length; i++){
				if (eleDOF[i] != -1){
					((Array1DVector) riGlobal).add(eleDOF[i], riElement.get(i));
				}
			}
			//System.out.println("riElement " + MatrixFormat.format(riElement));

		}

		return riGlobal;
	}
	
	public IVectorRO solve(){
		IMatrixRO K = this.assembleStiffnessMatrix();
		IVectorRO r = this.assembleExternalForces();
		LUDecomposition lud = new LUDecomposition(K);
		IVectorRO u = lud.solveFor(r);
		uGlobal = u;
		return u;
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
	
	public void selectPlasticDisplacements() {
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
	
	public void computeStrain(){
		for (Element e : this.listOfElements){
			((Quad) e).computeStrain();
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
		System.out.println();

	}

	public void setDisplacement() {
		for (Element e : this.listOfElements){
			e.setDisplacement();
		}
	}

	@Override
	public void addElement(Element e) {
		this.listOfElements.add(e);
	}

	@Override
	public void setSF(ShapeFunctions sf) {
		this.sf = sf;
	}
	
	public ShapeFunctions getSF(){
		return this.sf;
	}

	@Override
	public List<LineLoad> getLOLL() {
		return this.listOfLLoads;
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

	public void assignPlasticGaussPoints(int i) {
		for (Element e : this.listOfElements){
			e.assignPlasticGaussPoints(i);
		}
	}

}