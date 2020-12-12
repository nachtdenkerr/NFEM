package fem_2D;

import fem_abstract.Element;
import iceb.jnumerics.IMatrix;
import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVector;
import iceb.jnumerics.IVectorRO;
import iceb.jnumerics.Vector3D;

import java.util.ArrayList;
import java.util.List;

import fem_1D.GaussPoint;
import fem_1D.Node;
import linalg.*;
import shape_functions.LinearSF;
import shape_functions.ShapeFunctions;

public class Quad extends Element{
	
	private double thickness,nuy,eModulus;
	private String constitutive;
	private IMatrix C_matrix = new Array2DMatrix(3,3);
	
	private Node node1, node2, node3, node4;
	private int[] dofNumbers = new int[12];
	private int start,dofs;
	
	private IMatrixRO position = new Array2DMatrix(4,2);
	private List<Node> LoN = new ArrayList <Node> ();
	
	private LinearSF sf;
	private double alpha1, alpha2;
	
	private int no_GP;
	private List<GaussPoint> LoGP = new ArrayList <GaussPoint>();
	
	private IMatrixRO B = new Array2DMatrix (3,12);
	private IVector displacement = new Array1DVector(12);
	

	public Quad (double e, double nuy, String a, double h, Node n1, Node n2, Node n3, Node n4) {
		this.node1 = n1; this.node2 = n2; this.node3 = n3; this.node4 = n4;
		this.thickness = h; this.nuy = nuy; this.eModulus = e; this.constitutive = a;
		
		// create C_matrix
		if (a == "plane strain"){ // plain strain
			for (int i = 0; i <2; i++){
				((Array2DMatrix) this.C_matrix).set(i, i, 1-nuy);
			}
			this.C_matrix.set(0, 1, nuy);
			this.C_matrix.set(1 ,0, nuy);
			this.C_matrix.set(2, 2, 0.5 - nuy);
			C_matrix = (IMatrix) this.C_matrix.multiply(e/((1+nuy)*(1-2*nuy)));
		}
		if (a == "plane stress"){ // plane stress
			for (int i = 0; i <2; i++){
				((Array2DMatrix) this.C_matrix).set(i, i, 1);
			}
			((Array2DMatrix) this.C_matrix).set(0, 1, nuy);
			((Array2DMatrix) this.C_matrix).set(1 ,0, nuy);
			((Array2DMatrix) this.C_matrix).set(2, 2, 0.5 - 0.5*nuy);
			C_matrix = (IMatrix) this.C_matrix.multiply(e/((1-nuy*nuy)));
		}
		
		this.sf = new LinearSF("Quad");
		this.LoN.add(n1);
		this.LoN.add(n2);
		this.LoN.add(n3);
		this.LoN.add(n4); 
		for (int i = 0; i<4; i++){ // only consider positions in direction 1 and 2
			((Array2DMatrix) position).set(i, 0, LoN.get(i).getPosition().c1);
			((Array2DMatrix) position).set(i, 1, LoN.get(i).getPosition().c2);
		}
	}
	
	public IMatrixRO getC(){
		return this.C_matrix;
	}
	
	public void enumerateDOFs(){
		start = node1.enumerateDOFs(start);
		start = node2.enumerateDOFs(start); 
		start = node3.enumerateDOFs(start);
		start = node4.enumerateDOFs(start);
		this.dofs = start;
	}
	
	public int[] getDOFNumbers(){
		System.arraycopy(node1.getDOFNumbers(), 0, this.dofNumbers, 0, 3);
		System.arraycopy(node2.getDOFNumbers(), 0, this.dofNumbers, 3, 3);
		System.arraycopy(node3.getDOFNumbers(), 0, this.dofNumbers, 6, 3);
		System.arraycopy(node4.getDOFNumbers(), 0, this.dofNumbers, 9, 3);
		return dofNumbers;
	}
	
	public List<Node> getLON(){
		return this.LoN;
	}
	
	public List<GaussPoint> getLOGP(){
		return this.LoGP;
	}
	
	public IMatrixRO getJacobian(GaussPoint g){
		IMatrixRO jacobian = new Array2DMatrix (2,2);
		jacobian = this.sf.getdSF(g).multiply(this.position);
		return jacobian;
	}
	
	public void assignGaussPoints(int number_of_GaussPoints){
		this.no_GP = number_of_GaussPoints;
		double alpha1 = 0; double alpha2 = 0;
		double a = 1/Math.sqrt(3);
		if (this.no_GP == 1){
			GaussPoint g1 = new GaussPoint (0,0,0);
			alpha1 = 2; alpha2 = 1;
			this.LoGP.add(g1);
		}
		if (this.no_GP == 2){
			GaussPoint g1 = new GaussPoint (-a,0,0);
			GaussPoint g2 = new GaussPoint ( a,0,0);
			alpha1 = 1; alpha2 = 2;
			this.LoGP.add(g1); this.LoGP.add(g2);
		}
		if (this.no_GP == 4){
			GaussPoint g1 = new GaussPoint (-a,-a,0);
			GaussPoint g2 = new GaussPoint ( a,-a,0);
			GaussPoint g3 = new GaussPoint ( a, a,0);
			GaussPoint g4 = new GaussPoint (-a, a,0);
			alpha1 = 1; alpha2 = 1;
			this.LoGP.add(g1); this.LoGP.add(g2); this.LoGP.add(g3); this.LoGP.add(g4);
		}
		this.alpha1 = alpha1; this.alpha2 = alpha2;
		
		for (GaussPoint gp : LoGP){
			gp.getMaterial().setTangent(C_matrix);
		}
	}
	
	public void assignPlasticGaussPoints(int number_of_GaussPoints){
		this.no_GP = number_of_GaussPoints;
		double alpha1 = 0; double alpha2 = 0;
		double a = 1/Math.sqrt(3);
		if (this.no_GP == 1){
			GaussPoint g1 = new GaussPoint (0,0,0);
			alpha1 = 2; alpha2 = 1;
			this.LoGP.add(g1);
		}
		if (this.no_GP == 2){
			GaussPoint g1 = new GaussPoint (-a,0,0);
			GaussPoint g2 = new GaussPoint ( a,0,0);
			alpha1 = 1; alpha2 = 2;
			this.LoGP.add(g1); this.LoGP.add(g2);
		}
		if (this.no_GP == 4){
			GaussPoint g1 = new GaussPoint (-a,-a,0);
			GaussPoint g2 = new GaussPoint ( a,-a,0);
			GaussPoint g3 = new GaussPoint ( a, a,0);
			GaussPoint g4 = new GaussPoint (-a, a,0);
			alpha1 = 1; alpha2 = 1;
			this.LoGP.add(g1); this.LoGP.add(g2); this.LoGP.add(g3); this.LoGP.add(g4);
		}
		this.alpha1 = alpha1; this.alpha2 = alpha2;
		
		this.computeB_Plastic_Matrix();

		
//		gp.getMaterial().setTangent(C);
		
		/*for (GaussPoint gp : LoGP){
			double E = this.getEModulus();
			double K = E/(3*(1 - 2*this.getNuy()));
			double G = E/(2*(1 + this.getNuy()));
			
			IVectorRO I = new Array1DVector(1.0, 1.0, 0.0, 1.0);

			IMatrixRO II = new Array2DMatrix(4,4);
			((Array2DMatrix) II).set(0, 0, 1.0);
			((Array2DMatrix) II).set(1, 1, 1.0);
			((Array2DMatrix) II).set(2, 2, 0.5);
			((Array2DMatrix) II).set(3, 3, 1.0);
			
			IMatrixRO Ptens = new Array2DMatrix(4,4);
			for (int i = 0; i <4; i++){
				for (int b = 0; b <4; b++){
					((Array2DMatrix) Ptens).set(i,b,II.get(i, b) - 1.0/3.0*I.get(i)*I.get(b));
				}
			}
			
			IMatrixRO C = new Array2DMatrix(4,4);
			for (int m = 0; m <4; m++){
				for (int n = 0; n <4; n++){
					((Array2DMatrix) C).set(m,n,2*G*Ptens.get(m, n) + K*I.get(m)*I.get(n));
				}
			}
			for (int n =0; n < 3; n++){
				for (int m = n +1; m <4; m++){
					((Array2DMatrix) C).set(m,n, C.get(n, m));
				}
			}

			Material mat = new Material();
			gp.setMaterial(mat);
			gp.getMaterial().setTangent(C);
			this.computeB_Plastic_Matrix();
		}*/
	}
	
	public void computeBMatrix() {
		for (GaussPoint g : this.LoGP){
			IMatrixRO j = this.getJacobian(g);
			IMatrixRO j_inv = ((Array2DMatrix) j).inverse2x2();
			IMatrixRO dN_dX = j_inv.multiply(this.sf.getdSF(g));
			IMatrixRO B = new Array2DMatrix(3,this.getDOFNumbers().length);
			for (int i = 0; i < 12; i = i+3){
				//set B
				((Array2DMatrix) B).set(0, i, dN_dX.get(0, i/3));
				((Array2DMatrix) B).set(1, i+1, dN_dX.get(1, i/3));
				((Array2DMatrix) B).set(2, i, dN_dX.get(1, i/3));
				((Array2DMatrix) B).set(2, i+1, dN_dX.get(0, i/3));
			}
			g.setBMatrix(B);
		}
	}
	
	public IMatrixRO computeStiffnessMatrix(){ // before calling compute, we need to assign GPs
		this.enumerateDOFs();
		IMatrixRO kElement = new Array2DMatrix(12,12);
		this.computeBMatrix();
		
		for (GaussPoint g : this.LoGP){
			IMatrixRO j = this.getJacobian(g);
			double det_J = ((Array2DMatrix) j).determinant();
			IMatrixRO B = g.getBMatrix();
			
			//compute K  alpha1*alpha2*B*C*B*detJ*h
			IMatrixRO kGP = new Array2DMatrix(this.dofs,this.dofs);
			double al = this.alpha1*this.alpha2*det_J;
			kGP = (B.transpose().multiply((C_matrix.multiply(B)))).multiply(al*this.thickness);

			//add KgaussPoints to Kelement
			kElement = kElement.add(kGP);
		}
		return kElement;
	}
	
	public void computeB_Plastic_Matrix() {
		for (GaussPoint g : this.LoGP){
			IMatrixRO j = this.getJacobian(g);
			IMatrixRO j_inv = ((Array2DMatrix) j).inverse2x2();
			IMatrixRO dN_dX = j_inv.multiply(this.sf.getdSF(g));
			IMatrixRO B = new Array2DMatrix( 4, this.getDOFNumbers().length);
			for (int i = 0; i < 12; i = i+3){
				//set B
				((Array2DMatrix) B).set(0, i  , dN_dX.get(0, i/3));
				((Array2DMatrix) B).set(1, i+1, dN_dX.get(1, i/3));
				((Array2DMatrix) B).set(2, i  , dN_dX.get(1, i/3));
				((Array2DMatrix) B).set(2, i+1, dN_dX.get(0, i/3));
				((Array2DMatrix) B).set(3, i  , 0.0);
			}
			g.setBMatrix(B);
		}
	}
	
	public IMatrixRO computePlasticStiffnessMatrix(){ // before calling compute, we need to assign GPs
		this.enumerateDOFs();
		IMatrixRO kElement = new Array2DMatrix(12,12);
		
		for (GaussPoint g : this.LoGP){
			IMatrixRO j = this.getJacobian(g);
			double det_J = ((Array2DMatrix) j).determinant();
			
			IMatrixRO B = g.getBMatrix();
			IMatrixRO C = g.getMaterial().getTangent();

			//compute K  alpha1*alpha2*B*C*B*detJ*h
			IMatrixRO kGP = new Array2DMatrix(this.dofs,this.dofs);
			double al = this.alpha1*this.alpha2*det_J;
			kGP = (B.transpose().multiply((C.multiply(B)))).multiply(al*this.thickness);

			//add KgaussPoints to Kelement
			kElement = kElement.add(kGP);
		}
		return kElement;
	}
	
	public IVectorRO computeInternalForce () {
		IVectorRO r_int = new Array1DVector(12);
		//System.out.println("Sigma: ");

		for (GaussPoint g: this.LoGP){
			IMatrixRO j = this.getJacobian(g);
			double det_J = ((Array2DMatrix) j).determinant();
			
			IVectorRO sigma = g.getMaterial().getSigma();
			//System.out.println(MatrixFormat.format(sigma));
			IMatrixRO B = g.getBMatrix();
			IVectorRO rGP = B.transpose().multiply(sigma).multiply(det_J*this.getThickness());
			r_int = r_int.add(rGP);
//			System.out.println("sigma\n" + MatrixFormat.format(sigma));

		}
		//System.out.println();
		//System.out.println("B\n" + MatrixFormat.format(B));
		//System.out.println();

		return r_int;
	}
	
	public IVectorRO computeStrain(){
		IVectorRO eps = new Array1DVector(3);
		
		IVector u = new Array1DVector(12);
		for (Node n : LoN){
			Vector3D u0 = new Vector3D (0,0,0);
			if (n.getDisplacement() == null){
				n.setDisplacement(u0);
			}
		}
		for (int i = 0; i < 12; i = i+3){
			u.set(i , this.getLON().get(i/3).getDisplacement().getX1());
			u.set(i+1,this.getLON().get(i/3).getDisplacement().getX2());
			u.set(i+2,this.getLON().get(i/3).getDisplacement().getX3());
		}
		this.displacement = u;
		eps = this.B.multiply(u);
		return eps;
	}
	
	public double getEModulus(){
		return eModulus;
	}
	
	public double getNuy(){
		return nuy;
	}
	
	public double getThickness(){
		return thickness;
	}
	
	public String getConstitutive(){
		return constitutive;
	}
	
	public int getStart(){
		return start;
	}
	
	public IVectorRO getDisplacement(){
		return displacement;
	}
	
	public void setCMatrix(IMatrixRO c){
		this.C_matrix = (IMatrix) c;
	}

	public void setDisplacement() {
		IVector u = new Array1DVector(12);
		for (Node n : LoN){
			if (n.getDisplacement() == null){
				n.setDisplacement(new Vector3D (0.0, 0.0, 0.0));
			}
		}
		for (int i = 0; i < 12; i = i+3){
			u.set(i , this.getLON().get(i/3).getDisplacement().getX1());
			u.set(i+1,this.getLON().get(i/3).getDisplacement().getX2());
			u.set(i+2,this.getLON().get(i/3).getDisplacement().getX3());
		}
		this.displacement = u;
	}

	@Override
	public IMatrixRO computeStiffnessMatrix(ShapeFunctions sf) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public IVectorRO computeInternalForce(ShapeFunctions sf) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void computeBMatrix(ShapeFunctions sf) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double getForce() {
		// TODO Auto-generated method stub
		return 0;
	}
}