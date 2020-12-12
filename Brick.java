package fem_3D;

import fem_abstract.Element;
import iceb.jnumerics.IMatrix;
import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVector;
import iceb.jnumerics.IVectorRO;
import iceb.jnumerics.MatrixFormat;
import iceb.jnumerics.Vector3D;

import java.util.ArrayList;
import java.util.List;

import fem_1D.GaussPoint;
import fem_1D.Node;
import linalg.Array1DVector;
import linalg.Array2DMatrix;
import shape_functions.ShapeFunctions;

public class Brick extends Element{

	private double nuy,eModulus;
	private IMatrix C_matrix = new Array2DMatrix(6,6);
	
	private Node node1, node2, node3, node4, node5, node6, node7, node8;
	private List<Node> LoN = new ArrayList <Node> ();

	private int[] dofNumbers = new int[24];
	private int start,dofs;
	
	private double alpha1, alpha2, alpha3;
	
	private IMatrixRO position = new Array2DMatrix(8,3);
	private List<GaussPoint> LoGP = new ArrayList <GaussPoint>();
	
	//private IMatrixRO B = new Array2DMatrix (6, 24);
	private IVector displacement = new Array1DVector(24);

	
	public Brick(double e, double nuy, List<Node> list){
		this.nuy = nuy; this.eModulus = e;
		this.node1 = list.get(0); this.node2 = list.get(1); this.node3 = list.get(2); this.node4 = list.get(3);
		this.node5 = list.get(4); this.node6 = list.get(5); this.node7 = list.get(6); this.node8 = list.get(7);
		LoN.add(node1); LoN.add(node2); LoN.add(node3); LoN.add(node4);
		LoN.add(node5); LoN.add(node6); LoN.add(node7); LoN.add(node8);
		
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
		
		for (int i = 0; i<8; i++){
			((Array2DMatrix) position).set(i, 0, LoN.get(i).getPosition().c1);
			((Array2DMatrix) position).set(i, 1, LoN.get(i).getPosition().c2);
			((Array2DMatrix) position).set(i, 2, LoN.get(i).getPosition().c3);
		}
	}
	
	@Override
	public double getNuy() {
		return this.nuy;
	}
	
	
	public IMatrixRO getJacobian(GaussPoint g, ShapeFunctions sf){
		IMatrixRO jacobian = new Array2DMatrix (3,3);
		jacobian = sf.getdSF(g).multiply(this.position);
		
		return jacobian;
	}
	
	@Override
	public void enumerateDOFs() {
		start = node1.enumerateDOFs(start);
		start = node2.enumerateDOFs(start); 
		start = node3.enumerateDOFs(start);
		start = node4.enumerateDOFs(start);
		start = node5.enumerateDOFs(start);
		start = node6.enumerateDOFs(start); 
		start = node7.enumerateDOFs(start);
		start = node8.enumerateDOFs(start);
		this.dofs = start;
	}
	
	public int[] getDOFNumbers() {
		System.arraycopy(node1.getDOFNumbers(), 0, this.dofNumbers, 0, 3);
		System.arraycopy(node2.getDOFNumbers(), 0, this.dofNumbers, 3, 3);
		System.arraycopy(node3.getDOFNumbers(), 0, this.dofNumbers, 6, 3);
		System.arraycopy(node4.getDOFNumbers(), 0, this.dofNumbers, 9, 3);
		System.arraycopy(node5.getDOFNumbers(), 0, this.dofNumbers, 12, 3);
		System.arraycopy(node6.getDOFNumbers(), 0, this.dofNumbers, 15, 3);
		System.arraycopy(node7.getDOFNumbers(), 0, this.dofNumbers, 18, 3);
		System.arraycopy(node8.getDOFNumbers(), 0, this.dofNumbers, 21, 3);
		return dofNumbers;
	}
	
	public void assignGaussPoints(int ngp1, int ngp2, int ngp3, ShapeFunctions sf){
		
		double alpha1 = 0; double alpha2 = 0; double alpha3 = 0;
		if (ngp1 + ngp2 + ngp3 == 3){
			GaussPoint g1 = new GaussPoint (0,0,0);
			alpha1 = 2; alpha2 = 2; alpha3 = 2;
			this.LoGP.add(g1);
		}
		double a = 1/Math.sqrt(3);
		if (ngp1 + ngp2 + ngp3 == 4){
			if (ngp1 == 2){
				GaussPoint g1 = new GaussPoint (-a, 0, 0);
				GaussPoint g2 = new GaussPoint ( a, 0, 0);
				alpha1 = 1; alpha2 = 2; alpha3 = 2;
				this.LoGP.add(g1); this.LoGP.add(g2);
			}
			if (ngp2 == 2){
				GaussPoint g1 = new GaussPoint (0,-a, 0);
				GaussPoint g2 = new GaussPoint (0, a, 0);
				alpha1 = 2; alpha2 = 1; alpha3 = 2;
				this.LoGP.add(g1); this.LoGP.add(g2);
			}
			if (ngp3 == 2){
				GaussPoint g1 = new GaussPoint (0, 0,-a);
				GaussPoint g2 = new GaussPoint (0, 0, a);
				alpha1 = 2; alpha2 = 2; alpha3 = 1;
				this.LoGP.add(g1); this.LoGP.add(g2);
			}
		}
		
		if (ngp1 + ngp2 + ngp3 == 5){
			if (ngp1 == 1){
				GaussPoint g1 = new GaussPoint (0,-a,-a);
				GaussPoint g2 = new GaussPoint (0, a,-a);
				GaussPoint g3 = new GaussPoint (0, a, a);
				GaussPoint g4 = new GaussPoint (0,-a, a);
				alpha1 = 2; alpha2 = 1; alpha3 = 1;
				this.LoGP.add(g1); this.LoGP.add(g2); this.LoGP.add(g3); this.LoGP.add(g4);
			}
			
			if (ngp2 == 1){
				GaussPoint g1 = new GaussPoint (-a, 0,-a);
				GaussPoint g2 = new GaussPoint ( a, 0,-a);
				GaussPoint g3 = new GaussPoint ( a, 0, a);
				GaussPoint g4 = new GaussPoint (-a, 0, a);
				alpha1 = 1; alpha2 = 2; alpha3 = 1;
				this.LoGP.add(g1); this.LoGP.add(g2); this.LoGP.add(g3); this.LoGP.add(g4);
			}
			
			if (ngp3 == 1){
				GaussPoint g1 = new GaussPoint (-a,-a,0);
				GaussPoint g2 = new GaussPoint ( a,-a,0);
				GaussPoint g3 = new GaussPoint ( a, a,0);
				GaussPoint g4 = new GaussPoint (-a, a,0);
				alpha1 = 1; alpha2 = 1; alpha3 = 2;
				this.LoGP.add(g1); this.LoGP.add(g2); this.LoGP.add(g3); this.LoGP.add(g4);
			}
		}
		
		if (ngp1 + ngp2 + ngp3 == 6){
			GaussPoint g1 = new GaussPoint (-a,-a,-a);
			GaussPoint g2 = new GaussPoint ( a,-a,-a);
			GaussPoint g3 = new GaussPoint ( a, a,-a);
			GaussPoint g4 = new GaussPoint (-a, a,-a);
			GaussPoint g5 = new GaussPoint (-a,-a, a);
			GaussPoint g6 = new GaussPoint ( a,-a, a);
			GaussPoint g7 = new GaussPoint ( a, a, a);
			GaussPoint g8 = new GaussPoint (-a, a, a);
			alpha1 = 1; alpha2 = 1; alpha3 = 1;
			this.LoGP.add(g1); this.LoGP.add(g2); this.LoGP.add(g3); this.LoGP.add(g4);
			this.LoGP.add(g5); this.LoGP.add(g6); this.LoGP.add(g7); this.LoGP.add(g8);
		}
		
		this.alpha1 = alpha1;
		this.alpha2 = alpha2;
		this.alpha3 = alpha3;
		
	}
	
	@Override
	public void computeBMatrix(ShapeFunctions sf){
		this.enumerateDOFs();		
		for (GaussPoint g : this.LoGP){
			IMatrixRO j = this.getJacobian(g, sf);
			IMatrixRO j_inv = ((Array2DMatrix) j).inverse3x3();
			IMatrixRO dN_dX = j_inv.multiply(sf.getdSF(g));
			IMatrixRO B = new Array2DMatrix(6,this.getDOFNumbers().length);
			for (int i = 0; i < 24; i = i+3){
				//set B
				((Array2DMatrix) B).set(0, i,   dN_dX.get(0, i/3));
				((Array2DMatrix) B).set(1, i+1, dN_dX.get(1, i/3));
				((Array2DMatrix) B).set(2, i+2, dN_dX.get(2, i/3));
				
				((Array2DMatrix) B).set(3, i,   dN_dX.get(1, i/3));
				((Array2DMatrix) B).set(3, i+1, dN_dX.get(0, i/3));
				((Array2DMatrix) B).set(4, i+1, dN_dX.get(2, i/3));
				((Array2DMatrix) B).set(4, i+2, dN_dX.get(1, i/3));
				((Array2DMatrix) B).set(5, i,   dN_dX.get(2, i/3));
				((Array2DMatrix) B).set(5, i+2, dN_dX.get(0, i/3));
			}
			g.setBMatrix(B);
		}
	}
	
	@Override
	public IMatrixRO computeStiffnessMatrix(ShapeFunctions sf) {
		this.enumerateDOFs();
		this.computeBMatrix(sf);
		IMatrixRO kElement = new Array2DMatrix(24,24);

		for (GaussPoint g : this.LoGP){
			IMatrixRO j = this.getJacobian(g, sf);
			double det_J = ((Array2DMatrix) j).determinant();
			
			IMatrixRO B = g.getBMatrix();

			//compute K  alpha1*alpha2*B*C*B*detJ*h
			IMatrixRO kGP = new Array2DMatrix(this.dofs,this.dofs);
			double al = this.alpha1*this.alpha2*this.alpha3*det_J;
			kGP = (B.transpose().multiply((g.getMaterial().getTangent().multiply(B)))).multiply(al);

			//add KgaussPoints to Kelement
			kElement = kElement.add(kGP);
		}
		return kElement;
	}
	
	@Override
	public IVectorRO computeInternalForce (ShapeFunctions sf) {
		IVectorRO r_int = new Array1DVector(24);
		//System.out.println("Sigma: ");

		for (GaussPoint g: this.LoGP){
			IMatrixRO j = this.getJacobian(g, sf);
			double det_J = ((Array2DMatrix) j).determinant();
			
			IVectorRO sigma = g.getMaterial().getSigma();
			IMatrixRO B = g.getBMatrix();
			IVectorRO rGP = B.transpose().multiply(sigma).multiply(det_J);
			r_int = r_int.add(rGP);
		}
		//System.out.println();
		//System.out.println("r_int element " + MatrixFormat.format(r_int));
		//System.out.println();

		return r_int;
	}
	
	public void computeStrain(){
		IVector u = new Array1DVector(24);
		for (Node n : LoN){
			Vector3D u0 = new Vector3D (0,0,0);
			if (n.getDisplacement() == null){
				n.setDisplacement(u0);
			}
		}
		for (int i = 0; i < 24; i = i+3){
			u.set(i , this.getLON().get(i/3).getDisplacement().getX1());
			u.set(i+1,this.getLON().get(i/3).getDisplacement().getX2());
			u.set(i+2,this.getLON().get(i/3).getDisplacement().getX3());
		}
		this.displacement = u;
		for (GaussPoint gp : LoGP){
			IVectorRO eps = gp.getBMatrix().multiply(u);
			gp.getMaterial().setEpsElastic(eps);
		}
	}
	
	public void computeStress(){
		for (GaussPoint gp : LoGP){
			IVectorRO sigma = C_matrix.multiply(gp.getMaterial().getEpsElastic());
			System.out.println("Stress " + MatrixFormat.format(sigma));
		}
	}
	
	@Override
	public double getEModulus() {
		return eModulus;
	}
	
	public IMatrixRO getC() {
		return C_matrix;
	}
	
	@Override
	public List<Node> getLON() {
		return LoN;
	}
	
	@Override
	public int getStart() {
		return start;
	}

	@Override
	public List<GaussPoint> getLOGP() {
		return this.LoGP;
	}
	
	@Override
	public void setDisplacement(){
		IVector u = new Array1DVector(24);
		for (Node n : LoN){
			if (n.getDisplacement() == null){
				n.setDisplacement(new Vector3D (0,0,0));
			}
		}
		for (int i = 0; i < 24; i = i+3){
			u.set(i , this.getLON().get(i/3).getDisplacement().getX1());
			u.set(i+1,this.getLON().get(i/3).getDisplacement().getX2());
			u.set(i+2,this.getLON().get(i/3).getDisplacement().getX3());
		}
		this.displacement = u;
	}

	@Override
	public IVectorRO getDisplacement() {
		return displacement;
	}

	public void setCMatrix(IMatrixRO c){
		this.C_matrix = (IMatrix) c;
	}

	@Override
	public IMatrixRO computeStiffnessMatrix() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double getForce() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public IVectorRO computeInternalForce() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void assignPlasticGaussPoints(int i) {
		// TODO Auto-generated method stub
		
	}
}
