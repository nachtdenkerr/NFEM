package fem_1D;

import fem_abstract.Load;
import iceb.jnumerics.IVector;
import iceb.jnumerics.IVectorRO;
import inf.text.ArrayFormat;
import linalg.Array1DVector;

public class LineLoad extends Load{
	
	private double q;
	public Node n1, n2;
	private String dir;
	private double[] lineload;
	private IVectorRO direction;
	
	public LineLoad (double q, Node n1, Node n2, String a){
		//String is for up or down (out or into the element)
		this.q = q;
		this.n1 = n1;
		this.n2 = n2;
		this.dir = a;
	}
	
	public double [] computeLineLoad() {
		IVectorRO a = n1.getPosition().toVector();
		IVectorRO b = n2.getPosition().toVector();
		IVectorRO y = new Array1DVector(0, 1, 0);
		IVectorRO direction = b.subtract(a);
		this.direction = direction;
		
		double qlc2 = q*direction.normTwo()/2;
		direction = direction.multiply(1/direction.normTwo());

		IVector perpendicular = new Array1DVector(0,0,0);
		perpendicular.set(0,  direction.get(1));
		perpendicular.set(1, -direction.get(0));
		
		double dot = ((IVectorRO) perpendicular).dot(y);
		if ((dir == "up" && dot < 0) || (dir == "down" && dot > 0) ){
			perpendicular.set(0, -direction.get(1));
			perpendicular.set(1, direction.get(0));
		}
		
		double[] lineload = new double[3];
		lineload[0] = qlc2*perpendicular.get(0)/perpendicular.normTwo();
		lineload[1] = qlc2*perpendicular.get(1)/perpendicular.normTwo();
		lineload[2] = 0;
		this.lineload = lineload;
		return lineload;
	}
	
	public double getComponent (int c){
		return lineload[c];
	}
	
	public Node getNode (Node n){
		return n;
	}
	
	public IVectorRO getDirection (){
		return direction;
	}
	
	public void print(){
		System.out.println(ArrayFormat.format(this.lineload));
		System.out.println(ArrayFormat.format(this.lineload));
	}	

}
