package fem_1D;

import iceb.jnumerics.Vector3D;

public class Node {
	private int[] dofNumbers = new int[3];
	private Vector3D position;
	private Constraint bc = new Constraint(true, true, true);
	private Force fo;
	private Vector3D disp;
	private int counter;
	private boolean lineload = false;
	
	public Node(double x1, double x2, double x3){
		this.position =new Vector3D(x1,x2,x3);
	}
	
	public void setCounter(int start){
		counter = start;
	}
	
	public int getCounter(){
		return counter;
	}
	
	public void setConstraint (Constraint c){
		bc = c;
	}
	
	public void setConstraint (boolean a, boolean b, boolean c){
		bc.set(0, a);
		bc.set(1, b);
		bc.set(2, c);
	}
	
	public Constraint getConstraint(){
		return bc;
	}
	
	public void setForce(Force f){
		fo = f;
	}
	
	public Force getForce(){
		return fo;
	}
	
	public void setLineLoad(boolean a){
		lineload = a;
	}
	
	public boolean hasLineLoad(){
		return lineload;
	}
	
	public int enumerateDOFs (int start){
		for (int n = 0; n < 3; n++){
			if (this.bc.isFree(n) == true){
				dofNumbers[n] = start;
				start++;
			}
			else if(this.bc.isFree(n) == false){
				dofNumbers[n] = -1;
			}
		}
		return start;
	}
	
	public int[] getDOFNumbers(){
		this.enumerateDOFs(counter);
		return dofNumbers;
	}
	
	public Vector3D getPosition(){
		return position;
	}

	public void setDisplacement(Vector3D u){
		disp = u;
	}
	
	public Vector3D getDisplacement(){
		return disp;
	}
	
	public void print(){
		System.out.println(position);
	}
}
