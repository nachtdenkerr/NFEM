package fem_1D;

public class Constraint {
	private boolean[] free = new boolean[3];
	
	public Constraint(boolean u1, boolean u2, boolean u3){
		this.free[0] = u1;
		this.free[1] = u2;
		this.free[2] = u3;
	}
	
	public void set(int i, boolean a){
		this.free[i] = a;
	}
	
	public boolean isFree(int c){
		if (this.free[c] == true){
			return true; 
		}
		else return false;
	}
	
	public void print(){
		int n;
		for (n=0; n<3; n++){
			if (isFree(n) == true) System.out.printf("%-15s","free");
			else System.out.printf("%-15s","fixed");
		}
		System.out.println();
	}
}
