package fem_1D;

import iceb.jnumerics.IVectorRO;

public class Spring {
	
	private double k;
	private Node n1, n2;
	private IVectorRO dir;
	
	public Spring (double k, Node n, IVectorRO dir2){
		this.k = k;
		this.n1 = n;
		this.dir = dir2;
	}
	
	public Spring (double k, Node n1, Node n2){
		this.k = k;
		this.n1 = n1; this.n2 = n2;
	}
	
	public double getK() {
		return k;
	}
	
	public Node getNode1() {
		return n1;
	}
	
	public Node getNode2() {
		return n2;
	}
	
	public IVectorRO getdirection(){
		return dir;
	}
}
