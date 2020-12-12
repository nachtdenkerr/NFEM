package fem_1D;

import plasticity.PlasticModel;
import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVectorRO;
import iceb.jnumerics.Vector3D;

public class GaussPoint {
	private Vector3D position;
	private IMatrixRO B,Bnl;
	private Material mat;
	private PlasticModel pm;
	private IVectorRO STRAT;
	
	public GaussPoint (double x1, double x2, double x3){
		this.position =new Vector3D(x1,x2,x3);
	}
	
	public void setMaterial(Material mat){
		this.mat = mat;
	}
	
	public Material getMaterial(){
		return this.mat;
	}
	
	public boolean hasMaterial(){
		if (this.mat != null) return true;
		else return false;
	}
	
	public void setPlasticModel(PlasticModel pm){
		this.pm = pm;
	}
	
	public PlasticModel getPlasticModel(){
		return this.pm;
	}
		
	public Vector3D getPosition(){
		return position;
	}
	
	public void setBMatrix(IMatrixRO a){
		this.B = a;
	}
	
	public IMatrixRO getBMatrix(){
		return B;
	}
	
	public void setBnlMatrix(IMatrixRO a){
		this.Bnl = a;
	}
	
	public IMatrixRO getBnlMatrix() {
		return Bnl;
	}
	
	public void setSTRAT(IVectorRO STRAT){
		this.STRAT = STRAT;
	}
	
	public IVectorRO getSTRAT(){
		return this.STRAT;
	}
}
